#!/usr/bin/env python3
"""Watch matched processes and kill them if their RSS exceeds a memory fraction.

Default: if the sum of RSS of matched processes exceeds **40% of total RAM**
(``MemTotal``), send SIGTERM (then SIGKILL) to those processes.

Examples::

    # Watch GENIE syst workers (default match), 40% of total RAM, poll every 10s
    python3 mem_watchdog.py

    # Stricter, log to a file, run in tmux
    python3 mem_watchdog.py --threshold-pct 40 --interval 10 \\
        --match 'syst_genie_parallel|get_systematics_genie' \\
        --log /tmp/mem_watchdog.log

    # Dry-run (report only)
    python3 mem_watchdog.py --dry-run --once

Environment overrides (CLI wins when set):
  MEM_WATCHDOG_THRESHOLD_PCT   default 40
  MEM_WATCHDOG_INTERVAL_SEC    default 10
  MEM_WATCHDOG_MATCH           default regex (see --match help)
"""
from __future__ import annotations

import argparse
import os
import re
import signal
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_MATCH = (
    r"syst_genie_parallel|get_systematics_genie|run_syst_genie_chunked"
    r"|syst_multisim_(?:chunk|parallel|aggregate)"
    r"|syst_cosmics_(?:chunk|aggregate)"
    r"|syst_detvar_(?:chunk|aggregate)"
    r"|syst_histcounts_from_df"
    r"|run_chi2_subset_detector|run_chi2_track_subset"
    r"|chi2_track_subset"
)


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def log(msg: str, log_path: Optional[str] = None) -> None:
    line = f"[{_ts()}] {msg}"
    print(line, flush=True)
    if log_path:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def read_meminfo() -> Dict[str, int]:
    """Return selected /proc/meminfo fields in **bytes**."""
    out: Dict[str, int] = {}
    with open("/proc/meminfo", "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            key = parts[0].rstrip(":")
            # values are kB
            out[key] = int(parts[1]) * 1024
    return out


def basis_bytes(mem: Dict[str, int], basis: str) -> int:
    if basis == "total":
        return int(mem["MemTotal"])
    if basis == "available":
        # Fall back to MemFree+Buffers+Cached if MemAvailable missing
        if "MemAvailable" in mem:
            return int(mem["MemAvailable"])
        return int(mem.get("MemFree", 0) + mem.get("Buffers", 0) + mem.get("Cached", 0))
    raise ValueError(f"unknown basis {basis!r}")


def iter_user_pids(uid: int) -> Iterable[int]:
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            st = entry.joinpath("status").read_text(encoding="utf-8", errors="replace")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        for line in st.splitlines():
            if line.startswith("Uid:"):
                real_uid = int(line.split()[1])
                if real_uid == uid:
                    yield int(entry.name)
                break


def pid_cmdline(pid: int) -> str:
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return ""
    return raw.replace(b"\x00", b" ").decode("utf-8", errors="replace").strip()


def pid_rss_bytes(pid: int) -> Optional[int]:
    try:
        # /proc/pid/statm: size resident shared ... (pages)
        parts = Path(f"/proc/{pid}/statm").read_text(encoding="utf-8").split()
        resident_pages = int(parts[1])
        return resident_pages * os.sysconf("SC_PAGE_SIZE")
    except (FileNotFoundError, PermissionError, ProcessLookupError, IndexError, ValueError):
        return None


def collect_matched(
    uid: int,
    pattern: re.Pattern[str],
    exclude_self: bool = True,
) -> List[Tuple[int, int, str]]:
    """Return list of ``(pid, rss_bytes, cmdline)`` matching *pattern*."""
    self_pid = os.getpid()
    hits: List[Tuple[int, int, str]] = []
    for pid in iter_user_pids(uid):
        if exclude_self and pid == self_pid:
            continue
        cmd = pid_cmdline(pid)
        if not cmd or not pattern.search(cmd):
            continue
        rss = pid_rss_bytes(pid)
        if rss is None:
            continue
        hits.append((pid, rss, cmd))
    return hits


def fmt_bytes(n: int) -> str:
    x = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if x < 1024.0 or unit == "TiB":
            return f"{x:.1f} {unit}"
        x /= 1024.0
    return f"{n} B"


def terminate_pids(
    pids: Sequence[int],
    *,
    dry_run: bool,
    grace_sec: float,
    log_path: Optional[str],
) -> None:
    if not pids:
        return
    uniq = sorted(set(pids))
    log(f"Shutting down {len(uniq)} process(es): {uniq}", log_path)
    if dry_run:
        log("DRY-RUN: would send SIGTERM then SIGKILL", log_path)
        return

    for pid in uniq:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except PermissionError as exc:
            log(f"SIGTERM denied for pid={pid}: {exc}", log_path)

    deadline = time.time() + grace_sec
    while time.time() < deadline:
        alive = [p for p in uniq if Path(f"/proc/{p}").exists()]
        if not alive:
            log("All targeted processes exited after SIGTERM", log_path)
            return
        time.sleep(0.5)

    alive = [p for p in uniq if Path(f"/proc/{p}").exists()]
    if not alive:
        log("All targeted processes exited after SIGTERM", log_path)
        return

    log(f"SIGKILL remaining: {alive}", log_path)
    for pid in alive:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except PermissionError as exc:
            log(f"SIGKILL denied for pid={pid}: {exc}", log_path)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--threshold-pct",
        type=float,
        default=float(os.environ.get("MEM_WATCHDOG_THRESHOLD_PCT", "40")),
        help="Kill when matched RSS sum exceeds this %% of the memory basis (default: 40). "
        "Ignored if --threshold-gb is set.",
    )
    p.add_argument(
        "--threshold-gb",
        type=float,
        default=float(os.environ.get("MEM_WATCHDOG_THRESHOLD_GB", "0") or 0),
        help="Absolute RSS sum limit in GiB (overrides --threshold-pct when > 0).",
    )
    p.add_argument(
        "--basis",
        choices=("total", "available"),
        default="total",
        help="Denominator for the threshold: MemTotal (default) or MemAvailable.",
    )
    p.add_argument(
        "--interval",
        type=float,
        default=float(os.environ.get("MEM_WATCHDOG_INTERVAL_SEC", "10")),
        help="Polling interval in seconds (default: 10).",
    )
    p.add_argument(
        "--match",
        default=os.environ.get("MEM_WATCHDOG_MATCH", DEFAULT_MATCH),
        help="Regex matched against process cmdline (default: GENIE/multisim/cosmics/detvar drivers).",
    )
    p.add_argument(
        "--user",
        default=os.environ.get("USER") or os.environ.get("LOGNAME") or "",
        help="Only watch this user's processes (default: $USER).",
    )
    p.add_argument(
        "--grace-sec",
        type=float,
        default=15.0,
        help="Seconds to wait after SIGTERM before SIGKILL (default: 15).",
    )
    p.add_argument(
        "--log",
        default=None,
        help="Append status / kill events to this file.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Never kill; only log what would happen.",
    )
    p.add_argument(
        "--once",
        action="store_true",
        help="Check once and exit (exit code 2 if over threshold).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Only log when over threshold or when the matched set changes size.",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.threshold_gb and args.threshold_gb > 0:
        if args.threshold_gb <= 0:
            print("--threshold-gb must be > 0", file=sys.stderr)
            return 2
    elif args.threshold_pct <= 0 or args.threshold_pct > 100:
        print("--threshold-pct must be in (0, 100]", file=sys.stderr)
        return 2
    if args.interval <= 0:
        print("--interval must be > 0", file=sys.stderr)
        return 2

    try:
        import pwd

        uid = pwd.getpwnam(args.user).pw_uid if args.user else os.getuid()
    except KeyError:
        print(f"unknown user {args.user!r}", file=sys.stderr)
        return 2

    try:
        pattern = re.compile(args.match)
    except re.error as exc:
        print(f"invalid --match regex: {exc}", file=sys.stderr)
        return 2

    use_gb = bool(args.threshold_gb and args.threshold_gb > 0)
    if use_gb:
        log(
            f"mem_watchdog start threshold={args.threshold_gb:g} GiB (absolute) "
            f"interval={args.interval:g}s user={args.user or uid} dry_run={args.dry_run} "
            f"match={args.match!r}",
            args.log,
        )
    else:
        log(
            f"mem_watchdog start threshold={args.threshold_pct:g}% basis={args.basis} "
            f"interval={args.interval:g}s user={args.user or uid} dry_run={args.dry_run} "
            f"match={args.match!r}",
            args.log,
        )

    last_n = -1
    while True:
        mem = read_meminfo()
        den = basis_bytes(mem, args.basis)
        if use_gb:
            limit = int(args.threshold_gb * (1024**3))
            limit_label = f"{args.threshold_gb:g} GiB"
        else:
            limit = int(den * (args.threshold_pct / 100.0))
            limit_label = f"{args.threshold_pct:g}% of {args.basis}"
        hits = collect_matched(uid, pattern)
        rss_sum = sum(r for _p, r, _c in hits)
        frac = (100.0 * rss_sum / den) if den else 0.0

        status = (
            f"matched={len(hits)} rss={fmt_bytes(rss_sum)} "
            f"({frac:.1f}% of {args.basis} {fmt_bytes(den)}; "
            f"limit={limit_label} → {fmt_bytes(limit)})"
        )
        if not args.quiet or len(hits) != last_n or rss_sum > limit:
            log(status, args.log)
        last_n = len(hits)

        if hits and rss_sum > limit:
            log(
                f"THRESHOLD EXCEEDED — killing matched processes "
                f"(rss {fmt_bytes(rss_sum)} > {fmt_bytes(limit)})",
                args.log,
            )
            # Largest first so parents (often smaller) still get signal too
            pids = [p for p, _r, _c in sorted(hits, key=lambda t: -t[1])]
            terminate_pids(
                pids,
                dry_run=args.dry_run,
                grace_sec=args.grace_sec,
                log_path=args.log,
            )
            if args.once:
                return 2
            # After a kill, keep watching in case the launcher respawns workers
            time.sleep(args.interval)
            continue

        if args.once:
            return 0
        time.sleep(args.interval)


if __name__ == "__main__":
    sys.exit(main())
