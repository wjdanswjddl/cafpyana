#!/usr/bin/env python
"""Notebook / CLI launcher for Flux & G4 histcounts-from-DF processing.

Same campaign layout as ``run_syst_histcounts_from_df.sh``, but callable from a
Jupyter kernel so workers run on the connected machine (e.g. EAF).

Foreground (blocks, streams logs to the notebook)::

    from analysis_village.numucc_1p0pi.scripts.run_syst_histcounts_from_df import (
        run_histcounts_from_df,
    )
    h = run_histcounts_from_df(family="G4", workers=10, n_universe=1000)

Background (detached like ``nohup … &``; kernel stays free)::

    h = run_histcounts_from_df(
        family="G4", workers=10, n_universe=1000, background=True,
    )
    print(h.out_root, h.log_path, h.pid)

CLI (equivalent to the bash env-var driver)::

    python analysis_village/numucc_1p0pi/scripts/run_syst_histcounts_from_df.py \\
      --family G4 --workers 10 --n-universe 1000
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from os import path
from typing import List, Optional, Sequence

REPO_ROOT = path.dirname(
    path.dirname(path.dirname(path.dirname(path.abspath(__file__))))
)
PARALLEL_SCRIPT = path.join(
    REPO_ROOT,
    "analysis_village",
    "numucc_1p0pi",
    "scripts",
    "syst_histcounts_from_df_parallel.py",
)

DEFAULT_FLUX_GLOB = (
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_09_04_172234__sel_all-wgts_flux-corrected_updated/*.df"
)
DEFAULT_G4_GLOB = (
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_09_04_175940__sel_all-wgts_g4-corrected_updated/*.df"
)


@dataclass
class HistcountsFromDfHandle:
    """Result of a campaign launch."""

    out_root: str
    family: str
    workers: int
    n_universe: int
    background: bool
    log_path: str = ""
    pid: Optional[int] = None
    returncode: Optional[int] = None

    def poll(self) -> Optional[int]:
        """Return exit code if the background process has finished, else None."""
        if self.pid is None:
            return self.returncode
        try:
            os.kill(self.pid, 0)
        except ProcessLookupError:
            # Process gone — best-effort read of a sentinel if we wrote one later.
            return self.returncode if self.returncode is not None else 0
        except PermissionError:
            return None
        return None

    def is_running(self) -> bool:
        if not self.background or self.pid is None:
            return False
        try:
            os.kill(self.pid, 0)
            return True
        except OSError:
            return False


def default_out_root(
    stamp: Optional[str] = None,
    user: Optional[str] = None,
) -> str:
    stamp = stamp or datetime.now().strftime("%Y_%m_%d_%H%M%S")
    user = user or os.environ.get("USER", "munjung")
    return (
        f"/pnfs/sbnd/scratch/users/{user}/cafpyana_tmp/"
        f"syst_histcounts_from_df_{stamp}"
    )


def _normalize_families(family: str) -> List[str]:
    key = family.strip().lower()
    if key in ("flux",):
        return ["Flux"]
    if key in ("g4",):
        return ["G4"]
    if key in ("both", "all"):
        return ["Flux", "G4"]
    raise ValueError("family must be Flux, G4, or both (got %r)" % family)


def _thread_env(base: Optional[dict] = None) -> dict:
    env = dict(base or os.environ)
    env["PYTHONPATH"] = REPO_ROOT + (
        (os.pathsep + env["PYTHONPATH"]) if env.get("PYTHONPATH") else ""
    )
    for k in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        env.setdefault(k, "1")
    env.setdefault("MPLBACKEND", "Agg")
    return env


def _family_argv(
    *,
    fam: str,
    glob_pat: str,
    out_root: str,
    workers: int,
    n_universe: int,
    sample: str,
    max_files: int,
    include_slim: bool,
    flux_knob_groups: str,
    max_splits: int,
    hist_backend: str,
    skip_existing: bool,
) -> List[str]:
    out_dir = path.join(out_root, "dfs", "hist_mc_%s" % fam.lower())
    failed_log = path.join(out_root, "logs", "failed_%s.txt" % fam)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(path.join(out_root, "logs"), exist_ok=True)
    argv = [
        sys.executable,
        PARALLEL_SCRIPT,
        "--input-glob",
        glob_pat,
        "--out-dir",
        out_dir,
        "--family",
        fam,
        "--workers",
        str(int(workers)),
        "--n-universe",
        str(int(n_universe)),
        "--sample",
        sample,
        "--max-files",
        str(int(max_files)),
        "--flux-knob-groups",
        flux_knob_groups,
        "--max-splits",
        str(int(max_splits)),
        "--hist-backend",
        hist_backend,
        "--failed-log",
        failed_log,
        "--include-slim" if include_slim else "--no-include-slim",
        "--skip-existing" if skip_existing else "--no-skip-existing",
    ]
    return argv


def run_histcounts_from_df(
    family: str = "both",
    *,
    workers: int = 2,
    n_universe: int = 1000,
    out_root: Optional[str] = None,
    stamp: Optional[str] = None,
    flux_glob: str = DEFAULT_FLUX_GLOB,
    g4_glob: str = DEFAULT_G4_GLOB,
    sample: str = "mc",
    max_files: int = 0,
    include_slim: bool = True,
    flux_knob_groups: str = "all",
    max_splits: int = 0,
    hist_backend: str = "vector",
    skip_existing: bool = True,
    background: bool = False,
    log_path: Optional[str] = None,
    dry_run: bool = False,
) -> HistcountsFromDfHandle:
    """Run Flux and/or G4 histcounts-from-DF on the kernel host.

    Parameters
    ----------
    family
        ``\"Flux\"``, ``\"G4\"``, or ``\"both\"``.
    background
        If True, detach a process group (``nohup``-style) and return immediately.
        Progress goes to ``log_path`` (default under ``out_root/logs/`` or ``/tmp``).
    """
    families = _normalize_families(family)
    out_root = out_root or default_out_root(stamp=stamp)
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(path.join(out_root, "logs"), exist_ok=True)

    # Same breadcrumb as the bash driver.
    try:
        with open("/tmp/syst_histcounts_from_df_campaign.txt", "w") as fh:
            fh.write(out_root + "\n")
    except OSError:
        pass

    globs = {"Flux": flux_glob, "G4": g4_glob}
    cmd_blocks: List[List[str]] = []
    for fam in families:
        cmd_blocks.append(
            _family_argv(
                fam=fam,
                glob_pat=globs[fam],
                out_root=out_root,
                workers=workers,
                n_universe=n_universe,
                sample=sample,
                max_files=max_files,
                include_slim=include_slim,
                flux_knob_groups=flux_knob_groups,
                max_splits=max_splits,
                hist_backend=hist_backend,
                skip_existing=skip_existing,
            )
        )

    print("OUT_ROOT=%s" % out_root, flush=True)
    print(
        "family=%s workers=%d n_universe=%d background=%s"
        % (",".join(families), workers, n_universe, background),
        flush=True,
    )
    for argv in cmd_blocks:
        print("CMD:", " ".join(argv), flush=True)

    if dry_run:
        return HistcountsFromDfHandle(
            out_root=out_root,
            family=",".join(families),
            workers=workers,
            n_universe=n_universe,
            background=background,
        )

    env = _thread_env()

    if background:
        if not log_path:
            tag = "_".join(f.lower() for f in families)
            log_path = path.join(
                out_root, "logs", "run_%s_%s.log" % (tag, datetime.now().strftime("%H%M%S"))
            )
            # Prefer /tmp when out_root is slow pnfs for live tailing — still
            # keep a copy path under out_root for the campaign record.
            tmp_log = "/tmp/syst_histcounts_from_df_run_%s.log" % tag
            log_path = tmp_log
        # Wrapper script: run all family commands sequentially in one detached job.
        wrapper = path.join(out_root, "logs", "_run_wrapper.sh")
        with open(wrapper, "w") as fh:
            fh.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            fh.write("cd %s\n" % shlex.quote(REPO_ROOT))
            for argv in cmd_blocks:
                fh.write("echo '===== %s ====='\n" % argv[argv.index("--family") + 1])
                fh.write(" ".join(shlex.quote(a) for a in argv) + "\n")
            fh.write(
                "echo '===== DONE campaign=%s ====='\n" % shlex.quote(out_root)
            )
        os.chmod(wrapper, 0o755)
        log_fh = open(log_path, "w")
        proc = subprocess.Popen(
            ["bash", wrapper],
            cwd=REPO_ROOT,
            env=env,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )
        log_fh.close()
        print(
            "started background pid=%d  log=%s  out_root=%s"
            % (proc.pid, log_path, out_root),
            flush=True,
        )
        return HistcountsFromDfHandle(
            out_root=out_root,
            family=",".join(families),
            workers=workers,
            n_universe=n_universe,
            background=True,
            log_path=log_path,
            pid=proc.pid,
        )

    # Foreground: sequential family runs; inherit notebook stdout/stderr.
    rc = 0
    for argv in cmd_blocks:
        fam = argv[argv.index("--family") + 1]
        print("===== FAMILY=%s =====" % fam, flush=True)
        proc = subprocess.run(argv, cwd=REPO_ROOT, env=env)
        if proc.returncode != 0:
            rc = proc.returncode
            print(
                "[run_histcounts_from_df] family %s failed rc=%d" % (fam, proc.returncode),
                flush=True,
            )
            break
    print("===== DONE  campaign=%s =====" % out_root, flush=True)
    return HistcountsFromDfHandle(
        out_root=out_root,
        family=",".join(families),
        workers=workers,
        n_universe=n_universe,
        background=False,
        returncode=rc,
    )


def tail_log(handle: HistcountsFromDfHandle, n_lines: int = 40) -> str:
    """Return the last ``n_lines`` of a background run log (empty if missing)."""
    if not handle.log_path or not path.isfile(handle.log_path):
        return ""
    with open(handle.log_path) as fh:
        lines = fh.readlines()
    return "".join(lines[-n_lines:])


def parse_cli(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--family", default=os.environ.get("FAMILY", "both"))
    p.add_argument("--workers", type=int, default=int(os.environ.get("WORKERS", "2")))
    p.add_argument(
        "--n-universe",
        type=int,
        default=int(os.environ.get("N_UNIVERSE", "1000")),
    )
    p.add_argument("--out-root", default=os.environ.get("OUT_ROOT", ""))
    p.add_argument("--stamp", default=os.environ.get("STAMP", ""))
    p.add_argument("--flux-glob", default=os.environ.get("FLUX_GLOB", DEFAULT_FLUX_GLOB))
    p.add_argument("--g4-glob", default=os.environ.get("G4_GLOB", DEFAULT_G4_GLOB))
    p.add_argument("--sample", default=os.environ.get("SAMPLE", "mc"))
    p.add_argument("--max-files", type=int, default=int(os.environ.get("MAX_FILES", "0")))
    p.add_argument(
        "--include-slim",
        action=argparse.BooleanOptionalAction,
        default=os.environ.get("INCLUDE_SLIM", "1") not in ("0", "false", "False"),
    )
    p.add_argument("--flux-knob-groups", default="all")
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument("--hist-backend", default="vector", choices=("vector", "loop"))
    p.add_argument(
        "--skip-existing",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    p.add_argument(
        "--background",
        action="store_true",
        default=False,
        help="Detach like nohup (writes /tmp/syst_histcounts_from_df_run_*.log)",
    )
    p.add_argument("--log-path", default="")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_cli(argv)
    handle = run_histcounts_from_df(
        family=args.family,
        workers=args.workers,
        n_universe=args.n_universe,
        out_root=args.out_root or None,
        stamp=args.stamp or None,
        flux_glob=args.flux_glob,
        g4_glob=args.g4_glob,
        sample=args.sample,
        max_files=args.max_files,
        include_slim=args.include_slim,
        flux_knob_groups=args.flux_knob_groups,
        max_splits=args.max_splits,
        hist_backend=args.hist_backend,
        skip_existing=args.skip_existing,
        background=args.background,
        log_path=args.log_path or None,
        dry_run=args.dry_run,
    )
    if handle.background:
        # Keep CLI alive briefly so the child can start; then exit 0.
        time.sleep(0.2)
        print("background pid=%s log=%s" % (handle.pid, handle.log_path), flush=True)
        return 0
    return int(handle.returncode or 0)


if __name__ == "__main__":
    raise SystemExit(main())
