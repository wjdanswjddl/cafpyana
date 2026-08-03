#!/usr/bin/env python
"""Parallel dispatcher for the GENIE multisim chunk-map phase.

Replaces the serial ``for line in jobs; do python get_systematics_genie.py chunk-map ...; done``
loop in :mod:`run_syst_genie_chunked.sh`. Same design as
:mod:`syst_multisim_parallel`:

* Imports :mod:`get_systematics_genie` **once** in the master; workers fork
  from it (cheap on Linux).
* Calls :func:`get_systematics_genie.run_chunk_map` directly inside each
  worker (no subprocess relaunch per file).
* Writes the **same per-file pickle** (``genie__<GROUP>__<stem>.pkl``) the
  serial path would write, so:
    - ``chunk-merge`` is unchanged,
    - skip-existing semantics are preserved,
    - a bad file cannot corrupt a sibling (process-level isolation +
      atomic ``.tmp`` rename in the chunk writer).
* Optionally fans the per-group merge phase out across processes as well
  (``--merge-workers``); merges are CPU-bound on bigger groups.

The shell driver passes a single Python invocation through this module instead
of looping per (group, df) pair.
"""
from __future__ import annotations

# Cap BLAS / OpenMP thread fan-out BEFORE numpy / pandas / etc. import.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import glob
import multiprocessing as mp
import sys
import time
import traceback
from collections import Counter
from datetime import datetime
from os import path
from typing import List, Optional, Sequence, Tuple

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.dataset_locations import (
    GENIE_GROUP_KNOBS,
    GENIE_GROUP_ORDER,
    _genie_glob_map,
    iter_genie_chunk_map_tasks,
    sorted_glob,
)
from analysis_village.numucc_1p0pi.scripts import get_systematics_genie as genie_mod


# ---------------------------------------------------------------------------
# Job construction
# ---------------------------------------------------------------------------
def _ordered_group_tags(
    mc_df_stage: str,
    allowed: Optional[Sequence[str]] = None,
) -> List[str]:
    """``GENIE_GROUP_ORDER`` first, then any extra keys (sorted)."""
    gmap = _genie_glob_map(mc_df_stage)
    seen: set = set()
    out: List[str] = []
    for t in GENIE_GROUP_ORDER:
        if t not in gmap:
            continue
        if allowed is not None and t not in allowed:
            continue
        out.append(t)
        seen.add(t)
    for t in sorted(k for k in gmap if k not in seen):
        if allowed is not None and t not in allowed:
            continue
        out.append(t)
    return out


def _build_jobs(
    *,
    mc_df_stage: str,
    allowed: Optional[Sequence[str]],
    max_files_per_group: int,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, int, int]]]:
    """Return (jobs, per_group_stats).

    ``jobs`` = list of ``(genie_group, df_path)`` pairs ordered by group.
    ``per_group_stats`` = list of ``(group, glob_matches, jobs_queued)`` for log.
    """
    allowed_set = set(allowed) if allowed else None
    raw: List[Tuple[str, str]] = []
    for grp, p in iter_genie_chunk_map_tasks(mc_df_stage=mc_df_stage):
        if allowed_set is not None and grp not in allowed_set:
            continue
        raw.append((grp, p))
    if max_files_per_group > 0:
        per_grp = Counter()
        capped = []
        for grp, p in raw:
            if per_grp[grp] >= max_files_per_group:
                continue
            per_grp[grp] += 1
            capped.append((grp, p))
        jobs = capped
    else:
        jobs = raw
    gmap = _genie_glob_map(mc_df_stage)
    queued = Counter(g for g, _ in jobs)
    stats: List[Tuple[str, int, int]] = []
    for t in _ordered_group_tags(mc_df_stage):
        if allowed_set is not None and t not in allowed_set:
            stats.append((t, len(sorted_glob(gmap[t])), 0))
        else:
            stats.append((t, len(sorted_glob(gmap[t])), queued.get(t, 0)))
    return jobs, stats


def _expected_out_pkl(out_dir: str, group: str, df_file: str) -> str:
    stem = path.splitext(path.basename(df_file))[0]
    return path.join(out_dir, "genie__%s__%s.pkl" % (group, stem))


# ---------------------------------------------------------------------------
# Worker (chunk-map)
# ---------------------------------------------------------------------------
def _worker_chunk_map(job: dict) -> dict:
    """Run one (group, df_file) chunk-map job. Mirrors the CLI flags."""
    grp = job["group"]
    df_file = job["df_file"]
    out_dir = job["out_dir"]
    started = time.time()
    out_pkl = _expected_out_pkl(out_dir, grp, df_file)
    if path.exists(out_pkl):
        return {
            "ok": True,
            "status": "skipped",
            "group": grp,
            "df_file": df_file,
            "out_path": out_pkl,
            "elapsed": 0.0,
        }
    try:
        knobs = list(GENIE_GROUP_KNOBS[grp])
        if not knobs:
            raise SystemExit("[genie-chunk-map] no knobs for group %s" % grp)
        syst_names = [("mc", k) for k in knobs]
        var_configs = genie_mod.genie_all_var_configs(job["input_stage"])
        out_path = genie_mod.run_chunk_map(
            df_file=df_file,
            out_dir=out_dir,
            genie_group=grp,
            var_configs=var_configs,
            syst_names=syst_names,
            max_splits=int(job["max_splits"]),
            input_stage=job["input_stage"],
        )
        return {
            "ok": True,
            "status": "ok",
            "group": grp,
            "df_file": df_file,
            "out_path": out_path,
            "elapsed": time.time() - started,
        }
    except SystemExit as e:
        return {
            "ok": False,
            "status": "failed",
            "group": grp,
            "df_file": df_file,
            "err": "SystemExit: %s" % e,
            "elapsed": time.time() - started,
        }
    except BaseException as e:
        return {
            "ok": False,
            "status": "failed",
            "group": grp,
            "df_file": df_file,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=8)),
            "elapsed": time.time() - started,
        }


# ---------------------------------------------------------------------------
# Worker (chunk-merge) — independent per group
# ---------------------------------------------------------------------------
def _worker_chunk_merge(job: dict) -> dict:
    grp = job["group"]
    started = time.time()
    try:
        import logging
        logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
        var_configs = genie_mod.genie_all_var_configs(job["input_stage"])
        out_sub = path.join(job["merge_root"], grp)
        os.makedirs(out_sub, exist_ok=True)
        genie_mod.run_chunk_merge(
            chunks_dir=job["chunks_dir"],
            out_dir=out_sub,
            genie_group=grp,
            var_configs=var_configs,
            xsec_unit=float(job["xsec_unit"]),
            bkgd_subtract=True,
            save_figs=False,
            npz_path=None,
        )
        return {
            "ok": True,
            "group": grp,
            "out_dir": out_sub,
            "elapsed": time.time() - started,
        }
    except BaseException as e:
        return {
            "ok": False,
            "group": grp,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=8)),
            "elapsed": time.time() - started,
        }


# ---------------------------------------------------------------------------
# CLI / driver
# ---------------------------------------------------------------------------
def parse_cli(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--chunks-dir", required=True)
    p.add_argument("--merge-root", required=True)
    p.add_argument("--failed-log", required=True)
    p.add_argument(
        "--genie-groups",
        default="",
        help="Comma-separated subset of GENIE_GROUP_GLOBS keys (empty = all in active map).",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Per-group cap on map jobs (0 = no cap).",
    )
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument("--xsec-unit", type=float, default=1.0)
    p.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Map-phase worker count.",
    )
    p.add_argument(
        "--merge-workers",
        type=int,
        default=0,
        help="Merge-phase worker count (0 → min(N_groups, --workers)).",
    )
    p.add_argument("--skip-merge", action="store_true",
                   help="Run only the map phase; useful for grid-style fan-out + later merge.")
    return p.parse_args(argv)


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _failed_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main(argv: Optional[List[str]] = None) -> int:
    cli = parse_cli(argv)
    allowed = [s.strip() for s in cli.genie_groups.split(",") if s.strip()] or None
    if allowed:
        gmap = _genie_glob_map(cli.mc_df_stage)
        bad = set(allowed) - set(gmap.keys())
        if bad:
            print(
                "[genie-parallel] ERROR: unknown GENIE group(s) %r; valid: %s"
                % (sorted(bad), tuple(gmap.keys())),
                file=sys.stderr,
            )
            return 2

    os.makedirs(cli.chunks_dir, exist_ok=True)
    os.makedirs(cli.merge_root, exist_ok=True)
    os.makedirs(path.dirname(cli.failed_log) or ".", exist_ok=True)

    jobs, stats = _build_jobs(
        mc_df_stage=cli.mc_df_stage,
        allowed=allowed,
        max_files_per_group=int(cli.max_files or 0),
    )
    for g, ng, nq in stats:
        print(
            "[genie-parallel]   %s: glob_matches=%d  chunk_map_jobs_queued=%d" % (g, ng, nq)
        )
    total = len(jobs)
    print(
        "[genie-parallel] chunk-map queue: %d (group, .df) job(s) "
        "for MC_DF_STAGE=%s; workers=%d"
        % (total, cli.mc_df_stage, cli.workers)
    )

    map_failures = 0
    if total > 0:
        workers = max(1, min(int(cli.workers), total))
        payloads = [
            {
                "group": grp,
                "df_file": p,
                "out_dir": cli.chunks_dir,
                "input_stage": cli.mc_df_stage,
                "max_splits": int(cli.max_splits),
            }
            for (grp, p) in jobs
        ]
        done = 0
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=workers, maxtasksperchild=32) as pool, open(cli.failed_log, "a") as flog:
            for res in pool.imap_unordered(_worker_chunk_map, payloads, chunksize=1):
                done += 1
                grp = res["group"]
                df = res["df_file"]
                if res["ok"]:
                    tag = "(skip existing)" if res["status"] == "skipped" else "END"
                    print(
                        "[genie-parallel] progress map overall %d/%d  group=%s  %s %s "
                        "elapsed=%.1fs df_file=%s -> %s"
                        % (done, total, grp, tag, _now_iso(), res["elapsed"], df, res.get("out_path"))
                    )
                else:
                    map_failures += 1
                    print(
                        "[genie-parallel] progress map overall %d/%d  group=%s  FAILED %s "
                        "elapsed=%.1fs df_file=%s\n%s"
                        % (done, total, grp, _now_iso(), res["elapsed"], df, res["err"]),
                        file=sys.stderr,
                    )
                    flog.write("%s\t%s\t%s\t%s\n" % (_failed_ts(), cli.mc_df_stage, grp, df))
                    flog.flush()
        print(
            "[genie-parallel] map DONE total=%d failures=%d failed_log=%s"
            % (total, map_failures, cli.failed_log)
        )
    else:
        print("[genie-parallel] map queue empty; nothing to do")

    if cli.skip_merge:
        print("[genie-parallel] --skip-merge → chunks only at %s" % cli.chunks_dir)
        return 0 if (total == 0 or map_failures < total) else 1

    # Per-group merge — find groups that actually have pickles under chunks_dir.
    merge_targets: List[str] = []
    for grp in _ordered_group_tags(cli.mc_df_stage, allowed=allowed):
        if glob.glob(path.join(cli.chunks_dir, "genie__%s__*.pkl" % grp)):
            merge_targets.append(grp)
        else:
            print("[genie-parallel] skip chunk-merge (no genie__%s__*.pkl under %s)" % (grp, cli.chunks_dir))
    n_merge = len(merge_targets)
    print("[genie-parallel] chunk-merge queue: %d knob group(s)" % n_merge)

    if n_merge == 0:
        return 0 if (total == 0 or map_failures < total) else 1

    merge_workers = int(cli.merge_workers or 0) or min(n_merge, max(1, int(cli.workers)))
    merge_workers = max(1, min(merge_workers, n_merge))
    merge_payloads = [
        {
            "group": g,
            "chunks_dir": cli.chunks_dir,
            "merge_root": cli.merge_root,
            "input_stage": cli.mc_df_stage,
            "xsec_unit": float(cli.xsec_unit),
        }
        for g in merge_targets
    ]
    merge_failures = 0
    done = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=merge_workers) as pool:
        for res in pool.imap_unordered(_worker_chunk_merge, merge_payloads, chunksize=1):
            done += 1
            grp = res["group"]
            if res["ok"]:
                print(
                    "[genie-parallel] progress merge %d/%d  group=%s  END %s elapsed=%.1fs -> %s"
                    % (done, n_merge, grp, _now_iso(), res["elapsed"], res.get("out_dir"))
                )
            else:
                merge_failures += 1
                print(
                    "[genie-parallel] progress merge %d/%d  group=%s  FAILED %s elapsed=%.1fs\n%s"
                    % (done, n_merge, grp, _now_iso(), res["elapsed"], res["err"]),
                    file=sys.stderr,
                )

    print(
        "[genie-parallel] DONE map=%s merge=%s map_failures=%d merge_failures=%d"
        % (cli.chunks_dir, cli.merge_root, map_failures, merge_failures)
    )
    # Non-zero exit only when merge fully failed (matches multisim parallel convention).
    if merge_failures == n_merge and n_merge:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
