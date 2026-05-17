#!/usr/bin/env python
"""Parallel dispatcher for :mod:`syst_cc_joint_genie_chunk` (joint-bin GENIE rate map).

Same (GENIE group, ``.df``) queue as :mod:`syst_genie_parallel` chunk-map; workers call
:func:`syst_cc_joint_genie_chunk.run_with_args` so outputs are ``nu__joint_cc_genie__*.pkl``.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
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

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    GENIE_GROUP_ORDER,
    _genie_glob_map,
    iter_genie_chunk_map_tasks,
    sorted_glob,
)
from analysis_village.numucc_1p0pi.scripts import syst_cc_joint_genie_chunk as joint_genie_chunk_mod  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import joint_cc_genie_chunk_basename  # noqa: E402


def _ordered_group_tags(mc_df_stage: str, allowed: Optional[Sequence[str]] = None) -> List[str]:
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


def _worker(job: dict) -> dict:
    grp = job["group"]
    df_file = job["df_file"]
    started = time.time()
    stem = path.splitext(path.basename(df_file))[0]
    out_pkl = path.join(job["out_dir"], joint_cc_genie_chunk_basename(grp, stem))
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
        ns = argparse.Namespace(
            df_file=df_file,
            out_dir=job["out_dir"],
            genie_group=grp,
            input_stage=job["input_stage"],
            max_splits=int(job["max_splits"]),
            pairs=job.get("pairs"),
        )
        out_path, status = joint_genie_chunk_mod.run_with_args(ns, skip_existing=False)
        return {
            "ok": True,
            "status": status,
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


def parse_cli(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--chunks-dir", required=True, help="Output directory for nu__joint_cc_genie__*.pkl")
    p.add_argument("--failed-log", required=True)
    p.add_argument("--genie-groups", default="", help="Comma-separated subset (empty = all in active map).")
    p.add_argument("--max-files", type=int, default=0, help="Per-group cap on map jobs (0 = all).")
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--pairs", default=None, help="Optional CSV of preset pair slugs for chunk --pairs.")
    return p.parse_args(argv)


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _failed_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main(argv: Optional[List[str]] = None) -> int:
    cli = parse_cli(argv)
    if cli.mc_df_stage != "final":
        print(
            "[cc-joint-genie-parallel] ERROR: only --mc-df-stage final is supported.",
            file=sys.stderr,
        )
        return 2
    allowed = [s.strip() for s in cli.genie_groups.split(",") if s.strip()] or None
    jobs, stats = _build_jobs(
        mc_df_stage=cli.mc_df_stage,
        allowed=allowed,
        max_files_per_group=int(cli.max_files or 0),
    )
    os.makedirs(cli.chunks_dir, exist_ok=True)
    os.makedirs(path.dirname(cli.failed_log) or ".", exist_ok=True)

    for g, n_glob, n_q in stats:
        print("[cc-joint-genie-parallel] group=%s glob_files=%d queued_map_jobs=%d" % (g, n_glob, n_q))

    total = len(jobs)
    print("[cc-joint-genie-parallel] mc_df_stage=%s queue=%d workers=%d" % (cli.mc_df_stage, total, cli.workers))
    if total == 0:
        return 0

    workers = max(1, min(int(cli.workers), total))
    payloads = [
        {
            "group": g,
            "df_file": p,
            "out_dir": cli.chunks_dir,
            "input_stage": cli.mc_df_stage,
            "max_splits": int(cli.max_splits),
            "pairs": cli.pairs,
        }
        for (g, p) in jobs
    ]

    done = 0
    failures = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers, maxtasksperchild=16) as pool, open(cli.failed_log, "a") as flog:
        for res in pool.imap_unordered(_worker, payloads, chunksize=1):
            done += 1
            grp = res["group"]
            df = res["df_file"]
            if res["ok"]:
                tag = "(skip existing)" if res["status"] == "skipped" else "END"
                print(
                    "[cc-joint-genie-parallel] progress map overall %d/%d  group=%s  %s %s "
                    "elapsed=%.1fs df_file=%s -> %s"
                    % (done, total, grp, tag, _now_iso(), res["elapsed"], df, res.get("out_path"))
                )
            else:
                failures += 1
                print(
                    "[cc-joint-genie-parallel] progress map overall %d/%d  group=%s  FAILED %s "
                    "elapsed=%.1fs df_file=%s\n%s"
                    % (done, total, grp, _now_iso(), res["elapsed"], df, res["err"]),
                    file=sys.stderr,
                )
                flog.write("%s\t%s\t%s\t%s\n" % (_failed_ts(), cli.mc_df_stage, grp, df))
                flog.flush()
    print(
        "[cc-joint-genie-parallel] map DONE total=%d failures=%d failed_log=%s"
        % (total, failures, cli.failed_log)
    )
    if failures and failures == total:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
