#!/usr/bin/env python
"""Parallel dispatcher for :mod:`syst_histcounts_from_df_chunk`.

Processes a glob of sel_all weight ``.df`` files into histcounts HDFs with a
small worker pool (memory-heavy inputs — start with ``--workers 2``).
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
import glob
import multiprocessing as mp
import sys
import time
import traceback
from datetime import datetime
from os import path
from typing import List, Optional, Sequence

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.scripts import syst_histcounts_from_df_chunk as chunk_mod


def parse_cli(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input-glob",
        required=True,
        help="Glob of sel_all weight .df files (quote the glob)",
    )
    p.add_argument("--out-dir", required=True, help="Output directory for histcounts .df")
    p.add_argument("--family", required=True, choices=("Flux", "G4", "flux", "g4"))
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--n-universe", type=int, default=1000)
    p.add_argument("--include-slim", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--sample", default="mc", choices=("mc", "dirt"))
    p.add_argument("--flux-knob-groups", default="all")
    p.add_argument("--max-files", type=int, default=0, help="If >0, only first N files")
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument(
        "--hist-backend",
        default="vector",
        choices=("vector", "loop"),
        help="Universe histogram backend (default vector)",
    )
    p.add_argument("--failed-log", default="", help="Append failures here")
    p.add_argument(
        "--skip-existing",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip outputs that already exist (default: on)",
    )
    return p.parse_args(argv)


def _worker(job: dict) -> dict:
    started = time.time()
    df_file = job["df_file"]
    try:
        ns = argparse.Namespace(
            df_file=df_file,
            out_dir=job["out_dir"],
            family=job["family"],
            n_universe=job["n_universe"],
            include_slim=job["include_slim"],
            sample=job["sample"],
            flux_knob_groups=job["flux_knob_groups"],
            max_splits=job["max_splits"],
            out_name="",
            hist_backend=job.get("hist_backend", "vector"),
        )
        out_path, status = chunk_mod.run_with_args(ns, skip_existing=job["skip_existing"])
        return {
            "ok": True,
            "status": status,
            "df_file": df_file,
            "out_path": out_path,
            "elapsed": time.time() - started,
        }
    except SystemExit as e:
        return {
            "ok": False,
            "status": "failed",
            "df_file": df_file,
            "err": "SystemExit: %s" % e,
            "elapsed": time.time() - started,
        }
    except BaseException as e:
        return {
            "ok": False,
            "status": "failed",
            "df_file": df_file,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=8)),
            "elapsed": time.time() - started,
        }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_cli(argv)
    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit("[histcounts-from-df-parallel] no files match %r" % args.input_glob)
    if args.max_files > 0:
        files = files[: int(args.max_files)]

    os.makedirs(args.out_dir, exist_ok=True)
    jobs = [
        {
            "df_file": f,
            "out_dir": args.out_dir,
            "family": args.family,
            "n_universe": int(args.n_universe),
            "include_slim": bool(args.include_slim),
            "sample": args.sample,
            "flux_knob_groups": args.flux_knob_groups,
            "max_splits": int(args.max_splits),
            "hist_backend": str(args.hist_backend),
            "skip_existing": bool(args.skip_existing),
        }
        for f in files
    ]

    n_workers = max(1, int(args.workers))
    print(
        "[histcounts-from-df-parallel] %s  files=%d  workers=%d  family=%s  out=%s"
        % (datetime.now().isoformat(timespec="seconds"), len(jobs), n_workers, args.family, args.out_dir),
        flush=True,
    )

    ok = skip = fail = 0
    failed_lines: List[str] = []
    # fork: import chunk_mod in parent so workers inherit
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers, maxtasksperchild=8) as pool:
        for i, res in enumerate(pool.imap_unordered(_worker, jobs), start=1):
            if res["ok"]:
                if res["status"] == "skipped":
                    skip += 1
                else:
                    ok += 1
            else:
                fail += 1
                failed_lines.append("%s\t%s" % (res["df_file"], res.get("err", "")))
                print(
                    "[histcounts-from-df-parallel] FAIL %s: %s"
                    % (path.basename(res["df_file"]), res.get("err", "")[:200]),
                    flush=True,
                )
            if i % 10 == 0 or i == len(jobs):
                print(
                    "[histcounts-from-df-parallel] progress %d/%d  ok=%d skip=%d fail=%d"
                    % (i, len(jobs), ok, skip, fail),
                    flush=True,
                )

    if args.failed_log and failed_lines:
        with open(args.failed_log, "a") as fh:
            for line in failed_lines:
                fh.write(line.replace("\n", "\\n") + "\n")

    print(
        "[histcounts-from-df-parallel] done ok=%d skip=%d fail=%d" % (ok, skip, fail),
        flush=True,
    )
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
