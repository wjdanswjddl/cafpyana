#!/usr/bin/env python
"""Run event selection on one batched job (multiple ``.df`` files, ≤ size budget).

Processes each file sequentially through the notebook pipeline accumulators,
then writes one pickle compatible with ``event_selection_aggregate.py``.
"""
from __future__ import annotations

import argparse
import os
import sys
from os import path

os.environ.setdefault("MPLBACKEND", "Agg")

import warnings

import pandas as pd

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from analysis_village.numucc_1p0pi.event_selection_batch_core import run_batch_selection


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True, choices=["mc", "data", "intime", "offbeam", "dirt"])
    p.add_argument("--df_file", action="append", required=True, help="Input .df path (repeatable)")
    p.add_argument("--job_id", type=int, default=0)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--use-mc-genweight", action="store_true")
    p.add_argument("--mc-univ-syst", default="", help="MC only: Flux,G4,GENIE")
    p.add_argument("--trace", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    mc_univ = tuple(x.strip() for x in (args.mc_univ_syst or "").split(",") if x.strip())
    tag = f"batch_{args.job_id:04d}"
    out_path = path.join(args.out_dir, f"{args.sample}__{tag}.pkl")

    pipeline_trace = (lambda msg: print(msg, flush=True)) if args.trace else None

    print(
        f"[batch_map] sample={args.sample} job={tag} files={len(args.df_file)} → {out_path}",
        flush=True,
    )
    meta = run_batch_selection(
        args.sample,
        args.df_file,
        out_path,
        job_id=tag,
        use_mc_genweight=args.use_mc_genweight,
        mc_univ_syst_tags=mc_univ,
        pipeline_trace=pipeline_trace,
    )
    print(
        f"[batch_map] done  n_evt={meta['n_evt']}  pot={meta['chunk_pot']:.3e}  "
        f"per_file={[round(f['pot'], 3) for f in meta['per_file']]}",
        flush=True,
    )


if __name__ == "__main__":
    main()
