#!/usr/bin/env python3
"""Shard walk for matched CV (nominal sel_all Product A/B — no calo univs)."""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.scripts.reprocess_wiremod import _rss_gb
from analysis_village.numucc_1p0pi.syst_detvar_common import accumulate_matched_sel_all_products, log


def _merge_cv(acc: dict, chunk: dict) -> dict:
    if not acc:
        return {
            "hists_cut": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_cut"].items()},
            "hists_final": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_final"].items()},
            "pot": float(chunk["pot"]),
            "cut_var_names": list(chunk["cut_var_names"]),
            "final_var_names": list(chunk["final_var_names"]),
        }
    acc["pot"] = float(acc["pot"]) + float(chunk["pot"])
    for key in ("hists_cut", "hists_final"):
        for var, hist in chunk[key].items():
            acc[key][var] = np.asarray(acc[key].get(var, 0.0), dtype=float) + np.asarray(hist, dtype=float)
    return acc


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--files-pkl", required=True)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--shard-id", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--batch-size", type=int, default=25)
    p.add_argument("--rss-limit-gb", type=float, default=20.0)
    p.add_argument(
        "--mu-chi2mu-th",
        type=float,
        default=None,
        help="Override muon chi2 cut (default: MU_CHI2MU_TH=30 from selections.py)",
    )
    args = p.parse_args(argv)

    with open(args.files_pkl, "rb") as fh:
        all_files = list(pickle.load(fh))
    files = [f for i, f in enumerate(all_files) if i % args.n_shards == args.shard_id]
    pid_kw = {"mu_chi2mu_th": float(args.mu_chi2mu_th)} if args.mu_chi2mu_th is not None else None
    log(
        f"[CV s{args.shard_id}/{args.n_shards}] {len(files)}/{len(all_files)} files "
        f"mu_chi2mu_th={None if pid_kw is None else pid_kw['mu_chi2mu_th']}"
    )
    if not files:
        return 0

    final_defs = dc.build_final_var_defs()
    ck = Path(args.checkpoint)
    ck.parent.mkdir(parents=True, exist_ok=True)
    acc: dict = {}
    start_batch = 0
    if ck.is_file():
        with open(ck, "rb") as fh:
            state = pickle.load(fh)
        acc = state.get("acc") or {}
        start_batch = int(state.get("next_batch", 0))
        log(f"  resume next_batch={start_batch} pot={acc.get('pot', 0):.3e}")

    batch_size = max(int(args.batch_size), 1)
    n = len(files)
    n_batches = (n + batch_size - 1) // batch_size
    for bi in range(start_batch, n_batches):
        lo = bi * batch_size
        hi = min(n, lo + batch_size)
        batch = files[lo:hi]
        rss0 = _rss_gb()
        log(f"  batch {bi + 1}/{n_batches} files[{lo}:{hi}] rss={rss0:.2f} GiB")
        if rss0 > args.rss_limit_gb:
            raise RuntimeError(f"RSS {rss0:.2f} GiB exceeds limit")
        t0 = time.time()
        chunk = accumulate_matched_sel_all_products(
            batch,
            final_var_defs=final_defs,
            include_cut_stage=True,
            mu_p_candidate_kwargs=pid_kw,
        )
        acc = _merge_cv(acc, chunk)
        del chunk
        gc.collect()
        log(f"    pot_acc={acc['pot']:.3e} rss={_rss_gb():.2f} GiB dt={time.time() - t0:.1f}s")
        with open(ck, "wb") as fh:
            pickle.dump(
                {
                    "acc": acc,
                    "next_batch": bi + 1,
                    "n_files": n,
                    "mu_p_candidate_kwargs": dict(pid_kw or {}),
                },
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
    log(f"[CV s{args.shard_id}] done pot={acc.get('pot', 0):.3e} rss={_rss_gb():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
