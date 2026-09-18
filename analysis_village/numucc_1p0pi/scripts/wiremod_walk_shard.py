#!/usr/bin/env python3
"""Exclusive-shard WireMod hist walk (Product A/B accumulation).

Each worker owns ``files[i]`` where ``i % n_shards == shard_id`` (or an explicit
file list). Checkpoints after each batch; merge shards with
``wiremod_merge_walk_shards.py``.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.scripts.reprocess_wiremod import _merge_univ_products, _rss_gb
from analysis_village.numucc_1p0pi.syst_detvar_common import accumulate_wiremod_matched_products, log


def _load_drop_map(path: Optional[Path]) -> Optional[dict]:
    if path is None or not path.is_file():
        return None
    with open(path, "rb") as fh:
        payload = pickle.load(fh)
    if isinstance(payload, dict) and "drop_map" in payload:
        log(
            f"  drop_map {path.name}: files={payload.get('n_files_with_drops')} "
            f"drop_rows={payload.get('n_drop')}"
        )
        return payload["drop_map"]
    return payload


def _walk(
    files: Sequence[str],
    *,
    checkpoint: Path,
    batch_size: int,
    rss_limit_gb: float,
    final_var_defs,
    drop_map: Optional[dict] = None,
    mu_p_candidate_kwargs: Optional[dict] = None,
) -> dict:
    files = list(files)
    n = len(files)
    start_batch = 0
    acc: dict = {}
    if checkpoint.is_file():
        with open(checkpoint, "rb") as fh:
            state = pickle.load(fh)
        acc = state.get("acc") or {}
        start_batch = int(state.get("next_batch", 0))
        log(f"  resume {checkpoint.name}: next_batch={start_batch} pot={acc.get('pot', 0):.3e}")

    n_batches = max((n + batch_size - 1) // batch_size, 0)
    for bi in range(start_batch, n_batches):
        lo = bi * batch_size
        hi = min(n, lo + batch_size)
        batch = files[lo:hi]
        rss0 = _rss_gb()
        log(f"  batch {bi + 1}/{n_batches} files[{lo}:{hi}] rss={rss0:.2f} GiB")
        if rss0 > rss_limit_gb:
            raise RuntimeError(f"RSS {rss0:.2f} GiB exceeds limit {rss_limit_gb:.2f} GiB")
        t0 = time.time()
        chunk = accumulate_wiremod_matched_products(
            batch,
            final_var_defs=final_var_defs,
            include_cut_stage=True,
            drop_map=drop_map,
            mu_p_candidate_kwargs=mu_p_candidate_kwargs,
        )
        acc = _merge_univ_products(acc, chunk)
        del chunk
        gc.collect()
        rss1 = _rss_gb()
        log(f"    pot_acc={acc['pot']:.3e} rss={rss1:.2f} GiB dt={time.time() - t0:.1f}s")
        if rss1 > rss_limit_gb:
            raise RuntimeError(f"RSS {rss1:.2f} GiB exceeds limit after batch")
        with open(checkpoint, "wb") as fh:
            pickle.dump(
                {
                    "acc": acc,
                    "next_batch": bi + 1,
                    "n_files": n,
                    "files": files,
                    "mu_p_candidate_kwargs": dict(mu_p_candidate_kwargs or {}),
                },
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
    return acc


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--files-pkl", required=True, help="Pickle of list[str] input matched dfs")
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--shard-id", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--batch-size", type=int, default=25)
    p.add_argument("--rss-limit-gb", type=float, default=20.0)
    p.add_argument("--label", default="shard")
    p.add_argument(
        "--drop-map-pkl",
        default=None,
        help="Optional artkey drop map from dedupe_matched_artkeys.py (XTXW)",
    )
    p.add_argument(
        "--mu-chi2mu-th",
        type=float,
        default=None,
        help="Override muon chi2 cut (default: MU_CHI2MU_TH=30 from selections.py)",
    )
    args = p.parse_args(argv)

    with open(args.files_pkl, "rb") as fh:
        all_files = list(pickle.load(fh))
    if args.n_shards > 1:
        files = [f for i, f in enumerate(all_files) if i % args.n_shards == args.shard_id]
    else:
        files = all_files
    pid_kw = {"mu_chi2mu_th": float(args.mu_chi2mu_th)} if args.mu_chi2mu_th is not None else None
    log(
        f"[{args.label} s{args.shard_id}/{args.n_shards}] "
        f"{len(files)}/{len(all_files)} files batch={args.batch_size} "
        f"mu_chi2mu_th={None if pid_kw is None else pid_kw['mu_chi2mu_th']}"
    )
    if not files:
        log("nothing to do")
        return 0

    drop_map = _load_drop_map(Path(args.drop_map_pkl) if args.drop_map_pkl else None)
    final_defs = dc.build_final_var_defs()
    ck = Path(args.checkpoint)
    ck.parent.mkdir(parents=True, exist_ok=True)
    acc = _walk(
        files,
        checkpoint=ck,
        batch_size=max(int(args.batch_size), 1),
        rss_limit_gb=float(args.rss_limit_gb),
        final_var_defs=final_defs,
        drop_map=drop_map,
        mu_p_candidate_kwargs=pid_kw,
    )
    log(f"[{args.label} s{args.shard_id}] done pot={acc.get('pot', 0):.3e} rss={_rss_gb():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
