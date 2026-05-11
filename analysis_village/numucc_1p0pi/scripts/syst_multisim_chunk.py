#!/usr/bin/env python
"""Map phase: one MC ``.df`` file → pickle with summed Flux/G4/MCstat universe histograms.

HDF splits ``evt_0``, ``evt_1``, … are processed sequentially (low RAM).  Aggregation
(``syst_multisim_aggregate.py``) sums ``univ_events`` / ``cv_events`` across files then
runs ``get_covariance_matrix`` — additive across disjoint chunks.

Usage::
    python syst_multisim_chunk.py --df_file PATH.df --out_dir CHUNKS \\
        [--var-set final|intermediate|both] [--syst-names Flux,G4,...] \\
        [--n-universe 100] [--max-splits 0]

When ``--syst-names`` is a proper subset of MCstat/Flux/G4, output is
``nu__<names>__<stem>.pkl`` so chunks from different input dirs do not collide.
All three (default) keeps the legacy name ``nu__<stem>.pkl``.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from os import path

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x


os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.split_df_helpers import get_n_split

from analysis_village.numucc_1p0pi.syst_multisim_common import (
    NEUTRINO_SYST_ORDER,
    build_var_configs,
    drop_bad_g4_weights,
    syst_key_for_name,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates


def _parse_syst_names(spec: str | None):
    if not spec or not str(spec).strip():
        return tuple(NEUTRINO_SYST_ORDER)
    raw = [x.strip() for x in str(spec).split(",") if x.strip()]
    unk = set(raw) - set(NEUTRINO_SYST_ORDER)
    if unk:
        raise SystemExit("[multisim-chunk] unknown --syst-names entries: %s" % sorted(unk))
    # Preserve canonical NEUTRINO_SYST_ORDER for repeats reproducibility
    ordered = tuple(sn for sn in NEUTRINO_SYST_ORDER if sn in raw)
    if len(ordered) != len(raw):
        raise SystemExit("[multisim-chunk] duplicate syst name in --syst-names")
    return ordered


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--var-set", choices=("final", "intermediate", "both"), default="final")
    p.add_argument(
        "--syst-names",
        default=None,
        help="Comma-separated subset of MCstat,Flux,G4 (default: all). Use one name when that "
        "systematic's weights live in a separate directory of .df files.",
    )
    p.add_argument("--n-universe", type=int, default=100)
    p.add_argument(
        "--max-splits",
        type=int,
        default=0,
        help="Cap HDF splits processed (0 = all).",
    )
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    syst_names = _parse_syst_names(args.syst_names)
    syst_active = frozenset(syst_names)

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if args.max_splits <= 0 else min(args.max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[multisim-chunk] no splits")

    var_configs = build_var_configs(args.var_set)

    acc_syst = {sn: {} for sn in NEUTRINO_SYST_ORDER}

    def flush_evt(mc_evt_df: pd.DataFrame) -> None:
        if "G4" in syst_active:
            mc_evt_df = drop_bad_g4_weights(mc_evt_df, n_univ=args.n_universe)
        if len(mc_evt_df) == 0:
            return
        for sname in syst_names:
            sk = syst_key_for_name(sname)
            for vc in var_configs:
                try:
                    univ, cv = get_univ_rates(
                        cov_type="rate",
                        evtdf=mc_evt_df,
                        var_config=vc,
                        syst_name=sk,
                        n_univ=args.n_universe,
                        bkgd_subtract=True,
                    )
                except Exception as ex:
                    print(f"[multisim-chunk] skip var={vc.var_save_name} syst={sname}: {ex}")
                    continue
                slug = vc.var_save_name
                if slug not in acc_syst[sname]:
                    acc_syst[sname][slug] = {
                        "univ_events": np.array(univ, dtype=float),
                        "cv_events": np.array(cv, dtype=float),
                    }
                else:
                    acc_syst[sname][slug]["univ_events"] += univ
                    acc_syst[sname][slug]["cv_events"] += cv

    for i in tqdm(range(n_use), desc="HDF splits"):
        mc_evt_df = pd.read_hdf(args.df_file, key=f"evt_{i}")
        flush_evt(mc_evt_df)
        del mc_evt_df
        gc.collect()

    stem = path.splitext(path.basename(args.df_file))[0]
    if syst_names == NEUTRINO_SYST_ORDER:
        out_leaf = "nu__{}.pkl".format(stem)
    else:
        tag = "_".join(syst_names)
        out_leaf = "nu__{}__{}.pkl".format(tag, stem)
    out_path = path.join(args.out_dir, out_leaf)
    blob = {
        "kind": "multisim_syst_nu_chunk",
        "df_file": args.df_file,
        "syst_names_computed": list(syst_names),
        "var_set": args.var_set,
        "n_univ_requested": args.n_universe,
        "splits_processed": n_use,
        "syst": acc_syst,
    }
    with open(out_path, "wb") as f:
        pickle.dump(blob, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("[multisim-chunk] wrote", out_path)


if __name__ == "__main__":
    main()
