#!/usr/bin/env python
"""Map phase: one offbeam or intime ``.df`` → pickle with summed histograms + gate counts.

HDF splits ``evt_{i}`` / ``hdr_{i}`` are processed sequentially. Intime histograms are
**unweighted** (gates scaling is applied in ``syst_cosmics_aggregate.py`` after summing
gates across all chunk pickles), matching the global ``data_gates / mc_gates`` recipe in
``files_config.get_ana_dfs(option='cosmics_systs')``.

Usage::

    python syst_cosmics_chunk.py --sample offbeam --df_file PATH.df --out_dir CHUNKS
    python syst_cosmics_chunk.py --sample intime --df_file PATH.df --out_dir CHUNKS

Output: ``cosmics__<sample>__<stem>.pkl``. Input paths are usually taken from
``dataset_locations.iter_cosmics_chunk_df_paths``.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from os import path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x


os.environ.setdefault("MPLBACKEND", "Agg")
_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from pyanalib.split_df_helpers import get_n_split

from analysis_village.numucc_1p0pi.syst_cosmics_common import build_variable_configs


def _gate_sum_from_hdr(hdr: pd.DataFrame, sample: str) -> float:
    if hdr is None or len(hdr) == 0:
        return 0.0
    m = hdr["first_in_subrun"] == 1
    sub = hdr.loc[m]
    if len(sub) == 0:
        return 0.0
    if sample == "offbeam":
        return float(sub["noffbeambnb"].sum())
    return float(sub["ngenevt"].sum())


def _histogram_evt(
    evtdf: pd.DataFrame,
    var_config: Any,
    *,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    col = var_config.var_evt_reco_col
    bins = var_config.bins
    if var_config.var_save_name == "integrated":
        if weights is None:
            return np.array([float(len(evtdf))], dtype=np.float64)
        return np.array([float(np.sum(weights))], dtype=np.float64)
    w = weights if weights is not None else None
    h, _ = np.histogram(evtdf[col], bins=bins, weights=w)
    return h.astype(np.float64, copy=False)


def accumulate_file(
    df_file: str,
    sample: str,
    var_configs: Sequence[Any],
    max_splits: int = 0,
) -> Dict[str, Any]:
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise ValueError("no HDF splits in %s" % df_file)

    gates_total = 0.0
    hists: Dict[str, np.ndarray] = {}
    for vc in var_configs:
        if vc.var_save_name == "integrated":
            hists[vc.var_save_name] = np.zeros(1, dtype=np.float64)
        else:
            hists[vc.var_save_name] = np.zeros(len(vc.bins) - 1, dtype=np.float64)

    for i in tqdm(range(n_use), desc=f"{sample}:{path.basename(df_file)}", leave=False):
        hdr = pd.read_hdf(df_file, key=f"hdr_{i}")
        gates_total += _gate_sum_from_hdr(hdr, sample)
        del hdr

        evtdf = pd.read_hdf(df_file, key=f"evt_{i}")
        # intime: store raw counts here; aggregate applies data/mc gate scale once.
        w = None if sample == "offbeam" else None
        for vc in var_configs:
            hists[vc.var_save_name] += _histogram_evt(evtdf, vc, weights=w)
        del evtdf
        gc.collect()

    return {
        "kind": "cosmics_syst_chunk",
        "sample": sample,
        "df_file": df_file,
        "splits_processed": n_use,
        "gates": gates_total,
        "hists": hists,
        "var_save_names": [vc.var_save_name for vc in var_configs],
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True, choices=("offbeam", "intime"))
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--max-splits", type=int, default=0, help="Cap HDF splits (0 = all).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    var_configs = build_variable_configs(None)
    payload = accumulate_file(args.df_file, args.sample, var_configs, max_splits=args.max_splits)

    stem = path.splitext(path.basename(args.df_file))[0]
    out_path = path.join(args.out_dir, "cosmics__%s__%s.pkl" % (args.sample, stem))
    with open(out_path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(
        "[cosmics-chunk] wrote %s  gates=%.6e  splits=%d"
        % (out_path, payload["gates"], payload["splits_processed"]),
        flush=True,
    )


if __name__ == "__main__":
    main()
