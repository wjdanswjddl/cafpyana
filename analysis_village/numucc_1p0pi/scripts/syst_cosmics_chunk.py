#!/usr/bin/env python
"""Map phase: one offbeam/intime ``.df`` → pickle with summed histograms + gate counts.

Two input layouts are supported:

* ``--input-stage final``: ``.df`` already at final-selection level. Each HDF split
  has only ``evt_{i}`` / ``hdr_{i}`` columns for the final variables; we simply
  histogram those columns. Same behaviour as the original chunk script.

* ``--input-stage sel_all``: ``.df`` is a raw sel_all bundle (``evt_{i}`` / ``trk_{i}`` /
  ``hdr_{i}`` columns, no slice cuts yet). We re-run the full numuCC 1p0pi event
  selection (``event_selection_pipeline_def.build_pipeline``) on each split and
  histogram each stage's cut variable (e.g. ``nu_score``, ``n_trks``,
  ``track_score``, ``vtx_dist`` …) plus the final-selected variables.

In both layouts intime histograms are **unweighted** -- the aggregator applies the
global ``sum(offbeam gates) / sum(intime gates)`` scale once, matching
``files_config.get_ana_dfs(option='cosmics_systs')``.

Usage::

    python syst_cosmics_chunk.py --sample offbeam --df_file PATH.df --out_dir CHUNKS \\
        --input-stage final         # legacy: SELECTED_EVENTS_GLOBS-style inputs
    python syst_cosmics_chunk.py --sample intime  --df_file PATH.df --out_dir CHUNKS \\
        --input-stage sel_all       # new: EVENT_SELECTION_GLOBS-style inputs

Output: ``cosmics__<sample>__<stem>.pkl``. The pickle records ``input_stage`` so the
aggregator knows whether to expect ``hists`` (final) or ``stage_hists`` (sel_all).
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

import warnings
# turn off performance warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

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

from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.syst_cosmics_common import build_variable_configs
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    final_stage_var_configs,
    get_var_series,
    histogram_var,
    walk_pipeline,
)


# ---------------------------------------------------------------------------
# Gate accumulator (identical to legacy behaviour).
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Final-stage histogram (preserves the legacy "integrated == count" shortcut).
# ---------------------------------------------------------------------------
def _histogram_evt_final(evtdf: pd.DataFrame, var_config: Any) -> np.ndarray:
    if var_config.var_save_name == "integrated":
        return np.array([float(len(evtdf))], dtype=np.float64)
    try:
        x = multicol_get_series(evtdf, var_config.var_evt_reco_col)
    except KeyError:
        return np.zeros(len(var_config.bins) - 1, dtype=np.float64)
    return histogram_var(np.asarray(x, dtype=float), np.asarray(var_config.bins))


# ===========================================================================
# Final-stage input path: histogram already-selected variables from evt only.
# ===========================================================================
def accumulate_file_final(
    df_file: str,
    sample: str,
    var_configs: Sequence[Any],
    max_splits: int = 0,
    *,
    verbose: bool = False,
) -> Dict[str, Any]:
    def vlog(msg: str) -> None:
        if verbose:
            print("[cosmics-chunk]", msg, flush=True)

    vlog("[final] accumulate_file_final start sample=%s df_file=%s" % (sample, df_file))
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise ValueError("no HDF splits in %s" % df_file)
    vlog("[final] processing n_use=%d splits" % n_use)

    gates_total = 0.0
    hists: Dict[str, np.ndarray] = {}
    for vc in var_configs:
        n = 1 if vc.var_save_name == "integrated" else len(vc.bins) - 1
        hists[vc.var_save_name] = np.zeros(n, dtype=np.float64)

    for i in tqdm(range(n_use), desc=f"{sample}:{path.basename(df_file)}", leave=False):
        hdr = pd.read_hdf(df_file, key=f"hdr_{i}")
        gates_total += _gate_sum_from_hdr(hdr, sample)
        del hdr

        evtdf = pd.read_hdf(df_file, key=f"evt_{i}")
        for vc in var_configs:
            hists[vc.var_save_name] += _histogram_evt_final(evtdf, vc)
        del evtdf
        gc.collect()

    return {
        "kind": "cosmics_syst_chunk",
        "input_stage": "final",
        "sample": sample,
        "df_file": df_file,
        "splits_processed": n_use,
        "gates": gates_total,
        "hists": hists,
        "var_save_names": [vc.var_save_name for vc in var_configs],
    }


# ===========================================================================
# Sel_all input path: run the full selection pipeline, histogram cut + final vars.
# ===========================================================================
def _hdf_has_key(df_file: str, key: str) -> bool:
    try:
        with pd.HDFStore(df_file, mode="r") as store:
            keys = store.keys()
        return ("/" + key) in keys
    except Exception:
        return False


def _init_stage_hists(
    cut_specs: Sequence[Any],
    final_var_configs: Sequence[Any],
) -> Dict[str, Dict[str, np.ndarray]]:
    """``stage_hists[stage_key][var_save_name]`` -> zeroed 1D array."""
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for spec in cut_specs:
        nbin = len(spec.var_config.bins) - 1
        out.setdefault(spec.stage_key, {})[spec.var_config.var_save_name] = np.zeros(
            nbin, dtype=np.float64
        )
    for vc in final_var_configs:
        nbin = 1 if vc.var_save_name == "integrated" else len(vc.bins) - 1
        out.setdefault(FINAL_STAGE_KEY, {})[vc.var_save_name] = np.zeros(nbin, dtype=np.float64)
    return out


def accumulate_file_sel_all(
    df_file: str,
    sample: str,
    max_splits: int = 0,
    *,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Re-run the full pipeline on each split; histogram cut + final variables."""
    def vlog(msg: str) -> None:
        if verbose:
            print("[cosmics-chunk]", msg, flush=True)

    # Pipeline cuts only branch on evt/trk; intime/offbeam are non-MC so mcnu is
    # always absent. The pipeline tolerates missing mcnu.
    pipeline_sample = sample  # "intime" or "offbeam"

    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise ValueError("no HDF splits in %s" % df_file)
    vlog("[sel_all] processing n_use=%d splits (pipeline sample=%s)" % (n_use, pipeline_sample))

    cut_specs = list(CUT_STAGE_VAR_SPECS)
    final_var_configs = list(final_stage_var_configs())
    cut_by_stage: Dict[str, List[Any]] = {}
    for s in cut_specs:
        cut_by_stage.setdefault(s.stage_key, []).append(s)

    stage_hists = _init_stage_hists(cut_specs, final_var_configs)
    gates_total = 0.0

    for i in tqdm(range(n_use), desc=f"{sample}:{path.basename(df_file)}", leave=False):
        hdr = pd.read_hdf(df_file, key=f"hdr_{i}")
        gates_total += _gate_sum_from_hdr(hdr, sample)

        evt = pd.read_hdf(df_file, key=f"evt_{i}")
        trk = pd.read_hdf(df_file, key=f"trk_{i}")

        state: Dict[str, Any] = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}
        for stage_key, post_state in walk_pipeline(state, pipeline_sample,
                                                   trace=vlog if verbose and i == 0 else None):
            # Cut-variable histograms for this stage (if any).
            for spec in cut_by_stage.get(stage_key, ()):
                got = get_var_series(post_state, spec.var_config, spec.target)
                if got is None:
                    continue
                values, _idx = got
                stage_hists[stage_key][spec.var_config.var_save_name] += histogram_var(
                    values, np.asarray(spec.var_config.bins)
                )
            # Final-stage variables.
            if stage_key == FINAL_STAGE_KEY:
                evt_final = post_state.get("evt")
                if evt_final is not None and len(evt_final) > 0:
                    for vc in final_var_configs:
                        if vc.var_save_name == "integrated":
                            stage_hists[FINAL_STAGE_KEY][vc.var_save_name] += np.array(
                                [float(len(evt_final))], dtype=np.float64
                            )
                        else:
                            got = get_var_series(post_state, vc, "evt")
                            if got is None:
                                continue
                            values, _ = got
                            stage_hists[FINAL_STAGE_KEY][vc.var_save_name] += histogram_var(
                                values, np.asarray(vc.bins)
                            )

        del state, evt, trk, hdr
        gc.collect()

    var_save_names: List[str] = []
    for stage_key, var_map in stage_hists.items():
        for vsn in var_map.keys():
            if vsn not in var_save_names:
                var_save_names.append(vsn)

    return {
        "kind": "cosmics_syst_chunk",
        "input_stage": "sel_all",
        "sample": sample,
        "df_file": df_file,
        "splits_processed": n_use,
        "gates": gates_total,
        "stage_hists": stage_hists,
        "stage_var_save_names": {
            stage_key: list(var_map.keys()) for stage_key, var_map in stage_hists.items()
        },
        "var_save_names": var_save_names,
    }


# ===========================================================================
# CLI
# ===========================================================================
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True, choices=("offbeam", "intime"))
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument(
        "--input-stage",
        choices=("final", "sel_all"),
        default="final",
        help="``final``: histogram already-selected evt columns (legacy). "
             "``sel_all``: re-run the event selection pipeline on raw evt+trk+hdr "
             "and save histograms at every cut stage AND at the final stage.",
    )
    p.add_argument("--max-splits", type=int, default=0, help="Cap HDF splits (0 = all).")
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Log each HDF open stage (split key, row counts for first splits).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.verbose:
        print(
            "[cosmics-chunk] main start sample=%s input_stage=%s out_dir=%s max_splits=%d"
            % (args.sample, args.input_stage, args.out_dir, args.max_splits),
            flush=True,
        )

    if args.input_stage == "sel_all":
        payload = accumulate_file_sel_all(
            args.df_file,
            args.sample,
            max_splits=args.max_splits,
            verbose=args.verbose,
        )
    else:
        var_configs = build_variable_configs(None)
        payload = accumulate_file_final(
            args.df_file,
            args.sample,
            var_configs,
            max_splits=args.max_splits,
            verbose=args.verbose,
        )

    stem = path.splitext(path.basename(args.df_file))[0]
    out_path = path.join(args.out_dir, "cosmics__%s__%s.pkl" % (args.sample, stem))
    with open(out_path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(
        "[cosmics-chunk] wrote %s  input_stage=%s gates=%.6e  splits=%d"
        % (out_path, payload["input_stage"], payload["gates"], payload["splits_processed"]),
        flush=True,
    )


if __name__ == "__main__":
    main()
