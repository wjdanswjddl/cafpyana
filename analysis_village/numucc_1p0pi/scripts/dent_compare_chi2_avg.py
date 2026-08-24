#!/usr/bin/env python3
"""
DENT vs CV plots for plane-averaged chi2 (``avg_chi2``) at the ``2prong-vtxdist``
stage — the same quantity used by ``get_mu_p_candidate`` PID cuts.

Variables:
  * chi2_avg_mu  — avg of chi2_muon over I0/I1/I2 (zeros excluded)
  * chi2_avg_p   — avg of chi2_proton over I0/I1/I2 (zeros excluded)

Regions (same as ``dent_compare_regions.py``):
  full sample, 8 octants (x=0,y=0,z=250), halves (x=0 → E/W)

Outputs overlays + per-event (DENT−CV)/CV + DENT/CV, and caches paired
value arrays for replot / metric recalculation without re-reading dfs.

Example:
    python dent_compare_chi2_avg.py --max-files 3
"""

from __future__ import annotations

import argparse
import gc
import glob
import os
import sys
import warnings
from os import makedirs, path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
os.environ.setdefault("MPLBACKEND", "Agg")

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.event_selection_batch_core import (
    attach_intrinsic_weights,
    ensure_phi_and_kinematics_cols,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.syst_pipeline_walker import walk_pipeline
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from makedf.util import avg_chi2
from pyanalib.split_df_helpers_new import get_n_split

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.scripts import dent_compare_regions as reg

_OUT_BASE = (
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT-chi2avg"
)

STAGE_KEY = "2prong-vtxdist"
CHI2_AVG_MU_COL = ("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")
CHI2_AVG_P_COL = ("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")


def build_chi2_avg_var_defs() -> Dict[str, dict]:
    defs: Dict[str, dict] = {}
    for vc in (VariableConfig.chi2_avg_mu(), VariableConfig.chi2_avg_proton()):
        defs[vc.var_save_name] = {
            "label": vc.var_labels[0] if vc.var_labels else vc.var_save_name,
            "bins": np.asarray(vc.bins),
            "col": vc.var_evt_reco_col,
        }
    return defs


def _trk_chi2_avg(trk_block: pd.DataFrame, which: str) -> np.ndarray:
    """Return avg chi2 array for one trk1/trk2 block; compute on the fly if needed."""
    col = CHI2_AVG_MU_COL if which == "muon" else CHI2_AVG_P_COL
    try:
        return multicol_get_series(trk_block, col).to_numpy(dtype=float)
    except Exception:
        pass
    # Fallback: same as selections.get_mu_p_candidate / _attach_chi2_avgs
    name = "chi2_muon" if which == "muon" else "chi2_proton"
    return avg_chi2(trk_block, name).to_numpy(dtype=float)


def _extract_chi2_avg_pair(evt: pd.DataFrame) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Return {var: (v_trk1, v_trk2)} for chi2_avg_mu / chi2_avg_p."""
    return {
        "chi2_avg_mu": (
            _trk_chi2_avg(evt.trk1, "muon"),
            _trk_chi2_avg(evt.trk2, "muon"),
        ),
        "chi2_avg_p": (
            _trk_chi2_avg(evt.trk1, "proton"),
            _trk_chi2_avg(evt.trk2, "proton"),
        ),
    }


def process_sel_all_chi2_avg(
    df_file: str,
    *,
    accum: reg.RegionAccum,
    collect_keyed: bool,
) -> None:
    n_split = get_n_split(df_file)
    for i in range(n_split):
        split: Dict[str, Optional[pd.DataFrame]] = {}
        for key in ("evt", "trk", "hdr"):
            try:
                split[key] = pd.read_hdf(df_file, key=f"{key}_{i}")
            except Exception:
                split[key] = None
        hdr = split.get("hdr")
        evt = split.get("evt")
        trk = split.get("trk")
        if evt is None or len(evt) == 0:
            continue

        entry_table = dc._sel_all_entry_key_table(hdr, evt) if collect_keyed else None
        attach_intrinsic_weights(evt, trk, "mc", use_mc_genweight=False)
        evt, _ = ensure_phi_and_kinematics_cols(evt, trk, None)
        state = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}

        for stage_key, cur in walk_pipeline(state, sample="mc"):
            if stage_key != STAGE_KEY:
                continue
            cur_evt = cur.get("evt")
            if cur_evt is None or len(cur_evt) == 0:
                continue
            try:
                top = cur_evt.columns.get_level_values(0).unique()
            except Exception:
                continue
            if "trk1" not in top or "trk2" not in top:
                continue

            try:
                vx, vy, vz = reg._vertex_xyz(cur_evt)
                oct_ids = reg.octant_ids_from_xyz(vx, vy, vz)
                half_ids = reg.half_ids_from_x(vx)
            except Exception:
                n = len(cur_evt)
                oct_ids = np.full(n, -1, dtype=np.int8)
                half_ids = np.full(n, -1, dtype=np.int8)

            base_keys: Optional[List[Optional[dc.PairKey]]] = None
            if collect_keyed:
                base_keys = dc._evt_base_keys(cur_evt, entry_table)

            try:
                vals_by_var = _extract_chi2_avg_pair(cur_evt)
            except Exception as exc:
                print(f"  chi2_avg extract failed {path.basename(df_file)}: {exc}", flush=True)
                continue

            for var_name, (v1, v2) in vals_by_var.items():
                if var_name not in accum.var_defs:
                    continue
                if collect_keyed and base_keys is not None:
                    accum.fill_evt_var(
                        var_name, v1, oct_ids, half_ids, base_keys,
                        collect_keyed=True, trk_slot=0,
                    )
                    accum.fill_evt_var(
                        var_name, v2, oct_ids, half_ids, base_keys,
                        collect_keyed=True, trk_slot=1,
                    )
                else:
                    vals = np.concatenate([v1, v2])
                    oct_both = np.concatenate([oct_ids, oct_ids])
                    half_both = np.concatenate([half_ids, half_ids])
                    accum.fill_evt_var(
                        var_name, vals, oct_both, half_both, None,
                        collect_keyed=False,
                    )
            break  # only need this stage

        del split, state
        gc.collect()


def process_variation(
    variation: str,
    matched_dir: str,
    *,
    var_defs: Dict[str, dict],
    max_files: Optional[int],
    collect_keyed: bool = True,
) -> reg.RegionAccum:
    files = dc.list_matched_files(matched_dir, "sel_all")
    if not files:
        files = sorted(glob.glob(path.join(matched_dir, "*sel_all*.df")))
        files = [f for f in files if "_matched" not in path.basename(f)]
    if max_files is not None:
        files = files[:max_files]
    print(f"[{variation}] {len(files)} sel_all matched files", flush=True)

    accum = reg.RegionAccum(var_defs)
    for fpath in tqdm(files, desc=f"chi2avg {variation}"):
        process_sel_all_chi2_avg(
            fpath, accum=accum, collect_keyed=collect_keyed,
        )
    return accum


def parse_args():
    p = argparse.ArgumentParser(
        description="DENT CV vs DENT: chi2_avg (mu/p) at 2prong-vtxdist"
    )
    p.add_argument("--out-base", default=_OUT_BASE)
    p.add_argument("--cv-all-dir", default=dc.DEFAULT_DIRS["cv_all"])
    p.add_argument("--dent-all-dir", default=dc.DEFAULT_DIRS["dent_all"])
    p.add_argument("--max-files", type=int, default=None)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = path.join(args.out_base, "cache")
    fig_dir = path.join(args.out_base, "plots")
    makedirs(cache_dir, exist_ok=True)

    var_defs = build_chi2_avg_var_defs()
    print(
        f"Stage={STAGE_KEY}; vars={list(var_defs)}; "
        f"regions=all + {reg.N_OCT} octants + {reg.N_HALF} halves",
        flush=True,
    )

    accum_cv = process_variation(
        "cv", args.cv_all_dir, var_defs=var_defs, max_files=args.max_files,
    )
    accum_dent = process_variation(
        "dent", args.dent_all_dir, var_defs=var_defs, max_files=args.max_files,
    )
    reg.plot_and_package_stage(
        "sel_all_2prong-vtxdist",
        var_defs,
        accum_cv,
        accum_dent,
        fig_dir=fig_dir,
        cache_path=path.join(cache_dir, "dent_chi2_avg_regions.pkl"),
    )

    print(f"\nDone. Outputs under {args.out_base}", flush=True)
    print("  plots/all/, plots/octants/, plots/halves/", flush=True)
    print(
        "  cache/dent_chi2_avg_regions.pkl  "
        "(hists + paired cv/dent/frac/ratio per region)",
        flush=True,
    )
    print(
        "  Replot: python dent_replot_regions.py "
        f"--cache {path.join(cache_dir, 'dent_chi2_avg_regions.pkl')} "
        f"--fig-dir {fig_dir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
