#!/usr/bin/env python
"""Detector (calorimetry) unisim chunk processor.

This is the **map** phase for the WireMod + calo unisim systematic. The input
``.df`` files are the variation tables produced by configs like
``configs/numucc_1p0pi/sel_2prong-updatecalo.py``: each HDF5 split holds a
``evt_cv_<i>`` central value plus eight per-calo variations
``evt_<calovar>_{p,m}_<i>`` (calovar in {ccal, alpha, beta, R}). All evt frames
are already past the 2-prong selection (``trk1`` / ``trk2`` are saved on each
slice). The recalculated PID columns carry the ``_new`` suffix.

For every universe in a split, we re-run the rest of the selection
(2prong-contained -> trackscore -> vtxdist -> mu/p candidates -> kinematics ->
TKI), histogram each plot variable per stage, and accumulate per-bin POT-weighted
counts. One pickle is written per .df file and aggregated downstream by
``syst_detvar_aggregate.py``.

Usage
-----
    python syst_detvar_chunk.py \\
        --df_file PATH.df \\
        --wiremod_tag wiremod_yz \\
        --out_dir OUT_DIR

Pickles are saved as ``<wiremod_tag>__<basename>.pkl`` so the aggregator can
group by WireMod model.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
import warnings
from os import path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable=None, **kwargs):
        if iterable is None:

            class _Dummy:
                def __enter__(self):
                    return self

                def __exit__(self, *args):
                    pass

                def update(self, *args, **kwargs):
                    pass

            return _Dummy()
        return iterable


warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.makedf.selections import (
    cut_2prong_contained, cut_2prong_trackscore, cut_2prong_vtxdist,
    get_mu_p_candidate, cut_has_mu, cut_has_p,
    cut_mu_kinematics, cut_p_kinematics,
    NU_SCORE_TH, TRACKSCORE_TH, VTXDIST_TH,
    MU_CHI2MU_TH, MU_CHI2P_TH, MU_LEN_TH, QUAL_TH, P_CHI2P_TH, P_LEN_TH,
    MU_PLO_TH, MU_PHI_TH, P_PLO_TH, P_PHI_TH,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    FINAL_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.categories import DETECTOR
from analysis_village.numucc_1p0pi.selection_framework import (
    multicol_get_series, multicol_resolve_column_key,
)
from pyanalib.variable_calculator import add_reco_cc1p0pi_tki_evtdf
from pyanalib.split_df_helpers import get_n_split
from makedf.util import avg_chi2


# ---------------------------------------------------------------------------
# Universe layout
# ---------------------------------------------------------------------------
CALO_PARAMS: Tuple[str, ...] = ("ccal", "alpha", "beta", "R")
SHIFTS: Tuple[str, ...] = ("p", "m")
UNIVERSES: Tuple[str, ...] = ("cv",) + tuple(
    f"{c}_{s}" for c in CALO_PARAMS for s in SHIFTS
)


def evt_key_for_universe(univ: str, split_idx: int) -> str:
    """HDF5 key inside one split for a universe label (cv/ccal_p/...)."""
    return f"evt_{univ}_{split_idx}"


# ---------------------------------------------------------------------------
# Plot variable definitions
# ---------------------------------------------------------------------------
# For varied samples we want the chi2 averages computed from the **_new** columns
# (i.e. with updated calo). We store them under a chi2_*_avg slot that points to
# the recomputed values, and use dedicated VariableConfigs that consume that slot.
def _chi2_mu_new_var_config():
    base = VariableConfig.chi2_mu()
    # use the avg_new slot we attach below
    base.var_evt_reco_col = ('pfp', 'trk', 'chi2pid', 'avg', 'chi2_muon_new', '')
    return base


def _chi2_p_new_var_config():
    base = VariableConfig.chi2_proton()
    base.var_evt_reco_col = ('pfp', 'trk', 'chi2pid', 'avg', 'chi2_proton_new', '')
    return base


# Per-track plots taken on concat([trk1, trk2]).
PER_TRK_PLOTS: List[VariableConfig] = [
    VariableConfig.track_score(),
    VariableConfig.vtx_dist(),
    VariableConfig.trk_len(),
    VariableConfig.mcs_range_diff(),
    _chi2_mu_new_var_config(),
    _chi2_p_new_var_config(),
]

# Per-event plots taken on the surviving evt dataframe (core kinematics + Eν + merged finals).
PER_EVT_PLOTS: List[VariableConfig] = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
_ne_cfg = VariableConfig.neutrino_energy()
if _ne_cfg.var_save_name not in {c.var_save_name for c in PER_EVT_PLOTS}:
    PER_EVT_PLOTS.append(_ne_cfg)

_MORE_FINAL_EVT_VARS = FINAL_SELECTED_EVT_VARIABLE_CONFIGS
_seen_evt_save = {vc.var_save_name for vc in PER_EVT_PLOTS}
for _vc in _MORE_FINAL_EVT_VARS:
    if _vc.var_save_name not in _seen_evt_save:
        PER_EVT_PLOTS.append(_vc)
        _seen_evt_save.add(_vc.var_save_name)


# Stages run AFTER the 2prong saved state. Each entry is
# (stage_key, label, cut_fn(df) -> df, plot_specs).
# ``plot_specs`` is a list of ``("trk"|"evt", VariableConfig)`` pairs evaluated
# on the post-cut state.
@dataclass
class StageDef:
    key: str
    label: str
    cut: Optional[callable] = None  # df -> df
    plots: List[Tuple[str, VariableConfig]] = field(default_factory=list)


def _attach_chi2_avg_new(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the calo-aware chi2 averages on trk1 / trk2 columns.

    Mirrors the notebook block that fills
    ``("trkN", "pfp", "trk", "chi2pid", "avg", "chi2_{muon,proton}_new", "")``.
    """
    if df is None or len(df) == 0:
        return df
    # cut_2prong_vtxdist returns a slice; copy so the new col assignment doesn't
    # trip pandas' SettingWithCopy warnings (the result is not used elsewhere).
    df = df.copy()
    for trk_idx in (1, 2):
        trk_top = f"trk{trk_idx}"
        if trk_top not in df.columns.get_level_values(0).unique():
            continue
        for part in ("muon", "proton"):
            try:
                this_avg = avg_chi2(df[trk_top], f"chi2_{part}_new")
            except Exception:
                continue
            df.loc[:, (trk_top, "pfp", "trk", "chi2pid", "avg", f"chi2_{part}_new", "")] = this_avg
    return df


def build_stages() -> List[StageDef]:
    """Selection stages applied to each variation universe.

    Mirrors the cut sequence in
    ``analysis_village/numucc_1p0pi/event_selection_pipeline_def.build_pipeline``
    starting AFTER the 2prong cut. The PID candidate finder uses
    ``score_tag="_new"`` so the calo-varied chi2 columns drive mu/p choice.
    """
    stages: List[StageDef] = []

    # entry stage: data is already 2prong + trk1/trk2 attached
    stages.append(StageDef(
        key="2prong",
        label="At 2-prong (input)",
        cut=lambda df: df,
        plots=[],
    ))

    # contained
    stages.append(StageDef(
        key="2prong-contained",
        label="Both PFPs contained",
        cut=lambda df: cut_2prong_contained(df, det=DETECTOR),
        plots=[("trk", VariableConfig.track_score())],
    ))

    # track-score cut
    stages.append(StageDef(
        key="2prong-trackscore",
        label=f"Both tracks track-score > {TRACKSCORE_TH}",
        cut=lambda df: cut_2prong_trackscore(df, TRACKSCORE_TH),
        plots=[("trk", VariableConfig.vtx_dist())],
    ))

    # vtxdist cut + attach chi2 averages
    def _vtxdist_then_chi2(df):
        df = cut_2prong_vtxdist(df, VTXDIST_TH)
        df = _attach_chi2_avg_new(df)
        return df

    stages.append(StageDef(
        key="2prong-vtxdist",
        label=f"Both tracks vtxdist < {VTXDIST_TH} cm",
        cut=_vtxdist_then_chi2,
        plots=[
            ("trk", VariableConfig.trk_len()),
            ("trk", VariableConfig.mcs_range_diff()),
            ("trk", _chi2_mu_new_var_config()),
            ("trk", _chi2_p_new_var_config()),
        ],
    ))

    # mu candidate + kinematics
    def _muX_cut(df):
        df = get_mu_p_candidate(
            df,
            mu_chi2mu_th=MU_CHI2MU_TH, mu_chi2p_th=MU_CHI2P_TH,
            mu_len_th=MU_LEN_TH, qual_th=QUAL_TH,
            p_chi2mu_th=-1, p_chi2p_th=P_CHI2P_TH, p_len_th=P_LEN_TH,
            score_tag="_new",
        )
        df = cut_has_mu(df)
        df = cut_mu_kinematics(df, mu_Plo_th=MU_PLO_TH, mu_Phi_th=MU_PHI_TH)
        return df

    stages.append(StageDef(
        key="2prong-muX",
        label="One track muon-like",
        cut=_muX_cut,
        plots=[],
    ))

    # proton candidate + kinematics + TKI
    def _mup_cut(df):
        df = cut_has_p(df)
        df = cut_p_kinematics(df, p_Plo_th=P_PLO_TH, p_Phi_th=P_PHI_TH)
        df = add_reco_cc1p0pi_tki_evtdf(df)
        return df

    stages.append(StageDef(
        key="2prong-mup",
        label="The other proton-like (final)",
        cut=_mup_cut,
        plots=[("evt", vc) for vc in PER_EVT_PLOTS],
    ))

    return stages


# ---------------------------------------------------------------------------
# Per-(stage, universe, var_save_name) histogram accumulator
# ---------------------------------------------------------------------------
@dataclass
class UniHistAcc:
    """Per-bin POT-weighted counts and squared-weights, indexed by universe."""
    bins: np.ndarray
    universes: List[str]
    hist: np.ndarray = field(default=None)   # (n_univ, n_bin)
    err2: np.ndarray = field(default=None)   # (n_univ, n_bin)

    def __post_init__(self):
        if self.hist is None:
            self.hist = np.zeros((len(self.universes), len(self.bins) - 1))
            self.err2 = np.zeros_like(self.hist)

    def fill(self, univ: str, values: np.ndarray, weights: np.ndarray):
        if values is None or len(values) == 0:
            return
        if univ not in self.universes:
            return
        i = self.universes.index(univ)
        v = np.asarray(values, dtype=float)
        w = np.asarray(weights, dtype=float)
        # numpy.histogram(weights=w) propagates a single NaN to every bin -> sanitize
        w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
        # also clip values to bin range so under/overflows accumulate at the edge
        eps = (self.bins[-1] - self.bins[0]) * 1e-9
        v = np.clip(v, self.bins[0], self.bins[-1] - eps)
        finite = np.isfinite(v)
        if not finite.all():
            v = v[finite]
            w = w[finite]
        h, _ = np.histogram(v, bins=self.bins, weights=w)
        e2, _ = np.histogram(v, bins=self.bins, weights=np.square(w))
        self.hist[i] += h
        self.err2[i] += e2

    def __iadd__(self, other: "UniHistAcc"):
        assert np.array_equal(self.bins, other.bins)
        assert self.universes == other.universes
        self.hist += other.hist
        self.err2 += other.err2
        return self


# Key: (stage_key, var_save_name, "trk"|"evt")
HistDict = Dict[Tuple[str, str, str], UniHistAcc]


def _get_var_for_plot(df: pd.DataFrame, var_config: VariableConfig, target: str
                      ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Pluck the variable column and weights from ``df`` (or concat of trks)."""
    if df is None or len(df) == 0:
        return None
    if target == "trk":
        if "trk1" not in df.columns.get_level_values(0).unique():
            return None
        if "trk2" not in df.columns.get_level_values(0).unique():
            return None
        try:
            v1 = multicol_get_series(df.trk1, var_config.var_evt_reco_col)
            v2 = multicol_get_series(df.trk2, var_config.var_evt_reco_col)
        except KeyError:
            return None
        var = pd.concat([v1, v2])
    elif target == "evt":
        try:
            var = multicol_get_series(df, var_config.var_evt_reco_col)
        except KeyError:
            return None
    else:
        raise ValueError(target)
    var = np.asarray(var, dtype=float)
    weights = np.ones_like(var, dtype=float)
    return var, weights


# ---------------------------------------------------------------------------
# Main pipeline driver
# ---------------------------------------------------------------------------
def process_split(
    df_file: str,
    split_idx: int,
    stages: List[StageDef],
    histdict: HistDict,
    nevt_table: Dict[Tuple[str, str], int],
    verbose: bool = False,
    progress: bool = True,
):
    """Read one split of the HDF and accumulate histograms over all universes."""
    univ_iter = tqdm(
        UNIVERSES,
        desc=f"split {split_idx}",
        unit="univ",
        leave=False,
        disable=not progress,
    )
    for univ in univ_iter:
        key = evt_key_for_universe(univ, split_idx)
        try:
            df = pd.read_hdf(df_file, key=key)
        except KeyError:
            if verbose:
                print(f"  [warn] missing key {key} in {df_file}", flush=True)
            continue
        if df is None or len(df) == 0:
            continue

        # walk stages
        cur = df
        for stage in stages:
            if stage.cut is not None:
                try:
                    cur = stage.cut(cur)
                except Exception as e:
                    print(f"  [error] stage={stage.key} univ={univ}: {e}", flush=True)
                    cur = cur.iloc[0:0]
                    break
            # bookkeeping
            nevt_table[(stage.key, univ)] = nevt_table.get((stage.key, univ), 0) + int(len(cur))
            # plots
            for target, vc in stage.plots:
                got = _get_var_for_plot(cur, vc, target)
                if got is None:
                    continue
                values, weights = got
                k = (stage.key, vc.var_save_name, target)
                if k not in histdict:
                    histdict[k] = UniHistAcc(
                        bins=np.asarray(vc.bins).copy(),
                        universes=list(UNIVERSES),
                    )
                histdict[k].fill(univ, values, weights)

        del cur, df
    gc.collect()


# ---------------------------------------------------------------------------
# Exposure metadata
# ---------------------------------------------------------------------------
def read_chunk_pot(df_file: str, n_split: int) -> float:
    pot = 0.0
    for i in range(n_split):
        try:
            hp = pd.read_hdf(df_file, key=f"histpotdf_{i}")
            if "TotalPOT" in hp.columns:
                pot += float(hp["TotalPOT"].sum())
        except Exception:
            pass
    return pot


def read_chunk_genevts(df_file: str, n_split: int) -> float:
    n = 0.0
    for i in range(n_split):
        try:
            hg = pd.read_hdf(df_file, key=f"histgenevtdf_{i}")
            if "TotalGenEvents" in hg.columns:
                n += float(hg["TotalGenEvents"].sum())
        except Exception:
            pass
    return n


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--df_file", required=True)
    p.add_argument("--wiremod_tag", required=True,
                   help="Short label for the WireMod model (e.g. wiremod_yz, wiremod_xtxw)")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--max_splits", type=int, default=0,
                   help="Cap HDF5 splits processed (0 = all)")
    p.add_argument("--verbose", action="store_true")
    p.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm bars (clean logs / batch systems)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    _main_detvar(args)


def _main_detvar(args) -> None:
    os.makedirs(args.out_dir, exist_ok=True)

    n_split = int(get_n_split(args.df_file))
    n_use = min(args.max_splits, n_split) if args.max_splits > 0 else n_split
    print(f"[detvar-chunk] file={args.df_file} wiremod={args.wiremod_tag} "
          f"n_split={n_split} use={n_use}", flush=True)

    stages = build_stages()
    histdict: HistDict = {}
    nevt_table: Dict[Tuple[str, str], int] = {}

    stem = path.splitext(path.basename(args.df_file))[0]
    show_prog = not args.no_progress
    for i in tqdm(
        range(n_use),
        desc=f"{args.wiremod_tag}:{stem}",
        unit="split",
        disable=not show_prog,
    ):
        if args.verbose:
            print(f"  split {i + 1}/{n_use}", flush=True)
        process_split(
            args.df_file,
            i,
            stages,
            histdict,
            nevt_table,
            verbose=args.verbose,
            progress=show_prog,
        )

    chunk_pot = read_chunk_pot(args.df_file, n_use)
    chunk_genevts = read_chunk_genevts(args.df_file, n_use)

    out_payload = {
        "wiremod_tag": args.wiremod_tag,
        "df_file": args.df_file,
        "n_splits": n_use,
        "n_splits_in_file": n_split,
        "chunk_pot": chunk_pot,
        "chunk_genevts": chunk_genevts,
        "universes": list(UNIVERSES),
        "stages": [(s.key, s.label) for s in stages],
        "histdict": histdict,
        "nevt_table": nevt_table,
    }

    base = path.splitext(path.basename(args.df_file))[0]
    out_path = path.join(args.out_dir, f"{args.wiremod_tag}__{base}.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(out_payload, f)
    print(f"[detvar-chunk] wrote {out_path}  pot={chunk_pot:.3e}  "
          f"genevts={chunk_genevts:.3e}", flush=True)
    # sanity: print final-stage counts per universe
    final_key = stages[-1].key
    counts = {u: nevt_table.get((final_key, u), 0) for u in UNIVERSES}
    print(f"[detvar-chunk] final-stage counts (raw): {counts}", flush=True)


if __name__ == "__main__":
    main()
