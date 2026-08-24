#!/usr/bin/env python3
"""
DENT detector unisim comparison: matched CV vs DENT histograms, plots, and
selection efficiency / purity summary.

Unlike WireMod / SCE (matched at sel_mup / sel_2prong), DENT is evaluated from
``sel_all`` onward because the uncertainty affects early selection variables.

Outputs (parallel to WireMod / SCE):
  /exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT/
    cache/   — histogram pickles, summary CSV, per-event frac-diff hists
    plots/   — CV vs DENT overlays + per-event (DENT-CV)/CV distributions
    plots/sel_mup/ — final-selected variables (when sel_mup matched files exist)

Example (test on 3 files, using available DENT production if canonical dirs empty):
    python dent_compare.py --max-files 3 \\
        --dent-all-dir /pnfs/.../2026_08_18_120607__sel_all-mc-DENT \\
        --dent-mup-dir /pnfs/.../2026_08_18_120753__sel_mup-mc-DENT
"""

from __future__ import annotations

import argparse
import gc
import glob
import os
import pickle
import sys
import warnings
from dataclasses import dataclass, field
from os import makedirs, path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
os.environ.setdefault("MPLBACKEND", "Agg")

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.categories import DETECTOR, IsNuInFV_NumuCC_1p0pi
from analysis_village.numucc_1p0pi.event_selection_batch_core import (
    attach_intrinsic_weights,
    ensure_phi_and_kinematics_cols,
    hdr_chunk_pot,
)
from analysis_village.numucc_1p0pi.event_selection_batch_core import hdf_has_mcnu
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    get_var_series,
    histogram_var,
    walk_pipeline,
)
from analysis_village.numucc_1p0pi.makedf.selections import (
    MU_CHI2MU_TH,
    MU_CHI2P_TH,
    MU_LEN_TH,
    MU_PHI_TH,
    MU_PLO_TH,
    NU_SCORE_TH,
    P_PHI_TH,
    P_PLO_TH,
    QUAL_TH,
    TRACKSCORE_TH,
    VTXDIST_TH,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from pyanalib.split_df_helpers_new import get_n_split

# Cut markers for diagnostic overlays: [[x, direction], ...]
# direction 0 = keep left (<), 1 = keep right (>); matches event_selection_pipeline_def.
DIAGNOSTIC_CUT_VLINES: Dict[str, List[List[float]]] = {
    "nu_score": [[NU_SCORE_TH, 1]],
    "track_score": [[TRACKSCORE_TH, 1]],
    "vtx_dist": [[VTXDIST_TH, 0]],
    "trk_len": [[MU_LEN_TH, 1]],
    "mcs_range_diff": [[-QUAL_TH, 0], [QUAL_TH, 1]],
    "chi2_mu": [[MU_CHI2MU_TH, 0]],
    "chi2_p": [[MU_CHI2P_TH, 1]],
    "chi2_avg_mu": [[MU_CHI2MU_TH, 0]],
    "chi2_avg_p": [[MU_CHI2P_TH, 1]],
    "muon-p": [[MU_PLO_TH, 1], [MU_PHI_TH, 0]],
    "proton-p": [[P_PLO_TH, 1], [P_PHI_TH, 0]],
}

_DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"
_OUT_BASE = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT"

DEFAULT_DIRS = {
    "cv_all": f"{_DFS_ROOT}/2026_08_19_031254__sel_all-mc-CV",
    "dent_all": f"{_DFS_ROOT}/2026_08_18_120607__sel_all-mc-DENT",
    "cv_mup": f"{_DFS_ROOT}/2026_08_19_031423__sel_mup-mc-CV",
    "dent_mup": f"{_DFS_ROOT}/2026_08_18_120753__sel_mup-mc-DENT",
}

VARIATIONS = ("cv", "dent")
VARIATION_LABELS = {"cv": "CV", "dent": "DENT"}
VARIATION_COLORS = {"cv": "black", "dent": "C0"}
RATIO_PANEL_COLOR = "black"

# Per-event fractional difference: (x_DENT - x_CV) / x_CV for matched pairs.
# PairKey is (E, run, subrun, evt[, slc][, trk_slot]).
EventKey = Tuple[float, int, int, int]
PairKey = Tuple[Any, ...]
FRAC_DIFF_BINS = np.linspace(-1.5, 1.5, 61)
RATIO_BINS = np.linspace(0.0, 2.0, 61)
CV_ABS_EPS = 1e-9  # skip pairs with |CV| below this (avoids blow-ups)


# ---------------------------------------------------------------------------
# Histogram cache helpers
# ---------------------------------------------------------------------------
def save_hists(cache_path: str, payload: dict) -> None:
    makedirs(path.dirname(cache_path), exist_ok=True)
    safe = dict(payload)
    if "var_defs" in safe:
        safe["var_defs"] = {
            k: {"label": v.get("label", k), "bins": np.asarray(v["bins"])}
            for k, v in safe["var_defs"].items()
        }
    with open(cache_path, "wb") as fh:
        pickle.dump(safe, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved cache → {cache_path}", flush=True)


def load_hists(cache_path: str) -> dict:
    with open(cache_path, "rb") as fh:
        return pickle.load(fh)


def save_fig(fig, name: str, fig_dir: str, dpi: int = 150) -> None:
    makedirs(fig_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(path.join(fig_dir, f"{name}.{ext}"), dpi=dpi, bbox_inches="tight")


# ---------------------------------------------------------------------------
# sel_all cut-stage variable definitions (early selection)
# ---------------------------------------------------------------------------
def _per_trk_col(evt_df, col_tuple):
    out = []
    top_cols = evt_df.columns.get_level_values(0).unique()
    for trk in ("trk1", "trk2"):
        if trk not in top_cols:
            continue
        try:
            vals = multicol_get_series(evt_df[trk], col_tuple).to_numpy(dtype=float)
            out.append(vals)
        except Exception:
            pass
    return np.concatenate(out) if out else np.array([], dtype=float)


def _per_evt_col(evt_df, col_tuple):
    try:
        return multicol_get_series(evt_df, col_tuple).to_numpy(dtype=float)
    except Exception:
        return np.array([], dtype=float)


def build_sel_all_var_defs() -> Dict[str, dict]:
    defs: Dict[str, dict] = {}
    for spec in CUT_STAGE_VAR_SPECS:
        vc = spec.var_config
        name = vc.var_save_name
        if name in defs:
            continue
        col = vc.var_evt_reco_col
        if spec.target == "trk":
            extract = lambda df, col=col: _per_trk_col(df, col)
        else:
            extract = lambda df, col=col: _per_evt_col(df, col)
        defs[name] = {
            "label": vc.var_labels[0] if vc.var_labels else name,
            "bins": np.asarray(vc.bins),
            "extract": extract,
            "stage_key": spec.stage_key,
        }
    return defs


def build_mup_var_defs() -> Dict[str, dict]:
    from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
        CORE_SELECTED_EVT_VARIABLE_CONFIGS,
        with_final_selected_evt_variables,
    )

    defs: Dict[str, dict] = {}
    for vc in with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)):
        col = vc.var_evt_reco_col
        if vc.var_save_name == "integrated":
            defs["integrated"] = {
                "label": "Integrated",
                "bins": np.asarray(vc.bins),
                "extract": lambda df: np.full(len(df), 500.0, dtype=float),
            }
            continue
        defs[vc.var_save_name] = {
            "label": vc.var_labels[0] if vc.var_labels else vc.var_save_name,
            "bins": np.asarray(vc.bins),
            "extract": lambda df, col=col: _per_evt_col(df, col),
        }
    ne = VariableConfig.neutrino_energy()
    if ne.var_save_name not in defs:
        col = ne.var_evt_reco_col
        defs[ne.var_save_name] = {
            "label": ne.var_labels[0],
            "bins": np.asarray(ne.bins),
            "extract": lambda df, col=col: _per_evt_col(df, col),
        }
    return defs


# ---------------------------------------------------------------------------
# Pipeline metrics
# ---------------------------------------------------------------------------
@dataclass
class StageMetrics:
    n_evt_raw: int = 0
    n_signal_raw: int = 0
    n_evt_pot: float = 0.0
    n_signal_pot: float = 0.0
    n_truth_nu_raw: int = 0
    n_truth_nu_pot: float = 0.0


@dataclass
class SampleSummary:
    variation: str
    pot: float = 0.0
    n_matched_files: int = 0
    stages: Dict[str, StageMetrics] = field(default_factory=dict)
    final_stage_key: str = FINAL_STAGE_KEY

    def integrated_efficiency_pct(self) -> float:
        """Signal survival at final stage / truth signal at allreco (evt-based fallback)."""
        first = self.stages.get("allreco")
        last = self.stages.get(self.final_stage_key)
        if first is None or last is None:
            return float("nan")
        denom = first.n_signal_raw
        if denom <= 0:
            denom = first.n_truth_nu_raw
        if denom <= 0:
            return float("nan")
        return 100.0 * last.n_signal_raw / denom

    def integrated_purity_pct(self) -> float:
        last = self.stages.get(self.final_stage_key)
        if last is None or last.n_evt_raw <= 0:
            return float("nan")
        return 100.0 * last.n_signal_raw / last.n_evt_raw

    def mcnu_efficiency_pct(self) -> float:
        """Full efficiency when mcnu denominator is available."""
        last = self.stages.get(self.final_stage_key)
        first = self.stages.get("allreco")
        if last is None or first is None or first.n_truth_nu_raw <= 0:
            return float("nan")
        return 100.0 * last.n_signal_raw / first.n_truth_nu_raw


def _count_stage(state: dict, mcnu_df: Optional[pd.DataFrame]) -> StageMetrics:
    evt = state.get("evt")
    m = StageMetrics()
    if evt is None or len(evt) == 0:
        return m
    w = evt["pot_weight"].to_numpy(dtype=float) if "pot_weight" in evt.columns else np.ones(len(evt))
    sig = IsNuInFV_NumuCC_1p0pi(evt, detector=DETECTOR)
    sm = np.asarray(sig.values, dtype=bool) if hasattr(sig, "values") else np.asarray(sig, dtype=bool)
    m.n_evt_raw = int(len(evt))
    m.n_signal_raw = int(sm.sum())
    m.n_evt_pot = float(w.sum())
    m.n_signal_pot = float(w[sm].sum())
    if mcnu_df is not None and len(mcnu_df) > 0:
        nu_sig = IsNuInFV_NumuCC_1p0pi(mcnu_df, detector=DETECTOR)
        nu_sm = np.asarray(nu_sig.values, dtype=bool) if hasattr(nu_sig, "values") else np.asarray(nu_sig, dtype=bool)
        nu_w = (
            mcnu_df["pot_weight"].to_numpy(dtype=float)
            if "pot_weight" in mcnu_df.columns
            else np.ones(len(mcnu_df))
        )
        m.n_truth_nu_raw = int(nu_sm.sum())
        m.n_truth_nu_pot = float(nu_w[nu_sm].sum())
    return m


def source_file_pot(matched_df_file: str, n_split: int) -> float:
    """
    Read POT from the *original* (unmatched) source file that corresponds to a
    ``_matched.df`` file.  The original file carries ``histpotdf_{i}`` with a
    single column ``TotalPOT`` — the same value used by the rest of the
    framework.  The matched file only keeps ``evt_*`` / ``meta_*``.

    Reads via PyTables (block0_values) to avoid HDFStore overhead.
    """
    import tables as _tb

    orig = matched_df_file.replace("_matched.df", ".df")
    if not path.isfile(orig):
        return 0.0
    pot = 0.0
    try:
        with _tb.open_file(orig, "r") as h:
            for key in h.root._v_children:
                if not key.startswith("histpotdf_"):
                    continue
                try:
                    vals = h.get_node(f"/{key}/block0_values")[:]
                    pot += float(np.nansum(vals))
                except Exception:
                    pass
    except Exception:
        pass
    return pot


def _stage_specs_by_key() -> Dict[str, List[Tuple[str, Any, str]]]:
    out: Dict[str, List[Tuple[str, Any, str]]] = {}
    for spec in CUT_STAGE_VAR_SPECS:
        out.setdefault(spec.stage_key, []).append(
            (spec.var_config.var_save_name, spec.var_config, spec.target)
        )
    return out


# ---------------------------------------------------------------------------
# Per-event keyed values for fractional-difference plots
# ---------------------------------------------------------------------------
def _flatten_meta_event_cols(meta: pd.DataFrame) -> pd.DataFrame:
    """Return flat columns __ntuple, entry, E, run, subrun, evt from meta."""
    mr = meta.reset_index()

    def _col(name: str) -> pd.Series:
        if name in mr.columns:
            return mr[name]
        for c in mr.columns:
            if isinstance(c, tuple) and c[0] == name:
                return mr[c]
        raise KeyError(name)

    return pd.DataFrame(
        {
            "__ntuple": _col("__ntuple").to_numpy(),
            "entry": _col("entry").to_numpy(),
            "E": _col("E").astype(float).to_numpy(),
            "run": _col("run").astype(int).to_numpy(),
            "subrun": _col("subrun").astype(int).to_numpy(),
            "evt": _col("evt").astype(int).to_numpy(),
        }
    )


def _sel_all_entry_key_table(hdr: pd.DataFrame, evt: pd.DataFrame) -> pd.DataFrame:
    """Map (__ntuple, entry) → (E, run, subrun, evt) for one sel_all split."""
    from analysis_village.numucc_1p0pi.scripts.dent_match_common_events import (
        _event_energy_map,
    )

    if hdr is None or len(hdr) == 0 or evt is None or len(evt) == 0:
        return pd.DataFrame(columns=["__ntuple", "entry", "E", "run", "subrun", "evt"])
    E_map = _event_energy_map(evt)
    hdr_r = hdr.reset_index()
    m = hdr_r.merge(E_map.reset_index(), on=["__ntuple", "entry"], how="inner")
    m = m.dropna(subset=["E"])
    return pd.DataFrame(
        {
            "__ntuple": m["__ntuple"].to_numpy(),
            "entry": m["entry"].to_numpy(),
            "E": m["E"].astype(float).to_numpy(),
            "run": m["run"].astype(int).to_numpy(),
            "subrun": m["subrun"].astype(int).to_numpy(),
            "evt": m["evt"].astype(int).to_numpy(),
        }
    )


def _evt_base_keys(evt: pd.DataFrame, entry_table: pd.DataFrame) -> List[Optional[PairKey]]:
    """One PairKey per evt row: (E, run, subrun, evt[, slc])."""
    n = len(evt)
    if n == 0 or entry_table is None or len(entry_table) == 0:
        return [None] * n
    left = pd.DataFrame(
        {
            "__ntuple": evt.index.get_level_values("__ntuple").to_numpy(),
            "entry": evt.index.get_level_values("entry").to_numpy(),
            "row_i": np.arange(n, dtype=np.int64),
        }
    )
    has_slc = "rec.slc..index" in (evt.index.names or [])
    if has_slc:
        left["slc"] = evt.index.get_level_values("rec.slc..index").to_numpy()
    m = left.merge(entry_table, on=["__ntuple", "entry"], how="left")
    m = m.sort_values("row_i")
    keys: List[Optional[PairKey]] = [None] * n
    E = m["E"].to_numpy()
    run = m["run"].to_numpy()
    subrun = m["subrun"].to_numpy()
    ev = m["evt"].to_numpy()
    row_i = m["row_i"].to_numpy(dtype=np.int64)
    slc = m["slc"].to_numpy() if has_slc else None
    for j, i in enumerate(row_i):
        if not np.isfinite(E[j]):
            continue
        base: PairKey = (float(E[j]), int(run[j]), int(subrun[j]), int(ev[j]))
        if has_slc:
            base = base + (int(slc[j]),)
        keys[int(i)] = base
    return keys


def _store_keyed_value(
    keyed_maps: Dict[str, Dict[PairKey, float]],
    var_name: str,
    key: PairKey,
    value: float,
) -> None:
    if key is None or not np.isfinite(value):
        return
    keyed_maps.setdefault(var_name, {})[key] = float(value)


def _store_evt_var_keyed(
    keyed_maps: Dict[str, Dict[PairKey, float]],
    var_name: str,
    keys: Sequence[Optional[PairKey]],
    values: np.ndarray,
) -> None:
    for key, val in zip(keys, values):
        if key is None:
            continue
        _store_keyed_value(keyed_maps, var_name, key, float(val))


def _store_trk_var_keyed(
    keyed_maps: Dict[str, Dict[PairKey, float]],
    var_name: str,
    keys: Sequence[Optional[PairKey]],
    v1: np.ndarray,
    v2: np.ndarray,
) -> None:
    """Pair tracks by slot: trk1 ↔ trk1, trk2 ↔ trk2 via key + slot."""
    for key, a, b in zip(keys, v1, v2):
        if key is None:
            continue
        _store_keyed_value(keyed_maps, var_name, key + (0,), float(a))
        _store_keyed_value(keyed_maps, var_name, key + (1,), float(b))


def compute_frac_diffs(
    cv_map: Dict[PairKey, float],
    dent_map: Dict[PairKey, float],
    *,
    cv_eps: float = CV_ABS_EPS,
) -> np.ndarray:
    """Return (DENT - CV) / CV for keys present in both maps with |CV| > eps."""
    cv_a, dent_a = compute_paired_values(cv_map, dent_map, cv_eps=cv_eps)
    if len(cv_a) == 0:
        return np.array([], dtype=float)
    return (dent_a - cv_a) / cv_a


def compute_paired_values(
    cv_map: Dict[PairKey, float],
    dent_map: Dict[PairKey, float],
    *,
    cv_eps: float = CV_ABS_EPS,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return aligned (cv, dent) arrays for common keys with |CV| > eps."""
    if not cv_map or not dent_map:
        return np.array([], dtype=float), np.array([], dtype=float)
    cv_out: List[float] = []
    dent_out: List[float] = []
    for key, cv_v in cv_map.items():
        if key not in dent_map:
            continue
        if not np.isfinite(cv_v) or abs(cv_v) <= cv_eps:
            continue
        dent_v = dent_map[key]
        if not np.isfinite(dent_v):
            continue
        cv_out.append(float(cv_v))
        dent_out.append(float(dent_v))
    return np.asarray(cv_out, dtype=float), np.asarray(dent_out, dtype=float)


def compute_ratios_from_paired(
    cv_vals: np.ndarray,
    dent_vals: np.ndarray,
) -> np.ndarray:
    """Per-event DENT / CV (inputs already filtered for |CV| > eps)."""
    if len(cv_vals) == 0:
        return np.array([], dtype=float)
    return np.asarray(dent_vals, dtype=float) / np.asarray(cv_vals, dtype=float)


def process_sel_all_file(
    df_file: str,
    *,
    hists: Dict[str, np.ndarray],
    var_defs: Dict[str, dict],
    summary: SampleSummary,
    stage_specs: Dict[str, List[Tuple[str, Any, str]]],
    keyed_maps: Optional[Dict[str, Dict[PairKey, float]]] = None,
) -> float:
    """Walk the selection pipeline on one matched sel_all file."""
    chunk_pot = 0.0
    n_split = get_n_split(df_file)

    for i in range(n_split):
        split: Dict[str, Optional[pd.DataFrame]] = {}
        for key in ("evt", "trk", "hdr"):
            try:
                split[key] = pd.read_hdf(df_file, key=f"{key}_{i}")
            except Exception:
                split[key] = None

        hdr = split.get("hdr")
        chunk_pot += hdr_chunk_pot(hdr)
        evt = split.get("evt")
        trk = split.get("trk")
        if evt is None or len(evt) == 0:
            continue

        entry_table = (
            _sel_all_entry_key_table(hdr, evt) if keyed_maps is not None else None
        )

        attach_intrinsic_weights(evt, trk, "mc", use_mc_genweight=False)
        evt, _ = ensure_phi_and_kinematics_cols(evt, trk, None)
        state = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}

        for stage_key, cur in walk_pipeline(state, sample="mc"):
            sm = _count_stage(cur, None)
            if stage_key not in summary.stages:
                summary.stages[stage_key] = StageMetrics()
            acc = summary.stages[stage_key]
            acc.n_evt_raw += sm.n_evt_raw
            acc.n_signal_raw += sm.n_signal_raw
            acc.n_evt_pot += sm.n_evt_pot
            acc.n_signal_pot += sm.n_signal_pot

            cur_evt = cur.get("evt")
            base_keys: Optional[List[Optional[PairKey]]] = None
            if keyed_maps is not None and cur_evt is not None and len(cur_evt) > 0:
                base_keys = _evt_base_keys(cur_evt, entry_table)

            for var_name, vc, target in stage_specs.get(stage_key, []):
                if var_name not in hists:
                    continue
                got = get_var_series(cur, vc, target)
                if got is None:
                    continue
                vals, _ = got
                hists[var_name] += histogram_var(vals, var_defs[var_name]["bins"])

                if keyed_maps is None or base_keys is None or cur_evt is None:
                    continue
                if target == "evt":
                    _store_evt_var_keyed(keyed_maps, var_name, base_keys, vals)
                elif target == "trk":
                    try:
                        v1 = multicol_get_series(
                            cur_evt.trk1, vc.var_evt_reco_col
                        ).to_numpy(dtype=float)
                        v2 = multicol_get_series(
                            cur_evt.trk2, vc.var_evt_reco_col
                        ).to_numpy(dtype=float)
                    except Exception:
                        continue
                    _store_trk_var_keyed(keyed_maps, var_name, base_keys, v1, v2)

        del split, state
        gc.collect()

    return chunk_pot


def process_sel_mup_file(
    df_file: str,
    *,
    hists: Dict[str, np.ndarray],
    var_defs: Dict[str, dict],
    summary: SampleSummary,
    keyed_maps: Optional[Dict[str, Dict[PairKey, float]]] = None,
) -> float:
    """Histogram final-selected ``evt_cv`` and record final-stage eff / purity."""
    n_split = get_n_split(df_file)
    has_mcnu = hdf_has_mcnu(df_file)
    # POT comes from histpotdf_* in the original (unmatched) source file.
    chunk_pot = source_file_pot(df_file, n_split)

    for i in range(n_split):
        split: Dict[str, Optional[pd.DataFrame]] = {}
        evt = None
        for evt_key in ("evt_cv", "evt"):
            try:
                evt = pd.read_hdf(df_file, key=f"{evt_key}_{i}")
                break
            except Exception:
                continue
        meta = None
        try:
            meta = pd.read_hdf(df_file, key=f"meta_{i}")
        except Exception:
            meta = None
        for key in ("hdr", "mcnu"):
            try:
                split[key] = pd.read_hdf(df_file, key=f"{key}_{i}")
            except Exception:
                split[key] = None

        hdr = split.get("hdr")
        mcnu = split.get("mcnu") if has_mcnu else None
        if evt is None or len(evt) == 0:
            continue

        attach_intrinsic_weights(evt, None, "mc", use_mc_genweight=False)
        evt, mcnu = ensure_phi_and_kinematics_cols(evt, None, mcnu)
        if mcnu is not None and len(mcnu) > 0 and "pot_weight" not in mcnu.columns:
            mcnu = mcnu.copy()
            mcnu["pot_weight"] = np.ones(len(mcnu), dtype=float)

        if mcnu is not None and len(mcnu) > 0 and "allreco" not in summary.stages:
            nu_sig = IsNuInFV_NumuCC_1p0pi(mcnu, detector=DETECTOR)
            nu_sm = np.asarray(nu_sig.values, dtype=bool)
            nu_w = mcnu["pot_weight"].to_numpy(dtype=float)
            summary.stages["allreco"] = StageMetrics(
                n_truth_nu_raw=int(nu_sm.sum()),
                n_truth_nu_pot=float(nu_w[nu_sm].sum()),
            )

        sm = _count_stage({"evt": evt}, None)
        if FINAL_STAGE_KEY not in summary.stages:
            summary.stages[FINAL_STAGE_KEY] = StageMetrics()
        acc = summary.stages[FINAL_STAGE_KEY]
        acc.n_evt_raw += sm.n_evt_raw
        acc.n_signal_raw += sm.n_signal_raw
        acc.n_evt_pot += sm.n_evt_pot
        acc.n_signal_pot += sm.n_signal_pot

        base_keys: Optional[List[Optional[PairKey]]] = None
        if keyed_maps is not None and meta is not None and len(meta) > 0:
            mf = _flatten_meta_event_cols(meta)
            left = pd.DataFrame(
                {
                    "__ntuple": evt.index.get_level_values("__ntuple").to_numpy(),
                    "entry": evt.index.get_level_values("entry").to_numpy(),
                    "E": multicol_get_series(
                        evt, ("mc", "E", "", "", "", "")
                    ).to_numpy(dtype=np.float32),
                    "row_i": np.arange(len(evt), dtype=np.int64),
                }
            )
            mf = mf.copy()
            mf["E"] = mf["E"].astype(np.float32)
            m = left.merge(mf, on=["__ntuple", "entry", "E"], how="left")
            m = m.sort_values("row_i")
            base_keys = [None] * len(evt)
            E = m["E"].to_numpy()
            run = m["run"].to_numpy()
            subrun = m["subrun"].to_numpy()
            ev = m["evt"].to_numpy()
            row_i = m["row_i"].to_numpy(dtype=np.int64)
            for j, i in enumerate(row_i):
                if not np.isfinite(run[j]):
                    continue
                base_keys[int(i)] = (
                    float(E[j]),
                    int(run[j]),
                    int(subrun[j]),
                    int(ev[j]),
                )

        for var_name, cfg in var_defs.items():
            try:
                vals = cfg["extract"](evt)
            except Exception:
                continue
            if len(vals) == 0:
                continue
            hists[var_name] += histogram_var(vals, cfg["bins"])
            if keyed_maps is not None and base_keys is not None and var_name != "integrated":
                _store_evt_var_keyed(keyed_maps, var_name, base_keys, vals)

        del split
        gc.collect()

    return chunk_pot


def list_matched_files(search_dir: str, filename_str: str) -> List[str]:
    return sorted(glob.glob(path.join(search_dir, f"*{filename_str}*_matched.df")))


def run_matching(args) -> None:
    import dent_match_common_events as dent_match

    for fmt, cv_dir, dent_dir, fstr in (
        ("sel_all", args.cv_all_dir, args.dent_all_dir, "sel_all"),
        ("sel_mup", args.cv_mup_dir, args.dent_mup_dir, "sel_mup"),
    ):
        if not path.isdir(cv_dir):
            print(f"[match] skip {fmt}: CV dir missing: {cv_dir}", flush=True)
            continue
        if not path.isdir(dent_dir):
            print(f"[match] skip {fmt}: DENT dir missing: {dent_dir}", flush=True)
            continue
        print(f"\n[match] {fmt}: CV={cv_dir}", flush=True)
        print(f"[match] {fmt}: DENT={dent_dir}", flush=True)
        argv = [
            "--format", fmt,
            "--variation", "cv", cv_dir,
            "--variation", "dent", dent_dir,
            "--filename-str", fstr,
            "--phase", "all",
        ]
        if args.max_files is not None:
            argv.extend(["--max-files", str(args.max_files)])
        summary_csv = path.join(args.out_base, "cache", f"dent_match_summary-{fstr}.csv")
        argv.extend(["--summary-csv", summary_csv])
        dent_match.main(argv)


def process_variation(
    variation: str,
    matched_dir: str,
    filename_str: str,
    *,
    input_format: str,
    var_defs: Dict[str, dict],
    max_files: Optional[int],
    stage_specs: Optional[Dict[str, List[Tuple[str, Any, str]]]] = None,
    collect_keyed: bool = False,
) -> Tuple[Dict[str, np.ndarray], SampleSummary, float, Dict[str, Dict[PairKey, float]]]:
    files = list_matched_files(matched_dir, filename_str)
    if not files:
        files = sorted(glob.glob(path.join(matched_dir, f"*{filename_str}*.df")))
        files = [f for f in files if "_matched" not in path.basename(f)]
    if max_files is not None:
        files = files[:max_files]

    print(f"[{variation}] {len(files)} files in {matched_dir}", flush=True)
    hists = {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
    summary = SampleSummary(variation=variation, n_matched_files=len(files))
    total_pot = 0.0
    keyed_maps: Dict[str, Dict[PairKey, float]] = {} if collect_keyed else {}

    for fpath in tqdm(files, desc=f"process {variation}"):
        if input_format == "sel_all":
            total_pot += process_sel_all_file(
                fpath,
                hists=hists,
                var_defs=var_defs,
                summary=summary,
                stage_specs=stage_specs or {},
                keyed_maps=keyed_maps if collect_keyed else None,
            )
        else:
            total_pot += process_sel_mup_file(
                fpath,
                hists=hists,
                var_defs=var_defs,
                summary=summary,
                keyed_maps=keyed_maps if collect_keyed else None,
            )

    summary.pot = total_pot
    return hists, summary, total_pot, keyed_maps


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def draw_cut_vlines(ax, var_name: str, *, with_arrows: bool = True) -> None:
    """Draw selection cut lines (from ``makedf.selections``) on a diagnostic axis."""
    vlines = DIAGNOSTIC_CUT_VLINES.get(var_name)
    if not vlines:
        return
    ymin, ymax = ax.get_ylim()
    if not np.isfinite(ymax) or not np.isfinite(ymin) or ymax <= ymin:
        return
    xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
    yspan = ymax - ymin
    for v in vlines:
        xcut = float(v[0])
        if with_arrows:
            # Match event-selection style: line to 75% height + keep-side arrow.
            line_ymin = max(0.0, ymin)
            ax.vlines(
                x=xcut, ymin=line_ymin, ymax=ymin + 0.75 * yspan,
                color="red", linestyle="--", zorder=50,
            )
            if len(v) < 2:
                continue
            direction = int(v[1])
            dx = 0.18 * xspan
            arrow_kw = dict(
                width=0.01 * yspan,
                color="red",
                head_width=0.04 * yspan,
                head_length=0.03 * xspan,
                length_includes_head=True,
                clip_on=True,
                zorder=60,
            )
            if direction == 0:
                ax.arrow(xcut, ymin + 0.4 * yspan, -dx, 0, **arrow_kw)
            elif direction == 1:
                ax.arrow(xcut, ymin + 0.4 * yspan, dx, 0, **arrow_kw)
        else:
            ax.axvline(xcut, color="red", linestyle="--", zorder=50)


def _draw_var_comparison_axes(
    ax_top,
    ax_bot,
    var_name: str,
    cfg: dict,
    all_hists: Dict[str, Dict[str, np.ndarray]],
    *,
    reference: str = "cv",
    pot_scales: Optional[Dict[str, float]] = None,
    title: Optional[str] = None,
    show_legend: bool = True,
    show_ylabel: bool = True,
    show_xlabel: bool = True,
    legend_fontsize: Optional[float] = None,
    title_fontsize: Optional[float] = None,
) -> None:
    """Draw CV/DENT overlay + ratio onto an existing (top, bottom) axis pair."""
    bins = np.asarray(cfg["bins"])
    centers = 0.5 * (bins[:-1] + bins[1:])
    ref = all_hists[reference][var_name].astype(float)
    if pot_scales:
        ref = ref * pot_scales.get(reference, 1.0)

    for var in VARIATIONS:
        counts = all_hists[var][var_name].astype(float)
        if pot_scales:
            counts = counts * pot_scales.get(var, 1.0)
        ax_top.step(
            centers, counts, where="mid",
            label=VARIATION_LABELS[var], color=VARIATION_COLORS[var],
        )
        if var != reference and ref.sum() > 0:
            ratio = np.divide(counts, ref, out=np.zeros_like(counts), where=ref > 0)
            ax_bot.step(centers, ratio, where="mid", color=RATIO_PANEL_COLOR)

    if show_ylabel:
        ax_top.set_ylabel("Events")
        ax_bot.set_ylabel("DENT / CV")
    if show_legend:
        ax_top.legend(fontsize=legend_fontsize)
    if title:
        ax_top.set_title(title, fontsize=title_fontsize)
    ax_top.set_xlim(float(bins[0]), float(bins[-1]))
    ax_top.set_ylim(bottom=0)
    draw_cut_vlines(ax_top, var_name, with_arrows=True)
    ax_bot.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax_bot.set_ylim(0.85, 1.15)
    ax_bot.set_xlim(float(bins[0]), float(bins[-1]))
    draw_cut_vlines(ax_bot, var_name, with_arrows=False)
    if show_xlabel:
        ax_bot.set_xlabel(cfg["label"])


def plot_var_comparison(
    var_name: str,
    cfg: dict,
    all_hists: Dict[str, Dict[str, np.ndarray]],
    *,
    fig_dir: str,
    reference: str = "cv",
    pot_scales: Optional[Dict[str, float]] = None,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": [3, 1]})
    ax_top, ax_bot = axes
    _draw_var_comparison_axes(
        ax_top, ax_bot, var_name, cfg, all_hists,
        reference=reference, pot_scales=pot_scales,
    )
    fig.tight_layout()
    save_fig(fig, f"diagnostic_{var_name}", fig_dir)
    plt.close(fig)


def plot_var_comparison_panel(
    var_name: str,
    cfg: dict,
    region_hists: Sequence[Tuple[str, Dict[str, Dict[str, np.ndarray]]]],
    *,
    fig_dir: str,
    fig_name: str,
    nrows: int,
    ncols: int,
    figsize: Tuple[float, float],
    reference: str = "cv",
    pot_scales: Optional[Dict[str, float]] = None,
) -> None:
    """Wide multi-panel diagnostic (each cell = distribution + ratio).

    ``region_hists`` is an ordered list of ``(panel_title, {cv|dent: {var: hist}})``.
    """
    n = len(region_hists)
    if n == 0:
        return
    if nrows * ncols < n:
        raise ValueError(f"grid {nrows}x{ncols} too small for {n} panels")

    fig = plt.figure(figsize=figsize)
    outer = fig.add_gridspec(nrows, ncols, hspace=0.45, wspace=0.28)
    for i, (title, all_hists) in enumerate(region_hists):
        if var_name not in all_hists.get(reference, {}):
            continue
        r, c = divmod(i, ncols)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax_top = fig.add_subplot(inner[0])
        ax_bot = fig.add_subplot(inner[1], sharex=ax_top)
        show_ylab = c == 0
        _draw_var_comparison_axes(
            ax_top, ax_bot, var_name, cfg, all_hists,
            reference=reference, pot_scales=pot_scales,
            title=title,
            show_legend=(i == 0),
            show_ylabel=show_ylab,
            show_xlabel=True,
            legend_fontsize=8,
            title_fontsize=10,
        )
        ax_top.tick_params(labelbottom=False)
        if not show_ylab:
            ax_top.tick_params(labelleft=True)
            ax_bot.tick_params(labelleft=True)
    # hide unused grid cells
    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        fig.add_subplot(outer[r, c]).axis("off")

    # Nested gridspecs are incompatible with tight_layout; use constrained spacing.
    fig.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.12, wspace=0.28, hspace=0.45)
    save_fig(fig, fig_name, fig_dir)
    plt.close(fig)


def plot_frac_diff(
    var_name: str,
    cfg: dict,
    frac_vals: np.ndarray,
    *,
    fig_dir: str,
    bins: np.ndarray = FRAC_DIFF_BINS,
) -> np.ndarray:
    """Histogram of per-event (DENT - CV) / CV for matched pairs."""
    hist, _ = np.histogram(frac_vals, bins=bins)
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(centers, hist, where="mid", color=VARIATION_COLORS["dent"])
    ax.axvline(0.0, color="gray", ls="--", lw=0.8)
    ax.set_xlim(float(bins[0]), float(bins[-1]))
    ax.set_xlabel(rf"$({VARIATION_LABELS['dent']} - {VARIATION_LABELS['cv']}) / {VARIATION_LABELS['cv']}$"
                  f"  [{cfg['label']}]")
    ax.set_ylabel("Matched pairs")
    if len(frac_vals) > 0:
        stats = (
            f"N = {len(frac_vals)}\n"
            f"mean = {np.mean(frac_vals):.3g}\n"
            f"std = {np.std(frac_vals):.3g}\n"
            f"median = {np.median(frac_vals):.3g}"
        )
        ax.text(
            0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="0.8"),
        )
    fig.tight_layout()
    save_fig(fig, f"fracdiff_{var_name}", fig_dir)
    plt.close(fig)
    return hist


def plot_ratio_dist(
    var_name: str,
    cfg: dict,
    ratio_vals: np.ndarray,
    *,
    fig_dir: str,
    bins: np.ndarray = RATIO_BINS,
) -> np.ndarray:
    """Histogram of per-event DENT / CV for matched pairs."""
    hist, _ = np.histogram(ratio_vals, bins=bins)
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(centers, hist, where="mid", color=VARIATION_COLORS["dent"])
    ax.axvline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_xlim(float(bins[0]), float(bins[-1]))
    ax.set_xlabel(rf"${VARIATION_LABELS['dent']} / {VARIATION_LABELS['cv']}$  [{cfg['label']}]")
    ax.set_ylabel("Matched pairs")
    if len(ratio_vals) > 0:
        stats = (
            f"N = {len(ratio_vals)}\n"
            f"mean = {np.mean(ratio_vals):.3g}\n"
            f"std = {np.std(ratio_vals):.3g}\n"
            f"median = {np.median(ratio_vals):.3g}"
        )
        ax.text(
            0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="0.8"),
        )
    fig.tight_layout()
    save_fig(fig, f"ratio_{var_name}", fig_dir)
    plt.close(fig)
    return hist


def _stats_dict(arr: np.ndarray) -> dict:
    if len(arr) == 0:
        return {"n": 0, "mean": float("nan"), "std": float("nan"), "median": float("nan")}
    return {
        "n": int(len(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "median": float(np.median(arr)),
    }


def make_paired_diff_plots(
    var_defs: Dict[str, dict],
    keyed_cv: Dict[str, Dict[PairKey, float]],
    keyed_dent: Dict[str, Dict[PairKey, float]],
    *,
    fig_dir_frac: str,
    fig_dir_ratio: str,
    cache_path: Optional[str] = None,
) -> Dict[str, dict]:
    """
    Per-event (DENT-CV)/CV and DENT/CV plots, and persist paired values for replot.

    Cache layout (replot without re-reading dfs)::

        paired[var] = {cv, dent, frac, ratio}   # aligned 1d arrays
        hists_frac / hists_ratio / stats / bins
    """
    makedirs(fig_dir_frac, exist_ok=True)
    makedirs(fig_dir_ratio, exist_ok=True)
    payload: Dict[str, Any] = {
        "frac_bins": FRAC_DIFF_BINS,
        "ratio_bins": RATIO_BINS,
        "cv_abs_eps": CV_ABS_EPS,
        "var_defs": {
            k: {"label": v.get("label", k), "bins": np.asarray(v["bins"])}
            for k, v in var_defs.items()
            if k != "integrated"
        },
        "paired": {},
        "vars": {},
    }
    for var_name, cfg in var_defs.items():
        if var_name == "integrated":
            continue
        cv_map = keyed_cv.get(var_name, {})
        dent_map = keyed_dent.get(var_name, {})
        cv_a, dent_a = compute_paired_values(cv_map, dent_map)
        fracs = (dent_a - cv_a) / cv_a if len(cv_a) else np.array([], dtype=float)
        ratios = compute_ratios_from_paired(cv_a, dent_a)
        n_cv, n_dent = len(cv_map), len(dent_map)
        n_common = len(set(cv_map) & set(dent_map)) if cv_map and dent_map else 0
        print(
            f"[paired] {var_name}: cv={n_cv} dent={n_dent} "
            f"common_keys={n_common} usable={len(cv_a)}",
            flush=True,
        )
        hist_f = plot_frac_diff(var_name, cfg, fracs, fig_dir=fig_dir_frac)
        hist_r = plot_ratio_dist(var_name, cfg, ratios, fig_dir=fig_dir_ratio)
        payload["paired"][var_name] = {
            "cv": cv_a,
            "dent": dent_a,
            "frac": fracs,
            "ratio": ratios,
        }
        payload["vars"][var_name] = {
            "label": cfg.get("label", var_name),
            "n_cv": n_cv,
            "n_dent": n_dent,
            "n_common_keys": n_common,
            "n_paired": int(len(cv_a)),
            "frac_stats": _stats_dict(fracs),
            "ratio_stats": _stats_dict(ratios),
            "hist_frac": hist_f,
            "hist_ratio": hist_r,
        }
    if cache_path is not None:
        save_hists(cache_path, payload)
    return payload


def make_frac_diff_plots(
    var_defs: Dict[str, dict],
    keyed_cv: Dict[str, Dict[PairKey, float]],
    keyed_dent: Dict[str, Dict[PairKey, float]],
    *,
    fig_dir: str,
    cache_path: Optional[str] = None,
) -> Dict[str, dict]:
    """Backward-compatible wrapper: frac + ratio under ``fig_dir``/{fracdiff,ratio}. """
    return make_paired_diff_plots(
        var_defs,
        keyed_cv,
        keyed_dent,
        fig_dir_frac=fig_dir if fig_dir.endswith("fracdiff") else path.join(fig_dir, "fracdiff"),
        fig_dir_ratio=path.join(path.dirname(fig_dir.rstrip("/")), "ratio")
        if fig_dir.rstrip("/").endswith("fracdiff")
        else path.join(fig_dir, "ratio"),
        cache_path=cache_path,
    )


def write_summary_table(
    summaries: Dict[str, SampleSummary],
    out_csv: str,
) -> pd.DataFrame:
    rows = []
    for var, s in summaries.items():
        rows.append(
            {
                "sample": VARIATION_LABELS.get(var.replace("_mup", ""), var),
                "stage_input": "sel_mup" if var.endswith("_mup") else "sel_all",
                "variation": var,
                "n_files": s.n_matched_files,
                "pot": s.pot,
                "efficiency_pct_evt": round(s.integrated_efficiency_pct(), 3),
                "efficiency_pct_mcnu": round(s.mcnu_efficiency_pct(), 3),
                "purity_pct": round(s.integrated_purity_pct(), 3),
                "n_signal_final": s.stages.get(s.final_stage_key, StageMetrics()).n_signal_raw,
                "n_evt_final": s.stages.get(s.final_stage_key, StageMetrics()).n_evt_raw,
            }
        )
    df = pd.DataFrame(rows)
    makedirs(path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote summary table → {out_csv}", flush=True)
    print(df.to_string(index=False), flush=True)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="DENT CV vs DENT comparison workflow")
    p.add_argument("--out-base", default=_OUT_BASE)
    p.add_argument("--cv-all-dir", default=DEFAULT_DIRS["cv_all"])
    p.add_argument("--dent-all-dir", default=DEFAULT_DIRS["dent_all"])
    p.add_argument("--cv-mup-dir", default=DEFAULT_DIRS["cv_mup"])
    p.add_argument("--dent-mup-dir", default=DEFAULT_DIRS["dent_mup"])
    p.add_argument("--max-files", type=int, default=None, help="Cap files per variation (None = all; set to small N for tests)")
    p.add_argument("--skip-match", action="store_true", help="Skip event matching step")
    p.add_argument("--skip-plots", action="store_true")
    p.add_argument(
        "--skip-fracdiff",
        action="store_true",
        help="Skip per-event (DENT-CV)/CV distribution plots (overlay plots still made)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = path.join(args.out_base, "cache")
    fig_dir = path.join(args.out_base, "plots")
    fig_dir_mup = path.join(fig_dir, "sel_mup")
    fig_dir_frac = path.join(fig_dir, "fracdiff")
    fig_dir_ratio = path.join(fig_dir, "ratio")
    fig_dir_frac_mup = path.join(fig_dir_mup, "fracdiff")
    fig_dir_ratio_mup = path.join(fig_dir_mup, "ratio")
    makedirs(cache_dir, exist_ok=True)
    makedirs(fig_dir, exist_ok=True)

    if not args.skip_match:
        run_matching(args)

    collect_keyed = (not args.skip_plots) and (not args.skip_fracdiff)

    # ── sel_all early-selection histograms + pipeline metrics ──────────────
    sel_all_var_defs = build_sel_all_var_defs()
    stage_specs = _stage_specs_by_key()
    all_hists_all: Dict[str, Dict[str, np.ndarray]] = {}
    summaries: Dict[str, SampleSummary] = {}
    pots: Dict[str, float] = {}
    keyed_all: Dict[str, Dict[str, Dict[PairKey, float]]] = {}

    for var, dir_path in (("cv", args.cv_all_dir), ("dent", args.dent_all_dir)):
        if not path.isdir(dir_path):
            print(f"[sel_all] skip {var}: dir missing {dir_path}", flush=True)
            continue
        hists, summary, pot, keyed = process_variation(
            var,
            dir_path,
            "sel_all",
            input_format="sel_all",
            var_defs=sel_all_var_defs,
            max_files=args.max_files,
            stage_specs=stage_specs,
            collect_keyed=collect_keyed,
        )
        all_hists_all[var] = hists
        summaries[var] = summary
        pots[var] = pot
        keyed_all[var] = keyed

    if all_hists_all:
        pot_scales = None
        if pots.get("cv") and pots.get("dent") and pots["dent"] > 0:
            scale = pots["cv"] / pots["dent"]
            pot_scales = {"cv": 1.0, "dent": scale}
            print(f"POT scale DENT→CV: {scale:.4f}", flush=True)
        cache_all = path.join(cache_dir, "dent_sel_all_hists.pkl")
        save_hists(
            cache_all,
            {
                "var_defs": {k: {"label": v["label"], "bins": v["bins"]} for k, v in sel_all_var_defs.items()},
                "hists": all_hists_all,
                "pots": pots,
                "pot_scales": pot_scales,
                "summaries": {k: v.__dict__ for k, v in summaries.items()},
            },
        )

        if not args.skip_plots and "cv" in all_hists_all and "dent" in all_hists_all:
            for var_name, cfg in sel_all_var_defs.items():
                if var_name not in all_hists_all["cv"] or var_name not in all_hists_all["dent"]:
                    continue
                plot_var_comparison(
                    var_name, cfg, all_hists_all,
                    fig_dir=fig_dir, pot_scales=pot_scales,
                )
            if collect_keyed and "cv" in keyed_all and "dent" in keyed_all:
                make_paired_diff_plots(
                    sel_all_var_defs,
                    keyed_all["cv"],
                    keyed_all["dent"],
                    fig_dir_frac=fig_dir_frac,
                    fig_dir_ratio=fig_dir_ratio,
                    cache_path=path.join(cache_dir, "dent_sel_all_paired.pkl"),
                )

    # ── sel_mup final-selected variables ───────────────────────────────────
    mup_var_defs = build_mup_var_defs()
    all_hists_mup: Dict[str, Dict[str, np.ndarray]] = {}
    keyed_mup: Dict[str, Dict[str, Dict[PairKey, float]]] = {}

    for var, dir_path in (("cv", args.cv_mup_dir), ("dent", args.dent_mup_dir)):
        if not path.isdir(dir_path):
            print(f"[sel_mup] skip {var}: dir missing {dir_path}", flush=True)
            continue
        hists, summary_mup, pot, keyed = process_variation(
            var,
            dir_path,
            "sel_mup",
            input_format="sel_mup",
            var_defs=mup_var_defs,
            max_files=args.max_files,
            collect_keyed=collect_keyed,
        )
        all_hists_mup[var] = hists
        summaries[f"{var}_mup"] = summary_mup
        pots[f"{var}_mup"] = pot
        keyed_mup[var] = keyed

    if all_hists_mup:
        cache_mup = path.join(cache_dir, "dent_sel_mup_hists.pkl")
        save_hists(cache_mup, {"var_defs": mup_var_defs, "hists": all_hists_mup, "pots": pots})

        if not args.skip_plots and "cv" in all_hists_mup and "dent" in all_hists_mup:
            # Matched unisim events are 1:1 paired, so no POT rescaling is
            # needed (and would be wrong if the two productions happen to have
            # different numbers of subruns per file).  Use the sel_all CV POT
            # as the exposure label since that is the true shared exposure.
            pot_scales_mup = None
            for var_name, cfg in mup_var_defs.items():
                if var_name not in all_hists_mup["cv"] or var_name not in all_hists_mup["dent"]:
                    continue
                plot_var_comparison(
                    var_name, cfg, all_hists_mup,
                    fig_dir=fig_dir_mup, pot_scales=pot_scales_mup,
                )
            if collect_keyed and "cv" in keyed_mup and "dent" in keyed_mup:
                make_paired_diff_plots(
                    mup_var_defs,
                    keyed_mup["cv"],
                    keyed_mup["dent"],
                    fig_dir_frac=fig_dir_frac_mup,
                    fig_dir_ratio=fig_dir_ratio_mup,
                    cache_path=path.join(cache_dir, "dent_sel_mup_paired.pkl"),
                )

    summary_csv = path.join(cache_dir, "dent_efficiency_purity_summary.csv")
    if summaries:
        write_summary_table(summaries, summary_csv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
