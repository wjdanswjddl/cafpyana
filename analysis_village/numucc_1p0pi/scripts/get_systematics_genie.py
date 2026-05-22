#!/usr/bin/env python
"""GENIE multisim systematics: event-rate and cross-section covariances (chunked only).

Input ``.df`` paths and per-group directory globs live in
``analysis_village.numucc_1p0pi.dataset_locations`` (``GENIE_GROUP_GLOBS`` /
:func:`iter_genie_chunk_map_tasks`). Each **knob group** (CCQE, MEC, …) uses its own sample
directory; within a group, all knobs share the same files.

**Map** (HDF splits inside each ``.df`` file):

* ``chunk-map`` — sequential reads of ``evt_{i}`` / ``mcnu_{i}``; accumulates *additive*
  histograms for the **rate** path and xsec tensors (see module docstring in git history
  for the merge math).

* ``chunk-merge`` — sums ``genie__<GROUP>__*.pkl`` from ``chunk-map`` and builds covariances
  per knob, plus a **group-combined** block (independent-sum of per-knob **fractional**
  covariances via :func:`syst_multisim_common.combine_indep_knob_cov_packs`, same recipe as
  neutrino Flux/G4 aggregate) under key :data:`GENIE_MERGE_COMBINED_KEY` in the merge dict / NPZ.

Example::

    python get_systematics_genie.py chunk-map --df-file in.df --out-dir chunks_genie \\
        --genie-group CCQE --max-splits 0

    python get_systematics_genie.py chunk-merge --chunks-dir chunks_genie \\
        --genie-group CCQE --out-dir plots_genie --xsec-unit 1.0

Batch inputs::

    python3 -c "from analysis_village.numucc_1p0pi.dataset_locations import iter_genie_chunk_map_tasks
    for g,p in iter_genie_chunk_map_tasks(): print(g,p)"
"""
from __future__ import annotations

import argparse
import gc
import glob
from collections import defaultdict
import logging
import os
import sys
import pickle
from os import path
from typing import AbstractSet, Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x, **kwargs):
        return x

import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

sys.path.append(path.join(path.dirname(__file__), "..", "..", ".."))

from pyanalib.covariance import get_covariance_matrix  # noqa: E402
from pyanalib.split_df_helpers import get_n_split  # noqa: E402

from analysis_village.numucc_1p0pi.categories import (  # noqa: E402
    get_genie_category,
    get_topo_category,
    topology_list,
)
from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_ORDER, GENIE_GROUP_KNOBS  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    genie_univ_weight_series,
    get_clipped_evts,
    get_response_matrix,
    get_univ_rates,
    plot_heatmap,
    plot_univ_hists,
    signal_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import (  # noqa: E402
    INTEGRATED_VAR_SAVE_NAME,
    VariableConfig,
)
from analysis_village.numucc_1p0pi.evt_derived_kinematics import (  # noqa: E402
    ensure_derived_trk_kinematics_cols,
    ensure_mc_level_phi_mcnu,
)
from analysis_village.numucc_1p0pi.syst_multisim_common import combine_indep_knob_cov_packs  # noqa: E402
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (  # noqa: E402
    CUT_STAGE_RATE_ONLY_SLUGS,
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    walk_pipeline,
)
from pyanalib.variable_calculator import (  # noqa: E402
    add_mc_cc1p0pi_tki_mcnu,
    add_reco_cc1p0pi_tki_evtdf,
    add_truth_cc1p0pi_tki_evtdf,
)

logger = logging.getLogger(__name__)

SystName = Tuple[str, str]


def genie_final_var_configs() -> List[VariableConfig]:
    """Variable set used for GENIE xsec tensors in the legacy ``final`` layout."""
    return with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))


def genie_all_var_configs(input_stage: str) -> List[VariableConfig]:
    """All ``VariableConfig`` objects produced by ``chunk-map`` for this input layout.

    * ``final``: same as :func:`genie_final_var_configs` (unchanged behaviour).
    * ``sel_all``: cut-stage observables (rate-only downstream) plus the same final
      variables as ``final`` (rate + xsec), deduped by ``var_save_name``.
    """
    if input_stage == "final":
        return genie_final_var_configs()
    seen: set[str] = set()
    out: List[VariableConfig] = []
    for spec in CUT_STAGE_VAR_SPECS:
        sn = spec.var_config.var_save_name
        if sn not in seen:
            seen.add(sn)
            out.append(spec.var_config)
    for vc in genie_final_var_configs():
        if vc.var_save_name not in seen:
            seen.add(vc.var_save_name)
            out.append(vc)
    return out


def _align_evt_mcnu(evt_df: pd.DataFrame, mcnu_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Inner-join ``evt`` and ``mcnu`` on the index (same convention as split HDF producers)."""
    if evt_df is None or len(evt_df) == 0:
        return evt_df, mcnu_df.iloc[0:0]
    ix = evt_df.index.intersection(mcnu_df.index)
    return evt_df.loc[ix], mcnu_df.loc[ix]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_genie_dataframes(
    dfs: Mapping[str, Any],
    *,
    context: str = "",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(evt, mcnu)`` after verifying both frames exist."""
    missing = [k for k in ("evt", "mcnu") if k not in dfs]
    if missing:
        raise ValueError(
            f"{context} GENIE setup requires both 'evt' and 'mcnu' dataframes; "
            f"missing keys: {missing}. Present: {sorted(dfs.keys())}"
        )
    evt = dfs["evt"]
    mcnu = dfs["mcnu"] #.mc.copy()
    if not isinstance(evt, pd.DataFrame) or not isinstance(mcnu, pd.DataFrame):
        raise TypeError(f"{context} 'evt' and 'mcnu' must be pandas DataFrames.")
    if len(evt) == 0 and len(mcnu) == 0:
        logger.warning("%s empty evt and mcnu frames", context)
    return evt, mcnu


def validate_split_pair(evt_df: pd.DataFrame, mcnu_df: pd.DataFrame, split_idx: int) -> None:
    if len(evt_df) == 0 and len(mcnu_df) == 0:
        logger.debug("split %d: both evt and mcnu empty", split_idx)


def _attach_phi_degrees(mc_evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame) -> None:
    """Fill mu/p track phi in degrees (same recipe as the legacy monolithic driver).

    Event frames use top-level ``mu`` / ``p`` under ``pfp.trk``; ``mcnu`` frames (after
    :func:`_prefix_mcnu_columns`) nest the same under ``mc`` → ``mc.mu`` / ``mc.p``.
    """
    def _fill_trk_phi(df: pd.DataFrame, *, mcnu: bool) -> None:
        if df is None or len(df) == 0:
            return
        for pref in ("mu", "p"):
            if mcnu:
                head: Tuple[str, ...] = ("mc", pref)
                phi_col = head + ("phi", "")
                dir_x = head + ("dir", "x")
                dir_y = head + ("dir", "y")
                df.loc[:, phi_col] = np.degrees(np.arctan2(df.loc[:, dir_x], df.loc[:, dir_y]))

            else:
                head: Tuple[str, ...] = (pref,)
                phi_col = head + ("pfp", "trk", "phi", "", "", "")
                dir_x = head + ("pfp", "trk", "dir", "x", "", "")
                dir_y = head + ("pfp", "trk", "dir", "y", "", "")
                df.loc[:, phi_col] = np.degrees(np.arctan2(df.loc[:, dir_x], df.loc[:, dir_y]))

    _fill_trk_phi(mc_evt_df, mcnu=False)
    _fill_trk_phi(mc_nu_df, mcnu=True)


# ---------------------------------------------------------------------------
# Universe helpers
# ---------------------------------------------------------------------------


def _block_column_weight_leaf(col) -> str:
    """Leaf name for GENIE weight columns (univ_*, morph, ps*, ms*) under a knob block.

    After ``multicol_concat``, tuples are often padded with trailing ``''``; the meaningful
    leaf is the **first non-empty** segment (same idea as :mod:`makedf.getsyst`).
    """
    if not isinstance(col, tuple):
        return str(col)
    for x in col:
        if x != "" and x is not None:
            return str(x)
    return str(col[0])


def _iter_leaf_strings(block_cols: pd.Index) -> List[str]:
    return [_block_column_weight_leaf(c) for c in block_cols]


def _infer_multisim_n_univ(block_cols: pd.Index) -> int:
    max_i = -1
    for leaf_s in _iter_leaf_strings(block_cols):
        if not leaf_s.startswith("univ_"):
            continue
        try:
            max_i = max(max_i, int(leaf_s.split("_", 1)[1]))
        except ValueError:
            continue
    return max_i + 1


def _replace_first_nonempty_segment(col: Tuple, value: str) -> Tuple:
    parts = list(col)
    for i, x in enumerate(parts):
        if x != "" and x is not None:
            parts[i] = value
            return tuple(parts)
    parts[0] = value
    return tuple(parts)


def _ensure_univ0_from_leaf(df: pd.DataFrame, syst_key: Tuple[str, ...], src_leaf: str) -> None:
    """
    For non-multisim knobs, alias one leaf (e.g. 'ps1' or 'morph') into a synthetic 'univ_0'
    column under the same systematic block so downstream code can treat it like multisim.
    """
    block = df.loc[:, syst_key]
    src_rest = None
    for c in block.columns:
        if _block_column_weight_leaf(c) == src_leaf:
            src_rest = c if isinstance(c, tuple) else (c,)
            break
    if src_rest is None:
        # Nothing to do; caller will raise a clearer error.
        return

    dst_rest = _replace_first_nonempty_segment(tuple(src_rest), "univ_0")
    full_src = tuple(syst_key) + tuple(src_rest)
    full_dst = tuple(syst_key) + tuple(dst_rest)

    if full_dst in df.columns:
        return
    df.loc[:, full_dst] = df.loc[:, full_src]


def normalize_and_infer_n_univ(mc_evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame, syst_name: SystName) -> int:
    """
    Detect whether a knob is multisim vs multisigma vs unisim and normalize to a common
    interface: multisim keeps 'univ_*', while multisigma/unisim get a synthetic 'univ_0'.

    Rules (per user request):
    - multisigma: use 'ps1' as the one-universe unisim
    - unisim: use 'morph' as the one-universe unisim
    """
    key = tuple(syst_name)
    block_cols = mc_evt_df.loc[:, key].columns
    n = _infer_multisim_n_univ(block_cols)
    if n > 0:
        return n

    leaves = set(_iter_leaf_strings(block_cols))
    if "ps1" in leaves:
        _ensure_univ0_from_leaf(mc_evt_df, key, "ps1")
        _ensure_univ0_from_leaf(mc_nu_df, key, "ps1")
        return 1
    if "morph" in leaves:
        _ensure_univ0_from_leaf(mc_evt_df, key, "morph")
        _ensure_univ0_from_leaf(mc_nu_df, key, "morph")
        return 1

    raise ValueError(
        f"No univ_* columns under syst_name={syst_name!r}, and also no 'ps1' (multisigma) "
        f"or 'morph' (unisim) leaf found. Available leaves: {sorted(leaves)}"
    )


def _syst_plot_key(syst_name: SystName) -> SystName:
    return syst_name


def copy_matrix_pack(pack: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Deep copy of a ``{cov, cov_frac, corr}`` pack."""
    return {
        k: np.asarray(pack[k], dtype=np.float64).copy() for k in ("cov", "cov_frac", "corr")
    }


def xsec_pack_from_accumulators(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    n_univ: int,
    *,
    xsec_unit: float,
    bkgd_subtract: bool = True,
    plot: bool = False,
    save_fig: bool = False,
    save_fig_dir: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """Cross-section covariance for one knob (response-matrix path).

    For the integrated single-bin variable, holds the signal rate fixed at CV and varies
    efficiency, smearing, and background subtraction — do **not** copy the rate matrix.
    """
    nb = len(var_config.bins) - 1
    acc = _empty_xsec_tensor_acc(n_univ, nb)
    accumulate_xsec_path_chunk(mc_evt_df, mc_nu_df, var_config, syst_name, n_univ, acc)
    univ_xsec = finalize_xsec_univ_events(acc, xsec_unit=xsec_unit)
    cv_xsec = finalize_cv_sel_reco_xsec(
        acc, xsec_unit=xsec_unit, bkgd_subtract=bkgd_subtract
    )
    return covariance_bundle_univ_events(
        univ_xsec,
        cv_xsec,
        syst_name,
        var_config,
        "xsec",
        plot=plot,
        save_fig=save_fig,
        save_fig_dir=save_fig_dir,
    )


def sanitize_matrix_pack(
    ret: Mapping[str, np.ndarray], *, context: str = ""
) -> Dict[str, np.ndarray]:
    """Replace non-finite matrix elements with zero; log when that happens."""
    out: Dict[str, np.ndarray] = {}
    for k in ("cov", "cov_frac", "corr"):
        arr = np.asarray(ret[k], dtype=float)
        if np.any(~np.isfinite(arr)):
            logger.warning(
                "%s: non-finite values in %s matrix — zeroed (check empty MC bins / eff division)",
                context or "sanitize_matrix_pack",
                k,
            )
        out[k] = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return out


# ---------------------------------------------------------------------------
# XSEC-path accumulation (chunk-additive)
# ---------------------------------------------------------------------------


def _empty_xsec_tensor_acc(n_univ: int, nb: int) -> Dict[str, np.ndarray]:
    return {
        "nevts_allmc": np.zeros(nb, dtype=np.float64),
        "cv_sel_reco": np.zeros(nb, dtype=np.float64),
        "cv_allsel_reco": np.zeros(nb, dtype=np.float64),
        "reco_vs_true": np.zeros((n_univ, nb, nb), dtype=np.float64),
        "signal_allmc": np.zeros((n_univ, nb), dtype=np.float64),
        "signal_sel_truth": np.zeros((n_univ, nb), dtype=np.float64),
        "bg_cv": np.zeros(nb, dtype=np.float64),
        "bg_univ": np.zeros((n_univ, nb), dtype=np.float64),
    }


def accumulate_xsec_path_chunk(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    n_univ: int,
    acc: MutableMapping[str, np.ndarray],
) -> None:
    """Add contributions from one evt/mcnu chunk into ``acc`` (in-place)."""
    bins = var_config.bins
    nb = len(bins) - 1

    evtdf_signal = mc_evt_df[mc_evt_df.topo_categ == 1]
    nudf_signal = mc_nu_df[mc_nu_df.topo_categ == 1]
    evtdf_div_topo = [mc_evt_df[mc_evt_df.topo_categ == mode] for mode in topology_list]

    ret = signal_hists(mc_evt_df, mc_nu_df, var_config, return_data=True, plot=False)
    nevts_allmc = ret["nevts_allmc"]
    if nevts_allmc is None:
        raise ValueError(
            f"xsec path requires mcnu for var={var_config.var_save_name}; nevts_allmc is None"
        )

    acc["nevts_allmc"] += np.asarray(nevts_allmc, dtype=np.float64)
    acc["cv_sel_reco"] += np.asarray(ret["nevts_sel_reco"], dtype=np.float64)
    acc["cv_allsel_reco"] += np.asarray(ret["nevts_allsel_reco"], dtype=np.float64)

    wblock_evt = evtdf_signal[syst_name]
    wblock_nu = nudf_signal[syst_name]
    for uidx in range(n_univ):
        w_evt_univ = genie_univ_weight_series(wblock_evt, uidx)
        w_nu_univ = genie_univ_weight_series(wblock_nu, uidx)
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true, _, _ = np.histogram2d(
                ret["var_sel_truth"],
                ret["var_sel_reco"],
                weights=ret["wgt_sel_truth"] * w_evt_univ,
                bins=bins,
            )
        acc["reco_vs_true"][uidx] += reco_vs_true

        sam, _ = np.histogram(
            ret["var_allmc"],
            weights=ret["wgt_allmc"] * w_nu_univ,
            bins=bins,
        )
        sst, _ = np.histogram(
            ret["var_sel_truth"],
            weights=ret["wgt_sel_truth"] * w_evt_univ,
            bins=bins,
        )
        acc["signal_allmc"][uidx] += sam
        acc["signal_sel_truth"][uidx] += sst

    for this_evtdf in evtdf_div_topo[1:]:
        var, wgt = get_clipped_evts(
            this_evtdf,
            var_config.var_evt_reco_col,
            bins,
            var_save_name=var_config.var_save_name,
        )
        acc["bg_cv"] += np.histogram(var, bins=bins, weights=wgt)[0].astype(np.float64)
        wblock_bg = this_evtdf[syst_name]
        for uidx in range(n_univ):
            uw = genie_univ_weight_series(wblock_bg, uidx).copy()
            uw[np.isnan(uw)] = 1.0
            acc["bg_univ"][uidx] += np.histogram(var, bins=bins, weights=wgt * uw)[0].astype(
                np.float64
            )


def finalize_xsec_univ_events(
    acc: Mapping[str, np.ndarray],
    xsec_unit: float,
) -> np.ndarray:
    """Return ``univ_events`` of shape ``(n_univ, nbins)`` for covariance (xsec path)."""
    nevts_allmc = np.asarray(acc["nevts_allmc"], dtype=np.float64)
    nb = int(nevts_allmc.shape[0])
    n_univ = int(acc["reco_vs_true"].shape[0])
    scale = float(xsec_unit)

    rows: List[np.ndarray] = []
    for uidx in range(n_univ):
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true = acc["reco_vs_true"][uidx]

        # Match utils.get_univ_rates xsec branch: ratio of summed weighted histograms.
        eff = acc["signal_sel_truth"][uidx] / acc["signal_allmc"][uidx]
        response = get_response_matrix(reco_vs_true, eff)
        signal_univ = response @ nevts_allmc
        signal_univ = signal_univ + (acc["bg_univ"][uidx] - acc["bg_cv"])
        signal_univ *= scale
        rows.append(signal_univ)
    return np.asarray(rows, dtype=np.float64)


def finalize_cv_sel_reco_xsec(
    acc: Mapping[str, np.ndarray], xsec_unit: float, *, bkgd_subtract: bool = True
) -> np.ndarray:
    scale = float(xsec_unit)
    base = acc["cv_sel_reco"] if bkgd_subtract else acc["cv_allsel_reco"]
    return np.asarray(base, dtype=np.float64) * scale


XSEC_COMPONENTS = ("full", "efficiency", "smearing", "background", "signal")


def _xsec_cv_tensors(acc: Mapping[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """CV efficiency, truth→reco migration, and background yield (integrated: smear is 1×1)."""
    nevts_allmc = np.asarray(acc["nevts_allmc"], dtype=np.float64)
    nb = int(nevts_allmc.shape[0])
    eff_cv = np.asarray(acc["signal_sel_truth"][0], dtype=np.float64) / np.asarray(
        acc["signal_allmc"][0], dtype=np.float64
    )
    if nb == 1:
        reco_cv = np.array([[1.0]], dtype=np.float64)
    else:
        reco_cv = np.asarray(acc["reco_vs_true"][0], dtype=np.float64)
    bg_cv = np.asarray(acc["bg_cv"], dtype=np.float64)
    return eff_cv, reco_cv, bg_cv


def xsec_component_univ_events(
    acc: Mapping[str, np.ndarray],
    xsec_unit: float,
    component: str,
) -> np.ndarray:
    """Universe xsec vector varying only one part of the response-matrix path.

    * **efficiency** — ``eff(u)`` with CV smearing and no background term
    * **smearing** — ``R(reco(u), eff_CV)`` with no background term
    * **background** — CV signal ``R(reco_CV, eff_CV) @ N_gen`` plus ``bg(u)-bg_CV``
    * **signal** — efficiency + smearing, no background
    * **full** — same as :func:`finalize_xsec_univ_events`
    """
    if component not in XSEC_COMPONENTS:
        raise ValueError(f"component must be one of {XSEC_COMPONENTS}, got {component!r}")

    if component == "full":
        return finalize_xsec_univ_events(acc, xsec_unit=xsec_unit)

    nevts_allmc = np.asarray(acc["nevts_allmc"], dtype=np.float64)
    nb = int(nevts_allmc.shape[0])
    n_univ = int(acc["reco_vs_true"].shape[0])
    scale = float(xsec_unit)
    eff_cv, reco_cv, _bg_cv = _xsec_cv_tensors(acc)

    rows: List[np.ndarray] = []
    for uidx in range(n_univ):
        eff_u = np.asarray(acc["signal_sel_truth"][uidx], dtype=np.float64) / np.asarray(
            acc["signal_allmc"][uidx], dtype=np.float64
        )
        if nb == 1:
            reco_u = np.array([[1.0]], dtype=np.float64)
        else:
            reco_u = np.asarray(acc["reco_vs_true"][uidx], dtype=np.float64)
        bg_delta = np.asarray(acc["bg_univ"][uidx], dtype=np.float64) - np.asarray(
            acc["bg_cv"], dtype=np.float64
        )

        if component == "efficiency":
            response = get_response_matrix(reco_cv, eff_u)
            signal = response @ nevts_allmc
        elif component == "smearing":
            response = get_response_matrix(reco_u, eff_cv)
            signal = response @ nevts_allmc
        elif component == "background":
            response = get_response_matrix(reco_cv, eff_cv)
            signal = response @ nevts_allmc + bg_delta
        elif component == "signal":
            response = get_response_matrix(reco_u, eff_u)
            signal = response @ nevts_allmc
        else:
            raise AssertionError(component)

        rows.append(np.asarray(signal * scale, dtype=np.float64))
    return np.asarray(rows, dtype=np.float64)


def xsec_path_component_diagnostics(
    acc: Mapping[str, np.ndarray],
    var_config: VariableConfig,
    syst_name: SystName,
    *,
    xsec_unit: float = 1.0,
    bkgd_subtract: bool = True,
) -> Dict[str, Any]:
    """Per-component fractional covariance for one knob (integrated xsec decomposition)."""
    cv_events = finalize_cv_sel_reco_xsec(
        acc, xsec_unit=xsec_unit, bkgd_subtract=bkgd_subtract
    )
    out: Dict[str, Any] = {
        "knob": syst_name[1],
        "var": var_config.var_save_name,
        "n_univ": int(acc["reco_vs_true"].shape[0]),
        "n_bins": len(var_config.bin_centers),
        "cv_xsec": np.asarray(cv_events, dtype=np.float64),
        "nevts_allmc": np.asarray(acc["nevts_allmc"], dtype=np.float64).copy(),
    }
    eff_cv, reco_cv, bg_cv = _xsec_cv_tensors(acc)
    out["eff_cv"] = eff_cv
    out["reco_cv_sum"] = float(np.sum(reco_cv))
    out["bg_cv"] = bg_cv

    for comp in XSEC_COMPONENTS:
        univ = xsec_component_univ_events(acc, xsec_unit, comp)
        pack = sanitize_matrix_pack(
            get_covariance_matrix(univ, cv_events),
            context=f"{var_config.var_save_name}/{comp}/{syst_name[1]}",
        )
        frac_var = float(np.asarray(pack["cov_frac"]).flat[0])
        rel_pulls = (univ[:, 0] - cv_events[0]) / cv_events[0] if cv_events[0] else univ[:, 0] * 0.0
        out[comp] = {
            "univ_events": univ,
            "cov_frac": pack["cov_frac"],
            "frac_variance": frac_var,
            "unc_pct": 100.0 * np.sqrt(max(frac_var, 0.0)),
            "rel_pull_mean": float(np.mean(rel_pulls)),
            "rel_pull_rms": float(np.std(rel_pulls)),
            "rel_pull_min": float(np.min(rel_pulls)),
            "rel_pull_max": float(np.max(rel_pulls)),
        }
    return out


def combine_component_cov_fracs(
    diagnostics: Sequence[Mapping[str, Any]],
    component: str,
) -> np.ndarray:
    """Sum per-knob ``cov_frac`` for one component (independent-knob recipe)."""
    total: Optional[np.ndarray] = None
    for diag in diagnostics:
        block = np.asarray(diag[component]["cov_frac"], dtype=np.float64)
        total = block.copy() if total is None else total + block
    if total is None:
        raise ValueError("no diagnostics to combine")
    return total


def default_fsi_knob_other() -> str:
    """First pion FSI knob in the ``Other`` GENIE group (SBN v1 reweights)."""
    from makedf.geniesyst import other_genie_systematics

    for knob in other_genie_systematics:
        if knob.endswith("_pi") and "MFP" in knob:
            return knob
    for knob in other_genie_systematics:
        if knob.endswith("_pi"):
            return knob
    return other_genie_systematics[0]


def integrated_xsec_knob_probe_data(
    acc: Mapping[str, np.ndarray],
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    *,
    xsec_unit: float = 1.0,
    bkgd_subtract: bool = True,
) -> Dict[str, Any]:
    """Per-universe arrays for integrated xsec knob investigation plots.

    Returns CV + multisim pulls for truth-generated signal rate, background-subtracted
    selected rate (rate path), background, signal-only xsec path, efficiency, full xsec,
    and truth→reco smearing / response matrices.
    """
    n_univ = int(acc["reco_vs_true"].shape[0])
    nb = int(np.asarray(acc["nevts_allmc"]).shape[0])
    scale = float(xsec_unit)

    univ_sel_rate, cv_sel_rate = get_univ_rates(
        cov_type="rate",
        syst_type="GENIE",
        evtdf=mc_evt_df,
        nudf=mc_nu_df,
        var_config=var_config,
        syst_name=syst_name,
        n_univ=n_univ,
        bkgd_subtract=bkgd_subtract,
    )
    univ_sel_rate = np.asarray(univ_sel_rate, dtype=np.float64)
    cv_sel_rate = np.asarray(cv_sel_rate, dtype=np.float64)

    truth_gen = np.asarray(acc["signal_allmc"], dtype=np.float64)
    eff_cv, reco_cv, bg_cv = _xsec_cv_tensors(acc)
    eff_univ = np.asarray(acc["signal_sel_truth"], dtype=np.float64) / np.asarray(
        acc["signal_allmc"], dtype=np.float64
    )
    bg_univ = np.asarray(acc["bg_univ"], dtype=np.float64)

    signal_xsec = xsec_component_univ_events(acc, xsec_unit, "signal")
    full_xsec = finalize_xsec_univ_events(acc, xsec_unit=xsec_unit)
    cv_xsec = finalize_cv_sel_reco_xsec(acc, xsec_unit, bkgd_subtract=bkgd_subtract)

    reco_stack = np.asarray(acc["reco_vs_true"], dtype=np.float64)
    response_cv = get_response_matrix(reco_cv, eff_cv)
    response_univ = np.stack(
        [
            get_response_matrix(
                reco_stack[u] if nb > 1 else np.array([[1.0]]),
                eff_univ[u],
            )
            for u in range(n_univ)
        ],
        axis=0,
    )

    return {
        "knob": syst_name[1],
        "var": var_config.var_save_name,
        "n_univ": n_univ,
        "n_bins": nb,
        "nevts_allmc": np.asarray(acc["nevts_allmc"], dtype=np.float64).copy(),
        "truth_gen_rate": {"cv": truth_gen[0].copy(), "univ": truth_gen.copy()},
        "sel_rate": {"cv": cv_sel_rate.copy(), "univ": univ_sel_rate.copy()},
        "background": {"cv": bg_cv.copy(), "univ": bg_univ.copy()},
        "signal_xsec": {
            "cv": signal_xsec[0].copy(),
            "univ": signal_xsec.copy(),
        },
        "efficiency": {"cv": eff_cv.copy(), "univ": eff_univ.copy()},
        "integrated_xsec": {"cv": cv_xsec.copy(), "univ": full_xsec.copy()},
        "smearing": {
            "reco_vs_true_cv": reco_cv.copy(),
            "reco_vs_true_univ": reco_stack.copy(),
            "response_cv": response_cv.copy(),
            "response_univ": response_univ.copy(),
        },
    }


def print_xsec_component_table(
    diagnostics: Sequence[Mapping[str, Any]],
    *,
    title: str = "",
) -> None:
    """Print per-knob and summed component uncertainties [%]."""
    if title:
        print(title)
    comps = [c for c in XSEC_COMPONENTS if c != "full"]
    hdr = f"{'knob':40s}  " + "  ".join(f"{c:12s}" for c in ["full"] + comps)
    print(hdr)
    print("-" * len(hdr))
    for diag in diagnostics:
        row = f"{diag['knob'][:40]:40s}  "
        row += "  ".join(f"{diag[c]['unc_pct']:12.4f}" for c in ["full"] + comps)
        print(row)
    print("-" * len(hdr))
    row = f"{'SUM(indep)':40s}  "
    for c in ["full"] + comps:
        cf = combine_component_cov_fracs(diagnostics, c)
        row += f"{100.0 * np.sqrt(max(float(cf.flat[0]), 0.0)):12.4f}"
    print(row)


# ---------------------------------------------------------------------------
# Pickle blob layout for chunked GENIE
# ---------------------------------------------------------------------------

RATE_ACC_KEY = "rate_univ_cv"
XSEC_ACC_KEY = "xsec_accumulators"

# Top-level key in ``chunk-merge`` output dict / NPZ ``syst`` object: ``slug -> {rate, xsec?}``
# with the same matrix packs as per-knob entries, built by summing knob covariances as independent.
GENIE_MERGE_COMBINED_KEY = "__GENIE_group_combined__"


def accumulate_chunk_into_blob_root(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    blob_root: MutableMapping[str, Any],
    var_configs: Sequence[VariableConfig],
    syst_names: Sequence[SystName],
    *,
    bkgd_subtract: bool = True,
    skip_xsec_slugs: Optional[AbstractSet[str]] = None,
) -> None:
    validate_split_pair(mc_evt_df, mc_nu_df, -1)
    rate_blk = blob_root.setdefault(RATE_ACC_KEY, {})

    for syst_name in syst_names:
        knob = syst_name[1]
        n_univ = normalize_and_infer_n_univ(mc_evt_df, mc_nu_df, syst_name)

        for var_config in var_configs:
            slug = var_config.var_save_name

            # ---- rate path (identical to summing get_univ_rates per chunk)
            univ_r, cv_r = get_univ_rates(
                cov_type="rate",
                syst_type="GENIE",
                evtdf=mc_evt_df,
                nudf=mc_nu_df,
                var_config=var_config,
                syst_name=syst_name,
                n_univ=n_univ,
                bkgd_subtract=bkgd_subtract,
                plot=False,
            )
            slot_r = rate_blk.setdefault(knob, {}).setdefault(
                slug,
                {"univ": np.zeros_like(univ_r, dtype=np.float64), "cv": np.zeros_like(cv_r, dtype=np.float64)},
            )
            slot_r["univ"] += np.asarray(univ_r, dtype=np.float64)
            slot_r["cv"] += np.asarray(cv_r, dtype=np.float64)

            # ---- xsec tensors (response-matrix recipe; includes integrated single-bin)
            if skip_xsec_slugs is not None and slug in skip_xsec_slugs:
                continue
            xsec_blk = blob_root.setdefault(XSEC_ACC_KEY, {})
            knob_blk = xsec_blk.setdefault(knob, {})
            if slug not in knob_blk:
                nb = len(var_config.bins) - 1
                knob_blk[slug] = _empty_xsec_tensor_acc(n_univ, nb)
            accumulate_xsec_path_chunk(mc_evt_df, mc_nu_df, var_config, syst_name, n_univ, knob_blk[slug])


def merge_genie_chunk_pickles(paths: List[str]) -> Dict[str, Any]:
    merged: Optional[Dict[str, Any]] = None
    for fp in tqdm(paths, desc="merge GENIE chunks"):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        if merged is None:
            merged = {
                "kind": "genie_syst_merged",
                RATE_ACC_KEY: d[RATE_ACC_KEY],
                XSEC_ACC_KEY: d.get(XSEC_ACC_KEY, {}),
                "input_stage": d.get("input_stage", "final"),
                "meta": [d.get("meta", {})],
            }
            continue
        merged["meta"].append(d.get("meta", {}))
        if d.get("input_stage", "final") != merged.get("input_stage", "final"):
            logger.warning(
                "GENIE chunk input_stage mismatch: merged=%r vs chunk=%r (file=%s)",
                merged.get("input_stage"),
                d.get("input_stage"),
                fp,
            )

        # rate
        for knob, vars_d in d[RATE_ACC_KEY].items():
            for slug, pack in vars_d.items():
                mp = merged[RATE_ACC_KEY].setdefault(knob, {}).setdefault(
                    slug,
                    {"univ": pack["univ"].copy(), "cv": pack["cv"].copy()},
                )
                if mp["univ"].shape != pack["univ"].shape:
                    raise ValueError(f"rate shape mismatch {knob} {slug}")
                mp["univ"] += pack["univ"]
                mp["cv"] += pack["cv"]

        # xsec tensors
        for knob, vars_d in d.get(XSEC_ACC_KEY, {}).items():
            for slug, acc_new in vars_d.items():
                mp_acc = merged[XSEC_ACC_KEY].setdefault(knob, {}).setdefault(slug, {})
                if not mp_acc:
                    for k, arr in acc_new.items():
                        mp_acc[k] = np.array(arr, dtype=np.float64, copy=True)
                    mp_acc.setdefault(
                        "cv_allsel_reco", np.zeros_like(mp_acc["cv_sel_reco"], dtype=np.float64)
                    )
                else:
                    for k in acc_new:
                        mp_acc[k] += np.asarray(acc_new[k], dtype=np.float64)
    if merged is None:
        raise RuntimeError("no GENIE chunk pickles merged")
    return merged


# ---------------------------------------------------------------------------
# Plotting / covariance packaging (matches legacy get_systematics)
# ---------------------------------------------------------------------------


def covariance_bundle_univ_events(
    univ_events: np.ndarray,
    cv_events: np.ndarray,
    syst_name: SystName,
    var_config: VariableConfig,
    cov_tag: str,
    *,
    plot: bool = False,
    save_fig: bool = False,
    save_fig_dir: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    ret = sanitize_matrix_pack(
        get_covariance_matrix(univ_events, cv_events),
        context="%s/%s" % (var_config.var_save_name, cov_tag),
    )
    if save_fig and save_fig_dir:
        os.makedirs(save_fig_dir, exist_ok=True)
        sk = _syst_plot_key(syst_name)
        plot_univ_hists(
            univ_events,
            cv_events,
            sk,
            var_config,
            plot=plot,
            save_fig=True,
            save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-{syst_name[1]}_{cov_tag}-univ_hists"),
        )
        for matrix_type in ("cov", "cov_frac", "corr"):
            lab = matrix_type.replace("_", " ").title()
            plot_heatmap(
                ret[matrix_type],
                var_config.bins,
                plot_labels=[var_config.var_labels[1], var_config.var_labels[1], lab],
                plot=plot,
                save_fig=True,
                save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-{syst_name[1]}_{cov_tag}-{matrix_type}"),
            )
    return ret


def get_systematics(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    syst_type: str = "GENIE",
    plot: bool = False,
    save_fig: bool = False,
    save_fig_dir: Optional[str] = None,
    *,
    xsec_unit: float = 0.0,
):
    """Legacy helper — rate via :func:`utils.get_univ_rates`; xsec via response accumulators."""
    matrices: Dict[str, Any] = {}
    validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context="get_systematics: ")
    n_univ = normalize_and_infer_n_univ(mc_evt_df, mc_nu_df, syst_name)

    univ_rate, cv_rate = get_univ_rates(
        cov_type="rate",
        syst_type=syst_type,
        evtdf=mc_evt_df,
        nudf=mc_nu_df,
        var_config=var_config,
        syst_name=syst_name,
        n_univ=n_univ,
        xsec_unit=xsec_unit,
        plot=False,
    )
    matrices["rate"] = covariance_bundle_univ_events(
        univ_rate,
        cv_rate,
        syst_name,
        var_config,
        "rate",
        plot=plot,
        save_fig=save_fig,
        save_fig_dir=save_fig_dir,
    )

    if var_config.var_save_name == INTEGRATED_VAR_SAVE_NAME:
        matrices["xsec"] = xsec_pack_from_accumulators(
            mc_evt_df,
            mc_nu_df,
            var_config,
            syst_name,
            n_univ,
            xsec_unit=xsec_unit,
            plot=plot,
            save_fig=save_fig,
            save_fig_dir=save_fig_dir,
        )
    else:
        univ_xsec, cv_xsec = get_univ_rates(
            cov_type="xsec",
            syst_type=syst_type,
            evtdf=mc_evt_df,
            nudf=mc_nu_df,
            var_config=var_config,
            syst_name=syst_name,
            n_univ=n_univ,
            xsec_unit=xsec_unit,
            plot=False,
        )
        matrices["xsec"] = covariance_bundle_univ_events(
            univ_xsec,
            cv_xsec,
            syst_name,
            var_config,
            "xsec",
            plot=plot,
            save_fig=save_fig,
            save_fig_dir=save_fig_dir,
        )
    return matrices


# ---------------------------------------------------------------------------
# CLI chunk-map / chunk-merge
# ---------------------------------------------------------------------------


def _prefix_mcnu_columns(mc_nu_df: pd.DataFrame) -> None:
    """Ensure ``mcnu`` columns have a leading ``mc`` level (same as legacy chunk-map)."""
    if isinstance(mc_nu_df.columns, pd.MultiIndex):
        try:
            first_level = mc_nu_df.columns.get_level_values(0)
            need_prefix = not np.all(first_level == "mc")
        except Exception:
            need_prefix = True
        if need_prefix:
            mc_nu_df.columns = pd.MultiIndex.from_tuples(
                [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
            )


def _annotate_topo_genie_phi(evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame) -> None:
    # evt_df.loc[:, "topo_categ"] = get_topo_category(evt_df)
    # mc_nu_df.loc[:, "topo_categ"] = get_topo_category(mc_nu_df)
    # evt_df.loc[:, "genie_categ"] = get_genie_category(evt_df)
    # mc_nu_df.loc[:, "genie_categ"] = get_genie_category(mc_nu_df)
    _attach_phi_degrees(evt_df, mc_nu_df)


def run_chunk_map(
    df_file: str,
    out_dir: str,
    genie_group: str,
    var_configs: List[VariableConfig],
    syst_names: Sequence[SystName],
    n_univ_cap: int = 0,
    max_splits: int = 0,
    bkgd_subtract: bool = True,
    *,
    input_stage: str = "final",
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[genie-chunk-map] no HDF splits")

    blob_root: Dict[str, Any] = {
        "kind": "genie_syst_chunk",
        "input_stage": input_stage,
        "meta": {
            "df_file": df_file,
            "splits_processed": n_use,
            "genie_group": genie_group,
            "input_stage": input_stage,
        },
        RATE_ACC_KEY: {},
        XSEC_ACC_KEY: {},
    }

    cut_by_stage: Dict[str, List[Any]] = {}
    for spec in CUT_STAGE_VAR_SPECS:
        cut_by_stage.setdefault(spec.stage_key, []).append(spec)
    final_only_vcs = genie_final_var_configs()
    skip_xsec_cut = CUT_STAGE_RATE_ONLY_SLUGS

    for i in tqdm(range(n_use), desc="HDF splits"):
        if input_stage == "final":
            mc_evt_df = pd.read_hdf(df_file, key=f"evt_{i}")
            mc_nu_df = pd.read_hdf(df_file, key=f"mcnu_{i}")
            validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context=f"split {i}:")
            mc_evt_df = mc_evt_df.copy()
            mc_nu_df = mc_nu_df.copy()
            # _prefix_mcnu_columns(mc_nu_df)
            # mc_nu_df = ensure_mc_level_phi_mcnu(mc_nu_df)
            mc_evt_df = ensure_derived_trk_kinematics_cols(mc_evt_df)
            mc_evt_df = add_reco_cc1p0pi_tki_evtdf(mc_evt_df)
            mc_evt_df = add_truth_cc1p0pi_tki_evtdf(mc_evt_df)
            mc_nu_df = add_mc_cc1p0pi_tki_mcnu(mc_nu_df)
            _annotate_topo_genie_phi(mc_evt_df, mc_nu_df)
            if n_univ_cap > 0:
                pass
            accumulate_chunk_into_blob_root(
                mc_evt_df,
                mc_nu_df,
                blob_root,
                var_configs,
                syst_names,
                bkgd_subtract=bkgd_subtract,
                skip_xsec_slugs=None,
            )
            del mc_evt_df, mc_nu_df
            gc.collect()
            continue

        # ---- sel_all: raw evt / trk / hdr + mcnu; walk ``build_pipeline()`` like cosmics/multisim.
        evt = pd.read_hdf(df_file, key=f"evt_{i}")
        trk = pd.read_hdf(df_file, key=f"trk_{i}")
        try:
            hdr = pd.read_hdf(df_file, key=f"hdr_{i}")
        except Exception:
            hdr = None
        mcnu = pd.read_hdf(df_file, key=f"mcnu_{i}")
        validate_genie_dataframes({"evt": evt, "mcnu": mcnu}, context=f"split {i}:")
        evt = evt.copy()
        trk = trk.copy()
        mcnu = mcnu.copy()
        # _prefix_mcnu_columns(mcnu)
        # mcnu = ensure_mc_level_phi_mcnu(mcnu)
        _annotate_topo_genie_phi(evt, mcnu)
        mcnu_full = mcnu

        state0: Dict[str, Any] = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}
        for stage_key, post_state in walk_pipeline(state0, sample="mc"):
            post_evt = post_state.get("evt")
            if post_evt is None or len(post_evt) == 0:
                continue
            pe, pn = _align_evt_mcnu(post_evt, mcnu_full)
            if len(pe) == 0:
                continue
            for spec in cut_by_stage.get(stage_key, ()):
                accumulate_chunk_into_blob_root(
                    pe,
                    pn,
                    blob_root,
                    [spec.var_config],
                    syst_names,
                    bkgd_subtract=bkgd_subtract,
                    skip_xsec_slugs=skip_xsec_cut,
                )
            if stage_key == FINAL_STAGE_KEY:
                for vc in final_only_vcs:
                    accumulate_chunk_into_blob_root(
                        pe,
                        pn,
                        blob_root,
                        [vc],
                        syst_names,
                        bkgd_subtract=bkgd_subtract,
                        skip_xsec_slugs=None,
                    )

        del evt, trk, hdr, mcnu, mcnu_full, state0
        gc.collect()

    stem = path.splitext(path.basename(df_file))[0]
    out_path = path.join(out_dir, "genie__%s__%s.pkl" % (genie_group, stem))
    # Atomic write so concurrent dispatchers / re-runs never see partial files.
    tmp_path = out_path + ".tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(blob_root, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, out_path)
    logger.info("wrote %s", out_path)
    return out_path


def run_chunk_merge(
    chunks_dir: str,
    out_dir: str,
    genie_group: str,
    var_configs: List[VariableConfig],
    xsec_unit: float,
    *,
    bkgd_subtract: bool = True,
    save_figs: bool = False,
    npz_path: Optional[str] = None,
) -> Dict[str, Any]:
    paths = sorted(glob.glob(path.join(chunks_dir, "genie__%s__*.pkl" % genie_group)))
    if not paths:
        raise SystemExit(
            "[genie-chunk-merge] no genie__%s__*.pkl under %s" % (genie_group, chunks_dir)
        )
    vc_by = {v.var_save_name: v for v in var_configs}

    merged = merge_genie_chunk_pickles(paths)
    os.makedirs(out_dir, exist_ok=True)

    input_stage = merged.get("input_stage", "final")
    skip_xsec_out = CUT_STAGE_RATE_ONLY_SLUGS if input_stage == "sel_all" else frozenset()

    syst_dict_out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    rate_cov_packs_by_slug: Dict[str, List[Dict[str, np.ndarray]]] = defaultdict(list)
    cv_rate_ref_by_slug: Dict[str, np.ndarray] = {}
    xsec_cov_packs_by_slug: Dict[str, List[Dict[str, np.ndarray]]] = defaultdict(list)
    cv_xsec_ref_by_slug: Dict[str, np.ndarray] = {}

    for knob in tqdm(sorted(merged[RATE_ACC_KEY].keys()), desc="GENIE knobs"):
        syst_dict_out[knob] = {}
        syst_tuple: SystName = ("mc", knob)

        for slug, rate_pack in merged[RATE_ACC_KEY][knob].items():
            vc = vc_by.get(slug)
            if vc is None:
                continue
            univ_rate = np.asarray(rate_pack["univ"], dtype=np.float64)
            cv_rate = np.asarray(rate_pack["cv"], dtype=np.float64)
            rate_cov = covariance_bundle_univ_events(
                univ_rate,
                cv_rate,
                syst_tuple,
                vc,
                "rate",
                plot=False,
                save_fig=save_figs,
                save_fig_dir=out_dir if save_figs else None,
            )

            out_pack: Dict[str, Any] = {"rate": rate_cov}
            rate_cov_packs_by_slug[slug].append(rate_cov)
            if slug not in cv_rate_ref_by_slug:
                cv_rate_ref_by_slug[slug] = np.asarray(cv_rate, dtype=np.float64).copy()

            if slug not in skip_xsec_out:
                knob_x = merged[XSEC_ACC_KEY].get(knob, {})
                if slug not in knob_x:
                    logger.warning(
                        "[genie-chunk-merge] missing xsec accumulators for knob=%s slug=%s — skip xsec",
                        knob,
                        slug,
                    )
                else:
                    acc_x = knob_x[slug]
                    univ_xsec = finalize_xsec_univ_events(acc_x, xsec_unit=xsec_unit)
                    cv_xsec = finalize_cv_sel_reco_xsec(
                        acc_x, xsec_unit=xsec_unit, bkgd_subtract=bkgd_subtract
                    )
                    xsec_cov = covariance_bundle_univ_events(
                        univ_xsec,
                        cv_xsec,
                        syst_tuple,
                        vc,
                        "xsec",
                        plot=False,
                        save_fig=save_figs,
                        save_fig_dir=out_dir if save_figs else None,
                    )
                    out_pack["xsec"] = xsec_cov
                    xsec_cov_packs_by_slug[slug].append(xsec_cov)
                    if slug not in cv_xsec_ref_by_slug:
                        cv_xsec_ref_by_slug[slug] = np.asarray(cv_xsec, dtype=np.float64).copy()

            syst_dict_out[knob][slug] = out_pack

    combined_root: Dict[str, Dict[str, Any]] = {}
    for slug in sorted(rate_cov_packs_by_slug.keys()):
        rpacks = rate_cov_packs_by_slug[slug]
        if not rpacks:
            continue
        comb_rate = combine_indep_knob_cov_packs(rpacks, cv_rate_ref_by_slug[slug])
        merged_pack: Dict[str, Any] = {"rate": comb_rate}
        xpacks = xsec_cov_packs_by_slug.get(slug, [])
        if xpacks and slug in cv_xsec_ref_by_slug:
            merged_pack["xsec"] = combine_indep_knob_cov_packs(xpacks, cv_xsec_ref_by_slug[slug])
        combined_root[slug] = merged_pack
    if combined_root:
        syst_dict_out[GENIE_MERGE_COMBINED_KEY] = combined_root

    if npz_path:
        # Nested dict of covariance arrays: store as a 0-d object array (NumPy pickles contents).
        np.savez_compressed(npz_path, syst=np.array(syst_dict_out, dtype=object))
        logger.info("wrote %s (object array key 'syst')", npz_path)

    summ = path.join(out_dir, "genie_chunk_merge_summary_%s.txt" % genie_group)
    with open(summ, "w") as f:
        f.write("# merged chunk pickles\n")
        for p in paths:
            f.write("%s\n" % p)
    return syst_dict_out


def parse_chunk_cli(argv: Optional[Sequence[str]] = None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pm = sub.add_parser(
        "chunk-map",
        help="One .df file → genie__<GROUP>__stem.pkl accumulators for that group's knobs",
    )
    pm.add_argument("--df-file", required=True)
    pm.add_argument("--out-dir", required=True)
    pm.add_argument(
        "--genie-group",
        required=True,
        choices=list(GENIE_GROUP_ORDER),
        help="Knob group / sample layout (must match dataset_locations GENIE_GROUP_GLOBS).",
    )
    pm.add_argument(
        "--input-stage",
        choices=("final", "sel_all"),
        default="final",
        help="``final``: read ``evt``+``mcnu`` only (already selected). ``sel_all``: read "
        "``evt``+``trk``+``hdr``+``mcnu``, re-run ``build_pipeline()`` per split; cut-stage "
        "variables get **rate** systematics only (no xsec tensors). Final variables keep "
        "the existing GENIE xsec recipe unchanged.",
    )
    pm.add_argument("--max-splits", type=int, default=0, help="0 = all HDF splits")
    pm.add_argument(
        "--knobs",
        default=None,
        help="Comma-separated mc.* knob names (default: all knobs for this group from makedf.geniesyst)",
    )

    rg = sub.add_parser(
        "chunk-merge",
        help="Merge genie__<GROUP>__*.pkl for one group → covariance dict / NPZ "
        "(per-knob + ``%s`` group-combined rate/xsec packs)." % GENIE_MERGE_COMBINED_KEY,
    )
    rg.add_argument("--chunks-dir", required=True)
    rg.add_argument("--out-dir", required=True)
    rg.add_argument(
        "--genie-group",
        required=True,
        choices=list(GENIE_GROUP_ORDER),
        help="Must match the chunk-map --genie-group / pickle prefix.",
    )
    rg.add_argument(
        "--input-stage",
        default=None,
        help="Override ``final``/``sel_all`` layout (default: read from first chunk pickle).",
    )
    rg.add_argument("--xsec-unit", type=float, default=1.0)
    rg.add_argument("--save-figs", action="store_true")
    rg.add_argument(
        "--no-bkgd-subtract",
        action="store_true",
        help="Match get_univ_rates(..., bkgd_subtract=False) for rate+xsec CV baseline",
    )
    rg.add_argument("--npz", default=None)

    return p.parse_args(argv)


def main_cli_chunk(argv: Optional[Sequence[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_chunk_cli(argv)

    if args.cmd == "chunk-map":
        input_stage = getattr(args, "input_stage", "final")
        var_configs = genie_all_var_configs(input_stage)
        group = args.genie_group
        if args.knobs:
            knobs = [x.strip() for x in args.knobs.split(",") if x.strip()]
        else:
            knobs = list(GENIE_GROUP_KNOBS[group])
        if not knobs:
            raise SystemExit("[chunk-map] no knobs for group %s" % group)
        syst_names = [("mc", k) for k in knobs]
        run_chunk_map(
            args.df_file,
            args.out_dir,
            group,
            var_configs,
            syst_names,
            max_splits=args.max_splits,
            input_stage=input_stage,
        )

    elif args.cmd == "chunk-merge":
        probe_paths = sorted(glob.glob(path.join(args.chunks_dir, "genie__%s__*.pkl" % args.genie_group)))
        inferred = "final"
        if probe_paths:
            with open(probe_paths[0], "rb") as f:
                probe = pickle.load(f)
            inferred = probe.get("input_stage", "final")
            meta_g = probe.get("meta", {}).get("genie_group")
            if meta_g and meta_g != args.genie_group:
                logger.warning(
                    "chunk pickle meta genie_group=%r differs from CLI %r",
                    meta_g,
                    args.genie_group,
                )
        stage = args.input_stage or inferred
        if args.input_stage and args.input_stage != inferred:
            logger.warning(
                "[chunk-merge] --input-stage=%r overrides first-chunk value %r",
                args.input_stage,
                inferred,
            )
        var_configs = genie_all_var_configs(stage)
        run_chunk_merge(
            args.chunks_dir,
            args.out_dir,
            args.genie_group,
            var_configs,
            xsec_unit=args.xsec_unit,
            bkgd_subtract=not args.no_bkgd_subtract,
            save_figs=args.save_figs,
            npz_path=args.npz,
        )


if __name__ == "__main__":
    main_cli_chunk()
