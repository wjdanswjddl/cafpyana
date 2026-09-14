"""Shared GENIE covariance math (univ alias + rate/xsec path).

Canonical home for:
* ``normalize_and_infer_n_univ`` (multisim / ±σ / morph → ``univ_*``, optional ``/cv``)
* xsec tensor accumulate + finalize (``R_u @ N_gen^CV + Δbg``)

Used by ``scripts/get_systematics_genie.py``, ``syst_histcounts.py``, and notebooks.
Do not duplicate these recipes elsewhere.
"""
from __future__ import annotations

import os
from typing import Dict, List, Mapping, MutableMapping, Optional, Tuple

import numpy as np
import pandas as pd

from analysis_village.numucc_1p0pi.categories import topology_list
from analysis_village.numucc_1p0pi.utils import (
    genie_univ_weight_series,
    get_clipped_evts,
    get_response_matrix,
    signal_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

SystName = Tuple[str, str]


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


def _ensure_univ_from_leaf(
    df: pd.DataFrame,
    syst_key: Tuple[str, ...],
    src_leaf: str,
    dst_leaf: str,
    *,
    divide_by_cv: bool = False,
) -> None:
    """Alias a unisim leaf (``ps1``, ``ms1``, ``morph``) into ``univ_*`` for multisim-style loops.

    When ``divide_by_cv`` is set and a ``cv`` leaf exists, store ``src/cv`` so multisigma
    universes are relative to the knob central value (e.g. MINERvA Nature z-exp at ``b=0``).
    """
    block = df.loc[:, syst_key]
    src_rest = None
    cv_rest = None
    for c in block.columns:
        leaf = _block_column_weight_leaf(c)
        if leaf == src_leaf and src_rest is None:
            src_rest = c if isinstance(c, tuple) else (c,)
        if leaf == "cv" and cv_rest is None:
            cv_rest = c if isinstance(c, tuple) else (c,)
    if src_rest is None:
        return

    dst_rest = _replace_first_nonempty_segment(tuple(src_rest), dst_leaf)
    full_src = tuple(syst_key) + tuple(src_rest)
    full_dst = tuple(syst_key) + tuple(dst_rest)

    if full_dst in df.columns:
        return

    vals = np.asarray(df.loc[:, full_src], dtype=np.float64)
    vals = np.nan_to_num(vals, nan=1.0, posinf=1.0, neginf=0.0)
    if divide_by_cv and cv_rest is not None:
        full_cv = tuple(syst_key) + tuple(cv_rest)
        cv = np.asarray(df.loc[:, full_cv], dtype=np.float64)
        cv = np.nan_to_num(cv, nan=1.0, posinf=1.0, neginf=1.0)
        cv = np.where(cv == 0.0, 1.0, cv)
        vals = vals / cv
    vals = np.clip(vals, 0.0, None)
    df.loc[:, full_dst] = vals


def _multisigma_divide_by_cv_enabled() -> bool:
    """Whether multisigma ``ps|ms`` universes should be ``/ cv`` when a ``cv`` leaf exists.

    Controlled by ``GENIE_MULTISIGMA_DIVIDE_BY_CV`` (default ``1`` / true). Set to ``0`` /
    ``false`` / ``no`` to keep absolute CAF weights (e.g. side-by-side vs ``/cv``).
    """
    raw = os.environ.get("GENIE_MULTISIGMA_DIVIDE_BY_CV", "1").strip().lower()
    return raw not in ("0", "false", "no", "off")


def normalize_and_infer_n_univ(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: Optional[pd.DataFrame],
    syst_name: SystName,
    *,
    raise_if_missing: bool = True,
) -> int:
    """Detect multisim vs multisigma vs morph and map ±σ leaves to ``univ_*``.

    - **Multisim** (CAF type 0): existing ``univ_*`` columns.
    - **Multisigma**: ``ps1`` → ``univ_0``, ``ms1`` → ``univ_1`` when both exist; else ``ps1`` only.
      If a ``cv`` leaf is present and ``GENIE_MULTISIGMA_DIVIDE_BY_CV`` is enabled (default),
      universes are ``ps|ms / cv`` (uncertainty around knob CV).
    - **Morph unisim**: ``morph`` → ``univ_0``.

    ``mc_nu_df`` may be ``None`` (selection-stage rate-only jobs). When
    ``raise_if_missing`` is False, return ``0`` instead of raising if no leaves match.
    """
    key = tuple(syst_name)
    try:
        block_cols = mc_evt_df.loc[:, key].columns
    except Exception:
        if raise_if_missing:
            raise
        return 0
    n = _infer_multisim_n_univ(block_cols)
    if n > 0:
        return n

    leaves = set(_iter_leaf_strings(block_cols))
    div_cv = "cv" in leaves and _multisigma_divide_by_cv_enabled()

    def _alias_both(src: str, dst: str, *, divide_by_cv: bool = False) -> None:
        _ensure_univ_from_leaf(mc_evt_df, key, src, dst, divide_by_cv=divide_by_cv)
        if mc_nu_df is not None and len(mc_nu_df) > 0:
            try:
                _ensure_univ_from_leaf(mc_nu_df, key, src, dst, divide_by_cv=divide_by_cv)
            except Exception:
                pass

    if "ps1" in leaves and "ms1" in leaves:
        for src, dst in (("ps1", "univ_0"), ("ms1", "univ_1")):
            _alias_both(src, dst, divide_by_cv=div_cv)
        return 2
    if "ps1" in leaves:
        _alias_both("ps1", "univ_0", divide_by_cv=div_cv)
        return 1
    if "morph" in leaves:
        _alias_both("morph", "univ_0")
        return 1

    if not raise_if_missing:
        return 0
    raise ValueError(
        f"No univ_* columns under syst_name={syst_name!r}, and also no 'ps1' (multisigma) "
        f"or 'morph' (unisim) leaf found. Available leaves: {sorted(leaves)}"
    )


def ensure_univ_aliases(
    evt_df: pd.DataFrame,
    mc_nu_df: Optional[pd.DataFrame],
    syst_name: SystName,
) -> int:
    """Histcounts-friendly wrapper: never raises; returns 0 if knob missing."""
    return normalize_and_infer_n_univ(
        evt_df, mc_nu_df, syst_name, raise_if_missing=False
    )


# ---------------------------------------------------------------------------
# XSEC-path accumulation
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


# Histcounts / notebook aliases
empty_xsec_tensor_acc = _empty_xsec_tensor_acc
finalize_genie_xsec_univ = finalize_xsec_univ_events
finalize_genie_xsec_cv = finalize_cv_sel_reco_xsec
