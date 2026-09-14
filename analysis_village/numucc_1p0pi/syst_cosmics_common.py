"""Shared variable registry for cosmics chunk + aggregate + ``get_systematics_cosmics``."""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from pyanalib.covariance import cov_from_fraccov, corr_from_fraccov

from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig


def build_variable_configs(arg_vars: Sequence[str] | None) -> List[Any]:
    registry = {
        "integrated": VariableConfig.all_events,
        "vertex_x": VariableConfig.vertex_x,
        "vertex_y": VariableConfig.vertex_y,
        "vertex_z": VariableConfig.vertex_z,
        "muon-p": VariableConfig.muon_momentum,
        "muon-dir_z": VariableConfig.muon_direction,
        "muon-dir_x": VariableConfig.muon_direction_x,
        "muon-dir_y": VariableConfig.muon_direction_y,
        "muon-dir_phi": VariableConfig.muon_direction_phi,
        "proton-p": VariableConfig.proton_momentum,
        "proton-dir_z": VariableConfig.proton_direction,
        "proton-dir_x": VariableConfig.proton_direction_x,
        "proton-dir_y": VariableConfig.proton_direction_y,
        "proton-dir_phi": VariableConfig.proton_direction_phi,
        "muon-end_x": VariableConfig.muon_end_x,
        "muon-end_y": VariableConfig.muon_end_y,
        "muon-end_z": VariableConfig.muon_end_z,
        "opening_angle": VariableConfig.opening_angle,
        "tki-del_alpha": VariableConfig.tki_del_alpha,
        "tki-del_phi": VariableConfig.tki_del_phi,
        "tki-del_Tp": VariableConfig.tki_del_Tp,
        "tki-del_p": VariableConfig.tki_del_p,
        "tki-del_Tp_x": VariableConfig.tki_del_Tp_x,
        "tki-del_Tp_y": VariableConfig.tki_del_Tp_y,
    }
    if arg_vars:
        out = []
        for name in arg_vars:
            key = name.strip()
            if key not in registry:
                raise ValueError("Unknown variable key '%s'. Choices: %s" % (key, sorted(registry)))
            out.append(registry[key]())
        return out
    return with_final_selected_evt_variables(
        list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
        + [
            VariableConfig.vertex_x(),
            VariableConfig.vertex_y(),
            VariableConfig.vertex_z(),
            VariableConfig.muon_direction_x(),
            VariableConfig.muon_direction_y(),
            VariableConfig.proton_direction_x(),
            VariableConfig.proton_direction_y(),
            VariableConfig.opening_angle(),
        ]
    )


def frac_unc_from_cov_frac(cov_frac: np.ndarray) -> np.ndarray:
    """Per-bin fractional uncertainty sqrt(diag(cov_frac))."""
    return np.sqrt(np.maximum(np.diag(np.asarray(cov_frac, dtype=float)), 0.0))


def blown_up_frac_unc_mask(
    frac_unc: np.ndarray,
    *,
    blow_up_frac_unc_threshold: float = 1.0,
    cv_counts: np.ndarray | None = None,
    min_cv_count: float = 0.0,
) -> np.ndarray:
    """True where per-bin fractional uncertainty is unusable (limited stats / non-finite)."""
    u = np.asarray(frac_unc, dtype=float)
    bad = ~np.isfinite(u) | (u > float(blow_up_frac_unc_threshold))
    if cv_counts is not None and min_cv_count > 0:
        cv = np.asarray(cv_counts, dtype=float)
        bad |= cv < float(min_cv_count)
    return bad


def flat_uncorrelated_cov_frac(
    cov_frac: np.ndarray,
    *,
    blow_up_frac_unc_threshold: float = 1.0,
    cv_counts: np.ndarray | None = None,
    min_cv_count: float = 0.0,
) -> np.ndarray:
    """Diagonal fractional covariance: same uncertainty in every bin.

    Uses the largest sqrt(diag(cov_frac)) among bins that are not blown up.
    """
    c = np.asarray(cov_frac, dtype=float)
    n = c.shape[0]
    frac_unc = frac_unc_from_cov_frac(c)
    good = ~blown_up_frac_unc_mask(
        frac_unc,
        blow_up_frac_unc_threshold=blow_up_frac_unc_threshold,
        cv_counts=cv_counts,
        min_cv_count=min_cv_count,
    )
    if not np.any(good):
        good = np.isfinite(frac_unc)
    flat_unc = float(np.max(frac_unc[good])) if np.any(good) else 0.0
    return np.diag(np.full(n, flat_unc**2, dtype=float))


def apply_flat_cosmic_uncertainty(
    pay: MutableMapping[str, Any],
    *,
    blow_up_frac_unc_threshold: float = 1.0,
    min_cv_count: float = 0.0,
) -> Dict[str, Any]:
    """Replace cosmic ``cov_frac`` with a flat uncorrelated matrix (in-place on *pay*).

    Off-diagonal correlations from low-stat bins are dropped. Absolute ``cov`` and
    ``corr`` are updated consistently when ``cv_histogram`` is present.
    """
    cov_frac = np.asarray(pay["cov_frac"], dtype=float)
    cv = np.asarray(pay.get("cv_histogram", []), dtype=float)
    cv_counts = cv if cv.size == cov_frac.shape[0] else None

    new_cov_frac = flat_uncorrelated_cov_frac(
        cov_frac,
        blow_up_frac_unc_threshold=blow_up_frac_unc_threshold,
        cv_counts=cv_counts,
        min_cv_count=min_cv_count,
    )
    if cv_counts is not None:
        new_cov = cov_from_fraccov(new_cov_frac, cv_counts)
    else:
        old_diag = np.maximum(np.diag(cov_frac), 0.0)
        new_diag = np.diag(new_cov_frac)
        scale = np.ones_like(old_diag)
        nz = old_diag > 0
        scale[nz] = new_diag[nz] / old_diag[nz]
        new_cov = pay["cov"] * np.sqrt(np.outer(scale, scale))

    new_corr = corr_from_fraccov(new_cov_frac)
    np.fill_diagonal(new_corr, 1.0)

    pay["cov_frac"] = new_cov_frac
    pay["cov"] = np.asarray(new_cov, dtype=float)
    pay["corr"] = new_corr

    rate = pay.get("rate")
    if isinstance(rate, dict):
        rate["cov_frac"] = new_cov_frac
        rate["cov"] = pay["cov"]
        rate["corr"] = new_corr

    return dict(pay)


# ---------------------------------------------------------------------------
# Selected-rate (contamination-scaled) cosmic uncertainty
# ---------------------------------------------------------------------------


def topo_cosmic_contamination_fraction(
    mc_df: pd.DataFrame,
    var_config: Any,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bin weighted cosmic fraction in selected MC (``get_topo_category`` cut 0).

    Returns ``(frac, n_cosmic, n_total)``.
    """
    from analysis_village.numucc_1p0pi.categories import get_topo_category
    from analysis_village.numucc_1p0pi.utils import get_clipped_evts

    cut_cosmic, *_ = get_topo_category(mc_df, ret_cuts=True)
    mc_cosmic = mc_df.loc[cut_cosmic]

    if getattr(var_config, "var_save_name", "") == "integrated":
        if "pot_weight" in mc_df.columns:
            n_total = np.array([float(mc_df["pot_weight"].sum())])
            n_cosmic = np.array([float(mc_cosmic["pot_weight"].sum())])
        else:
            n_total = np.array([float(len(mc_df))])
            n_cosmic = np.array([float(len(mc_cosmic))])
    else:
        var_all, w_all = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)
        var_cos, w_cos = get_clipped_evts(mc_cosmic, var_config.var_evt_reco_col, var_config.bins)
        n_total, _ = np.histogram(var_all, bins=var_config.bins, weights=w_all)
        n_cosmic, _ = np.histogram(var_cos, bins=var_config.bins, weights=w_cos)

    frac = np.where(n_total > 0, n_cosmic / n_total, 0.0)
    return frac.astype(float), n_cosmic.astype(float), n_total.astype(float)


def scale_cov_frac_by_contamination(cov_frac: np.ndarray, contam_frac: np.ndarray) -> np.ndarray:
    """``cov_selected = cov_template * outer(f, f)``."""
    f = np.asarray(contam_frac, dtype=float).reshape(-1)
    c = np.asarray(cov_frac, dtype=float)
    return c * np.outer(f, f)


def selected_rate_uncertainty_from_cosmics(
    cosmic_pay: Mapping[str, Any],
    contam_frac: np.ndarray,
) -> Dict[str, Any]:
    """Propagate cosmic-template fractional uncertainty onto the selected event rate."""
    rate = cosmic_pay.get("rate") if isinstance(cosmic_pay.get("rate"), dict) else cosmic_pay
    cov_template = np.asarray(rate["cov_frac"], dtype=float)
    cov_selected = scale_cov_frac_by_contamination(cov_template, contam_frac)
    frac_unc_template = frac_unc_from_cov_frac(cov_template)
    f = np.asarray(contam_frac, dtype=float).reshape(-1)
    frac_unc_selected = frac_unc_template * f
    return {
        "cov_frac": cov_selected,
        "frac_unc": frac_unc_selected,
        "frac_unc_template": frac_unc_template,
        "contamination_fraction": f,
    }


def selected_rate_cell_from_cosmics(
    cosmic_pay: Mapping[str, Any],
    mc_df: pd.DataFrame,
    var_config: Any,
) -> Dict[str, Any]:
    """Build the NPZ ``SelectedRate`` cell used by the summary notebook."""
    contam_frac, n_cosmic, n_total = topo_cosmic_contamination_fraction(mc_df, var_config)
    sel_pay = selected_rate_uncertainty_from_cosmics(cosmic_pay, contam_frac)
    return {
        "rate": {"cov_frac": sel_pay["cov_frac"]},
        "contamination_fraction": sel_pay["contamination_fraction"],
        "n_cosmic": n_cosmic,
        "n_total": n_total,
    }


def attach_selected_rate_to_syst_dict(
    syst_dict: MutableMapping[str, MutableMapping[str, Any]],
    mc_df: pd.DataFrame,
    var_configs: Optional[Sequence[Any]] = None,
) -> None:
    """In-place: add ``SelectedRate`` under each variable that has a ``Cosmics`` pack.

    Summary notebooks read ``cell["SelectedRate"]["rate"]["cov_frac"]``. Without this
    step, aggregate NPZs only carry the raw cosmic template under ``Cosmics``.
    """
    if mc_df is None or len(mc_df) == 0:
        return
    vcs = list(var_configs) if var_configs is not None else []
    vsn_to_vc = {vc.var_save_name: vc for vc in vcs}
    for slug, cell in syst_dict.items():
        if not isinstance(cell, dict) or "Cosmics" not in cell:
            continue
        vc = vsn_to_vc.get(slug)
        if vc is None:
            # Best-effort: try integrated-style if bins unknown — skip
            continue
        cell["SelectedRate"] = selected_rate_cell_from_cosmics(cell["Cosmics"], mc_df, vc)
