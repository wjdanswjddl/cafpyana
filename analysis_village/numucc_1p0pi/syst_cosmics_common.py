"""Shared variable registry for cosmics chunk + aggregate + ``get_systematics_cosmics``."""
from __future__ import annotations

from typing import Any, Dict, List, MutableMapping, Sequence

import numpy as np

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
