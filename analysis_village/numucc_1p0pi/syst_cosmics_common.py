"""Shared cosmics registry and the offbeam/intime covariance used by the chunked aggregate."""
from __future__ import annotations

import logging
from os import path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from pyanalib.covariance import corr_from_fraccov, cov_from_fraccov, get_covariance_matrix

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
# Product A cosmic template: |offbeam − intime| (+ NN fill + smooth)
# ---------------------------------------------------------------------------


def _moving_average_odd(x: np.ndarray, window: int) -> np.ndarray:
    """Odd-length moving average with edge padding (shape-preserving smoother)."""
    y = np.asarray(x, dtype=float).copy()
    n = y.size
    if n == 0:
        return y
    w = int(window)
    if w < 1:
        return y
    if w % 2 == 0:
        w += 1
    if w == 1 or w > n:
        return y
    pad = w // 2
    yp = np.pad(y, pad, mode="edge")
    kernel = np.ones(w, dtype=float) / float(w)
    return np.convolve(yp, kernel, mode="valid")


def _fill_zero_denominator_frac_unc(
    frac_unc: np.ndarray,
    *,
    denom: np.ndarray,
    numer_abs: np.ndarray,
) -> np.ndarray:
    """Where denom==0 but |difference|>0, copy frac unc from the nearest valid bin."""
    u = np.asarray(frac_unc, dtype=float).copy()
    d = np.asarray(denom, dtype=float)
    num = np.asarray(numer_abs, dtype=float)
    valid = np.isfinite(d) & (d > 0.0) & np.isfinite(u)
    need = (
        np.isfinite(d)
        & (d <= 0.0)
        & np.isfinite(num)
        & (num > 0.0)
    )
    valid_idx = np.flatnonzero(valid)
    if valid_idx.size == 0:
        return u
    for i in np.flatnonzero(need):
        j = valid_idx[np.argmin(np.abs(valid_idx - i))]
        u[i] = u[j]
    return u


def product_a_cosmic_template_frac_unc_components(
    h_offbeam: np.ndarray,
    h_intime: np.ndarray,
    *,
    smooth_window: int = 5,
) -> Dict[str, np.ndarray]:
    """Return raw / intermediate / final pieces of the Product A cosmic template unc.

    Keys
    ----
    h_offbeam, h_intime
        Input histograms (1-D).
    delta
        ``|intime − offbeam|`` (unsmoothed).
    delta_smooth
        Moving-average of ``delta``.
    u_raw
        ``delta / offbeam`` (no smooth, no NN fill); 0 where offbeam==0.
    u_from_delta_smooth
        ``delta_smooth / offbeam`` before NN fill / final smooth.
    u
        Final fractional uncertainty (NN fill + smooth), same as
        :func:`product_a_cosmic_template_frac_unc`.
    """
    h_off = np.asarray(h_offbeam, dtype=float).ravel()
    h_in = np.asarray(h_intime, dtype=float).ravel()
    if h_off.shape != h_in.shape:
        raise ValueError(
            "offbeam/intime histogram shape mismatch: %s vs %s"
            % (h_off.shape, h_in.shape)
        )
    delta = np.abs(h_in - h_off)
    delta_s = _moving_average_odd(delta, smooth_window)

    u_raw = np.zeros_like(h_off, dtype=float)
    u_mid = np.zeros_like(h_off, dtype=float)
    valid = h_off > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        u_raw[valid] = delta[valid] / h_off[valid]
        u_mid[valid] = delta_s[valid] / h_off[valid]
    u_raw = np.nan_to_num(u_raw, nan=0.0, posinf=0.0, neginf=0.0)
    u_mid = np.nan_to_num(u_mid, nan=0.0, posinf=0.0, neginf=0.0)
    u = _fill_zero_denominator_frac_unc(u_mid, denom=h_off, numer_abs=delta)
    u = _moving_average_odd(u, smooth_window)
    u = np.maximum(u, 0.0)
    return {
        "h_offbeam": h_off,
        "h_intime": h_in,
        "delta": delta,
        "delta_smooth": delta_s,
        "u_raw": u_raw,
        "u_from_delta_smooth": u_mid,
        "u": u,
    }


def product_a_cosmic_template_frac_unc(
    h_offbeam: np.ndarray,
    h_intime: np.ndarray,
    *,
    smooth_window: int = 5,
) -> np.ndarray:
    """Per-bin cosmic *template* fractional uncertainty for Product A.

    Recipe (not a unisim / not flatten):
      1. ``delta = |intime − offbeam|`` (gate-scaled intime vs offbeam CV)
      2. Smooth ``delta`` with an odd moving average (damp low-stat spikes, keep shape)
      3. ``u = delta_smooth / offbeam`` where offbeam > 0
      4. If offbeam == 0 but delta > 0, copy ``u`` from the nearest valid bin
      5. Light smooth on ``u`` with the same window

    Returns 1-D ``u`` (fractional). Product A builds an uncorrelated envelope
    ``cov_frac = diag(u**2)``.
    """
    return product_a_cosmic_template_frac_unc_components(
        h_offbeam, h_intime, smooth_window=smooth_window
    )["u"]


def product_a_cosmic_template_cov_frac(
    h_offbeam: np.ndarray,
    h_intime: np.ndarray,
    *,
    smooth_window: int = 5,
) -> np.ndarray:
    """Uncorrelated fractional cov from :func:`product_a_cosmic_template_frac_unc`.

    Per-bin envelope ``|intime − offbeam| / offbeam`` with no bin–bin correlation:
    ``cov_frac = diag(u**2)``.
    """
    u = product_a_cosmic_template_frac_unc(
        h_offbeam, h_intime, smooth_window=smooth_window
    )
    return np.diag(np.square(u))


def product_a_cosmic_unisim_cov_frac(
    h_offbeam: np.ndarray,
    h_intime: np.ndarray,
) -> np.ndarray:
    """Rank-1 fractional cov: CV = offbeam, single universe = intime.

    ``cov_frac_ij = f_i f_j`` with ``f = (intime − offbeam) / max(offbeam, eps)``.
    Preserves coherent shape shifts (bin–bin correlations).
    """
    from pyanalib.covariance import get_covariance_matrix

    h_off = np.asarray(h_offbeam, dtype=float).ravel()
    h_in = np.asarray(h_intime, dtype=float).ravel()
    if h_off.shape != h_in.shape:
        raise ValueError(
            "offbeam/intime histogram shape mismatch: %s vs %s"
            % (h_off.shape, h_in.shape)
        )
    ret = get_covariance_matrix(np.stack([h_in], axis=0), h_off.copy())
    return np.asarray(ret["cov_frac"], dtype=float)


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
        try:
            cell["SelectedRate"] = selected_rate_cell_from_cosmics(cell["Cosmics"], mc_df, vc)
        except Exception:
            # Cut-stage slugs often lack the reco column on a sel_mup table.
            continue


def _sanitize_matrix_pack(pack: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    cov = np.nan_to_num(np.asarray(pack["cov"], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    cov_frac = np.nan_to_num(np.asarray(pack["cov_frac"], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    corr = np.asarray(pack["corr"], dtype=float)
    d = np.sqrt(np.maximum(np.diag(cov), 0.0))
    outer = np.outer(d, d)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(outer > 0, cov / outer, 0.0)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


def univ_stack_and_cv(h_off: np.ndarray, h_in: np.ndarray, mode: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return (univ_events, cv_events) for ``get_covariance_matrix``."""
    mode = mode.strip().lower()
    h_off = np.asarray(h_off, dtype=float)
    h_in = np.asarray(h_in, dtype=float)
    if mode == "offbeam":
        # CV = data-driven cosmic estimate; intime MC is the single alternate shape
        return np.stack([h_in], axis=0), h_off.copy()
    if mode == "intime":
        return np.stack([h_off], axis=0), h_in.copy()
    if mode == "mean":
        return np.stack([h_off, h_in], axis=0), 0.5 * (h_off + h_in)
    raise ValueError(f"Unknown --cv-mode {mode!r}; use mean|intime|offbeam")


def _plot_matrices(
    ret: Mapping[str, np.ndarray],
    var_config: Any,
    save_fig_dir: str,
    save_plots: bool,
    prefix: str,
) -> None:
    from analysis_village.numucc_1p0pi.utils import plot_heatmap

    labels = {
        "cov": "Covariance",
        "cov_frac": "Fractional covariance",
        "corr": "Correlation",
    }
    for matrix_type in ("cov", "cov_frac", "corr"):
        save_fig_name = path.join(save_fig_dir, f"{prefix}-{matrix_type}")
        cmap = "coolwarm" if matrix_type == "corr" else "viridis"
        plot_heatmap(
            ret[matrix_type],
            var_config.bins,
            plot_labels=[var_config.var_labels[1], var_config.var_labels[1], labels[matrix_type]],
            plot=False,
            cmap=cmap,
            approval="",
            save_fig=save_plots,
            save_name=save_fig_name,
        )


def process_variable_cosmics_from_histograms(
    h_off: np.ndarray,
    h_in: np.ndarray,
    var_config: Any,
    cv_mode: str,
    save_fig_dir: str,
    save_plots: bool,
    flat_uncertainty: bool = True,
    blow_up_frac_unc_threshold: float = 1.0,
) -> Dict[str, Any]:
    """Covariance + optional plots from precomputed offbeam / intime histograms."""
    h_off = np.asarray(h_off, dtype=float)
    h_in = np.asarray(h_in, dtype=float)
    univ_events, cv_events = univ_stack_and_cv(h_off, h_in, cv_mode)

    if save_plots:
        import matplotlib.pyplot as plt

        from analysis_village.numucc_1p0pi.utils import dpi, fig_ext

        plt.figure(figsize=(8, 6))
        mode_l = cv_mode.strip().lower()
        if mode_l == "offbeam":
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                linewidth=2.0,
                label="CV (offbeam data)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                label="Unisim variation (intime MC)",
            )
        elif mode_l == "intime":
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                linewidth=2.0,
                label="CV (intime MC)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                label="Unisim variation (offbeam data)",
            )
        else:
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                label="Offbeam data",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                label="Intime MC (scaled)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=cv_events,
                bins=var_config.bins,
                histtype="step",
                color="tab:blue",
                linestyle="--",
                label="CV (mean)",
            )
        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.legend(frameon=False)
        plt.savefig(
            path.join(save_fig_dir, f"{var_config.var_save_name}-cosmics-offbeam-vs-intime{fig_ext}"),
            bbox_inches="tight",
            dpi=dpi,
        )
        plt.close()

    ret = _sanitize_matrix_pack(get_covariance_matrix(univ_events, cv_events))

    if save_plots:
        from analysis_village.numucc_1p0pi.utils import plot_univ_hists

        plot_univ_hists(
            univ_events,
            cv_events,
            "Cosmics",
            var_config,
            plot=False,
            ax_titles=[
                "",
                "Events / Bin",
                f"Cosmics ({cv_mode}: CV vs variation, n_univ={univ_events.shape[0]})",
            ],
            save_fig=True,
            save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-cosmics_univ"),
        )
        _plot_matrices(
            ret,
            var_config,
            save_fig_dir,
            True,
            f"{var_config.var_save_name}-cosmics",
        )

    pay = {
        "cov": ret["cov"],
        "cov_frac": ret["cov_frac"],
        "corr": ret["corr"],
        "rate": ret,
        "univ_offbeam": h_off,
        "univ_intime": h_in,
        "cv_mode": cv_mode,
        "cv_histogram": cv_events,
        "n_univ": int(univ_events.shape[0]),
    }
    if flat_uncertainty:
        apply_flat_cosmic_uncertainty(
            pay, blow_up_frac_unc_threshold=blow_up_frac_unc_threshold
        )
    return pay


def save_cosmics_npz(syst_dict_by_var: Mapping[str, Any], out_path: str) -> None:
    np.savez_compressed(out_path, **syst_dict_by_var)
    logging.info("Wrote %s", out_path)
