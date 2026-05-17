#!/usr/bin/env python3
"""Data-driven conditional Gaussian constraint test (muon kinematics -> proton kinematics).

Implements the multivariate Gaussian conditional update used in model-validation studies;
methodology is aligned with the MicroBooNE data-driven model-validation prescription
(Abratenko *et al.*, "Data-driven model validation for neutrino-nucleus cross section
measurements", `arXiv:2411.03280`).

Notation (same event / same POT normalization as selected-event plots):
  * Y: constraining observables (one **or more** distributions, e.g. ``muon_p`` plus
    ``cos(theta_mu)`` bin counts concatenated).
  * X: constrained observables (e.g. ``proton_p`` bin counts).

Given joint (Gaussian) prior on stacked counts with mean ``mu = [mu_X, mu_Y]`` and covariance
``Σ`` with blocks ``Σ_XX, Σ_XY, Σ_YX, Σ_YY``, conditioning on **observed** counts ``n_Y``
(data histogram, possibly a concatenation of several constraining channels) should use an
effective ``Σ^eff_YY = Σ_YY + diag(Var_\\mathrm{data}(n_Y))`` where the added diagonal is
**only** the **data** counting (Poisson) variance ``\\approx n_Y`` evaluated from **observed**
bin contents — not Poisson noise from nominal MC ``\\mu_Y`` and not an MC finite-sample
diagonal. With that replacement, conditioning yields::

    mu_X|Y  = mu_X + Σ_XY @ (Σ^eff_YY)^{-1} @ (n_Y - mu_Y)
    Σ_XX|Y = Σ_XX - Σ_XY @ (Σ^eff_YY)^{-1} @ Σ_YX

The same formulas apply whether ``Y`` is a single distribution or a **concatenation** of
multiple constraining distributions ``Y = [Y_1; Y_2; …]``; the joint ``Σ`` simply has a
larger ``Y`` block with internal ``Y_i × Y_j`` cross-correlations.

Goodness-of-fit on proton data n_X (data--MC comparison) uses the combined covariance
``Σ^eff = Σ_MC + diag(V_\\mathrm{data})`` with ``V_\\mathrm{data}`` from observed bin counts
(same gamma/Poisson treatment as :func:`overlay_hists` in ``utils.py``)::

    chi2 = (n_X - mu)^T @ (Σ^eff)^{-1} @ (n_X - mu)

For post-fit, ``mu = mu_X|Y`` and ``Σ^eff = Σ_XX|Y + diag(V_\\mathrm{data})``; for pre-fit,
``mu = mu_X`` and ``Σ^eff = Σ_XX + diag(V_\\mathrm{data})``.

Optional decorrelation (principal axes of ``Σ_XX|Y``)::

    Σ_XX|Y = Q @ diag(lambda) @ Q^T   (eigh)
    Delta' = Q^T @ (n_X - mu_X|Y)
    epsilon_i = Delta'_i / sqrt(lambda_i)

Joint covariance model
----------------------
**Default (no ``--syst-disk-cc-root``):** per-variable fractional covariances from
:func:`analysis_village.numucc_1p0pi.utils.get_syst_unc` give ``Σ_XX`` and ``Σ_YY``.
Disk files do not encode multivariate ``Σ_XY`` across different ``VariableConfig`` objects;
this path therefore adds a **leading** cross block from same-universe multisim histograms on
``mc_df`` plus an MC-stat-style event sum (see :func:`mc_cross_cov_xy` and
:func:`multisim_cross_cov_xy`). This legacy path is **only** available for a single ``Y``.

**Recommended (``--syst-disk-cc-root``):** load **joint** covariance from ``syst_disk_CC``:
``JointMCstat/``, ``JointFlux/``, ``JointG4/`` (per-category ``joint_*_combined.npz``; summed by
:mod:`cc_joint_cov`) and, when present, ``JointGenie/joint_genie_combined.npz`` (GENIE
**rate** reweight universes). For a **single Y**, :func:`cc_joint_cov.build_joint_covariance_abs`
assembles ``[X; Y]``; for a **multi-Y** stack the new
:func:`cc_joint_cov.build_joint_multi_covariance_abs` walks every required pair
``(X, Y_i)`` and ``(Y_i, Y_j)``, re-orienting each pair NPZ to ``[A; B]``, then augments the
block diagonal with marginal sources (GENIE omitted when joint GENIE is loaded so the diagonal
GENIE term is not double-counted).

Usage
-----
Set ``NUMUCC_SYST_DISK_ROOT`` / ``--syst-disk-root`` for marginal syst files. For joint CC files,
set ``NUMUCC_SYST_DISK_CC_ROOT`` or ``--syst-disk-cc-root`` (producers:
``run_syst_cc_joint_multisim_chunked.sh``, ``run_syst_cc_joint_genie_chunked.sh``).

Single-variable constraint: ``--kinematic-pair muon_p__proton_costheta``.
Multi-variable constraint: ``--constrain-with muon_p,muon_costheta --target proton_p`` —
``Y`` is the concatenation of both muon distributions and the script reads all three joint
NPZs (proton_p × muon_p, proton_p × muon_costheta, muon_p × muon_costheta).

Data and MC load like ``selected_events.ipynb`` via
:func:`analysis_village.numucc_1p0pi.files_config.get_ana_dfs`.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

# Repository root (…/cafpyana)
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from analysis_village.numucc_1p0pi.categories import get_topo_category  # noqa: E402
from analysis_village.numucc_1p0pi.files_config import get_ana_dfs  # noqa: E402
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_clipped_evts,
    get_syst_unc as get_syst_unc_disk,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
from pyanalib.covariance import cov_from_fraccov  # noqa: E402
from pyanalib.stat_helpers import return_data_stat_err  # noqa: E402


def _safe_pinv(a: np.ndarray, rcond: float | None = None) -> np.ndarray:
    if rcond is None:
        rcond = 1e-10 * max(a.shape)
    return np.linalg.pinv(a, rcond=rcond, hermitian=True)


def _symmetrize(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def data_stat_covariance_diagonal(n_data: np.ndarray) -> np.ndarray:
    """Per-bin data counting variance (diagonal), matching ``overlay_hists`` in ``utils.py``."""
    n = np.asarray(n_data, dtype=float).reshape(-1)
    d_lo, d_hi = return_data_stat_err(n)
    sigma = 0.5 * (np.asarray(d_lo, dtype=float) + np.asarray(d_hi, dtype=float))
    return np.maximum(sigma**2, 1e-30)


def sigma_with_data_stat(sigma_mc: np.ndarray, n_data: np.ndarray) -> np.ndarray:
    """``Σ_MC + diag(V_data)`` for data--MC χ² (symmetrized)."""
    sigma_s = _symmetrize(np.asarray(sigma_mc, dtype=float))
    v = data_stat_covariance_diagonal(n_data)
    if sigma_s.shape[0] != len(v):
        raise ValueError(
            "sigma_mc shape %s incompatible with n_data length %d" % (sigma_s.shape, len(v))
        )
    return _symmetrize(sigma_s + np.diag(v))


def data_poisson_variance_diagonal(n_Y_data: np.ndarray) -> np.ndarray:
    """Per-bin Poisson variance ``\\approx n`` from **observed data** counts ``n_Y_data`` only.

    This is not derived from nominal MC ``mu_Y`` and does not add MC statistical uncertainty
    from simulation; pass the same **data** histogram used as ``n_Y`` in conditioning.
    """
    n = np.asarray(n_Y_data, dtype=float).reshape(-1)
    return np.maximum(n, 1e-30)


def sigma_yy_with_data_poisson_on_y(sigma_YY: np.ndarray, n_Y_data: np.ndarray) -> np.ndarray:
    """``Sigma^eff_YY = Sigma_YY + diag(V_\\mathrm{data})`` with ``V_\\mathrm{data}`` from **data** only.

    ``sigma_YY`` is the MC + systematic covariance on ``Y`` (disk / multisim). ``n_Y_data``
    must be the **data** histogram in the same bins (e.g. :func:`data_histogram` on the data
    dataframe). The diagonal add-on is the Gaussian/Poisson counting variance of the **data**
    observations (:func:`data_poisson_variance_diagonal`), never ``mu_Y`` and never an MC-only
    Poisson substitute.

    Augment before :func:`conditional_gaussian_update` so ``n_Y`` is not treated as a noiseless
    draw from ``Sigma_YY`` alone (avoids Kalman gain too large and ``Sigma_{XX|Y}`` too small).
    """
    yy = np.asarray(sigma_YY, dtype=float)
    v = data_poisson_variance_diagonal(n_Y_data)
    if yy.shape[0] != len(v):
        raise ValueError(
            "sigma_YY shape %s incompatible with n_Y_data length %d" % (yy.shape, len(v))
        )
    return _symmetrize(yy + np.diag(v))


# Backward-compatible alias (older notebook / callers)
sigma_yy_with_y_observation_noise = sigma_yy_with_data_poisson_on_y


def topology_total_mc(mc_df: pd.DataFrame, intime_df: pd.DataFrame | None, var_config) -> np.ndarray:
    """Stacked MC + intime cosmics prediction (same stacking as ``overlay_hists`` topology mode)."""
    vardf, _ = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)
    cuts = get_topo_category(mc_df, ret_cuts=True)
    # Boolean masks can yield pandas Series slices; np.histogram needs ndarray — coerce every row.
    var_categ = [np.asarray(vardf[i], dtype=float) for i in cuts]
    weights_categ = [list(mc_df.loc[cuts[i], "pot_weight"]) for i in range(len(cuts))]
    if intime_df is not None:
        v_int, _ = get_clipped_evts(intime_df, var_config.var_evt_reco_col, var_config.bins)
        var_categ[0] = np.concatenate(
            [np.asarray(v_int, dtype=float), np.asarray(var_categ[0], dtype=float)]
        )
        weights_categ[0] = list(intime_df["pot_weight"]) + list(weights_categ[0])
    hists = []
    for v, w in zip(var_categ, weights_categ):
        hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
        hists.append(hist_vals)
    return np.sum(hists, axis=0)


def data_histogram(data_df: pd.DataFrame, var_config) -> np.ndarray:
    vardf, _ = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
    w = data_df["pot_weight"] if "pot_weight" in data_df.columns else np.ones(len(data_df))
    vals, _ = np.histogram(vardf, bins=var_config.bins, weights=w)
    return vals


def mc_cross_cov_xy(
    mc_df: pd.DataFrame,
    intime_df: pd.DataFrame | None,
    var_X,
    var_Y,
) -> np.ndarray:
    """**MC-stat only** cross-block ``Cov(hist X, hist Y)`` from shared weighted events.

    For each row, accumulates ``w_e^2`` into ``out[ix, iy]`` where ``ix`` / ``iy`` are the
    X / Y bin indices. This captures only the finite-MC piece of the cross block, **not**
    the systematic-dominated correlation between X- and Y-bin counts. To make the
    conditional constraint meaningful you almost always want :func:`multisim_cross_cov_xy`
    on top of this — see the docstring there.
    """
    n_x = len(var_X.bin_centers)
    n_y = len(var_Y.bin_centers)
    out = np.zeros((n_x, n_y))

    def _accumulate(df: pd.DataFrame) -> None:
        vx, w_evt = get_clipped_evts(df, var_X.var_evt_reco_col, var_X.bins)
        vy, _ = get_clipped_evts(df, var_Y.var_evt_reco_col, var_Y.bins)
        w_evt = np.nan_to_num(np.asarray(w_evt, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        ix = np.clip(np.searchsorted(var_X.bins, vx, side="right") - 1, 0, n_x - 1)
        iy = np.clip(np.searchsorted(var_Y.bins, vy, side="right") - 1, 0, n_y - 1)
        ok = np.isfinite(vx) & np.isfinite(vy) & np.isfinite(w_evt)
        np.add.at(out, (ix[ok], iy[ok]), w_evt[ok] ** 2)

    _accumulate(mc_df)
    if intime_df is not None:
        _accumulate(intime_df)
    return out


# ---------------------------------------------------------------------------
# Universe-based systematic cross covariance (Σ_XY) from same-universe weights.
# ---------------------------------------------------------------------------
def _discover_multisim_knob_blocks(
    mc_df: pd.DataFrame,
    syst_names: Sequence[str] = ("Flux", "G4", "MCstat"),
    knob_probe_cap: int = 4096,
    univ_probe_cap: int = 512,
) -> list[tuple[str, str | None, int]]:
    """Find ``(syst, knob_or_None, n_univ)`` triples whose ``univ_*`` columns live in ``mc_df``.

    Two layouts are supported per category, matching the multisim chunk producers:

    * **Bundled**: ``(mc, <syst>, univ_i)`` — e.g. ``("mc", "MCstat", "univ_0")``.
    * **Knobs**:   ``(mc, <knob>, univ_i)`` for each knob name from ``g4_systematics`` /
      ``bnbsyst.regen_systematics``.

    Returns the union, deduplicated, sorted; missing categories are silently skipped so
    the caller can report which ones were and were not used.
    """
    if not isinstance(mc_df.columns, pd.MultiIndex):
        return []
    found: list[tuple[str, str | None, int]] = []

    def _count_univ_under(prefix: tuple[str, ...]) -> int:
        n = 0
        for i in range(univ_probe_cap):
            probe = prefix + (f"univ_{i}",)
            # Match either the full padded tuple or the unpadded probe.
            try:
                if probe in mc_df.columns:
                    n += 1
                    continue
            except TypeError:
                pass
            # Try padded:
            try:
                padded = probe + ("",) * (mc_df.columns.nlevels - len(probe))
                if padded in mc_df.columns:
                    n += 1
                    continue
            except Exception:
                pass
            break
        return n

    for sn in syst_names:
        # Bundled layout: (mc, sn, univ_i)
        n_bundled = _count_univ_under(("mc", sn))
        if n_bundled > 0:
            found.append((sn, None, n_bundled))
        # Knob layout: (mc, <knob>, univ_i)
        try:
            if sn == "G4":
                from analysis_village.numucc_1p0pi.syst_multisim_common import g4_mc_knob_names

                knobs = g4_mc_knob_names()
            elif sn == "Flux":
                from analysis_village.numucc_1p0pi.syst_multisim_common import flux_mc_knob_names

                knobs = flux_mc_knob_names("all")
            else:
                knobs = ()
        except Exception:
            knobs = ()
        for knob in knobs:
            n_k = _count_univ_under(("mc", knob))
            if n_k > 0:
                found.append((sn, knob, n_k))
    return found


def _multisim_univ_hists(
    mc_df: pd.DataFrame,
    var_cfg,
    syst: str,
    knob: str | None,
    n_univ: int,
    pot_scale: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-universe and CV weighted histograms of one variable for one (syst, knob) source.

    ``mc_df`` carries the per-row ``pot_weight``; each universe weight is multiplied in.
    Returns ``(univ_hist[n_univ, n_bin], cv_hist[n_bin])`` already in the same POT
    normalization as ``topology_total_mc`` (no extra scaling done here).
    """
    vals, w_pot = get_clipped_evts(mc_df, var_cfg.var_evt_reco_col, var_cfg.bins)
    w_pot = np.nan_to_num(np.asarray(w_pot, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    bins = np.asarray(var_cfg.bins)
    if pot_scale is not None:
        w_pot = w_pot * float(pot_scale)

    cv_hist, _ = np.histogram(vals, bins=bins, weights=w_pot)
    univ_hist = np.zeros((n_univ, len(bins) - 1), dtype=float)
    prefix = ("mc", syst) if knob is None else ("mc", knob)
    nlev = mc_df.columns.nlevels if isinstance(mc_df.columns, pd.MultiIndex) else 1
    for u in range(n_univ):
        probe = prefix + (f"univ_{u}",)
        if isinstance(mc_df.columns, pd.MultiIndex) and probe not in mc_df.columns:
            probe = probe + ("",) * max(0, nlev - len(probe))
        try:
            w_u = np.asarray(mc_df[probe], dtype=float).reshape(-1)
        except Exception:
            continue
        w_u = np.nan_to_num(w_u, nan=1.0, posinf=1.0, neginf=1.0)
        w = w_pot * w_u
        h, _ = np.histogram(vals, bins=bins, weights=w)
        univ_hist[u] = h
    return univ_hist, cv_hist


def multisim_cross_cov_xy(
    mc_df: pd.DataFrame,
    var_X,
    var_Y,
    syst_names: Sequence[str] = ("Flux", "G4", "MCstat"),
    n_univ_cap: int = 100,
    verbose: bool = True,
) -> tuple[np.ndarray, dict]:
    """Systematic Σ_XY from per-universe weighted histograms of X and Y on the **same** MC events.

    For each multisim source ``s`` (Flux / G4 / MCstat — flat bundled or knob-nested),
    builds per-universe histograms ``n^X(u)`` and ``n^Y(u)``, then accumulates::

        Σ_XY^s_ij = (1/N_univ^s) Σ_u (n_i^X(u) - cv_i^X) (n_j^Y(u) - cv_j^Y)

    The total Σ_XY is the sum across discovered sources (independent multisim families).
    A meta dict reports which sources were used and with how many universes — read it to
    confirm that the systematics in your Σ_XX / Σ_YY are also represented in Σ_XY.

    **Important caveats / mismatches to be aware of**

    * Within-channel Σ_XX / Σ_YY come from disk and include **MCstat + Flux + G4 + GENIE +
      Cosmics + Detector** (see :func:`analysis_village.numucc_1p0pi.utils.get_syst_unc`).
      This function only reconstructs the multisim subset present in ``mc_df`` columns.
      GENIE / Detector / Cosmics will not contribute to Σ_XY unless their universes are
      available on ``mc_df``.
    * The disk Σ_XX / Σ_YY may have been built with **background subtraction** (signal
      topology only). Here we histogram every event in ``mc_df`` with ``pot_weight`` —
      matching the **inclusive** stacked prediction used by ``topology_total_mc``.
    * Add :func:`mc_cross_cov_xy` if you also want the finite-MC stat coupling.
    """
    n_x = len(var_X.bin_centers)
    n_y = len(var_Y.bin_centers)
    blocks = _discover_multisim_knob_blocks(mc_df, syst_names=tuple(syst_names))
    meta: dict = {
        "discovered": [],
        "skipped_no_univ": [],
        "shape": (n_x, n_y),
    }
    if not blocks:
        meta["warning"] = (
            "No multisim universe columns found on mc_df under (mc, *, univ_i). Σ_XY from this "
            "function will be all zero, so the conditional update Σ_XX|Y will be ~ Σ_XX and "
            "looks unchanged in pre/post plots. The selected_events bundle written by "
            "selected_events.py does not retain universe weights; reload mc_df with a config "
            "that keeps them, e.g. get_ana_dfs(option='systs', "
            "systs_mc_df_tag='-sel_all-wgts', systs_chunk_tags=generate_tags('ah')[1:]) and "
            "apply the same final selection. Alternatively expand syst_names to whatever knob "
            "blocks your mc_df actually has."
        )
        if verbose:
            print("[multisim_cross_cov_xy] " + meta["warning"])
        return np.zeros((n_x, n_y), dtype=float), meta

    sigma_xy_total = np.zeros((n_x, n_y), dtype=float)
    for syst, knob, n_present in blocks:
        n_use = min(int(n_univ_cap), int(n_present))
        if n_use < 2:
            meta["skipped_no_univ"].append((syst, knob, n_present))
            continue
        univ_x, cv_x = _multisim_univ_hists(mc_df, var_X, syst, knob, n_use, pot_scale=None)
        univ_y, cv_y = _multisim_univ_hists(mc_df, var_Y, syst, knob, n_use, pot_scale=None)
        # Deviation matrices (n_univ, n_bin)
        dx = univ_x - cv_x[None, :]
        dy = univ_y - cv_y[None, :]
        # MC-style covariance estimate: 1/N_univ (NOT 1/(N-1)), matches get_covariance_matrix.
        block = (dx.T @ dy) / float(n_use)
        sigma_xy_total += block
        meta["discovered"].append({
            "syst": syst,
            "knob": knob,
            "n_universes_used": int(n_use),
            "frobenius_block": float(np.linalg.norm(block, ord="fro")),
            "max_abs_block": float(np.max(np.abs(block))),
        })
        if verbose:
            label = syst if knob is None else "%s/%s" % (syst, knob)
            print(
                "[multisim_cross_cov_xy] %-24s n_univ=%-4d ||Σ_XY^s||_F=%.3e max|.|=%.3e"
                % (label, n_use, meta["discovered"][-1]["frobenius_block"], meta["discovered"][-1]["max_abs_block"])
            )

    meta["frobenius_total"] = float(np.linalg.norm(sigma_xy_total, ord="fro"))
    meta["max_abs_total"] = float(np.max(np.abs(sigma_xy_total)))
    if verbose:
        print(
            "[multisim_cross_cov_xy] TOTAL ||Σ_XY||_F=%.3e max|.|=%.3e from %d source(s)"
            % (meta["frobenius_total"], meta["max_abs_total"], len(meta["discovered"]))
        )
    return sigma_xy_total, meta


def constraint_diagnostics(
    sigma_XX: np.ndarray,
    sigma_XY: np.ndarray,
    sigma_YY: np.ndarray,
    mu_X: np.ndarray,
    mu_Y: np.ndarray,
    n_Y: np.ndarray,
    sigma_XX_c: np.ndarray | None = None,
    mu_X_c: np.ndarray | None = None,
    pinv_rcond: float | None = None,
) -> dict:
    """Numerical summary used by the notebook to verify the constraint is actually doing something.

    Returns a dict with magnitudes / ratios; small values across the board mean Σ_XY is
    too weak (e.g. only finite-MC stat) and the constraint is a no-op.
    """
    if sigma_XX_c is None or mu_X_c is None:
        mu_X_c, sigma_XX_c, _ = conditional_gaussian_update(
            mu_X, mu_Y, n_Y, sigma_XX, sigma_XY, sigma_YY, pinv_rcond=pinv_rcond
        )
    diag_pre = np.diag(sigma_XX)
    diag_post = np.diag(sigma_XX_c)
    # Avoid /0 in the relative-change report
    safe_pre = np.where(np.abs(diag_pre) > 0, diag_pre, np.nan)
    rel_dia_change = (diag_post - diag_pre) / safe_pre
    sigma_x_pre = np.sqrt(np.maximum(diag_pre, 0.0))
    sigma_x_post = np.sqrt(np.maximum(diag_post, 0.0))
    sigma_y = np.sqrt(np.maximum(np.diag(sigma_YY), 0.0))
    # Bin-pair correlation strength embedded in the cross block (capped at 1 in magnitude)
    denom = np.maximum(np.outer(sigma_x_pre, sigma_y), 1e-30)
    rho_xy = sigma_XY / denom
    mu_shift = mu_X_c - mu_X
    return {
        "fro_sigma_XX": float(np.linalg.norm(sigma_XX, ord="fro")),
        "fro_sigma_YY": float(np.linalg.norm(sigma_YY, ord="fro")),
        "fro_sigma_XY": float(np.linalg.norm(sigma_XY, ord="fro")),
        "fro_sigma_XX_c": float(np.linalg.norm(sigma_XX_c, ord="fro")),
        "fro_delta_sigma_XX": float(np.linalg.norm(sigma_XX_c - sigma_XX, ord="fro")),
        "ratio_fro_delta_over_pre": float(
            np.linalg.norm(sigma_XX_c - sigma_XX, ord="fro")
            / max(np.linalg.norm(sigma_XX, ord="fro"), 1e-30)
        ),
        "max_abs_rho_XY": float(np.max(np.abs(np.where(np.isfinite(rho_xy), rho_xy, 0.0)))),
        "median_abs_rho_XY": float(
            np.median(np.abs(np.where(np.isfinite(rho_xy), rho_xy, 0.0)))
        ),
        "rel_diag_change": rel_dia_change.tolist(),
        "max_abs_rel_diag_change": float(np.nanmax(np.abs(rel_dia_change))),
        "frac_sigma_x_pre_over_mu": (sigma_x_pre / np.maximum(np.abs(mu_X), 1e-30)).tolist(),
        "frac_sigma_x_post_over_mu": (
            sigma_x_post / np.maximum(np.abs(mu_X_c), 1e-30)
        ).tolist(),
        "frac_sigma_y_over_mu": (sigma_y / np.maximum(np.abs(mu_Y), 1e-30)).tolist(),
        "mu_X": mu_X.tolist(),
        "mu_X_c": mu_X_c.tolist(),
        "mu_shift": mu_shift.tolist(),
        "mu_shift_over_sigma_x_pre": (
            mu_shift / np.maximum(sigma_x_pre, 1e-30)
        ).tolist(),
        "n_Y_minus_mu_Y_over_sigma_y": (
            (n_Y - mu_Y) / np.maximum(sigma_y, 1e-30)
        ).tolist(),
    }


def slice_data_time_batch(
    data_evt_df: pd.DataFrame,
    data_hdr_df: pd.DataFrame,
    n_splits: int,
    batch_index: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match ``selected_events.py`` exposure batching on sorted hdr rows."""
    _sorted = data_hdr_df.sort_values(["run", "evt"], kind="mergesort")
    splits = [_sorted.iloc[idx] for idx in np.array_split(np.arange(len(_sorted)), n_splits)]
    hdr = splits[batch_index]
    evt_idxs = data_evt_df.index
    common = data_evt_df.reset_index(level=[2]).index.intersection(hdr.index)
    sliced = (
        data_evt_df.reset_index(level=[2]).loc[common].reset_index().set_index(["__ntuple", "entry", "rec.slc..index"])
    )
    return sliced, hdr


def conditional_gaussian_update(
    mu_X: np.ndarray,
    mu_Y: np.ndarray,
    n_Y: np.ndarray,
    sigma_XX: np.ndarray,
    sigma_XY: np.ndarray,
    sigma_YY: np.ndarray,
    pinv_rcond: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return mu_X_cond, sigma_XX_cond, K where mu_X_cond = mu_X + K @ (n_Y - mu_Y)."""
    sigma_YX = sigma_XY.T
    inv_yy = _safe_pinv(sigma_YY, rcond=pinv_rcond)
    k = sigma_XY @ inv_yy
    mu_c = mu_X + k @ (n_Y - mu_Y)
    sigma_c = _symmetrize(sigma_XX - k @ sigma_YX)
    return mu_c, sigma_c, k


def data_mc_chi2_ndof(
    n: np.ndarray,
    mu_pred: np.ndarray,
    sigma_pred: np.ndarray,
    pinv_rcond: float | None = None,
    *,
    include_data_stat: bool = True,
) -> tuple[float, int, float]:
    """Gaussian χ² for data--MC: ``(n - μ)ᵀ (Σ_MC + V_data)⁻¹ (n - μ)``.

    When ``include_data_stat`` is true (default), adds per-bin data counting variance on
    the diagonal (same treatment as :func:`analysis_village.numucc_1p0pi.utils.overlay_hists`).

    Returns ``(chi2, ndof, chi2/ndof)`` with ``ndof = len(n)`` (full bin vector).
    """
    delta = np.asarray(n, dtype=float).reshape(-1) - np.asarray(mu_pred, dtype=float).reshape(-1)
    sigma_s = (
        sigma_with_data_stat(sigma_pred, n)
        if include_data_stat
        else _symmetrize(np.asarray(sigma_pred, dtype=float))
    )
    inv = _safe_pinv(sigma_s, rcond=pinv_rcond)
    chi2 = float(delta @ inv @ delta)
    ndof = int(len(delta))
    return chi2, ndof, chi2 / max(ndof, 1)


def chi2_and_pull(
    n_X: np.ndarray,
    mu_X_c: np.ndarray,
    sigma_XX_c: np.ndarray,
    pinv_rcond: float | None = None,
    *,
    include_data_stat: bool = True,
) -> tuple[float, float, int, np.ndarray]:
    sigma_c = (
        sigma_with_data_stat(sigma_XX_c, n_X)
        if include_data_stat
        else _symmetrize(np.asarray(sigma_XX_c, dtype=float))
    )
    inv = _safe_pinv(sigma_c, rcond=pinv_rcond)
    delta = np.asarray(n_X, dtype=float).reshape(-1) - np.asarray(mu_X_c, dtype=float).reshape(-1)
    chi2 = float(delta @ inv @ delta)
    ndof = int(len(n_X))
    pval = float(chi2_dist.sf(chi2, df=max(ndof, 1)))
    # Pull vector in bin space using same metric as chi2 (not independent per bin)
    pull = inv @ delta
    return chi2, pval, ndof, pull


def eigh_decomposition_tensions(
    sigma_XX_c: np.ndarray,
    delta_X: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Eigenvalues (ascending), Q columns = eigenvectors, epsilon_i in uncorrelated basis."""
    lam, q = np.linalg.eigh(_symmetrize(sigma_XX_c))
    # eigh returns ascending order; use all modes with lambda > 0 for tension
    delta_prime = q.T @ delta_X
    eps = np.zeros_like(lam)
    mask = lam > 1e-18
    eps[mask] = delta_prime[mask] / np.sqrt(lam[mask])
    return lam, q, eps


def plot_chi2_ndof_comparison(
    chi2_prefit_ndof: float,
    chi2_postfit_ndof: float,
    save_path: Path,
    *,
    ndof: int | None = None,
) -> None:
    """Horizontal bar chart of proton-channel χ²/Ndof pre- vs post-Gaussian constraint."""
    fig, ax = plt.subplots(figsize=(6.8, 2.9))
    y = np.arange(2)
    vals = [chi2_prefit_ndof, chi2_postfit_ndof]
    colors = ["steelblue", "darkorange"]
    ax.barh(y, vals, height=0.52, color=colors, edgecolor="0.25", linewidth=0.6)
    ax.set_yticks(y, [r"Pre-fit (data vs MC)", r"Post-fit (data vs MC$|Y$)"])
    ax.set_xlabel(r"$\chi^2 / N_{\mathrm{dof}}$ (proton channel)")
    ndof_s = (" (ndof=%d)" % ndof) if ndof is not None else ""
    ax.set_title(r"Data--MC goodness of fit, proton bins" + ndof_s, fontsize=11)
    xmax = max(vals + [1e-6])
    ax.set_xlim(0.0, xmax * 1.12)
    for yi, v in zip(y, vals):
        ax.text(v + 0.02 * xmax, yi, "%.2f" % v, va="center", fontsize=10)
    ax.grid(True, axis="x", alpha=0.28)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches="tight", dpi=200)
    plt.close(fig)


def _plot_proton_panel(
    var_X,
    n_X: np.ndarray,
    mu_X: np.ndarray,
    mu_X_c: np.ndarray,
    sigma_XX: np.ndarray,
    sigma_XX_c: np.ndarray,
    pot_label: str,
    save_path: Path,
    pinv_rcond: float | None = None,
) -> Path:
    centers = var_X.bin_centers
    bins = np.asarray(var_X.bins, dtype=float)
    widths = np.diff(bins)

    _, ndof_pre, r_pre = data_mc_chi2_ndof(n_X, mu_X, sigma_XX, pinv_rcond=pinv_rcond)
    _, ndof_post, r_post = data_mc_chi2_ndof(n_X, mu_X_c, sigma_XX_c, pinv_rcond=pinv_rcond)

    d_lo, d_hi = return_data_stat_err(n_X)
    unc_diag = np.sqrt(np.maximum(np.diag(sigma_XX), 0.0))
    con_diag = np.sqrt(np.maximum(np.diag(sigma_XX_c), 0.0))

    fig, (ax, ax_r) = plt.subplots(
        2,
        1,
        figsize=(8.5, 8.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0]},
    )
    ax.bar(
        centers,
        mu_X,
        width=widths,
        facecolor="steelblue",
        edgecolor=None,
        alpha=0.35,
        label="MC (unconstr.)",
    )
    ax.errorbar(centers, n_X, yerr=[d_lo, d_hi], fmt="ko", capsize=3, label="Data")
    ax.step(
        bins,
        np.append(mu_X_c, mu_X_c[-1]),
        where="post",
        color="darkorange",
        linewidth=2.2,
        label="MC (constr.)",
    )
    ax.bar(
        centers,
        2 * unc_diag,
        width=widths,
        bottom=mu_X - unc_diag,
        facecolor="none",
        edgecolor="gray",
        hatch="xxx",
        linewidth=0.0,
        label=r"$\pm 1\sigma$ unc. (unconstr.)",
    )
    ax.bar(
        centers,
        2 * con_diag,
        width=widths,
        bottom=mu_X_c - con_diag,
        facecolor="none",
        edgecolor="darkorange",
        hatch="xxx",
        linewidth=0.0,
        label=r"$\pm 1\sigma$ unc. (constr.)",
    )
    ylabel = pot_label
    if "POT=" in pot_label:
        ylabel = pot_label.split("(POT=")[0].rstrip()
    ax.set_ylabel(ylabel)
    # χ² lowest, SBND Internal above that, legend upper-right on top (no frames).
    ax.text(
        0.98,
        0.68,
        (
            rf"$\chi^2/\mathrm{{ndof}}={r_pre:.2f}$ (unconstr.)"
            "\n"
            rf"$\chi^2/\mathrm{{ndof}}={r_post:.2f}$ (constr.)"
        ),
        transform=ax.transAxes,
        fontsize=12,
        ha="right",
        va="top",
    )
    ax.text(
        0.98,
        0.82,
        r"$\mathbf{SBND}$ Internal",
        transform=ax.transAxes,
        fontsize=18,
        ha="right",
        va="top",
        color="rosybrown",
    )
    ax.legend(
        fontsize=12,
        ncol=2,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
    )
    ax.set_xlim(float(bins[0]), float(bins[-1]))

    r_unc = n_X / np.maximum(mu_X, 1e-12)
    r_con = n_X / np.maximum(mu_X_c, 1e-12)
    r_err_lo = d_lo / np.maximum(mu_X, 1e-12)
    r_err_hi = d_hi / np.maximum(mu_X, 1e-12)
    pred_band_unc = unc_diag / np.maximum(mu_X, 1e-12)
    pred_band_c = con_diag / np.maximum(mu_X_c, 1e-12)

    ax_r.axhline(1.0, color="k", linestyle="--", linewidth=0.8)
    ax_r.fill_between(
        centers,
        1.0 - pred_band_unc,
        1.0 + pred_band_unc,
        step="mid",
        facecolor="none",
        edgecolor="gray",
        hatch="xxx",
        linewidth=0.0,
        label=r"MC $\pm 1\sigma$ (unconstr.)",
    )
    ax_r.fill_between(
        centers,
        1.0 - pred_band_c,
        1.0 + pred_band_c,
        step="mid",
        facecolor="none",
        edgecolor="darkorange",
        hatch="xxx",
        linewidth=0.0,
        label=r"MC $\pm 1\sigma$ (constr.)",
    )
    ax_r.errorbar(centers, r_unc, yerr=[r_err_lo, r_err_hi], fmt="o", color="black", capsize=2, markersize=3, label="Data / MC (unconstr.)")
    ax_r.plot(centers, r_con, "s", color="darkorange", markersize=4, label="Data / MC (constr.)")
    ax_r.set_ylabel("Data / Pred")
    ax_r.set_xlabel(var_X.var_labels[1])
    ax_r.set_ylim(0.0, 2.0)
    ax_r.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches="tight", dpi=200)
    plt.close(fig)

    chi2_bar_path = save_path.with_name(
        save_path.stem.replace("_conditional", "") + "_chi2_ndof" + save_path.suffix
    )
    if chi2_bar_path == save_path:
        chi2_bar_path = save_path.with_name(save_path.stem + "_chi2_ndof" + save_path.suffix)
    plot_chi2_ndof_comparison(r_pre, r_post, chi2_bar_path, ndof=ndof_pre)
    return chi2_bar_path


VARIABLE_PRESET_FACTORIES: dict[str, "callable"] = {
    "muon_p": VariableConfig.muon_momentum,
    "muon_costheta": VariableConfig.muon_direction,
    "proton_p": VariableConfig.proton_momentum,
    "proton_costheta": VariableConfig.proton_direction,
}


def _resolve_variable(name: str) -> VariableConfig:
    key = str(name).strip()
    if key not in VARIABLE_PRESET_FACTORIES:
        raise SystemExit(
            "Unknown variable %r; choose one of: %s"
            % (key, ", ".join(sorted(VARIABLE_PRESET_FACTORIES)))
        )
    return VARIABLE_PRESET_FACTORIES[key]()


def _parse_var_pair(arg: str) -> tuple[VariableConfig, VariableConfig]:
    """Legacy single-Y preset: returns ``(var_Y, var_X)`` as in the original script."""
    presets = {
        "muon_p__proton_p": (VariableConfig.muon_momentum(), VariableConfig.proton_momentum()),
        "muon_costheta__proton_costheta": (
            VariableConfig.muon_direction(),
            VariableConfig.proton_direction(),
        ),
        "muon_p__proton_costheta": (VariableConfig.muon_momentum(), VariableConfig.proton_direction()),
        "muon_costheta__proton_p": (VariableConfig.muon_direction(), VariableConfig.proton_momentum()),
    }
    if arg not in presets:
        raise SystemExit(
            "Unknown --kinematic-pair %r; choose one of: %s" % (arg, ", ".join(sorted(presets)))
        )
    return presets[arg]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--kinematic-pair",
        default=None,
        choices=[
            "muon_p__proton_p",
            "muon_costheta__proton_costheta",
            "muon_p__proton_costheta",
            "muon_costheta__proton_p",
        ],
        help=(
            "Legacy single-Y preset (slug ``Y__X``). First variable is Y (constraint channel), "
            "second is X (constrained channel). Mutually exclusive with --target / --constrain-with."
        ),
    )
    p.add_argument(
        "--target",
        default=None,
        help=(
            "Multi-Y mode: name of the constrained variable X (e.g. ``proton_p``). "
            "Requires --constrain-with and --syst-disk-cc-root."
        ),
    )
    p.add_argument(
        "--constrain-with",
        default=None,
        help=(
            "Multi-Y mode: comma-separated list of constraining variables Y_1,Y_2,… "
            "(e.g. ``muon_p,muon_costheta``). Requires --target and --syst-disk-cc-root."
        ),
    )
    p.add_argument(
        "--syst-disk-root",
        default=None,
        help="Override NUMUCC_SYST_DISK_ROOT for :func:`utils.get_syst_unc` (marginal syst_disk).",
    )
    p.add_argument(
        "--syst-disk-cc-root",
        default=None,
        help="Root of syst_disk_CC (joint ``JointMCstat/``, ``JointFlux/``, ``JointG4/`` NPZs and optional ``JointGenie/*.npz``). "
        "When set, Σ is built from those joint files plus marginal diagonal augmentation; "
        "otherwise marginal-only covariances and multisim-on-mc_df cross cov are used.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for figures and summary JSON (default: cwd / conditional_constraint_out).",
    )
    p.add_argument(
        "--n-time-splits",
        type=int,
        default=1,
        help="If >1, slice data hdr into time batches like selected_events.py (default 1 = full data).",
    )
    p.add_argument(
        "--exposure-batch-index",
        type=int,
        default=0,
        help="Which time batch to use when --n-time-splits > 1.",
    )
    p.add_argument(
        "--pinv-rcond",
        type=float,
        default=None,
        help="rcond for pinv of Sigma_YY and Sigma_XX|Y (default: numpy pinv default scaled by size).",
    )
    p.add_argument(
        "--no-y-poisson-diag",
        action="store_true",
        help=(
            "Do not add **data** Poisson diag(n_Y) to Sigma_YY before conditioning "
            "(n_Y = data histogram; matches legacy behavior)."
        ),
    )
    args = p.parse_args()

    multi_mode = bool(args.target) or bool(args.constrain_with)
    if multi_mode:
        if not (args.target and args.constrain_with):
            raise SystemExit(
                "Multi-Y mode requires both --target and --constrain-with "
                "(e.g. --target proton_p --constrain-with muon_p,muon_costheta)."
            )
        if args.kinematic_pair:
            raise SystemExit(
                "Cannot combine --kinematic-pair with --target / --constrain-with; pick one mode."
            )
        if not args.syst_disk_cc_root and not os.environ.get("NUMUCC_SYST_DISK_CC_ROOT"):
            raise SystemExit(
                "Multi-Y conditional constraint needs joint syst_disk_CC NPZs "
                "(set --syst-disk-cc-root or NUMUCC_SYST_DISK_CC_ROOT)."
            )
        var_X = _resolve_variable(args.target)
        var_Ys = tuple(
            _resolve_variable(tok) for tok in str(args.constrain_with).split(",") if tok.strip()
        )
        if len(var_Ys) == 0:
            raise SystemExit("--constrain-with must list at least one variable")
        kinematic_label = "%s|%s" % (var_X.var_save_name, ",".join(v.var_save_name for v in var_Ys))
    else:
        if not args.kinematic_pair:
            args.kinematic_pair = "muon_p__proton_costheta"
        var_Y, var_X = _parse_var_pair(args.kinematic_pair)
        var_Ys = (var_Y,)
        kinematic_label = args.kinematic_pair

    out_dir = args.output_dir or (Path.cwd() / "conditional_constraint_out")
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.syst_disk_root:
        os.environ["NUMUCC_SYST_DISK_ROOT"] = os.path.abspath(args.syst_disk_root)
    if args.syst_disk_cc_root:
        os.environ["NUMUCC_SYST_DISK_CC_ROOT"] = os.path.abspath(args.syst_disk_cc_root)

    dfs = get_ana_dfs(option="selected_events")
    mc_df = dfs["mc"]
    data_evt = dfs["data"]
    intime_df = dfs["intime"]
    data_hdr = dfs["data_hdr"]

    mc_df = mc_df.copy()
    data_evt = data_evt.copy()
    intime_df = intime_df.copy()
    mc_df.loc[mc_df.mc.iscc.isna(), ("mc", "iscc")] = 999
    data_evt["mc", "iscc"] = 999

    if args.n_time_splits > 1:
        data_evt, data_hdr = slice_data_time_batch(data_evt, data_hdr, args.n_time_splits, args.exposure_batch_index)
        data_tot_pot = data_hdr["pot"].sum()
        mc_tot_pot = dfs["mc_hdr"]["pot"].sum()
        mc_pot_scale = data_tot_pot / mc_tot_pot
        mc_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_df))
        intime_gates = dfs["intime_hdr"][dfs["intime_hdr"]["first_in_subrun"] == 1]["noffbeambnb"].sum()
        data_gates = data_hdr.nbnbinfo.sum()
        f_beam = 0.0753
        scale_intime = (1 - f_beam) * data_gates / intime_gates
        intime_df["pot_weight"] = scale_intime * np.ones(len(intime_df))

    mu_X = topology_total_mc(mc_df, intime_df, var_X)
    n_X = data_histogram(data_evt, var_X)
    mu_Ys = tuple(topology_total_mc(mc_df, intime_df, vy) for vy in var_Ys)
    n_Ys = tuple(data_histogram(data_evt, vy) for vy in var_Ys)
    mu_Y = np.concatenate(mu_Ys)
    n_Y = np.concatenate(n_Ys)

    marg_root = args.syst_disk_root or os.environ.get("NUMUCC_SYST_DISK_ROOT")

    if args.syst_disk_cc_root or os.environ.get("NUMUCC_SYST_DISK_CC_ROOT"):
        from analysis_village.numucc_1p0pi.cc_joint_cov import (  # noqa: PLC0415
            build_joint_multi_covariance_abs,
            split_joint_multi_sigma,
        )

        sigma_joint = build_joint_multi_covariance_abs(
            var_X,
            var_Ys,
            mu_X,
            mu_Ys,
            syst_cc_root=args.syst_disk_cc_root,
            syst_marginal_root=marg_root,
            marginal_syst_components=None,
            marginal_genie_cov_frac_key="genie_rate",
        )
        sigma_XX, sigma_XY, sigma_YY, _slices_Y = split_joint_multi_sigma(
            sigma_joint, len(mu_X), [len(m) for m in mu_Ys]
        )
        sigma_XX = _symmetrize(sigma_XX)
        sigma_YY = _symmetrize(sigma_YY)
    else:
        if len(var_Ys) != 1:
            raise SystemExit(
                "Multi-Y constraint requires --syst-disk-cc-root (joint CC NPZs); the legacy "
                "mc_df cross-covariance path only supports a single Y variable."
            )
        var_Y_single = var_Ys[0]
        mu_Y_single = mu_Ys[0]
        _, frac_Y = get_syst_unc_disk(
            var_Y_single, syst_disk_root=args.syst_disk_root, genie_cov_frac_key="genie_rate"
        )
        _, frac_X = get_syst_unc_disk(
            var_X, syst_disk_root=args.syst_disk_root, genie_cov_frac_key="genie_rate"
        )

        sigma_YY = _symmetrize(cov_from_fraccov(frac_Y, mu_Y_single))
        sigma_XX = _symmetrize(cov_from_fraccov(frac_X, mu_X))
        sigma_XY = mc_cross_cov_xy(mc_df, intime_df, var_X, var_Y_single)

        sigma_joint = np.block([[sigma_XX, sigma_XY], [sigma_XY.T, sigma_YY]])

    sigma_YY_for_update = (
        sigma_yy_with_data_poisson_on_y(sigma_YY, n_Y) if not args.no_y_poisson_diag else sigma_YY
    )
    mu_X_c, sigma_XX_c, k_gain = conditional_gaussian_update(
        mu_X, mu_Y, n_Y, sigma_XX, sigma_XY, sigma_YY_for_update, pinv_rcond=args.pinv_rcond
    )
    chi2_prefit, _, chi2_prefit_ndof = data_mc_chi2_ndof(
        n_X, mu_X, sigma_XX, pinv_rcond=args.pinv_rcond
    )
    chi2, pval, ndof, pull = chi2_and_pull(n_X, mu_X_c, sigma_XX_c, pinv_rcond=args.pinv_rcond)
    chi2_postfit_ndof = chi2 / max(ndof, 1)
    lam, q, eps = eigh_decomposition_tensions(sigma_XX_c, n_X - mu_X_c)

    summary = {
        "kinematic_pair": kinematic_label,
        "mode": "multi_Y" if multi_mode else "single_Y",
        "joint_covariance_source": (
            "syst_disk_CC_joint_plus_marginal_diag"
            if (args.syst_disk_cc_root or os.environ.get("NUMUCC_SYST_DISK_CC_ROOT"))
            else "marginal_disk_plus_mc_multisim_cross"
        ),
        "add_data_poisson_diag_sigma_yy": (not args.no_y_poisson_diag),
        "syst_disk_cc_root": os.path.abspath(args.syst_disk_cc_root) if args.syst_disk_cc_root else None,
        "var_Y_list": [vy.var_save_name for vy in var_Ys],
        "var_X": var_X.var_save_name,
        "n_bins_Y_per_var": [int(m.size) for m in mu_Ys],
        "mu_Y": mu_Y.tolist(),
        "mu_X": mu_X.tolist(),
        "n_Y": n_Y.tolist(),
        "n_X": n_X.tolist(),
        "mu_X_conditional": mu_X_c.tolist(),
        "chi2_prefit_proton": chi2_prefit,
        "chi2_prefit_per_ndof": chi2_prefit_ndof,
        "chi2_postfit_per_ndof": chi2_postfit_ndof,
        "chi2_proton_vs_conditional_mc": chi2,
        "p_value": pval,
        "ndof": ndof,
        "eigenvalues_Sigma_XX_cond": lam.tolist(),
        "epsilon_principal": eps.tolist(),
        "frobenius_joint_sigma": float(np.linalg.norm(sigma_joint, ord="fro")),
    }
    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    pot_label = dfs.get("pot_label", "Events / bin")
    chi2_path = _plot_proton_panel(
        var_X,
        n_X,
        mu_X,
        mu_X_c,
        sigma_XX,
        sigma_XX_c,
        pot_label,
        out_dir / ("proton_%s_conditional.png" % var_X.var_save_name),
        pinv_rcond=args.pinv_rcond,
    )

    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "mode",
                    "chi2_prefit_per_ndof",
                    "chi2_postfit_per_ndof",
                    "chi2_proton_vs_conditional_mc",
                    "p_value",
                    "ndof",
                    "var_X",
                    "var_Y_list",
                )
            },
            indent=2,
        )
    )
    print("Wrote", out_dir / "summary.json")
    print("Wrote", out_dir / ("proton_%s_conditional.png" % var_X.var_save_name))
    print("Wrote", chi2_path)


if __name__ == "__main__":
    main()
