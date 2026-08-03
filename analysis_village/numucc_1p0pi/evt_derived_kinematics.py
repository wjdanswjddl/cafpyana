"""Derived per-event / per-track columns expected by :class:`variable_configs.VariableConfig`.

Used by chunked syst drivers when reading pre-stored ``evt`` / ``mcnu`` tables that may omit
reco ``phi``, opening-angle helpers, truth-level ``(..., truth, p, phi, )``, or
``(mc, mu|p|trk1|trk2, phi, ...)`` on prefixed ``mcnu`` frames.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import pad_column_name

from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key


def ensure_derived_trk_kinematics_cols(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Add reco + truth kinematic columns used by :class:`~variable_configs.VariableConfig`.

    Reco (CAF-style): ``theta_mu_p``, ``(mu|p).pfp.trk.phi`` from ``...dir.{x,y,z}``;
    truth: ``mc_theta_mu_p`` from ``mu``/``p`` ``...truth.p.dir.*`` (same ``arccos(dot)`` in rad);
    ``(mu|p).pfp.trk.truth.p.phi`` from ``...truth.p.dir.x/y`` in degrees (same ``arctan2`` as reco).
    ``opening_angle`` bins are radians on ``[0, \\pi]`` for both reco and truth openers.
    """
    if evtdf is None or len(evtdf) == 0:
        return evtdf
    if not isinstance(evtdf.columns, pd.MultiIndex):
        return evtdf

    def _col(df: pd.DataFrame, *parts: str) -> pd.Series | None:
        key = multicol_resolve_column_key(df, parts)
        if key is None:
            return None
        try:
            return df.loc[:, key]
        except Exception:
            return None

    df = evtdf
    add_theta = multicol_resolve_column_key(df, ("theta_mu_p", "", "", "", "", "", "")) is None
    add_mu_phi = multicol_resolve_column_key(df, ("mu", "pfp", "trk", "phi", "", "", "")) is None
    add_p_phi = multicol_resolve_column_key(df, ("p", "pfp", "trk", "phi", "", "", "")) is None
    add_mc_theta = multicol_resolve_column_key(df, ("mc_theta_mu_p", "", "", "", "", "", "")) is None
    add_mu_t_phi = multicol_resolve_column_key(df, ("mu", "pfp", "trk", "truth", "p", "phi", "")) is None
    add_p_t_phi = multicol_resolve_column_key(df, ("p", "pfp", "trk", "truth", "p", "phi", "")) is None

    dirs_ok = all(
        _col(df, "mu", "pfp", "trk", "dir", ax, "", "") is not None
        and _col(df, "p", "pfp", "trk", "dir", ax, "", "") is not None
        for ax in ("x", "y", "z")
    )
    mu_xy_ok = _col(df, "mu", "pfp", "trk", "dir", "x", "", "") is not None and _col(
        df, "mu", "pfp", "trk", "dir", "y", "", ""
    ) is not None
    p_xy_ok = _col(df, "p", "pfp", "trk", "dir", "x", "", "") is not None and _col(
        df, "p", "pfp", "trk", "dir", "y", "", ""
    ) is not None

    truth_dirs_ok = all(
        _col(df, "mu", "pfp", "trk", "truth", "p", "dir", ax) is not None
        and _col(df, "p", "pfp", "trk", "truth", "p", "dir", ax) is not None
        for ax in ("x", "y", "z")
    )
    mu_truth_xy_ok = _col(df, "mu", "pfp", "trk", "truth", "p", "dir", "x") is not None and _col(
        df, "mu", "pfp", "trk", "truth", "p", "dir", "y"
    ) is not None
    p_truth_xy_ok = _col(df, "p", "pfp", "trk", "truth", "p", "dir", "x") is not None and _col(
        df, "p", "pfp", "trk", "truth", "p", "dir", "y"
    ) is not None

    if not (
        (add_theta and dirs_ok)
        or (add_mu_phi and mu_xy_ok)
        or (add_p_phi and p_xy_ok)
        or (add_mc_theta and truth_dirs_ok)
        or (add_mu_t_phi and mu_truth_xy_ok)
        or (add_p_t_phi and p_truth_xy_ok)
    ):
        return df

    out = df.copy()
    if add_theta and dirs_ok:
        mx = np.asarray(_col(out, "mu", "pfp", "trk", "dir", "x", "", ""), dtype=float)
        my = np.asarray(_col(out, "mu", "pfp", "trk", "dir", "y", "", ""), dtype=float)
        mz = np.asarray(_col(out, "mu", "pfp", "trk", "dir", "z", "", ""), dtype=float)
        px = np.asarray(_col(out, "p", "pfp", "trk", "dir", "x", "", ""), dtype=float)
        py = np.asarray(_col(out, "p", "pfp", "trk", "dir", "y", "", ""), dtype=float)
        pz = np.asarray(_col(out, "p", "pfp", "trk", "dir", "z", "", ""), dtype=float)
        dot = mx * px + my * py + mz * pz
        dot = np.clip(dot, -1.0, 1.0)
        out.loc[:, pad_column_name(("theta_mu_p", "", "", "", "", "", ""), out)] = np.arccos(dot)
    if add_mu_phi and mu_xy_ok:
        mux = np.asarray(_col(out, "mu", "pfp", "trk", "dir", "x", "", ""), dtype=float)
        muy = np.asarray(_col(out, "mu", "pfp", "trk", "dir", "y", "", ""), dtype=float)
        out.loc[:, pad_column_name(("mu", "pfp", "trk", "phi", "", "", ""), out)] = np.degrees(
            np.arctan2(mux, muy)
        )
    if add_p_phi and p_xy_ok:
        px = np.asarray(_col(out, "p", "pfp", "trk", "dir", "x", "", ""), dtype=float)
        py = np.asarray(_col(out, "p", "pfp", "trk", "dir", "y", "", ""), dtype=float)
        out.loc[:, pad_column_name(("p", "pfp", "trk", "phi", "", "", ""), out)] = np.degrees(
            np.arctan2(px, py)
        )
    if add_mc_theta and truth_dirs_ok:
        mx = np.asarray(_col(out, "mu", "pfp", "trk", "truth", "p", "dir", "x"), dtype=float)
        my = np.asarray(_col(out, "mu", "pfp", "trk", "truth", "p", "dir", "y"), dtype=float)
        mz = np.asarray(_col(out, "mu", "pfp", "trk", "truth", "p", "dir", "z"), dtype=float)
        px = np.asarray(_col(out, "p", "pfp", "trk", "truth", "p", "dir", "x"), dtype=float)
        py = np.asarray(_col(out, "p", "pfp", "trk", "truth", "p", "dir", "y"), dtype=float)
        pz = np.asarray(_col(out, "p", "pfp", "trk", "truth", "p", "dir", "z"), dtype=float)
        dot = mx * px + my * py + mz * pz
        dot = np.clip(dot, -1.0, 1.0)
        out.loc[:, pad_column_name(("mc_theta_mu_p", "", "", "", "", "", ""), out)] = np.arccos(dot)
    if add_mu_t_phi and mu_truth_xy_ok:
        mux = np.asarray(_col(out, "mu", "pfp", "trk", "truth", "p", "dir", "x"), dtype=float)
        muy = np.asarray(_col(out, "mu", "pfp", "trk", "truth", "p", "dir", "y"), dtype=float)
        out.loc[:, pad_column_name(("mu", "pfp", "trk", "truth", "p", "phi", ""), out)] = np.degrees(
            np.arctan2(mux, muy)
        )
    if add_p_t_phi and p_truth_xy_ok:
        px = np.asarray(_col(out, "p", "pfp", "trk", "truth", "p", "dir", "x"), dtype=float)
        py = np.asarray(_col(out, "p", "pfp", "trk", "truth", "p", "dir", "y"), dtype=float)
        out.loc[:, pad_column_name(("p", "pfp", "trk", "truth", "p", "phi", ""), out)] = np.degrees(
            np.arctan2(px, py)
        )
    return out


def _mcnu_series(df: pd.DataFrame, parts: tuple) -> pd.Series | None:
    key = multicol_resolve_column_key(df, parts)
    if key is None:
        return None
    try:
        return df.loc[:, key]
    except Exception:
        return None


def _mcnu_has_mc_branch_phi(df: pd.DataFrame, branch: str) -> bool:
    for probe in (
        ("mc", branch, "phi", "", "", "", ""),
        ("mc", branch, "phi", ""),
        ("mc", branch, "phi"),
    ):
        if multicol_resolve_column_key(df, probe) is not None:
            return True
    return False


def ensure_mc_level_phi_mcnu(mc_nu_df: pd.DataFrame) -> pd.DataFrame:
    """Add ``mc.<branch>.phi`` (degrees) from ``mc.<branch>.dir.{x,y}`` when missing.

    ``VariableConfig`` xsec nu columns use e.g. ``('mc', 'mu', 'phi', '', '', '', '')``.
    Call **after** :func:`get_systematics_genie._prefix_mcnu_columns` so the ``mc`` group exists.
    Same ``arctan2(dir.x, dir.y)`` convention as reco phi in the GENIE driver.
    """
    if mc_nu_df is None or len(mc_nu_df) == 0:
        return mc_nu_df
    if not isinstance(mc_nu_df.columns, pd.MultiIndex):
        return mc_nu_df

    df = mc_nu_df
    modified = False
    out = df

    for branch in ("mu", "p", "trk1", "trk2"):
        if _mcnu_has_mc_branch_phi(df, branch):
            continue
        mux = _mcnu_series(df, ("mc", branch, "dir", "x"))
        muy = _mcnu_series(df, ("mc", branch, "dir", "y"))
        if mux is None or muy is None:
            continue
        if not modified:
            out = df.copy()
            modified = True
        phi_col = pad_column_name(("mc", branch, "phi"), out)
        out.loc[:, phi_col] = np.degrees(
            np.arctan2(np.asarray(mux, dtype=float), np.asarray(muy, dtype=float))
        )
    return out if modified else df
