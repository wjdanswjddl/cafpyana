#!/usr/bin/env python
"""Map phase: one MC ``.df`` file → pickle with summed Flux/G4/MCstat universe histograms.

Two input layouts are supported:

* ``--input-stage final``: ``.df`` is at the final-selection level; each HDF split has
  one ``evt_{i}`` table on which we run :func:`utils.get_univ_rates` (signal+
  background-subtracted rate covariance) per (syst, var).  For **G4**, the default is ``--g4-mode knobs``: weights under ``(mc, <knob>, univ_i)``
  from ``makedf.g4syst.g4_systematics``.  For **Flux**, the default is ``--flux-mode knobs``:
  weights under ``(mc, <knob>, univ_i)`` from ``makedf.bnbsyst`` (see ``--flux-knob-groups``).
  Use ``--g4-mode bundled`` / ``--flux-mode bundled`` only when the HDF has consolidated
  ``mc.G4`` / ``mc.Flux`` blocks.

* ``--input-stage sel_all``: ``.df`` is a sel_all bundle (raw ``evt`` / ``trk`` /
  ``hdr`` with the multi-universe weight columns attached). We re-run the full
  numuCC 1p0pi event-selection pipeline
  (``event_selection_pipeline_def.build_pipeline``) on each split and at every
  cut stage AND at the final stage compute per-universe weighted rate
  histograms vs an unweighted CV. This yields systematic covariances for
  cut-stage variables (nu_score, n_trks, track_score, vtx_dist, …) in addition
  to the final-selected ones.

Aggregation (``syst_multisim_aggregate.py``) sums ``univ_events`` / ``cv_events``
across files then runs ``get_covariance_matrix`` — additive across disjoint
chunks regardless of input stage.

Usage::

    python syst_multisim_chunk.py --df_file PATH.df --out_dir CHUNKS \\
        --input-stage final \\
        [--var-set final|intermediate|both] [--syst-names Flux,G4,...] \\
        [--n-universe 100] [--max-splits 0]

When ``--syst-names`` is a proper subset of MCstat/Flux/G4, output is
``nu__<names>__<stem>.pkl`` so chunks from different input dirs do not collide.
All three (default) keeps the legacy name ``nu__<stem>.pkl``.
``run_syst_multisim_chunked.sh`` sets ``--out_dir`` to ``multisim_syst-chunked-*/chunks/{Combined,MCstat}``,
``g4_syst-chunked-*/chunks``, or ``flux_syst-chunked-*/chunks`` by job type; aggregate merges all chunk roots.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from os import path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x

# turn off performance warnings
import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.pandas_helpers import pad_column_name
from pyanalib.split_df_helpers import get_n_split

from analysis_village.numucc_1p0pi.categories import get_topo_category
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    NEUTRINO_SYST_ORDER,
    build_var_configs,
    drop_bad_flux_knob_weights,
    drop_bad_g4_knob_weights,
    drop_bad_g4_weights,
    flux_mc_knob_names,
    g4_mc_knob_names,
    syst_acc_bucket_nonempty,
    syst_key_for_name,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    final_stage_var_configs,
    get_var_series,
    histogram_var,
    walk_pipeline,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates


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


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------
def _parse_syst_names(spec: str | None):
    if not spec or not str(spec).strip():
        return tuple(NEUTRINO_SYST_ORDER)
    raw = [x.strip() for x in str(spec).split(",") if x.strip()]
    unk = set(raw) - set(NEUTRINO_SYST_ORDER)
    if unk:
        raise SystemExit("[multisim-chunk] unknown --syst-names entries: %s" % sorted(unk))
    ordered = tuple(sn for sn in NEUTRINO_SYST_ORDER if sn in raw)
    if len(ordered) != len(raw):
        raise SystemExit("[multisim-chunk] duplicate syst name in --syst-names")
    return ordered


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument(
        "--input-stage",
        choices=("final", "sel_all"),
        default="final",
        help="``final``: legacy. ``sel_all``: re-run the event selection pipeline "
        "on raw evt+trk+hdr and accumulate per-universe rate histograms at every "
        "cut stage AND at the final stage.",
    )
    p.add_argument("--var-set", choices=("final", "intermediate", "both"), default="final")
    p.add_argument(
        "--syst-names",
        default=None,
        help="Comma-separated subset of MCstat,Flux,G4 (default: all). Use one name when that "
        "systematic's weights live in a separate directory of .df files.",
    )
    p.add_argument("--n-universe", type=int, default=100)
    p.add_argument(
        "--max-splits",
        type=int,
        default=0,
        help="Cap HDF splits processed (0 = all).",
    )
    p.add_argument(
        "--g4-mode",
        choices=("knobs", "bundled"),
        default="knobs",
        help="``knobs``: loop ``makedf.g4syst.g4_systematics`` → ``(mc, knob, univ_i)``. "
        "``bundled``: single ``(mc, G4, univ_i)`` / ``mc.G4`` block.",
    )
    p.add_argument(
        "--flux-mode",
        choices=("knobs", "bundled"),
        default="knobs",
        help="``knobs``: loop ``makedf.bnbsyst`` flux knobs → ``(mc, knob, univ_i)``. "
        "``bundled``: single ``(mc, Flux, univ_i)`` / ``mc.Flux`` block.",
    )
    p.add_argument(
        "--flux-knob-groups",
        default="all",
        help="With --flux-mode knobs: ``all`` → bnbsyst.regen_systematics; else e.g. "
        "beam,hadron,xsec (keys of BNB_FLUX_GROUPS).",
    )
    return p.parse_args()


# ===========================================================================
# Final-stage input path (legacy: signal+bkgd-subtracted rate via get_univ_rates).
# ===========================================================================
def _accumulate_final(args, syst_names) -> Dict[str, Any]:
    syst_active = frozenset(syst_names)
    var_configs = build_var_configs(args.var_set)
    acc_syst: Dict[str, Any] = {sn: {} for sn in NEUTRINO_SYST_ORDER}
    g4_knobs = g4_mc_knob_names() if getattr(args, "g4_mode", "knobs") == "knobs" else ()
    flux_knobs = (
        flux_mc_knob_names(args.flux_knob_groups)
        if getattr(args, "flux_mode", "knobs") == "knobs"
        else ()
    )

    def flush_evt(mc_evt_df: pd.DataFrame) -> None:
        # get_univ_rates expects evtdf.topo_categ for signal vs background topology (utils.py).
        if "topo_categ" not in mc_evt_df.columns:
            mc_evt_df = mc_evt_df.copy()
            mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
        mc_evt_df = ensure_derived_trk_kinematics_cols(mc_evt_df)
        if "G4" in syst_active:
            if args.g4_mode == "knobs":
                mc_evt_df = drop_bad_g4_knob_weights(mc_evt_df, knobs=g4_knobs, n_univ=args.n_universe)
            else:
                mc_evt_df = drop_bad_g4_weights(mc_evt_df, n_univ=args.n_universe)
        if "Flux" in syst_active and args.flux_mode == "knobs" and flux_knobs:
            mc_evt_df = drop_bad_flux_knob_weights(mc_evt_df, knobs=flux_knobs, n_univ=args.n_universe)
        if len(mc_evt_df) == 0:
            return
        for sname in syst_names:
            if sname == "Flux" and args.flux_mode == "knobs":
                if flux_knobs:
                    _accum_univ_rates_for_mc_knobs(
                        mc_evt_df, var_configs, flux_knobs, acc_syst["Flux"], "Flux", args.n_universe
                    )
                continue
            if sname == "G4" and args.g4_mode == "knobs":
                if g4_knobs:
                    _accum_univ_rates_for_mc_knobs(
                        mc_evt_df, var_configs, g4_knobs, acc_syst["G4"], "G4", args.n_universe
                    )
                continue
            sk = syst_key_for_name(sname)
            n_u = min(int(args.n_universe), _count_univ_columns(mc_evt_df, sk))
            if n_u <= 0:
                continue
            for vc in var_configs:
                try:
                    univ, cv = get_univ_rates(
                        cov_type="rate",
                        evtdf=mc_evt_df,
                        var_config=vc,
                        syst_name=sk,
                        n_univ=n_u,
                        bkgd_subtract=True,
                    )
                except Exception as ex:
                    print(f"[multisim-chunk] skip var={vc.var_save_name} syst={sname}: {ex}")
                    continue
                slug = vc.var_save_name
                if slug not in acc_syst[sname]:
                    acc_syst[sname][slug] = {
                        "univ_events": np.array(univ, dtype=float),
                        "cv_events": np.array(cv, dtype=float),
                    }
                else:
                    acc_syst[sname][slug]["univ_events"] += univ
                    acc_syst[sname][slug]["cv_events"] += cv

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if args.max_splits <= 0 else min(args.max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[multisim-chunk] no splits")
    for i in tqdm(range(n_use), desc="HDF splits"):
        mc_evt_df = pd.read_hdf(args.df_file, key=f"evt_{i}")
        flush_evt(mc_evt_df)
        del mc_evt_df
        gc.collect()
    return acc_syst


# ===========================================================================
# Sel_all input path: re-run pipeline, accumulate per-universe rates at every stage.
# ===========================================================================
def _count_univ_columns(evt_df: pd.DataFrame, syst_col_key: Any, cap: int = 512) -> int:
    """Count consecutive ``univ_0``, ``univ_1``, … columns present for this systematic."""
    n = 0
    for i in range(cap):
        if isinstance(syst_col_key, tuple):
            probe = syst_col_key + (f"univ_{i}",)
        else:
            probe = (syst_col_key, f"univ_{i}")
        if multicol_resolve_column_key(evt_df, probe) is None:
            break
        n += 1
    return n


def _univ_weight_matrix(
    evt_df: pd.DataFrame,
    syst_col_key: Any,
    n_univ: int,
) -> np.ndarray | None:
    """Return ``(n_univ_present, n_evt)`` weight matrix for one systematic, or None.

    Stops at the first missing ``univ_i`` (files may store fewer universes than
    ``--n-universe``). This matches :func:`selection_framework.mc_univ_weight_matrix`.
    """
    if evt_df is None or len(evt_df) == 0:
        return None
    rows = []
    for i in range(n_univ):
        if isinstance(syst_col_key, tuple):
            probe = syst_col_key + (f"univ_{i}",)
        else:
            probe = (syst_col_key, f"univ_{i}")
        key = multicol_resolve_column_key(evt_df, probe)
        if key is None:
            break
        try:
            w = np.asarray(evt_df[key], dtype=float).reshape(-1)
        except Exception:
            break
        rows.append(np.nan_to_num(w, nan=1.0, posinf=1.0, neginf=1.0))
    if not rows:
        return None
    return np.vstack(rows)  # shape (<= n_univ, n_evt)


def _accum_univ_rates_for_mc_knobs(
    mc_evt_df: pd.DataFrame,
    var_configs: list,
    knobs: tuple[str, ...],
    acc_cat: dict,
    log_tag: str,
    n_univ_requested: int,
) -> None:
    """Sum ``get_univ_rates`` histograms per ``(mc, knob, univ_i)`` into ``acc_cat[knob][var]``."""
    for knob in knobs:
        sk = ("mc", knob)
        n_u = min(int(n_univ_requested), _count_univ_columns(mc_evt_df, sk))
        if n_u <= 0:
            continue
        knob_acc = acc_cat.setdefault(knob, {})
        for vc in var_configs:
            try:
                univ, cv = get_univ_rates(
                    cov_type="rate",
                    evtdf=mc_evt_df,
                    var_config=vc,
                    syst_name=sk,
                    n_univ=n_u,
                    bkgd_subtract=True,
                )
            except Exception as ex:
                print(f"[multisim-chunk] skip var={vc.var_save_name} {log_tag} knob={knob}: {ex}")
                continue
            slug = vc.var_save_name
            if slug not in knob_acc:
                knob_acc[slug] = {
                    "univ_events": np.array(univ, dtype=float),
                    "cv_events": np.array(cv, dtype=float),
                }
            else:
                knob_acc[slug]["univ_events"] += univ
                knob_acc[slug]["cv_events"] += cv


def _stage_var_list(stage_key: str, final_var_configs: Sequence[Any]) -> List[Tuple[Any, str]]:
    """Variables to histogram at this stage: list of (var_config, target)."""
    out: List[Tuple[Any, str]] = []
    for spec in CUT_STAGE_VAR_SPECS:
        if spec.stage_key == stage_key:
            out.append((spec.var_config, spec.target))
    if stage_key == FINAL_STAGE_KEY:
        for vc in final_var_configs:
            out.append((vc, "evt"))
    return out


def _init_acc_sel_all(syst_names, final_var_configs) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """``acc[sn][var_save_name] = {univ_events: (n_univ, n_bin), cv_events: (n_bin,), stage_key}``."""
    acc: Dict[str, Dict[str, Dict[str, Any]]] = {sn: {} for sn in NEUTRINO_SYST_ORDER}
    # Seed shapes lazily on first split (per-bin shape depends on var_config).
    return acc


def _seed_or_add(acc_entry: Dict[str, Dict[str, Any]], var_save_name: str,
                 univ_hist: np.ndarray, cv_hist: np.ndarray, stage_key: str,
                 n_univ: int) -> None:
    if var_save_name not in acc_entry:
        acc_entry[var_save_name] = {
            "univ_events": np.array(univ_hist, dtype=float, copy=True),
            "cv_events": np.array(cv_hist, dtype=float, copy=True),
            "stage_key": stage_key,
        }
    else:
        entry = acc_entry[var_save_name]
        entry["univ_events"] += univ_hist
        entry["cv_events"] += cv_hist


def _wmats_for_mc_knobs(post_evt: pd.DataFrame, knobs: tuple[str, ...], n_univ: int) -> Dict[str, np.ndarray]:
    """``knob -> (n_univ, n_evt)`` weight matrix for each ``(mc, knob, univ_i)`` present on ``post_evt``."""
    out: Dict[str, np.ndarray] = {}
    for knob in knobs:
        sk = ("mc", knob)
        wmat = _univ_weight_matrix(post_evt, sk, n_univ)
        if wmat is not None:
            out[knob] = wmat
    return out


def _seed_sel_all_univ_hists(
    acc: dict,
    vc,
    stage_key: str,
    bins: np.ndarray,
    values: np.ndarray,
    idx: np.ndarray,
    post_evt: pd.DataFrame,
    wmats: Dict[str, np.ndarray],
    wmats_flux: Dict[str, np.ndarray],
    wmats_g4: Dict[str, np.ndarray],
    *,
    integrated: bool,
) -> None:
    """One plot row: fill ``acc`` from bundled ``wmats`` plus optional Flux/G4 knob maps."""
    if integrated:
        cv_h = np.array([float(len(post_evt))], dtype=np.float64)
        for sname, wmat in wmats.items():
            n_row = int(wmat.shape[0])
            univ_h = np.sum(wmat, axis=1).reshape(n_row, 1).astype(np.float64)
            _seed_or_add(acc[sname], vc.var_save_name, univ_h, cv_h, stage_key, n_row)
        for cat_key, wmap in (("Flux", wmats_flux), ("G4", wmats_g4)):
            for knob, wmat in wmap.items():
                n_row = int(wmat.shape[0])
                univ_h = np.sum(wmat, axis=1).reshape(n_row, 1).astype(np.float64)
                _seed_or_add(acc[cat_key].setdefault(knob, {}), vc.var_save_name, univ_h, cv_h, stage_key, n_row)
        return
    cv_h = histogram_var(values, bins)
    for sname, wmat in wmats.items():
        n_row = int(wmat.shape[0])
        univ_h = np.zeros((n_row, len(bins) - 1), dtype=np.float64)
        for u in range(n_row):
            w_per_value = wmat[u][idx]
            univ_h[u] = histogram_var(values, bins, weights=w_per_value)
        _seed_or_add(acc[sname], vc.var_save_name, univ_h, cv_h, stage_key, n_row)
    for cat_key, wmap in (("Flux", wmats_flux), ("G4", wmats_g4)):
        for knob, wmat in wmap.items():
            n_row = int(wmat.shape[0])
            univ_h = np.zeros((n_row, len(bins) - 1), dtype=np.float64)
            for u in range(n_row):
                w_per_value = wmat[u][idx]
                univ_h[u] = histogram_var(values, bins, weights=w_per_value)
            _seed_or_add(acc[cat_key].setdefault(knob, {}), vc.var_save_name, univ_h, cv_h, stage_key, n_row)


def _accumulate_sel_all(args, syst_names) -> Dict[str, Any]:
    syst_active = frozenset(syst_names)
    final_var_configs = list(final_stage_var_configs())
    g4_knobs = g4_mc_knob_names() if getattr(args, "g4_mode", "knobs") == "knobs" else ()
    flux_knobs = (
        flux_mc_knob_names(args.flux_knob_groups)
        if getattr(args, "flux_mode", "knobs") == "knobs"
        else ()
    )

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if args.max_splits <= 0 else min(args.max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[multisim-chunk] no splits")

    acc = _init_acc_sel_all(syst_names, final_var_configs)
    n_univ = int(args.n_universe)

    for i in tqdm(range(n_use), desc="HDF splits"):
        evt = pd.read_hdf(args.df_file, key=f"evt_{i}")
        try:
            trk = pd.read_hdf(args.df_file, key=f"trk_{i}")
        except KeyError:
            print(f"[multisim-chunk] no trk_{i} in {args.df_file}; sel_all requires evt+trk+hdr")
            continue
        try:
            hdr = pd.read_hdf(args.df_file, key=f"hdr_{i}")
        except KeyError:
            hdr = None

        if "G4" in syst_active:
            if args.g4_mode == "knobs":
                evt = drop_bad_g4_knob_weights(evt, knobs=g4_knobs, n_univ=n_univ)
            else:
                evt = drop_bad_g4_weights(evt, n_univ=n_univ)
        if "Flux" in syst_active and args.flux_mode == "knobs" and flux_knobs:
            evt = drop_bad_flux_knob_weights(evt, knobs=flux_knobs, n_univ=n_univ)
        if len(evt) == 0:
            del evt, trk, hdr
            continue

        state: Dict[str, Any] = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}
        for stage_key, post_state in walk_pipeline(state, sample="mc"):
            plots = _stage_var_list(stage_key, final_var_configs)
            if not plots:
                continue
            post_evt = post_state.get("evt")
            if post_evt is None or len(post_evt) == 0:
                continue
            post_evt = ensure_derived_trk_kinematics_cols(post_evt)
            post_state["evt"] = post_evt

            wmats: Dict[str, np.ndarray] = {}
            for sname in syst_names:
                if sname == "Flux" and args.flux_mode == "knobs":
                    continue
                if sname == "G4" and args.g4_mode == "knobs":
                    continue
                key = syst_key_for_name(sname)
                wmat = _univ_weight_matrix(post_evt, key, n_univ)
                if wmat is not None:
                    wmats[sname] = wmat
            wmats_flux_knob = (
                _wmats_for_mc_knobs(post_evt, flux_knobs, n_univ)
                if ("Flux" in syst_active and args.flux_mode == "knobs")
                else {}
            )
            wmats_g4_knob = (
                _wmats_for_mc_knobs(post_evt, g4_knobs, n_univ)
                if ("G4" in syst_active and args.g4_mode == "knobs")
                else {}
            )

            for vc, target in plots:
                got = get_var_series(post_state, vc, target)
                if got is None:
                    continue
                values, idx = got
                bins = np.asarray(vc.bins)
                _seed_sel_all_univ_hists(
                    acc,
                    vc,
                    stage_key,
                    bins,
                    values,
                    idx,
                    post_evt,
                    wmats,
                    wmats_flux_knob,
                    wmats_g4_knob,
                    integrated=(vc.var_save_name == "integrated"),
                )

        del state, evt, trk, hdr
        gc.collect()

    return acc


# ===========================================================================
# Driver
# ===========================================================================
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    syst_names = _parse_syst_names(args.syst_names)
    if "G4" in syst_names and args.g4_mode == "knobs" and not g4_mc_knob_names():
        raise SystemExit(
            "[multisim-chunk] ERROR: --g4-mode knobs but makedf.g4syst.g4_systematics is empty."
        )
    if "Flux" in syst_names and args.flux_mode == "knobs":
        try:
            _fk = flux_mc_knob_names(args.flux_knob_groups)
        except ValueError as ex:
            raise SystemExit("[multisim-chunk] ERROR: %s" % ex)
        if not _fk:
            raise SystemExit(
                "[multisim-chunk] ERROR: --flux-mode knobs produced an empty knob list "
                "(check --flux-knob-groups / makedf.bnbsyst)."
            )

    if args.input_stage == "sel_all":
        acc_syst = _accumulate_sel_all(args, syst_names)
        var_set_label = "sel_all"
    else:
        acc_syst = _accumulate_final(args, syst_names)
        var_set_label = args.var_set

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if args.max_splits <= 0 else min(args.max_splits, n_keys)

    stem = path.splitext(path.basename(args.df_file))[0]
    if syst_names == NEUTRINO_SYST_ORDER:
        out_leaf = "nu__{}.pkl".format(stem)
    else:
        tag = "_".join(syst_names)
        out_leaf = "nu__{}__{}.pkl".format(tag, stem)
    out_path = path.join(args.out_dir, out_leaf)
    empty_syst = [sn for sn in syst_names if not syst_acc_bucket_nonempty(sn, acc_syst.get(sn))]
    if empty_syst:
        raise SystemExit(
            "[multisim-chunk] ERROR: no per-variable histograms were accumulated for requested "
            "systematic(s) %s (check HDF weight columns, --n-universe vs stored universes, and "
            "map-phase logs for skipped variables). Refusing to write an empty chunk pickle."
            % (empty_syst,)
        )
    blob = {
        "kind": "multisim_syst_nu_chunk",
        "input_stage": args.input_stage,
        "df_file": args.df_file,
        "syst_names_computed": list(syst_names),
        "g4_mode": getattr(args, "g4_mode", "knobs"),
        "flux_mode": getattr(args, "flux_mode", "knobs"),
        "flux_knob_groups": getattr(args, "flux_knob_groups", "all"),
        "var_set": var_set_label,
        "n_univ_requested": args.n_universe,
        "splits_processed": n_use,
        "syst": acc_syst,
    }
    with open(out_path, "wb") as f:
        pickle.dump(blob, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("[multisim-chunk] wrote", out_path, "input_stage=", args.input_stage)


if __name__ == "__main__":
    main()
