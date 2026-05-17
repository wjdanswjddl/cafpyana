#!/usr/bin/env python
"""Aggregate detector-variation chunk pickles, plot envelopes, build unisim covariances.

Inputs are pickles produced by ``syst_detvar_chunk.py`` in ``--in_dir``.
Each pickle holds per-(stage, variable, universe) histograms for a single .df file
of one WireMod model. Filenames follow ``<wiremod_tag>__<basename>.pkl``.

For every (variable, stage, WireMod model) we

    1. sum histograms across chunks (additive),
    2. plot selection-stage variables as separate figures per calorimetry parameter
       (CV ± shifts per WireMod; area-normalized main panel; ratio-to-CV panel; no titles),
    3. plot **final stage** distributions per WireMod as CV plus max(|+|,|−|) unisim
       per calo (area-normalized density main panel; ratio panel; legend describes the max-shift convention),
    4. for the **final stage**, build a unisim universe per calo parameter:
       ``n_unisim[i] = n_cv[i] + sign(d) * max(|n_p[i]-n_cv[i]|, |n_m[i]-n_cv[i]|)``
       where the sign comes from whichever shift wins,
    5. compute (cov, cov_frac, corr) for each unisim universe via
       ``pyanalib.covariance.get_covariance_matrix``,
    6. combine the four per-calo **fractional** covariances (sum ``cov_frac``), rebuild absolute
       ``cov`` from the CV spectrum, then combine WireMod models the same way into the global
       detector covariance,
    7. save ``detector_syst_dict.npz`` whose ``detector`` key plugs into
       ``analysis_village.numucc_1p0pi.utils.get_syst_unc`` (same layout as the
       chunked event-selection workflow). Also writes ``detector_by_wiremod``:
       ``var_save_name -> {wiremod_tag -> pack}``, mirroring Flux/G4
       ``*_by_knob`` nesting while keeping legacy top-level ``detector-<tag>``
       dicts for backward compatibility,
    8. save ``detector_syst_selection_dict.npz`` with the same fields for each
       selection-stage histogram, keyed by ``stage__tgt__var_save_name``, plus
       ``detector_by_wiremod`` for the same composite keys.

Usage
-----
    python syst_detvar_aggregate.py \\
        --in_dir CHUNKS_DIR \\
        --syst-disk-root SYST_DISK_ROOT

All outputs go under ``<SYST_DISK_ROOT>/Detector/`` (see ``syst_disk_layout``). If ``--syst-disk-root``
is omitted, ``NUMUCC_SYST_DISK_ROOT`` must be set.

Input globs / default work dirs for chunked drivers live in
``analysis_village.numucc_1p0pi.dataset_locations``.
"""
from __future__ import annotations

import argparse
import glob
import os
import pickle
import sys
from itertools import product
from os import path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable=None, **kwargs):
        if iterable is None:

            class _Dummy:
                def __enter__(self):
                    return self

                def __exit__(self, *args, **kwargs):
                    pass

                def update(self, *args, **kwargs):
                    pass

            return _Dummy()
        return iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.scripts.syst_detvar_chunk import (
    UniHistAcc,
    UNIVERSES,
    CALO_PARAMS,
    PER_EVT_PLOTS,
    PER_TRK_PLOTS,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.syst_disk_layout import SUB_DETECTOR, SYST_DISK_ENV
from pyanalib.covariance import cov_from_fraccov, get_covariance_matrix


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
# Color = WireMod model. Map fixed so YZ / XThetaXW always look the same.
WIREMOD_COLORS = {
    "wiremod_yz":   "C0",
    "wiremod_xtxw": "C1",
    "wiremod_xz":   "C2",  # spare for future variants
}
# Line style = which calo shift (p / m). Marker = which calo parameter.
SHIFT_STYLES = {"p": "--", "m": ":"}
CALO_LABELS = {
    "ccal":  r"$C_{cal}$",
    "alpha": r"$\alpha$",
    "beta":  r"$\beta$",
    "R":     r"$R$",
}

# Distinct colors for max-envelope lines (final-stage plots, one per calo).
CALO_COLORS_FINAL = {"ccal": "C0", "alpha": "C1", "beta": "C2", "R": "C3"}


def wiremod_display_title(tag: str) -> str:
    """Short title label for the WireMod variation type (for frac-unc / cov titles)."""
    t = tag.strip().replace("_", " ")
    tl = t.lower()
    if tl.startswith("wiremod"):
        rest = t[len("wiremod") :].strip().upper()
        return "WireMod " + rest if rest else "WireMod"
    return t.title()


def color_for_wiremod(tag: str) -> str:
    if tag in WIREMOD_COLORS:
        return WIREMOD_COLORS[tag]
    # stable fallback color
    n = len(WIREMOD_COLORS)
    return f"C{(n + abs(hash(tag)) % 10) % 10}"


# ---------------------------------------------------------------------------
# Aggregation utilities
# ---------------------------------------------------------------------------
def collect_chunks(in_dir: str) -> Dict[str, List[str]]:
    """Group pickles by wiremod_tag (the prefix before ``__``)."""
    out: Dict[str, List[str]] = {}
    for fp in sorted(glob.glob(path.join(in_dir, "*__*.pkl"))):
        base = path.basename(fp)
        tag = base.split("__", 1)[0]
        out.setdefault(tag, []).append(fp)
    for tag, files in out.items():
        print(f"[detvar-agg] wiremod_tag={tag}  {len(files)} pickle(s)")
    return out


def build_detector_by_wiremod_block(
    detector_root: Dict[str, Dict[str, dict]],
    wiremod_tags: Sequence[str],
) -> Dict[str, Dict[str, dict]]:
    """``primary_key -> {wiremod_tag -> pack}`` for keys under ``detector_root['detector']``.

    Same logical layout as multisim ``flux_by_knob`` / ``G4_by_knob`` before those
    are folded into per-variable NPZ cells: one nested map for all WireMod types
    from ``WIREMOD_DIRS`` instead of scanning ``detector-<tag>`` top-level keys.
    """
    out: Dict[str, Dict[str, dict]] = {}
    base = detector_root.get("detector")
    if not base:
        return out
    for pk in base.keys():
        nested: Dict[str, dict] = {}
        for tag in wiremod_tags:
            ent = detector_root.get(f"detector-{tag}", {}).get(pk)
            if ent is not None:
                nested[tag] = ent
        if nested:
            out[pk] = nested
    return out


def aggregate_one_wiremod(
    chunk_files: List[str],
    wiremod_tag: str = "",
    progress: bool = True,
) -> dict:
    """Sum the contents of multiple chunk pickles for one WireMod tag."""
    if not chunk_files:
        raise ValueError("no chunk files")

    merge_desc = f"merge {wiremod_tag}" if wiremod_tag else "merge chunks"
    out: Optional[dict] = None
    for cf in tqdm(
        chunk_files,
        desc=merge_desc,
        unit="pkl",
        disable=not progress,
    ):
        with open(cf, "rb") as f:
            d = pickle.load(f)
        if out is None:
            out = {
                "wiremod_tag": d["wiremod_tag"],
                "stages": d["stages"],
                "universes": list(d["universes"]),
                "chunk_pot": float(d.get("chunk_pot", 0.0)),
                "chunk_genevts": float(d.get("chunk_genevts", 0.0)),
                "histdict": {k: v for k, v in d["histdict"].items()},
                "nevt_table": dict(d.get("nevt_table", {})),
                "n_files": 1,
            }
            continue
        assert d["wiremod_tag"] == out["wiremod_tag"]
        out["chunk_pot"] += float(d.get("chunk_pot", 0.0))
        out["chunk_genevts"] += float(d.get("chunk_genevts", 0.0))
        out["n_files"] += 1
        for k, hd in d["histdict"].items():
            if k in out["histdict"]:
                out["histdict"][k] += hd
            else:
                out["histdict"][k] = hd
        for k, v in d.get("nevt_table", {}).items():
            out["nevt_table"][k] = out["nevt_table"].get(k, 0) + int(v)
    return out


# ---------------------------------------------------------------------------
# Plotting: distributions + ratio-to-CV bottom panel (no titles on distributions)
# ---------------------------------------------------------------------------
def _area_normalize_stairs(h: np.ndarray, bin_edges: np.ndarray) -> np.ndarray:
    """Scale per-bin heights so ``sum(h * Δx) == 1`` (probability density on bin centers/steps).

    Matches ``numpy.histogram(..., density=True)`` convention for histogram/stairs plots.
    """
    h = np.asarray(h, dtype=float)
    edges = np.asarray(bin_edges, dtype=float)
    w = np.diff(edges)
    if h.shape != w.shape:
        raise ValueError("hist length must match len(bin_edges) - 1")
    area = float(np.dot(h, w))
    if not np.isfinite(area) or area <= 1e-30:
        return np.zeros_like(h)
    return h / area


def _ratio_to_cv(num: np.ndarray, cv: np.ndarray) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    cv = np.asarray(cv, dtype=float)
    out = np.ones_like(num, dtype=float)
    mask = cv > 1e-12
    np.divide(num, cv, out=out, where=mask)
    out[~mask] = np.nan
    return out


def _make_distribution_ratio_axes(figsize=(7.5, 5.8)):
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.15], hspace=0.09)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    plt.setp(ax.get_xticklabels(), visible=False)
    axr.axhline(1.0, color="0.45", linestyle=":", linewidth=1.0, zorder=1)
    axr.set_ylabel("Ratio to CV")
    axr.set_ylim(0.5, 1.5)
    return fig, ax, axr


def plot_selection_var_per_calo(
    calo: str,
    bin_edges: np.ndarray,
    per_wiremod: Dict[str, Tuple[np.ndarray, List[str]]],
    cv_per_wiremod: Dict[str, np.ndarray],
    plot_xlabel: str,
    save_path: str,
    show_fig: bool = False,
):
    """Selection-variable plot: one figure per calorimetry parameter (CV + ± shifts only).

    Panels: main distribution area-normalized (density, ∫ = 1), bottom ratio of
    normalized shift curves to each WireMod's normalized CV.
    """
    fig, ax, axr = _make_distribution_ratio_axes()

    calo_p = f"{calo}_p"
    calo_m = f"{calo}_m"
    for wm, cv_h in cv_per_wiremod.items():
        wm_color = color_for_wiremod(wm)
        cv_lab = "CV (%s)" % wiremod_display_title(wm)
        ax.stairs(
            _area_normalize_stairs(cv_h, bin_edges),
            bin_edges,
            color=wm_color,
            linewidth=2,
            label=cv_lab,
            zorder=5,
        )

    for wm, (hists, universes) in per_wiremod.items():
        color = color_for_wiremod(wm)
        cv_h = cv_per_wiremod[wm]
        if calo_p not in universes or calo_m not in universes:
            continue
        ip = universes.index(calo_p)
        im = universes.index(calo_m)
        hp = hists[ip]
        hm = hists[im]
        cv_n = _area_normalize_stairs(cv_h, bin_edges)
        hp_n = _area_normalize_stairs(hp, bin_edges)
        hm_n = _area_normalize_stairs(hm, bin_edges)
        wm_title = wiremod_display_title(wm)
        ax.stairs(
            hp_n,
            bin_edges,
            color=color,
            linestyle=SHIFT_STYLES["p"],
            linewidth=2,
            alpha=0.9,
            label=f"{wm_title} " + r"$+1\sigma$",
            zorder=3,
        )
        ax.stairs(
            hm_n,
            bin_edges,
            color=color,
            linestyle=SHIFT_STYLES["m"],
            linewidth=2,
            alpha=0.9,
            label=f"{wm_title} " + r"$-1\sigma$",
            zorder=3,
        )
        axr.stairs(_ratio_to_cv(hp_n, cv_n), bin_edges, color=color, linestyle=SHIFT_STYLES["p"], linewidth=2)
        axr.stairs(_ratio_to_cv(hm_n, cv_n), bin_edges, color=color, linestyle=SHIFT_STYLES["m"], linewidth=2)

    ax.set_xlim(bin_edges[0], bin_edges[-1])
    ax.set_xlabel("")
    ax.set_ylabel("Probability density")
    ax.legend(fontsize=7, frameon=False, loc="best", ncol=1)
    axr.set_xlim(bin_edges[0], bin_edges[-1])
    axr.set_xlabel(plot_xlabel)
    fig.align_ylabels([ax, axr])
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)


def plot_final_max_per_calo(
    bin_edges: np.ndarray,
    hist_univ: np.ndarray,
    universes: List[str],
    plot_xlabel: str,
    save_path: str,
    show_fig: bool = False,
    wiremod_tag: str = "",
):
    """Final selection: CV + one curve per calo = max(|+|,|−|) unisim envelope.

    Legend notes the max-shift convention. Main panel is area-normalized density;
    ratio panel compares normalized envelope to normalized CV. CV is labeled with
    the WireMod tag and drawn in that model's color when ``wiremod_tag`` is set.
    """
    fig, ax, axr = _make_distribution_ratio_axes()

    cv_idx = universes.index("cv")
    n_cv = hist_univ[cv_idx]
    n_cv_n = _area_normalize_stairs(n_cv, bin_edges)

    cv_color = color_for_wiremod(wiremod_tag) if wiremod_tag else "black"
    cv_label = ("CV (%s)" % wiremod_display_title(wiremod_tag)) if wiremod_tag else "CV"
    ax.stairs(n_cv_n, bin_edges, color=cv_color, linewidth=2, label=cv_label, zorder=5)

    for calo in CALO_PARAMS:
        p_u = f"{calo}_p"
        m_u = f"{calo}_m"
        if p_u not in universes or m_u not in universes:
            continue
        n_p = hist_univ[universes.index(p_u)]
        n_m = hist_univ[universes.index(m_u)]
        n_max = unisim_per_calo(n_cv, n_p, n_m)
        n_max_n = _area_normalize_stairs(n_max, bin_edges)
        clr = CALO_COLORS_FINAL.get(calo, "gray")
        lbl = CALO_LABELS.get(calo, calo) + r" $\mathrm{max}(|+\sigma|, |-\sigma|)$"
        ax.stairs(n_max_n, bin_edges, color=clr, linewidth=2, label=lbl, zorder=4)
        axr.stairs(_ratio_to_cv(n_max_n, n_cv_n), bin_edges, color=clr, linewidth=2)

    ax.set_xlim(bin_edges[0], bin_edges[-1])
    ax.set_xlabel("")
    ax.set_ylabel("Probability density")
    ax.legend(fontsize=8, frameon=False, loc="best")
    axr.set_xlim(bin_edges[0], bin_edges[-1])
    axr.set_xlabel(plot_xlabel)
    fig.align_ylabels([ax, axr])
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Unisim covariance
# ---------------------------------------------------------------------------
def unisim_per_calo(
    n_cv: np.ndarray,
    n_p: np.ndarray,
    n_m: np.ndarray,
) -> np.ndarray:
    """Return n_unisim such that the per-bin shift is the larger-magnitude of (p, m).

    Mirrors the notebook: per bin, pick whichever of (p - cv, m - cv) is larger
    in absolute value, retain its sign, add to CV.
    """
    n_cv = np.asarray(n_cv, dtype=float)
    diffs_p = np.asarray(n_p, dtype=float) - n_cv
    diffs_m = np.asarray(n_m, dtype=float) - n_cv
    pick_p = np.abs(diffs_p) >= np.abs(diffs_m)
    chosen = np.where(pick_p, diffs_p, diffs_m)
    return n_cv + chosen


def cov_pack_for_universe(n_cv: np.ndarray, n_var: np.ndarray) -> dict:
    """Wrap get_covariance_matrix; sanitize NaNs from zero-CV bins."""
    ret = get_covariance_matrix(np.array([n_var], dtype=float),
                                np.asarray(n_cv, dtype=float))
    cov = np.nan_to_num(ret["cov"], nan=0.0, posinf=0.0, neginf=0.0)
    cov_frac = np.nan_to_num(ret["cov_frac"], nan=0.0, posinf=0.0, neginf=0.0)
    corr = np.nan_to_num(ret["corr"], nan=0.0, posinf=0.0, neginf=0.0)
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


def unisim_wiremod_cov_bundle_from_histacc(
    hd: UniHistAcc,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, dict], np.ndarray, np.ndarray]:
    """Per-WireMod summed fractional/calo covariance from one aggregated histogram."""
    cv_idx = hd.universes.index("cv")
    n_cv = hd.hist[cv_idx]
    n_bin = hd.hist.shape[1]
    bins = hd.bins.copy()
    cov_frac_wm = np.zeros((n_bin, n_bin))
    per_calo: Dict[str, dict] = {}
    for calo in CALO_PARAMS:
        p_u = f"{calo}_p"
        m_u = f"{calo}_m"
        if p_u not in hd.universes or m_u not in hd.universes:
            continue
        n_p = hd.hist[hd.universes.index(p_u)]
        n_m = hd.hist[hd.universes.index(m_u)]
        n_var = unisim_per_calo(n_cv, n_p, n_m)
        pack = cov_pack_for_universe(n_cv, n_var)
        pack["n_cv"] = n_cv.copy()
        pack["n_p"] = n_p.copy()
        pack["n_m"] = n_m.copy()
        pack["n_unisim"] = n_var
        pack["bins"] = bins.copy()
        per_calo[calo] = pack
        cov_frac_wm += pack["cov_frac"]
    cov_wm = cov_from_fraccov(cov_frac_wm, n_cv)
    return cov_frac_wm, cov_wm, per_calo, n_cv.copy(), bins


def selection_cov_npz_key(stage_key: str, tgt: str, vsn: str) -> str:
    """Composite key for selection-stage covariances in ``detector_syst_selection_dict.npz``."""
    return f"{stage_key}__{tgt}__{vsn}"


def plot_frac_unc_cov_wiremod(
    tag: str,
    bins: np.ndarray,
    xlab: str,
    vsn: str,
    per_calo: Dict[str, dict],
    cov_frac_wm: np.ndarray,
    plot_dir: str,
    plot_slug: Optional[str],
    show_fig: bool,
):
    """Fractional-uncertainty curve + frac-cov heatmap for one WireMod tag.

    ``plot_slug`` None uses legacy ``*_tag__vsn.png`` names; otherwise
    ``*_tag__{plot_slug}__vsn.png`` (selection plots).
    """
    tail = f"{plot_slug}__{vsn}.png" if plot_slug else f"{vsn}.png"
    plot_unisim_frac_unc(
        bins,
        xlab,
        per_calo,
        title=wiremod_display_title(tag),
        save_path=path.join(plot_dir, f"frac_unc__{tag}__{tail}"),
        show_fig=show_fig,
    )
    plot_cov_heatmap(
        bins,
        cov_frac_wm,
        title=wiremod_display_title(tag),
        save_path=path.join(plot_dir, f"cov_frac__{tag}__{tail}"),
        cbar_label="Frac. Cov.",
        show_fig=show_fig,
    )


def plot_combined_detector_frac_unc(
    bins: np.ndarray,
    xlab: str,
    vsn: str,
    cov_frac_total: np.ndarray,
    per_wiremod_packs: Dict[str, dict],
    plot_dir: str,
    plot_slug: Optional[str],
):
    tail = f"{plot_slug}__{vsn}.png" if plot_slug else f"{vsn}.png"
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.stairs(np.sqrt(np.diag(cov_frac_total)) * 100.0, bins,
              label="Total (WM + calo, quad.)", color="black", linewidth=2)
    for tag, pack in per_wiremod_packs.items():
        ax.stairs(np.sqrt(np.diag(pack["cov_frac"])) * 100.0, bins,
                  label=wiremod_display_title(tag),
                  color=color_for_wiremod(tag), linewidth=2)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(xlab)
    ax.set_ylabel("Fractional uncertainty (%)")
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(path.join(plot_dir, f"frac_unc__detector__{tail}"),
                dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Final-stage helpers
# ---------------------------------------------------------------------------
def _final_stage_key(stages: List[Tuple[str, str]]) -> str:
    return stages[-1][0]


def _hists_at_stage(agg: dict, stage_key: str) -> Dict[Tuple[str, str], UniHistAcc]:
    """Return histdict entries for the given stage_key, keyed by (var_save_name, target)."""
    out: Dict[Tuple[str, str], UniHistAcc] = {}
    for (sk, vsn, tgt), hd in agg["histdict"].items():
        if sk == stage_key:
            out[(vsn, tgt)] = hd
    return out


def plot_unisim_frac_unc(
    bin_edges: np.ndarray,
    var_label: str,
    per_calo: Dict[str, dict],
    title: str,
    save_path: str,
    show_fig: bool = False,
):
    """Per-calo fractional uncertainty curves for one variable (one WireMod model).

    Uncertainty is plotted in percent; ``title`` should name the WireMod variation only.
    """
    n_bins = max(0, len(bin_edges) - 1)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    total_sq = np.zeros(n_bins)
    for calo, pack in per_calo.items():
        frac = np.sqrt(np.diag(pack["cov_frac"])) * 100.0
        ax.stairs(frac, bin_edges, label=CALO_LABELS.get(calo, calo), linewidth=2)
        total_sq = total_sq + (frac / 100.0) ** 2
    ax.stairs(np.sqrt(total_sq) * 100.0, bin_edges, label="Total", color="black", linewidth=2)
    ax.set_xlim(bin_edges[0], bin_edges[-1])
    ax.set_xlabel(var_label)
    ax.set_ylabel("Fractional uncertainty (%)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)


def plot_cov_heatmap(
    bin_edges: np.ndarray,
    cov: np.ndarray,
    title: str,
    save_path: str,
    cbar_label: str = "Cov.",
    show_fig: bool = False,
):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    extent = [bin_edges[0], bin_edges[-1], bin_edges[0], bin_edges[-1]]
    im = ax.imshow(cov, origin="lower", aspect="auto", extent=extent, cmap="RdBu_r",
                   vmin=-np.nanmax(np.abs(cov)) if cov.size else None,
                   vmax=np.nanmax(np.abs(cov)) if cov.size else None)
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Variable label / bins lookup
# ---------------------------------------------------------------------------
def build_var_lookup() -> Dict[Tuple[str, str], VariableConfig]:
    """Map (var_save_name, target) -> VariableConfig for axis labels.

    Uses the same configs as ``syst_detvar_chunk`` (including chi2 *_new*).
    """
    lut: Dict[Tuple[str, str], VariableConfig] = {}
    for vc in PER_EVT_PLOTS:
        lut[(vc.var_save_name, "evt")] = vc
    for vc in PER_TRK_PLOTS:
        lut[(vc.var_save_name, "trk")] = vc
    return lut


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in_dir", required=True, help="Dir with chunk pickles")
    p.add_argument(
        "--syst-disk-root",
        dest="syst_disk_root",
        default=None,
        help=(
            "Syst disk layout root; all outputs under <root>/%s/. "
            "If omitted, uses environment variable %s." % (SUB_DETECTOR, SYST_DISK_ENV)
        ),
    )
    p.add_argument("--show_fig", action="store_true")
    p.add_argument("--syst_save_name", default="detector_syst_dict.npz",
                   help="npz file consumed by utils.get_syst_unc()")
    p.add_argument(
        "--selection_syst_save_name",
        default="detector_syst_selection_dict.npz",
        help="npz for unisim cov on selection-cut variables (composite keys stage__tgt__var)",
    )
    p.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm bars (clean logs / batch systems)",
    )
    args = p.parse_args()
    root = args.syst_disk_root or os.environ.get(SYST_DISK_ENV)
    if not root:
        p.error(
            "Pass --syst-disk-root or set %s (outputs always go to <root>/%s/)."
            % (SYST_DISK_ENV, SUB_DETECTOR)
        )
    root_abs = os.path.abspath(os.path.expanduser(root.rstrip("/")))
    args.syst_disk_root = root_abs
    args.out_dir = os.path.join(root_abs, SUB_DETECTOR)
    return args


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    plot_dir = path.join(args.out_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    show_prog = not args.no_progress

    chunk_groups = collect_chunks(args.in_dir)
    if not chunk_groups:
        print("[detvar-agg] no chunks found; nothing to do")
        return

    aggs: Dict[str, dict] = {}
    for tag, files in chunk_groups.items():
        print(f"[detvar-agg] aggregating wiremod={tag}  {len(files)} files")
        aggs[tag] = aggregate_one_wiremod(files, wiremod_tag=tag, progress=show_prog)

    var_lookup = build_var_lookup()

    # ------------------------------------------------------------------
    # 1) Envelope plots: per stage / variable, all WireMod + all calo vars
    # ------------------------------------------------------------------
    # Union of stage_keys / (var_save_name, target) across WireMod models
    stage_keys: List[str] = []
    stage_labels: Dict[str, str] = {}
    var_targets: List[Tuple[str, str]] = []
    for tag, agg in aggs.items():
        for sk, lab in agg["stages"]:
            if sk not in stage_keys:
                stage_keys.append(sk)
                stage_labels[sk] = lab
        for (sk, vsn, tgt) in agg["histdict"].keys():
            if (vsn, tgt) not in var_targets:
                var_targets.append((vsn, tgt))

    # Deepest stage (final selection) — needed for selection vs final envelope styles
    all_final_keys = [agg["stages"][-1][0] for agg in aggs.values()]
    if len(set(all_final_keys)) > 1:
        print(f"[detvar-agg] WARN: WireMod tags disagree on final stage: {all_final_keys}; "
              f"using the one from the first aggregator")
    final_key = all_final_keys[0]
    print(f"[detvar-agg] final stage for envelopes / unisim covariance: {final_key}")

    env_total = len(stage_keys) * len(var_targets)
    if CALO_PARAMS:
        env_total_sel = sum(
            1 for sk in stage_keys if sk != final_key for _ in var_targets
        ) * len(CALO_PARAMS)
        env_total_fin = sum(
            1 for sk in stage_keys if sk == final_key for _ in var_targets
        ) * len(aggs)
        env_total = env_total_sel + env_total_fin
    print(f"[detvar-agg] envelope plots: up to ~{env_total} figure(s)")
    for stage_key, (vsn, tgt) in tqdm(
        list(product(stage_keys, var_targets)),
        desc="envelope plots",
        unit="fig",
        disable=not show_prog,
    ):
        per_wiremod = {}
        cv_per_wiremod = {}
        bins = None
        for tag, agg in aggs.items():
            k = (stage_key, vsn, tgt)
            if k not in agg["histdict"]:
                continue
            hd: UniHistAcc = agg["histdict"][k]
            if bins is None:
                bins = hd.bins
            per_wiremod[tag] = (hd.hist.copy(), list(hd.universes))
            cv_idx = hd.universes.index("cv")
            cv_per_wiremod[tag] = hd.hist[cv_idx].copy()
        if not per_wiremod:
            continue
        vc = var_lookup.get((vsn, tgt))
        xlab = vc.var_labels[0] if vc is not None else vsn

        if stage_key == final_key:
            # Final stage: max(|+|,|−|) unisim envelope per calo (one file per WireMod tag).
            for tag, agg in aggs.items():
                k = (stage_key, vsn, tgt)
                if k not in agg["histdict"]:
                    continue
                hd = agg["histdict"][k]
                save_p = path.join(
                    plot_dir,
                    f"envelope__{stage_key}__{tgt}__{vsn}__max__{tag}.png",
                )
                plot_final_max_per_calo(
                    hd.bins,
                    hd.hist,
                    list(hd.universes),
                    xlab,
                    save_p,
                    show_fig=args.show_fig,
                    wiremod_tag=tag,
                )
        else:
            # Selection-stage variables: separate figure per calorimetry variation.
            for calo in CALO_PARAMS:
                save_p = path.join(
                    plot_dir,
                    f"envelope__{stage_key}__{tgt}__{vsn}__calo_{calo}.png",
                )
                plot_selection_var_per_calo(
                    calo,
                    bins,
                    per_wiremod,
                    cv_per_wiremod,
                    xlab,
                    save_p,
                    show_fig=args.show_fig,
                )

    # ------------------------------------------------------------------
    # 2) Unisim covariance + frac unc + heatmap on FINAL stage variables (evt only)
    # ------------------------------------------------------------------

    # detector_dict[..."detector"] is the one consumed by utils.get_syst_unc.
    # We also keep per-WireMod and per-WireMod x per-calo breakdowns for plotting.
    detector_dict: Dict[str, Dict[str, dict]] = {
        "detector": {},
    }
    for tag in aggs.keys():
        detector_dict[f"detector-{tag}"] = {}

    breakdown_dict: Dict[str, dict] = {}  # detailed nested breakdown

    final_evt_tasks = [
        (vsn, tgt)
        for vsn, tgt in var_targets
        if tgt == "evt"
        and any(((final_key, vsn, tgt) in agg["histdict"]) for agg in aggs.values())
    ]
    for vsn, tgt in tqdm(
        final_evt_tasks,
        desc="unisim cov + plots (final)",
        unit="var",
        disable=not show_prog,
    ):
        n_bin = None
        cov_frac_total = None
        per_wiremod_packs: Dict[str, dict] = {}

        for tag, agg in aggs.items():
            key = (final_key, vsn, tgt)
            if key not in agg["histdict"]:
                continue
            hd = agg["histdict"][key]
            if n_bin is None:
                n_bin = hd.hist.shape[1]
                cov_frac_total = np.zeros((n_bin, n_bin))
            if hd.hist.shape[1] != n_bin:
                print(f"[detvar-agg] WARN: bin mismatch for {vsn} ({tag}); skipping")
                continue

            cov_frac_wm, cov_wm, per_calo, n_cv, bins_wm = (
                unisim_wiremod_cov_bundle_from_histacc(hd)
            )

            per_wiremod_packs[tag] = {
                "per_calo": per_calo,
                "cov_frac": cov_frac_wm,
                "cov": cov_wm,
                "n_cv": n_cv,
                "bins": bins_wm,
            }
            cov_frac_total += cov_frac_wm

            detector_dict[f"detector-{tag}"][vsn] = {
                "cov_frac": cov_frac_wm,
                "cov": cov_wm,
                "corr": _safe_corr_from_cov(cov_wm),
                "n_cv": n_cv,
                "per_calo": {c: per_calo[c] for c in per_calo},
                "bins": bins_wm,
            }

            vc = var_lookup.get((vsn, tgt))
            xlab = vc.var_labels[0] if vc is not None else vsn
            plot_frac_unc_cov_wiremod(
                tag,
                bins_wm,
                xlab,
                vsn,
                per_calo,
                cov_frac_wm,
                plot_dir,
                None,
                args.show_fig,
            )

        if cov_frac_total is None:
            continue
        any_wm = next(iter(per_wiremod_packs.values()))
        cov_total = cov_from_fraccov(cov_frac_total, any_wm["n_cv"])
        detector_dict["detector"][vsn] = {
            "cov_frac": cov_frac_total,
            "cov": cov_total,
            "corr": _safe_corr_from_cov(cov_total),
            "n_cv": any_wm["n_cv"],
            "bins": any_wm["bins"],
        }
        breakdown_dict[vsn] = {
            "bins": any_wm["bins"],
            "n_cv": any_wm["n_cv"],
            "per_wiremod": per_wiremod_packs,
            "combined_cov_frac": cov_frac_total,
            "combined_cov": cov_total,
        }

        vc = var_lookup.get((vsn, tgt))
        xlab = vc.var_labels[0] if vc is not None else vsn
        plot_combined_detector_frac_unc(
            any_wm["bins"],
            xlab,
            vsn,
            cov_frac_total,
            per_wiremod_packs,
            plot_dir,
            None,
        )

    wiremod_tags_ordered = list(aggs.keys())
    det_by_wm = build_detector_by_wiremod_block(detector_dict, wiremod_tags_ordered)
    if det_by_wm:
        detector_dict["detector_by_wiremod"] = det_by_wm

    # ------------------------------------------------------------------
    # 2b) Same unisim covariance treatment for SELECTION-cut stages (all tgt)
    # ------------------------------------------------------------------
    selection_detector_dict: Dict[str, Dict[str, dict]] = {"detector": {}}
    for tag in aggs.keys():
        selection_detector_dict[f"detector-{tag}"] = {}
    selection_breakdown_dict: Dict[str, dict] = {}

    selection_tasks = [
        (sk, vsn, tgt)
        for sk in stage_keys
        if sk != final_key
        for vsn, tgt in var_targets
        if any((sk, vsn, tgt) in agg["histdict"] for agg in aggs.values())
    ]
    for stage_key, vsn, tgt in tqdm(
        selection_tasks,
        desc="unisim cov + plots (selection)",
        unit="panel",
        disable=not show_prog,
    ):
        cov_key = selection_cov_npz_key(stage_key, tgt, vsn)
        plot_slug = f"sel__{stage_key}__{tgt}"

        n_bin = None
        cov_frac_total = None
        per_wiremod_packs: Dict[str, dict] = {}

        for tag, agg in aggs.items():
            key = (stage_key, vsn, tgt)
            if key not in agg["histdict"]:
                continue
            hd = agg["histdict"][key]
            if n_bin is None:
                n_bin = hd.hist.shape[1]
                cov_frac_total = np.zeros((n_bin, n_bin))
            if hd.hist.shape[1] != n_bin:
                print(f"[detvar-agg] WARN: bin mismatch {cov_key} ({tag}); skipping")
                continue

            cov_frac_wm, cov_wm, per_calo, n_cv, bins_wm = (
                unisim_wiremod_cov_bundle_from_histacc(hd)
            )

            per_wiremod_packs[tag] = {
                "per_calo": per_calo,
                "cov_frac": cov_frac_wm,
                "cov": cov_wm,
                "n_cv": n_cv,
                "bins": bins_wm,
            }
            cov_frac_total += cov_frac_wm

            selection_detector_dict[f"detector-{tag}"][cov_key] = {
                "cov_frac": cov_frac_wm,
                "cov": cov_wm,
                "corr": _safe_corr_from_cov(cov_wm),
                "n_cv": n_cv,
                "per_calo": {c: per_calo[c] for c in per_calo},
                "bins": bins_wm,
                "stage_key": stage_key,
                "tgt": tgt,
                "var_save_name": vsn,
            }

            vc = var_lookup.get((vsn, tgt))
            xlab = vc.var_labels[0] if vc is not None else vsn
            plot_frac_unc_cov_wiremod(
                tag,
                bins_wm,
                xlab,
                vsn,
                per_calo,
                cov_frac_wm,
                plot_dir,
                plot_slug,
                args.show_fig,
            )

        if cov_frac_total is None:
            continue
        any_wm = next(iter(per_wiremod_packs.values()))
        cov_total = cov_from_fraccov(cov_frac_total, any_wm["n_cv"])
        selection_detector_dict["detector"][cov_key] = {
            "cov_frac": cov_frac_total,
            "cov": cov_total,
            "corr": _safe_corr_from_cov(cov_total),
            "n_cv": any_wm["n_cv"],
            "bins": any_wm["bins"],
            "stage_key": stage_key,
            "tgt": tgt,
            "var_save_name": vsn,
        }
        selection_breakdown_dict[cov_key] = {
            "stage_key": stage_key,
            "tgt": tgt,
            "var_save_name": vsn,
            "bins": any_wm["bins"],
            "n_cv": any_wm["n_cv"],
            "per_wiremod": per_wiremod_packs,
            "combined_cov_frac": cov_frac_total,
            "combined_cov": cov_total,
        }

        vc = var_lookup.get((vsn, tgt))
        xlab = vc.var_labels[0] if vc is not None else vsn
        plot_combined_detector_frac_unc(
            any_wm["bins"],
            xlab,
            vsn,
            cov_frac_total,
            per_wiremod_packs,
            plot_dir,
            plot_slug,
        )

    sel_by_wm = build_detector_by_wiremod_block(
        selection_detector_dict, wiremod_tags_ordered
    )
    if sel_by_wm:
        selection_detector_dict["detector_by_wiremod"] = sel_by_wm

    # ------------------------------------------------------------------
    # 3) Save outputs
    # ------------------------------------------------------------------
    # (a) npz that plugs into utils.get_syst_unc:
    #     dict(detector_syst)['detector'].item()[var_save_name]['cov_frac']
    #     Per-WireMod breakdown (Flux/G4 knob-style nesting):
    #     dict(detector_syst)['detector_by_wiremod'].item()[var_save_name][tag].
    syst_save_path = path.join(args.out_dir, args.syst_save_name)
    np.savez(syst_save_path, **detector_dict)
    print(f"[detvar-agg] wrote {syst_save_path}  keys={list(detector_dict.keys())}")
    print(f"[detvar-agg]   variables in detector key: "
          f"{sorted(detector_dict['detector'].keys())}")
    if "detector_by_wiremod" in detector_dict:
        print(
            "[detvar-agg]   detector_by_wiremod: %d variable(s) with per-WireMod breakdown"
            % (len(detector_dict["detector_by_wiremod"]),)
        )

    selection_npz_path = path.join(args.out_dir, args.selection_syst_save_name)
    np.savez(selection_npz_path, **selection_detector_dict)
    print(f"[detvar-agg] wrote {selection_npz_path}")
    print(f"[detvar-agg]   selection detector keys: "
          f"{sorted(selection_detector_dict['detector'].keys())}")
    if "detector_by_wiremod" in selection_detector_dict:
        print(
            "[detvar-agg]   selection detector_by_wiremod: %d composite key(s)"
            % (len(selection_detector_dict["detector_by_wiremod"]),)
        )

    # (b) full per-WireMod + per-calo breakdown pickle (richer than the npz)
    breakdown_path = path.join(args.out_dir, "detector_unisim_breakdown.pkl")
    with open(breakdown_path, "wb") as f:
        pickle.dump({
            "wiremod_tags": list(aggs.keys()),
            "calo_params": list(CALO_PARAMS),
            "universes":  list(UNIVERSES),
            "stages":     stage_keys,
            "stage_labels": stage_labels,
            "final_stage": final_key,
            "breakdown":   breakdown_dict,
            "selection_breakdown": selection_breakdown_dict,
            "exposure": {tag: {"chunk_pot": agg["chunk_pot"],
                                "chunk_genevts": agg["chunk_genevts"],
                                "n_files": agg["n_files"]}
                          for tag, agg in aggs.items()},
        }, f)
    print(f"[detvar-agg] wrote {breakdown_path}")

    # (c) human-readable text summary of CV / unisim event counts and frac unc
    summary_lines = []
    summary_lines.append("# Detector unisim summary\n")
    summary_lines.append(f"# Final stage: {final_key}\n")
    for vsn in sorted(detector_dict["detector"].keys()):
        e = detector_dict["detector"][vsn]
        bins = e["bins"]
        n_cv = e["n_cv"]
        diag = np.sqrt(np.diag(e["cov_frac"]))
        summary_lines.append(f"\n## {vsn}\n")
        summary_lines.append(f"   bins: {np.array2string(bins, precision=3, max_line_width=140)}\n")
        summary_lines.append(f"   n_cv: {np.array2string(n_cv, precision=3, max_line_width=140)}\n")
        summary_lines.append(f"   frac_unc(detector total): "
                             f"{np.array2string(diag, precision=4, max_line_width=140)}\n")
        for tag, agg in aggs.items():
            wm_entry = detector_dict.get(f"detector-{tag}", {}).get(vsn)
            if wm_entry is None:
                continue
            wm_diag = np.sqrt(np.diag(wm_entry["cov_frac"]))
            summary_lines.append(f"   - {tag}: frac_unc = "
                                 f"{np.array2string(wm_diag, precision=4, max_line_width=140)}\n")
    summary_path = path.join(args.out_dir, "detector_unisim_summary.txt")
    with open(summary_path, "w") as f:
        f.writelines(summary_lines)
    print(f"[detvar-agg] wrote {summary_path}")

    sel_summary_lines = []
    sel_summary_lines.append("# Detector unisim summary — selection-cut variables\n")
    sel_summary_lines.append(f"# Composite keys: stage__tgt__var_save_name\n")
    for cov_key in sorted(selection_detector_dict["detector"].keys()):
        e = selection_detector_dict["detector"][cov_key]
        bins_k = e["bins"]
        n_cv_k = e["n_cv"]
        diag_k = np.sqrt(np.diag(e["cov_frac"]))
        sel_summary_lines.append(f"\n## {cov_key}\n")
        sel_summary_lines.append(
            f"   bins: {np.array2string(bins_k, precision=3, max_line_width=140)}\n")
        sel_summary_lines.append(
            f"   n_cv: {np.array2string(n_cv_k, precision=3, max_line_width=140)}\n")
        sel_summary_lines.append(
            f"   frac_unc(detector total): "
            f"{np.array2string(diag_k, precision=4, max_line_width=140)}\n")
        for tag in aggs.keys():
            wm_entry = selection_detector_dict.get(f"detector-{tag}", {}).get(cov_key)
            if wm_entry is None:
                continue
            wm_diag = np.sqrt(np.diag(wm_entry["cov_frac"]))
            sel_summary_lines.append(
                f"   - {tag}: frac_unc = "
                f"{np.array2string(wm_diag, precision=4, max_line_width=140)}\n")
    selection_summary_path = path.join(args.out_dir, "detector_unisim_selection_summary.txt")
    with open(selection_summary_path, "w") as f:
        f.writelines(sel_summary_lines)
    print(f"[detvar-agg] wrote {selection_summary_path}")

    print(f"[detvar-agg] DONE -> plots in {plot_dir}")


def _safe_corr_from_cov(cov: np.ndarray) -> np.ndarray:
    """Element-wise correlation matrix, with safe handling of zero-diag bins."""
    cov = np.asarray(cov, dtype=float)
    n = cov.shape[0]
    out = np.zeros_like(cov)
    diag = np.diag(cov)
    valid = diag > 0
    if not np.any(valid):
        return out
    sd = np.sqrt(diag)
    for i in range(n):
        for j in range(n):
            denom = sd[i] * sd[j]
            if denom > 0:
                out[i, j] = cov[i, j] / denom
    return out


if __name__ == "__main__":
    main()
