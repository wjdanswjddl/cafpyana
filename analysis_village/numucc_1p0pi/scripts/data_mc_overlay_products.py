#!/usr/bin/env python3
"""Data–MC comparison overlays for PRL Product A (sel_all cut-stage) and Product B (sel_mup).

Product B
  - Style / counts cache: same as ``selected_xsec_overlay.py`` / notebook.
  - Data: ``sel_mup`` (not sel_all). After beam-quality cuts, ``n_evt_good`` must be
    **12,804** (used to pin the correct χ²μ / FV campaign among variants).
  - Cosmics: gate-scaled OffBeamLight — track-PDG as separate ``Intime Cosmics``;
    topology/genie/genie_sb folded into one ``Cosmics`` (MC + offbeam)
    (``cosmic_estimate="offbeam"``).
  - Systematics: PRL Product B consumer tree CategorySummary
    (``dataset_locations.PRL_PRODUCT_B_DIR`` = ``productB_sel_mup``
    since 2026-09-29: GENIE ``GENIE_slim_v3`` + MEC→May, DENT = rolling 80% w=3 + Gauss σ=1,
    detector smear26),
    ``syst_kind="rate"`` (cosmics already contamination-scaled ``SelectedRate``).
    Raw (unsmoothed) DENT is ``productB_sel_mup__dent_raw``. Former overlay GENIE
    (``GENIE_slim_both`` / FSI v1×v3) is kept under
    ``data_mc_overlays/productB_sel_mup__FSI_v1v3/`` (deprecated).

Product A
  - Prefer live-PRL batched counts under ``event_selection-batched-live-PRL``;
    aggregate + save ``merged_histdata.pkl`` (May ``…-20260525`` is outdated /
    pre-DQ and is not used).
  - Systematics: PRL ``productA_sel_all`` combined on the fly (no CategorySummary).
    Exclude MCstat; GENIE ``genie_rate``; **nominal Detector** is nested
    WireMod+DENT (geometry vs matched CV; calo/efield on WM cv; YZ/XTXW
    distinct) consumed as uncorrelated ``diag(u**2)`` by default
    (``--product-a-detector-mode diag``); ``unisim_full`` keeps on-disk
    unisim correlations (test). Cosmics nominal is coherent offbeam→intime
    unisim (rank-1), × contamination. Uncorrelated ``|Δ|`` / ``diag(u**2)``
    cosmics is a test via ``--product-a-cosmics-mode absdiff_diag``.
    Overlay χ² uses the full absolute covariance (not diagonal-only).
  - Track-PDG: keep MC Other; stack offbeam as ``Intime Cosmics``.
    Topology/genie/genie_sb: fold offbeam into MC cosmics as one ``Cosmics``.
  - Flat POT / Ntargets: fully correlated (same as CategorySummary).

Outputs under::

    /exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/data_mc_overlays/
      productB_sel_mup/   nominal overlays (real files;
                          GENIE_slim_v3 + MEC May, DENT rolling 80% w=3 + Gauss σ=1, smear26 Detector)
      productB_sel_mup__dent_raw/   deprecated raw DENT (symlink -> productB_sel_mup__genie_FSIv3_MEC_May)
      productB_sel_mup__FSI_v1v3/   deprecated former GENIE (GENIE_slim_both / FSI v1×v3)
      productA_sel_all/   nominal overlays (nested Detector + diag; cosmics unisim_rank1)
      productA_sel_all__cosmics_absdiff_diag/          test: uncorrelated cosmics
      productA_sel_all__detector_legacy_3knob/         test: former 3-knob max-envelope
      productA_sel_all__detector_unisim3/              test: 3-knob full unisim cov
      productA_sel_all__detector_wiremod10/            test: flat 10 WireMod + DENT
      productA_sel_all__detector_wiremod_nested/       test: nested + Product A unisim_full
      productA_sel_all__detector_wiremod_nested_diag/  test alias of nominal Product A
      productB_sel_mup__detector_*/                    matching Product B tests

Python::

    /exp/sbnd/app/users/munjung/xsec/freeze/cafpyana/envs/venv_py310_cafpyana/bin/python \\
      analysis_village/numucc_1p0pi/scripts/data_mc_overlay_products.py [--product A|B|both]
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import pickle
import sys
import warnings
from datetime import datetime, timezone
from os import makedirs, path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

_REPO = "/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana"
_SCRIPTS = path.join(_REPO, "analysis_village/numucc_1p0pi/scripts")
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from analysis_village.numucc_1p0pi.beam_quality import (  # noqa: E402
    EXPECTED_N_EVT_GOOD,
)
from analysis_village.numucc_1p0pi.constants import MC_POT_FIX  # noqa: E402
from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.categories import topology_labels  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cosmics_common import (  # noqa: E402
    product_a_cosmic_template_cov_frac,
    product_a_cosmic_unisim_cov_frac,
    scale_cov_frac_by_contamination,
)
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    CAT_COSMICS,
    CAT_DETECTOR,
    CAT_FLUX,
    CAT_G4,
    CAT_GENIE_RATE,
    CAT_MCSTAT,
    CAT_NTARGETS,
    CAT_POT,
    assemble_total_rate_cov_frac,
    load_category_syst_summary,
    rebase_fraccov_signal_to_total,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import SYST_DISK_ENV  # noqa: E402
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_syst_unc,
    strip_pot_from_ylabel,
    format_pot_corner_text,
)
from analysis_village.numucc_1p0pi.selection_framework import OverlayHistData  # noqa: E402
from analysis_village.numucc_1p0pi.selected_xsec_overlay_hist import (  # noqa: E402
    build_overlay_histdata_map,
    histdata_pkl_path,
    load_overlay_counts,
    plot_overlay_counts_map,
    save_overlay_counts,
)

import importlib.util as _ilu  # noqa: E402

_sxo_spec = _ilu.spec_from_file_location(
    "selected_xsec_overlay",
    path.join(_SCRIPTS, "selected_xsec_overlay.py"),
)
assert _sxo_spec is not None and _sxo_spec.loader is not None
sxo = _ilu.module_from_spec(_sxo_spec)
_sxo_spec.loader.exec_module(sxo)

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
OUTPUT_ROOT = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/data_mc_overlays"
)
PRODUCT_B_OUT = path.join(OUTPUT_ROOT, "productB_sel_mup")
PRODUCT_A_OUT = path.join(OUTPUT_ROOT, "productA_sel_all")

DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

# Exact BQ match (n_evt_good == 12804): FV-fix sel_mup, not chi2fix-real (12802)
# or the cut-campaign χ²μ variants (11k–13k).
PRODUCT_B_DATA_DIR = path.join(DFS_ROOT, "2026_09_01_064250__sel_mup-data-1e20-fvfix")
PRODUCT_B_DATA_FN = "sel_mup-data-1e20-fvfix"
PRODUCT_B_MC_DIR = path.join(DFS_ROOT, "2026_09_01_063924__sel_mup-mc-fvfix")
PRODUCT_B_MC_FN = "sel_mup-mc-fvfix"
PRODUCT_B_OFFBEAM_DIR = path.join(DFS_ROOT, "2026_09_01_140130__sel_mup-data-OffBeamLight")
PRODUCT_B_OFFBEAM_FN = "sel_mup-data-OffBeamLight"
EXPECTED_DATA_N = EXPECTED_N_EVT_GOOD
# Fraction of OffBeam gates that coincide with BNB (exclude from scale)
OFFBEAM_COINCIDENT_FRAC = 0.08
COSMIC_ESTIMATE = "offbeam"

PRODUCT_A_BATCHES_CHI2MCS20 = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL-chi2mcs20/batches"
)
PRODUCT_A_BATCHES_CHI2AVG = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL-chi2avg/batches"
)
PRODUCT_A_BATCHES_LEGACY_I2 = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL/batches"
)


def _resolve_product_a_batches() -> str:
    """Prefer chi2mcs20 remake, then chi2avg, else legacy I2 live-PRL."""
    import glob as _glob

    n_mcs = len(_glob.glob(path.join(PRODUCT_A_BATCHES_CHI2MCS20, "*__batch_*.pkl")))
    if n_mcs >= 180:
        return PRODUCT_A_BATCHES_CHI2MCS20
    n_avg = len(_glob.glob(path.join(PRODUCT_A_BATCHES_CHI2AVG, "*__batch_*.pkl")))
    if n_avg >= 180:
        return PRODUCT_A_BATCHES_CHI2AVG
    return PRODUCT_A_BATCHES_LEGACY_I2


PRODUCT_A_BATCHES = PRODUCT_A_BATCHES_CHI2MCS20  # default target; resolved at run time


PRODUCT_A_SYST_ROOT = str(prl_syst_disk_root("A"))
PRODUCT_B_SYST_ROOT = str(prl_syst_disk_root("B"))

COUNTS_REPORT_NAME = "counts_report.npz"
COUNTS_MANIFEST_NAME = "counts_report_manifest.json"

# Topology fill order from get_topo_category(ret_cuts=True): index 0 = cosmic,
# last index = signal (νμ CC 1p0π). Display labels are reversed for stacking.
_TOPO_COSMIC_IDX = 0
_TOPO_SIGNAL_IDX = -1  # last layer


# ===========================================================================
# Shared: reportable bin counts (data + MC signal / background)
# ===========================================================================


def _topo_signal_bkg_from_mc_hist(mc_hist: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split stacked topology MC into signal, ν-background, and MC-cosmic layers.

    Topology fill order (``get_topo_category(ret_cuts=True)``): index 0 = cosmic,
    last index = signal. Returns ``(signal, nu_bkg, mc_cosmic)``.
    """
    mh = np.asarray(mc_hist, dtype=float)
    if mh.ndim != 2 or mh.shape[0] == 0:
        z = np.zeros(mh.shape[-1] if mh.ndim else 0, dtype=float)
        return z, z, z
    sig = mh[_TOPO_SIGNAL_IDX].copy()
    mc_cosmic = mh[_TOPO_COSMIC_IDX].copy()
    nu_bkg = mh.sum(axis=0) - sig - mc_cosmic
    return sig, nu_bkg, mc_cosmic


def _overlay_signal_and_total(
    hd: OverlayHistData,
    *,
    cosmic_estimate: str = "offbeam",
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Per-bin ``(n_signal, n_total)`` for rebasing signal-CV rate fracs.

    Signal is the topology signal layer (last MC stack row). Total matches the
    overlay stack used for ``frac × total_mc``: MC + data-driven cosmics + dirt.
    Returns ``None`` when the hist is not a topology-style stack.
    """
    if hd is None or getattr(hd, "mc_hist", None) is None:
        return None
    mc = np.asarray(hd.mc_hist, dtype=float)
    bt = str(getattr(hd, "breakdown_type", "") or "")
    n_topo = len(topology_labels)
    if mc.ndim != 2 or mc.shape[0] < 2:
        return None
    if bt not in ("topology", "genie", "genie_sb") and mc.shape[0] != n_topo:
        return None
    sig, _nu_bkg, _mc_cosmic = _topo_signal_bkg_from_mc_hist(mc)
    n = int(sig.shape[0])
    prefer_offbeam = bt == "pdg" or str(cosmic_estimate).lower() == "offbeam"
    offbeam = (
        np.asarray(hd.offbeam_hist, dtype=float).ravel()
        if getattr(hd, "has_offbeam", False) and hd.offbeam_hist is not None
        else np.zeros(n, dtype=float)
    )
    intime = (
        np.asarray(hd.intime_hist, dtype=float).ravel()
        if getattr(hd, "has_intime", False) and hd.intime_hist is not None
        else np.zeros(n, dtype=float)
    )
    if prefer_offbeam and np.any(offbeam):
        cosmic_est = offbeam
    elif np.any(intime):
        cosmic_est = intime
    elif np.any(offbeam):
        cosmic_est = offbeam
    else:
        cosmic_est = np.zeros(n, dtype=float)
    dirt = np.zeros(n, dtype=float)
    dirt_cat = getattr(hd, "dirt_cat_hist", None)
    if (
        bt == "pdg"
        and dirt_cat is not None
        and np.any(np.asarray(dirt_cat, dtype=float))
    ):
        dirt = np.asarray(dirt_cat, dtype=float).sum(axis=0)
    elif getattr(hd, "has_dirt", False) and hd.dirt_hist is not None:
        dirt = np.asarray(hd.dirt_hist, dtype=float).ravel()
    n_total = mc.sum(axis=0).astype(float) + cosmic_est + dirt
    return np.asarray(sig, dtype=float), np.asarray(n_total, dtype=float)


def export_counts_report(
    out_dir: str,
    *,
    product: str,
    histdata_items: Sequence[Tuple[str, OverlayHistData]],
    extra_meta: Optional[Mapping[str, Any]] = None,
    cosmic_estimate: str = "offbeam",
) -> str:
    """Write NPZ + JSON manifest of bin-by-bin data / MC signal / MC background.

    Background for unfolding matches the overlay stack: keep MC cosmics and
    *add* the data-driven estimate —
    ``mc_background = nu_bkg + mc_cosmic + offbeam (+ dirt)``.
    """
    makedirs(out_dir, exist_ok=True)
    arrays: Dict[str, np.ndarray] = {}
    rows = []
    prefer_offbeam = str(cosmic_estimate).lower() == "offbeam"
    for slug, hd in histdata_items:
        if not getattr(hd, "has_data", False) and hd.data_hist is None:
            continue
        data = np.asarray(hd.data_hist, dtype=float).ravel()
        bins = np.asarray(hd.bins, dtype=float)
        if getattr(hd, "has_mc", False) and hd.mc_hist is not None:
            sig, nu_bkg, mc_cosmic = _topo_signal_bkg_from_mc_hist(hd.mc_hist)
            mc_tot_raw = np.asarray(hd.mc_hist, dtype=float).sum(axis=0)
        else:
            sig = nu_bkg = mc_cosmic = mc_tot_raw = np.zeros_like(data)

        intime = (
            np.asarray(hd.intime_hist, dtype=float).ravel()
            if getattr(hd, "has_intime", False) and hd.intime_hist is not None
            else np.zeros_like(data)
        )
        offbeam = (
            np.asarray(hd.offbeam_hist, dtype=float).ravel()
            if getattr(hd, "has_offbeam", False) and hd.offbeam_hist is not None
            else np.zeros_like(data)
        )
        dirt = (
            np.asarray(hd.dirt_hist, dtype=float).ravel()
            if getattr(hd, "has_dirt", False) and hd.dirt_hist is not None
            else np.zeros_like(data)
        )

        if prefer_offbeam and np.any(offbeam):
            cosmic_est = offbeam
        elif np.any(intime):
            cosmic_est = intime
        elif np.any(offbeam):
            cosmic_est = offbeam
        else:
            cosmic_est = np.zeros_like(data)

        # Unfold / overlay background: ν non-signal + MC cosmics + data-driven (+ dirt)
        bkg = nu_bkg + mc_cosmic + cosmic_est + dirt
        # Stack total matches overlay (MC kept + Intime Cosmics + dirt)
        mc_tot = sig + bkg
        other = intime + dirt + offbeam  # raw alternate samples (diagnostics)

        key = slug.replace("/", "_")
        arrays[f"{key}__bins"] = bins
        arrays[f"{key}__data"] = data
        arrays[f"{key}__mc_signal"] = sig
        arrays[f"{key}__mc_background"] = bkg
        arrays[f"{key}__mc_nu_background"] = nu_bkg
        arrays[f"{key}__mc_cosmic"] = mc_cosmic
        arrays[f"{key}__cosmic_estimate"] = cosmic_est
        arrays[f"{key}__mc_total"] = mc_tot
        arrays[f"{key}__mc_total_raw"] = mc_tot_raw
        arrays[f"{key}__mc_other"] = other
        arrays[f"{key}__offbeam"] = offbeam
        arrays[f"{key}__intime"] = intime
        rows.append(
            {
                "slug": slug,
                "var_save_name": getattr(hd, "var_save_name", ""),
                "breakdown_type": getattr(hd, "breakdown_type", ""),
                "n_bins": int(len(data)),
                "data_sum": float(np.nansum(data)),
                "mc_signal_sum": float(np.nansum(sig)),
                "mc_background_sum": float(np.nansum(bkg)),
                "mc_nu_background_sum": float(np.nansum(nu_bkg)),
                "cosmic_estimate_sum": float(np.nansum(cosmic_est)),
                "mc_total_sum": float(np.nansum(mc_tot)),
                "mc_other_sum": float(np.nansum(other)),
                "offbeam_sum": float(np.nansum(offbeam)),
                "topology_labels": list(topology_labels),
                "signal_layer_index": len(topology_labels) - 1,
                "cosmic_layer_index": 0,
                "cosmic_estimate": cosmic_estimate,
            }
        )

    npz_path = path.join(out_dir, COUNTS_REPORT_NAME)
    np.savez_compressed(npz_path, **arrays)
    manifest = {
        "schema": "data_mc_overlay_counts_v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "product": product,
        "npz_path": npz_path,
        "n_variables": len(rows),
        "variables": rows,
        "cosmic_estimate": cosmic_estimate,
        "extra": dict(extra_meta or {}),
    }
    man_path = path.join(out_dir, COUNTS_MANIFEST_NAME)
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[counts] wrote {npz_path} ({len(rows)} vars) + {man_path}", flush=True)
    return npz_path


# ===========================================================================
# Product B — sel_mup overlays (selected_xsec_overlay style)
# ===========================================================================


def _product_b_syst_cov(
    var_config,
    *,
    syst_disk_root: Optional[str] = None,
    hd: Optional[OverlayHistData] = None,
    cosmic_estimate: str = COSMIC_ESTIMATE,
):
    """Product B overlay rate fractional covariance.

    Prefers reassembling CategorySummary ``categories`` with an overlay-time
    rebase of flux/g4/mcstat/genie_rate (signal CV → total-selected CV) so old
    NPZs with unrebased ``total_rate`` still plot correctly. Detector / cosmics
    SelectedRate / pot / ntargets stay unrebased. Falls back to disk
    ``get_syst_unc`` (with the same rebase) when CategorySummary bins disagree
    with ``var_config``.
    """
    root = syst_disk_root or PRODUCT_B_SYST_ROOT
    centers = getattr(var_config, "bin_centers", None)
    n = int(len(centers)) if centers is not None else 0
    n_sig = n_tot = None
    counts = _overlay_signal_and_total(hd, cosmic_estimate=cosmic_estimate) if hd is not None else None
    if counts is not None:
        n_sig, n_tot = counts
        if n > 0 and (n_sig.shape[0] != n or n_tot.shape[0] != n):
            print(
                f"  [B syst] overlay counts shape mismatch for "
                f"{var_config.var_save_name}: sig={n_sig.shape} tot={n_tot.shape} "
                f"vs n={n}; skip rebase",
                flush=True,
            )
            n_sig = n_tot = None

    # --- Prefer CategorySummary categories (resilient to unrebased total_rate) ---
    try:
        summary = load_category_syst_summary(syst_disk_root=root)
        vsn = (
            getattr(var_config, "category_syst_var_save_name", None)
            or var_config.var_save_name
        )
        pack = summary["by_var"].get(vsn)
        if pack is not None:
            cats = pack.get("categories") or {}
            cat_covs = {
                k: np.asarray(cats[k]["cov_frac"], dtype=float)
                for k in (
                    CAT_FLUX,
                    CAT_G4,
                    CAT_MCSTAT,
                    CAT_DETECTOR,
                    CAT_COSMICS,
                    CAT_GENIE_RATE,
                    CAT_POT,
                    CAT_NTARGETS,
                )
                if k in cats and isinstance(cats[k], dict) and "cov_frac" in cats[k]
            }
            if cat_covs:
                rebase = n_sig is not None and n_tot is not None
                if rebase:
                    cov = assemble_total_rate_cov_frac(
                        cat_covs,
                        n_signal=n_sig,
                        n_total=n_tot,
                        rebase=True,
                    )
                else:
                    # Prefer on-disk total_rate (rebased by
                    # rebase_category_summary_total_rate.py) over an unrebased
                    # category sum when overlay counts are unavailable.
                    tr = pack.get("total_rate") or {}
                    cov = tr.get("cov_frac")
                    if cov is None:
                        print(
                            f"  [B syst] {vsn}: no overlay n_sig/n_tot and no "
                            "total_rate; using unrebased category sum",
                            flush=True,
                        )
                        cov = assemble_total_rate_cov_frac(
                            cat_covs,
                            n_signal=None,
                            n_total=None,
                            rebase=False,
                        )
                    else:
                        print(
                            f"  [B syst] {vsn}: no overlay n_sig/n_tot; "
                            "using CategorySummary total_rate",
                            flush=True,
                        )
                if cov is not None and np.any(np.isfinite(cov)):
                    cov = np.asarray(cov, dtype=float)
                    if n > 0 and cov.shape == (n, n):
                        return _sanitize_cov_frac(cov)
                    print(
                        f"  [B syst] CategorySummary shape {cov.shape} != ({n},{n}) for "
                        f"{vsn}; assembling from disk (skip mismatches)",
                        flush=True,
                    )
    except Exception as ex:
        print(f"  [B syst] CategorySummary miss {var_config.var_save_name}: {ex}", flush=True)

    # --- Disk fallback: rebase signal-CV components, leave flats as-is ---
    try:
        _unc, cov_sig = get_syst_unc(
            var_config,
            syst_disk_root=root,
            syst_components=("flux", "g4", "mcstat", "genie"),
            genie_cov_frac_key="genie_rate",
            skip_missing_vars=True,
            plot=False,
        )
        cov_sig = np.asarray(cov_sig if cov_sig is not None else 0.0, dtype=float)
        if cov_sig.ndim != 2:
            cov_sig = np.zeros((n, n), dtype=float)
        if n_sig is not None and n_tot is not None and cov_sig.shape == (n, n):
            cov_sig = rebase_fraccov_signal_to_total(cov_sig, n_sig, n_tot)

        _unc, cov_rest = get_syst_unc(
            var_config,
            syst_disk_root=root,
            syst_components=("detector", "cosmics", "pot", "ntargets"),
            skip_missing_vars=True,
            plot=False,
        )
        cov_rest = np.asarray(cov_rest if cov_rest is not None else 0.0, dtype=float)
        if cov_rest.ndim != 2:
            cov_rest = np.zeros((n, n), dtype=float)

        cov = _sanitize_cov_frac(cov_sig + cov_rest)
        if n > 0 and cov.shape != (n, n):
            print(
                f"  [B syst] disk cov shape {cov.shape} != ({n},{n}) for "
                f"{var_config.var_save_name}; skip",
                flush=True,
            )
            return None
        if not np.any(np.isfinite(cov)) or not np.any(cov):
            return None
        return cov
    except Exception as ex:
        print(f"  [B syst] skip {var_config.var_save_name}: {ex}", flush=True)
        return None


def _scale_overlay_mc(hd, factor: float) -> None:
    """Multiply neutrino-MC histograms by ``factor`` (weights and weight^2)."""
    if hd is None or getattr(hd, "mc_hist", None) is None:
        return
    hd.mc_hist = np.asarray(hd.mc_hist, dtype=float) * factor
    if getattr(hd, "mc_err2", None) is not None:
        hd.mc_err2 = np.asarray(hd.mc_err2, dtype=float) * (factor ** 2)
    univ = getattr(hd, "mc_univ_hist", None)
    if univ:
        for key, arr in list(univ.items()):
            univ[key] = np.asarray(arr, dtype=float) * factor


def _scale_merged_mc(merged: dict, factor: float) -> None:
    for hd in merged.get("histdata", {}).values():
        _scale_overlay_mc(hd, factor)
    for by_bt in merged.get("bar", {}).values():
        for bb in by_bt.values():
            bb.mc_counts = np.asarray(bb.mc_counts, dtype=float) * factor


def _refuse_in_place_mcpotfix(out: str, production: str) -> None:
    if path.abspath(out) == path.abspath(production):
        raise RuntimeError(
            f"MC POT fix must be written to a parallel directory, not {production}"
        )


def _scale_offbeam_to_data(
    offbeam_evt: pd.DataFrame,
    offbeam_hdr: pd.DataFrame,
    data_hdr: pd.DataFrame,
    *,
    f_coincident: float = OFFBEAM_COINCIDENT_FRAC,
) -> float:
    """Attach gate-scaled ``pot_weight``; return the scale factor."""
    data_gates = float(data_hdr["nbnbinfo"].sum())
    ob_hdr = offbeam_hdr
    if "first_in_subrun" in ob_hdr.columns:
        ob_hdr = ob_hdr[ob_hdr["first_in_subrun"] == 1]
    offbeam_gates = float(ob_hdr["noffbeambnb"].sum())
    if offbeam_gates <= 0:
        raise RuntimeError("OffBeam gates sum to 0; check hdr.noffbeambnb")
    scale = (1.0 - float(f_coincident)) * data_gates / offbeam_gates
    offbeam_evt["pot_weight"] = scale * np.ones(len(offbeam_evt))
    offbeam_evt["gates_weight"] = offbeam_evt["pot_weight"].copy()
    print(
        f"  offbeam scale: data_gates={data_gates:.3e}  "
        f"offbeam_gates={offbeam_gates:.3e}  f={f_coincident}  scale={scale:.4f}",
        flush=True,
    )
    return scale


def run_product_b(
    *,
    force_rebuild: bool = False,
    out_dir: Optional[str] = None,
    syst_disk_root: Optional[str] = None,
    mc_pot_fix: Optional[float] = None,
) -> None:
    out = out_dir or PRODUCT_B_OUT
    if mc_pot_fix:
        _refuse_in_place_mcpotfix(out, PRODUCT_B_OUT)
    syst_root = syst_disk_root or PRODUCT_B_SYST_ROOT
    print("\n" + "=" * 72 + "\nProduct B (sel_mup)\n" + "=" * 72, flush=True)
    print(f"  out_dir={out}", flush=True)
    print(f"  syst_disk_root={syst_root}", flush=True)
    os.environ[SYST_DISK_ENV] = syst_root
    makedirs(out, exist_ok=True)

    plot_set = {
        "tag": "productB_sel_mup",
        "output_dir": out,
        "mc_dir": PRODUCT_B_MC_DIR,
        "mc_filename_str": PRODUCT_B_MC_FN,
        "data_dir": PRODUCT_B_DATA_DIR,
        "data_filename_str": PRODUCT_B_DATA_FN,
        "offbeam_dir": PRODUCT_B_OFFBEAM_DIR,
        "offbeam_filename_str": PRODUCT_B_OFFBEAM_FN,
        "cosmic_estimate": COSMIC_ESTIMATE,
    }

    # Align sxo module knobs with this driver
    sxo.APPLY_BEAM_QUALITY = True
    sxo.LOAD_SYST = True
    sxo.SYST_DISK_ROOT = syst_root
    sxo.FORCE_REBUILD_COUNTS = force_rebuild
    sxo.SAVE_FIG = True
    sxo.PLOT = False
    sxo.APPROVAL = ""
    sxo.TEXTCHI2 = True
    # Node often sits ~70% used; allow Product B DF fill without aborting.
    sxo.MEMORY_LIMIT_FRAC = max(float(getattr(sxo, "MEMORY_LIMIT_FRAC", 0.6)), 0.90)

    # Prefer *this* out_dir counts when present (vertex×N / test trees). Fall back
    # to production Product B cache only when out_dir has no pickle yet.
    payload = None
    if not force_rebuild:
        payload = load_overlay_counts(out)
        if payload is None and out != PRODUCT_B_OUT:
            payload = load_overlay_counts(PRODUCT_B_OUT)
            cache_pkl = histdata_pkl_path(PRODUCT_B_OUT)
        else:
            cache_pkl = histdata_pkl_path(out)
    else:
        cache_pkl = histdata_pkl_path(out)
    bq_n = None
    if payload is not None:
        src_pkl = cache_pkl if path.isfile(cache_pkl) else histdata_pkl_path(out)
        print(f"  replot from counts: {src_pkl}", flush=True)
        raw_pot_label = payload.get("pot_label") or "Events / Bin"
        pot_text = format_pot_corner_text(raw_pot_label)
        pot_label = strip_pot_from_ylabel(raw_pot_label) or "Events / Bin"
        histdata_map = payload["histdata"]
        bq_n = (payload.get("plot_set") or {}).get("n_evt_good")
        if out != PRODUCT_B_OUT and path.isfile(cache_pkl):
            with open(path.join(out, "overlay_histdata_source.json"), "w") as f:
                json.dump(
                    {
                        "overlay_histdata_pkl": cache_pkl,
                        "syst_disk_root": syst_root,
                        "note": "counts reused; overlays regenerated only",
                    },
                    f,
                    indent=2,
                )
    else:
        print("  filling counts from dataframes...", flush=True)
        mc_evt, mc_hdr, n_mc_loaded, n_mc_total = sxo.load_mc_sample(
            plot_set["mc_dir"], plot_set["mc_filename_str"]
        )
        if n_mc_loaded < n_mc_total:
            print(
                f"  MC subsample: {n_mc_loaded}/{n_mc_total} files",
                flush=True,
            )
        data_evt, data_hdr = sxo.load_data_sample(
            plot_set["data_dir"], plot_set["data_filename_str"]
        )
        bq_n = int(len(data_evt))
        if bq_n != EXPECTED_DATA_N:
            raise RuntimeError(
                f"Product B data after beam quality is {bq_n}, expected "
                f"{EXPECTED_DATA_N}. Check data dir {PRODUCT_B_DATA_DIR}."
            )
        print(f"  beam-quality check OK: n_evt_good={bq_n}", flush=True)
        pot_label_raw = sxo.setup_pot_weights(mc_evt, mc_hdr, data_evt, data_hdr)
        plot_set["mc_pot_fix_factor"] = float(MC_POT_FIX)
        plot_set["mc_pot_fix_note"] = (
            f"Recorded neutrino-MC POT is high by {MC_POT_FIX}; "
            f"neutrino-MC pot_weight includes ×{MC_POT_FIX}."
        )
        pot_text = format_pot_corner_text(pot_label_raw)
        pot_label = strip_pot_from_ylabel(pot_label_raw) or "Events / Bin"

        print("  loading OffBeamLight for cosmic estimate...", flush=True)
        offbeam_evt, offbeam_hdr = sxo.load_sample(
            PRODUCT_B_OFFBEAM_DIR,
            PRODUCT_B_OFFBEAM_FN,
            label="offbeam",
        )
        _scale_offbeam_to_data(offbeam_evt, offbeam_hdr, data_hdr)

        # Derive φ (and related kinematics) if CAFs omitted precomputed phi.
        from analysis_village.numucc_1p0pi.evt_derived_kinematics import (
            ensure_derived_trk_kinematics_cols,
        )

        mc_evt = ensure_derived_trk_kinematics_cols(mc_evt)
        data_evt = ensure_derived_trk_kinematics_cols(data_evt)
        offbeam_evt = ensure_derived_trk_kinematics_cols(offbeam_evt)

        histdata_map = build_overlay_histdata_map(
            sxo.VAR_CONFIGS,
            sxo.BREAKDOWN_TYPES,
            mc_df=mc_evt,
            data_df=data_evt,
            offbeam_df=offbeam_evt,
        )
        plot_set["n_evt_good"] = bq_n
        pkl = save_overlay_counts(
            out,
            histdata_map,
            # Keep POT in pickle so corner text survives ylabel stripping on replot.
            pot_label=pot_label_raw,
            plot_set=plot_set,
            var_save_names=[vc.var_save_name for vc in sxo.VAR_CONFIGS],
            breakdown_types=sxo.BREAKDOWN_TYPES,
        )
        print(f"  wrote counts -> {pkl}", flush=True)
        del mc_evt, mc_hdr, data_evt, data_hdr, offbeam_evt, offbeam_hdr
        gc.collect()
        payload = {
            "pot_label": pot_label_raw,
            "plot_set": dict(plot_set),
            "histdata": histdata_map,
        }

    if mc_pot_fix:
        saved_ps = dict((payload or {}).get("plot_set") or {})
        already = saved_ps.get("mc_pot_fix_factor")
        if already:
            print(
                f"  MC POT fix already applied (factor={already}); not scaling again",
                flush=True,
            )
        else:
            sample_key = next(
                (k for k, h in histdata_map.items() if getattr(h, "has_mc", False)),
                None,
            )
            before = (
                float(np.sum(histdata_map[sample_key].mc_hist))
                if sample_key is not None
                else float("nan")
            )
            for hd in histdata_map.values():
                _scale_overlay_mc(hd, float(mc_pot_fix))
            after = (
                float(np.sum(histdata_map[sample_key].mc_hist))
                if sample_key is not None
                else float("nan")
            )
            print(
                f"  MC POT fix ×{mc_pot_fix}: example {sample_key} "
                f"sum(mc) {before:.6g} -> {after:.6g}",
                flush=True,
            )
            saved_ps["mc_pot_fix_factor"] = float(mc_pot_fix)
            saved_ps["mc_pot_fix_note"] = (
                f"Recorded neutrino-MC POT is high by {MC_POT_FIX}; "
                f"neutrino-MC histograms multiplied by {MC_POT_FIX}. "
                "Dirt and cosmics unchanged."
            )
            if bq_n is not None:
                saved_ps["n_evt_good"] = bq_n
            raw_for_save = (payload or {}).get("pot_label") or pot_label
            pkl = save_overlay_counts(
                out,
                histdata_map,
                pot_label=raw_for_save,
                plot_set=saved_ps,
                var_save_names=[vc.var_save_name for vc in sxo.VAR_CONFIGS],
                breakdown_types=sxo.BREAKDOWN_TYPES,
            )
            print(f"  wrote MC-POT-fixed counts -> {pkl}", flush=True)

    if bq_n is not None and int(bq_n) != EXPECTED_DATA_N:
        print(
            f"  WARN: cached counts report n_evt_good={bq_n} "
            f"(expected {EXPECTED_DATA_N}); rebuild with --force-rebuild",
            flush=True,
        )

    def _syst(vc):
        hd = histdata_map.get((vc.var_save_name, "topology"))
        if hd is None:
            for bt in ("genie", "genie_sb", "pdg"):
                hd = histdata_map.get((vc.var_save_name, bt))
                if hd is not None:
                    break
        return _product_b_syst_cov(
            vc,
            syst_disk_root=syst_root,
            hd=hd,
            cosmic_estimate=COSMIC_ESTIMATE,
        )

    plot_overlay_counts_map(
        histdata_map,
        sxo.VAR_CONFIGS,
        sxo.BREAKDOWN_TYPES,
        pot_label=pot_label,
        out_dir=out,
        get_syst=_syst,
        ax_ylim_ratio=sxo.AX_YLIM_RATIO,
        ratio=sxo.RATIO,
        textloc=sxo.TEXTLOC,
        approval="",
        pot_text=pot_text,
        save_fig=True,
        plot=False,
        textchi2=True,
        cosmic_estimate=COSMIC_ESTIMATE,
    )

    # Prefer topology breakdown for the reportable signal/bkg split
    items = []
    for vc in sxo.VAR_CONFIGS:
        key = (vc.var_save_name, "topology")
        if key in histdata_map:
            items.append((vc.var_save_name, histdata_map[key]))
    export_counts_report(
        out,
        product="B",
        histdata_items=items,
        cosmic_estimate=COSMIC_ESTIMATE,
        extra_meta={
            "data_dir": PRODUCT_B_DATA_DIR,
            "mc_dir": PRODUCT_B_MC_DIR,
            "offbeam_dir": PRODUCT_B_OFFBEAM_DIR,
            "syst_disk_root": syst_root,
            "syst_kind": "rate",
            "n_evt_good": bq_n,
            "expected_n_evt_good": EXPECTED_DATA_N,
            "cosmic_estimate": COSMIC_ESTIMATE,
            "offbeam_coincident_frac": OFFBEAM_COINCIDENT_FRAC,
            "mc_pot_fix_factor": (
                (payload or {}).get("plot_set", {}).get("mc_pot_fix_factor")
                or mc_pot_fix
            ),
        },
    )


# ===========================================================================
# Product A — sel_all cut-stage from live-PRL batches
# ===========================================================================


def _product_a_syst_slug(
    var_save_name: str,
    stage_key: str,
    *,
    name_suffix: Optional[str] = None,
) -> Optional[str]:
    """Map pipeline PlotSpec names onto Product A PRL syst keys.

    χ² overlays use plane-averaged columns (``VariableConfig.chi2_mu/proton``)
    but syst packs are stored under ``chi2_avg_*__at_<stage>``. ``not_mu``
    PlotSpecs must use the dedicated subset packs — never the all-track avg.
    """
    # Suffix variants before the all-track direct map (mcs_range_diff is in
    # ``direct``; len50 must not fall through to the all-track pack).
    if name_suffix == "len50" and var_save_name == "mcs_range_diff":
        return "mcs_range_diff_len50"

    # Direct overlap (cut-flow observables)
    direct = {
        "nu_score",
        "n_trks",
        "track_score",
        "trk_len",
        "vtx_dist",
        "mcs_range_diff",
    }
    if var_save_name in direct:
        return var_save_name
    # Stage-tagged χ² averages at muon-ID stages
    if stage_key in ("2prong-vtxdist", "2prong-muX", "2prong-mup"):
        if name_suffix == "not_mu":
            if var_save_name == "chi2_mu":
                return f"chi2_avg_mu_not_mu__at_{stage_key}"
            if var_save_name == "chi2_p":
                return f"chi2_avg_p_not_mu__at_{stage_key}"
            return None
        if name_suffix == "len50":
            # Subset packs only exist at vtxdist (see CHI2_TRACK_SUBSET_SLUGS).
            if stage_key != "2prong-vtxdist":
                return None
            if var_save_name == "chi2_mu":
                return "chi2_avg_mu_len50__at_2prong-vtxdist"
            if var_save_name == "chi2_p":
                return "chi2_avg_p_len50__at_2prong-vtxdist"
            return None
        if name_suffix == "len50_qual":
            if stage_key != "2prong-vtxdist":
                return None
            if var_save_name == "chi2_mu":
                return "chi2_avg_mu_len50_qual__at_2prong-vtxdist"
            if var_save_name == "chi2_p":
                return "chi2_avg_p_len50_qual__at_2prong-vtxdist"
            return None
        if var_save_name in ("chi2_mu", "chi2_avg_mu"):
            return f"chi2_avg_mu__at_{stage_key}"
        if var_save_name in ("chi2_p", "chi2_avg_p"):
            return f"chi2_avg_p__at_{stage_key}"
    return None

def _contam_frac_from_histdata(
    hd: OverlayHistData, *, cosmic_estimate: str = "offbeam"
) -> np.ndarray:
    """Per-bin cosmic contamination for scaling the cosmic-template frac-cov.

    Stack total = full MC + data-driven offbeam/intime (+ dirt).

    Numerator (the cosmic contribution the template uncertainty multiplies):
      - topology / genie / genie_sb: MC cosmics (layer 0) + data-driven sample
      - track-PDG: data-driven ``Intime Cosmics`` only (no MC cosmic layer)

    Fractional template cov is then scaled by ``outer(f, f)``.
    """
    n = len(hd.bins) - 1
    if hd.mc_hist is None:
        return np.zeros(n, dtype=float)

    mc = np.asarray(hd.mc_hist, dtype=float)
    bt = getattr(hd, "breakdown_type", "")
    prefer_offbeam = bt == "pdg" or str(cosmic_estimate).lower() == "offbeam"
    offbeam = (
        np.asarray(hd.offbeam_hist, dtype=float)
        if getattr(hd, "has_offbeam", False) and hd.offbeam_hist is not None
        else np.zeros(n, dtype=float)
    )
    intime = (
        np.asarray(hd.intime_hist, dtype=float)
        if getattr(hd, "has_intime", False) and hd.intime_hist is not None
        else np.zeros(n, dtype=float)
    )
    if prefer_offbeam and np.any(offbeam):
        data_driven = offbeam
    elif np.any(intime):
        data_driven = intime
    elif np.any(offbeam):
        data_driven = offbeam
    else:
        data_driven = np.zeros(n, dtype=float)

    dirt = np.zeros(n, dtype=float)
    dirt_cat = getattr(hd, "dirt_cat_hist", None)
    if (
        bt == "pdg"
        and dirt_cat is not None
        and np.any(np.asarray(dirt_cat, dtype=float))
    ):
        dirt = np.asarray(dirt_cat, dtype=float).sum(axis=0)
    elif getattr(hd, "has_dirt", False) and hd.dirt_hist is not None:
        dirt = np.asarray(hd.dirt_hist, dtype=float)

    total = mc.sum(axis=0) + data_driven + dirt

    # Cosmic contribution the template unc applies to (matches stacked Cosmics).
    if bt in ("topology", "genie", "genie_sb") and mc.shape[0] > _TOPO_COSMIC_IDX:
        cosmic_contrib = np.asarray(mc[_TOPO_COSMIC_IDX], dtype=float) + data_driven
    else:
        cosmic_contrib = data_driven

    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(total > 0, cosmic_contrib / total, 0.0)
    return np.nan_to_num(frac, nan=0.0, posinf=0.0, neginf=0.0)


def _sanitize_cov_frac(cov: np.ndarray, *, max_diag: float = 1e10) -> np.ndarray:
    """Zero rows/cols with non-finite or placeholder-huge diagonal (e.g. 1e24).

    Only strip unphysical placeholders — large-but-finite detector fractional
    uncertainties in sparse bins must remain visible on the overlay band.
    """
    c = np.asarray(cov, dtype=float).copy()
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        return np.zeros_like(c, dtype=float)
    diag = np.diag(c).copy()
    bad = ~np.isfinite(diag) | (diag > max_diag) | (diag < 0)
    if np.any(bad):
        c[bad, :] = 0.0
        c[:, bad] = 0.0
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    return c


class _VarProxy:
    """Minimal stand-in so get_syst_unc can look up a remapped var_save_name."""

    def __init__(self, var_save_name: str, bins: np.ndarray, bin_centers: np.ndarray):
        self.var_save_name = var_save_name
        self.bins = bins
        self.bin_centers = bin_centers
        self.var_labels = ["", var_save_name, ""]


# Disk-backed + flat terms for Product A overlays (except detector + cosmics).
# Nominal Detector (nested WireMod+DENT) is loaded separately (diag by default;
# unisim_full keeps on-disk correlations for tests). Cosmics are rebuilt from
# stored offbeam/intime hists. MCstat excluded; overlay Poisson MC-stat is
# added on the hatch/χ² path. Flux/G4/GENIE are loaded separately and rebased
# from signal CV → total-selected CV before combining with pot/ntargets.
_PRODUCT_A_SIGNAL_CV_COMPONENTS = ("flux", "g4", "genie")
_PRODUCT_A_FLAT_COMPONENTS = ("pot", "ntargets")


def _product_a_detector_cov_frac(
    proxy: _VarProxy,
    *,
    mode: str = "diag",
    syst_disk_root: Optional[str] = None,
) -> np.ndarray:
    """Detector frac-cov for Product A.

    ``mode``:
      - ``diag`` (nominal): keep nested Detector magnitude, drop bin–bin
        correlation (``diag(u**2)``).
      - ``unisim_full`` (test): keep the on-disk unisim / summed-knob correlations.
    """
    n = len(proxy.bin_centers)
    root = syst_disk_root or PRODUCT_A_SYST_ROOT
    mode = str(mode).lower()
    try:
        _unc, det = get_syst_unc(
            proxy,
            syst_disk_root=root,
            syst_components=("detector",),
            skip_missing_vars=True,
            plot=False,
        )
    except Exception as ex:
        print(f"  [A syst] detector skip {proxy.var_save_name}: {ex}", flush=True)
        return np.zeros((n, n), dtype=float)
    det = np.asarray(det, dtype=float)
    if det.shape != (n, n):
        return np.zeros((n, n), dtype=float)
    if mode in ("diag", "absdiff_diag", "diagonal"):
        return _sanitize_cov_frac(np.diag(np.diag(det)))
    if mode in ("unisim_full", "unisim", "full", "rank1"):
        return _sanitize_cov_frac(det)
    raise ValueError(
        f"Unknown Product A detector mode {mode!r}; use diag or unisim_full"
    )

def _product_a_cosmics_cov_frac(
    slug: str,
    hd: OverlayHistData,
    *,
    mode: str = "unisim_rank1",
    syst_disk_root: Optional[str] = None,
) -> Optional[np.ndarray]:
    """Contamination-scaled Product A cosmics fractional covariance.

    ``mode``:
      - ``unisim_rank1`` (nominal): CV=offbeam, univ=intime (coherent rank-1)
      - ``absdiff_diag`` (test): smoothed ``|Δ|/offbeam`` as uncorrelated ``diag(u**2)``
    """
    from analysis_village.numucc_1p0pi.syst_disk_layout import syst_disk_paths

    mode = str(mode).lower()
    root = syst_disk_root or PRODUCT_A_SYST_ROOT
    cpath = syst_disk_paths(root)["cosmics"]
    blob = dict(np.load(cpath, allow_pickle=True))
    if slug not in blob:
        return None
    cell = blob[slug].item()["Cosmics"]
    h_off = np.asarray(cell["univ_offbeam"], dtype=float).ravel()
    h_in = np.asarray(cell["univ_intime"], dtype=float).ravel()
    if mode in ("unisim_rank1", "unisim", "rank1"):
        raw = product_a_cosmic_unisim_cov_frac(h_off, h_in)
    elif mode in ("absdiff_diag", "absdiff", "diag"):
        raw = product_a_cosmic_template_cov_frac(h_off, h_in)
    else:
        raise ValueError(
            f"Unknown Product A cosmics mode {mode!r}; "
            "use absdiff_diag or unisim_rank1"
        )
    contam = _contam_frac_from_histdata(hd, cosmic_estimate="offbeam")
    if contam.shape[0] != raw.shape[0]:
        raise ValueError(
            f"cosmics bin mismatch for {slug}: contam {contam.shape} vs cov {raw.shape}"
        )
    return _sanitize_cov_frac(scale_cov_frac_by_contamination(raw, contam))


def product_a_rate_cov(
    var_config,
    stage_key: str,
    hd: OverlayHistData,
    *,
    name_suffix: Optional[str] = None,
    cosmics_mode: str = "unisim_rank1",
    detector_mode: str = "diag",
    syst_disk_root: Optional[str] = None,
) -> Optional[np.ndarray]:
    """PRL Product A rate fractional covariance for overlay bands / χ².

    Assembles finalized syst-disk components via :func:`get_syst_unc`
    (``genie_rate``, fully correlated POT/Ntargets; no MCstat, no cosmics pack),
    rebases flux/g4/genie from signal CV onto total-selected CV (topology signal
    layer + stack total from ``hd``), adds detector (see ``detector_mode``;
    nominal ``diag``), then adds Product A cosmics (see ``cosmics_mode``;
    nominal ``unisim_rank1``), scaled by cosmic contamination.
    """
    slug = _product_a_syst_slug(
        var_config.var_save_name, stage_key, name_suffix=name_suffix
    )
    if slug is None:
        return None

    root = syst_disk_root or PRODUCT_A_SYST_ROOT
    # Prefer OverlayHistData bins so rebinned test trees (e.g. χ² 12/20)
    # match syst packs; VariableConfig may still carry production fine bins.
    bins = np.asarray(hd.bins if hd is not None else var_config.bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    n = len(centers)
    proxy = _VarProxy(slug, bins, centers)

    # Signal-CV sources (bkgd_subtract=True): rebase before combining with flats.
    try:
        _unc, cov_sig = get_syst_unc(
            proxy,
            syst_disk_root=root,
            syst_components=_PRODUCT_A_SIGNAL_CV_COMPONENTS,
            genie_cov_frac_key="genie_rate",
            skip_missing_vars=True,
            plot=False,
        )
    except Exception as ex:
        print(f"  [A syst] base skip {slug}: {ex}", flush=True)
        cov_sig = np.zeros((n, n), dtype=float)

    cov_sig = np.asarray(cov_sig if cov_sig is not None else 0.0, dtype=float)
    if cov_sig.shape != (n, n):
        cov_sig = np.zeros((n, n), dtype=float)
    cov_sig = _sanitize_cov_frac(cov_sig)

    counts = _overlay_signal_and_total(hd, cosmic_estimate="offbeam")
    if counts is not None:
        n_sig, n_tot = counts
        if n_sig.shape[0] == n and n_tot.shape[0] == n:
            cov_sig = rebase_fraccov_signal_to_total(cov_sig, n_sig, n_tot)
        else:
            print(
                f"  [A syst] {slug}: overlay counts len "
                f"sig={n_sig.shape[0]} tot={n_tot.shape[0]} != {n}; skip rebase",
                flush=True,
            )
    else:
        print(
            f"  [A syst] {slug}: no topology signal layer on OverlayHistData; "
            "flux/g4/genie left unrebased",
            flush=True,
        )

    try:
        _unc, cov_flat = get_syst_unc(
            proxy,
            syst_disk_root=root,
            syst_components=_PRODUCT_A_FLAT_COMPONENTS,
            skip_missing_vars=True,
            plot=False,
        )
    except Exception as ex:
        print(f"  [A syst] flat skip {slug}: {ex}", flush=True)
        cov_flat = np.zeros((n, n), dtype=float)
    cov_flat = np.asarray(cov_flat if cov_flat is not None else 0.0, dtype=float)
    if cov_flat.shape != (n, n):
        cov_flat = np.zeros((n, n), dtype=float)

    cov = _sanitize_cov_frac(cov_sig + cov_flat)
    cov = cov + _product_a_detector_cov_frac(
        proxy, mode=detector_mode, syst_disk_root=root
    )

    try:
        cosmics = _product_a_cosmics_cov_frac(
            slug, hd, mode=cosmics_mode, syst_disk_root=root
        )
        if cosmics is not None:
            cov = cov + cosmics
    except Exception as ex:
        print(f"  [A syst] cosmics skip {slug}: {ex}", flush=True)

    cov = _sanitize_cov_frac(cov)
    if not np.any(cov):
        return None
    return cov


def run_product_a(
    *,
    force_rebuild: bool = False,
    out_dir: Optional[str] = None,
    cosmics_mode: str = "unisim_rank1",
    detector_mode: str = "diag",
    syst_disk_root: Optional[str] = None,
    merged_pkl: Optional[str] = None,
    mc_pot_fix: Optional[float] = None,
) -> None:
    out = out_dir or PRODUCT_A_OUT
    if mc_pot_fix:
        _refuse_in_place_mcpotfix(out, PRODUCT_A_OUT)
    syst_root = syst_disk_root or PRODUCT_A_SYST_ROOT
    print("\n" + "=" * 72 + "\nProduct A (sel_all / cut-stage)\n" + "=" * 72, flush=True)
    print(f"  out_dir={out}", flush=True)
    print(f"  cosmics_mode={cosmics_mode}", flush=True)
    print(f"  detector_mode={detector_mode}", flush=True)
    print(f"  syst_disk_root={syst_root}", flush=True)
    os.environ[SYST_DISK_ENV] = syst_root
    makedirs(out, exist_ok=True)

    batches_dir = _resolve_product_a_batches()
    print(f"  batches_dir={batches_dir}", flush=True)

    # Load aggregate helpers from the reduce script
    import importlib.util

    agg_path = path.join(_SCRIPTS, "event_selection_aggregate.py")
    spec = importlib.util.spec_from_file_location("event_selection_aggregate", agg_path)
    assert spec is not None and spec.loader is not None
    agg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(agg)

    # Prefer explicit / production cache for counts so test dirs only replot.
    cache_pkl = merged_pkl or path.join(PRODUCT_A_OUT, "merged_histdata.pkl")
    local_pkl = path.join(out, "merged_histdata.pkl")
    if (not force_rebuild) and path.isfile(cache_pkl):
        print(f"  loading cached merged counts: {cache_pkl}", flush=True)
        with open(cache_pkl, "rb") as f:
            merged_payload = pickle.load(f)
        merged = merged_payload["merged"]
        pot_str = merged_payload["pot_str"]
        data_pot = merged_payload.get("data_pot")
        if cache_pkl != local_pkl:
            # Lightweight pointer so the test dir records provenance.
            meta_path = path.join(out, "merged_histdata_source.json")
            with open(meta_path, "w") as f:
                json.dump(
                    {
                        "merged_histdata_pkl": cache_pkl,
                        "cosmics_mode": cosmics_mode,
                        "detector_mode": detector_mode,
                        "syst_disk_root": syst_root,
                        "note": "counts reused; overlays regenerated only",
                    },
                    f,
                    indent=2,
                )
    elif (not force_rebuild) and path.isfile(local_pkl):
        print(f"  loading cached merged counts: {local_pkl}", flush=True)
        with open(local_pkl, "rb") as f:
            merged_payload = pickle.load(f)
        merged = merged_payload["merged"]
        pot_str = merged_payload["pot_str"]
        data_pot = merged_payload.get("data_pot")
    else:
        if not path.isdir(batches_dir):
            raise FileNotFoundError(
                f"Product A live-PRL batches missing: {batches_dir}. "
                "Re-run event_selection_batched map phase / chi2avg remake."
            )
        print(f"  aggregating batches from {batches_dir}", flush=True)
        chunk_groups = agg.collect_chunks(batches_dir)
        samples = {}
        for s, files in chunk_groups.items():
            if not files:
                continue
            print(f"  aggregating {len(files)} batches for sample={s}", flush=True)
            samples[s] = agg.aggregate_chunk_files(files)
        if not samples:
            raise RuntimeError(f"No batch pickles under {batches_dir}")
        merged = agg.merge_samples(samples)
        n_hd_fixed, n_bar_fixed = agg.sanitize_merged_histdata_finite(merged)
        if n_hd_fixed or n_bar_fixed:
            print(
                f"  sanitized NaN/inf: {n_hd_fixed} hist keys, {n_bar_fixed} bar rows",
                flush=True,
            )
        totals = agg.accumulate_exposure_totals(chunk_groups)
        exposure_scales = agg.apply_global_exposure_scales(
            merged, totals, f_offbeam_coincident=0.08
        )
        _scale_merged_mc(merged, MC_POT_FIX)
        if exposure_scales.get("scale_mc") is not None:
            exposure_scales["scale_mc"] = float(exposure_scales["scale_mc"]) * MC_POT_FIX
        print(f"  exposure scales: {exposure_scales} (MC includes MC_POT_FIX={MC_POT_FIX})", flush=True)
        data_pot = totals.data_pot if totals.data_pot > 0 else 1.0
        pot_str = agg.get_pot_str(data_pot)
        merged_payload = {
            "merged": merged,
            "data_pot": data_pot,
            "pot_str": pot_str,
            "exposure_totals": totals,
            "exposure_scales": exposure_scales,
            "batches_dir": batches_dir,
            "note": "chi2 from VariableConfig.chi2_mu/proton (plane avg)",
            "cosmics_mode": cosmics_mode,
            "mc_pot_fix_factor": float(MC_POT_FIX),
        }
        with open(local_pkl, "wb") as f:
            pickle.dump(merged_payload, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"  wrote {local_pkl}", flush=True)

    if mc_pot_fix:
        already = merged_payload.get("mc_pot_fix_factor")
        if already:
            print(
                f"  MC POT fix already applied (factor={already}); not scaling again",
                flush=True,
            )
        else:
            sample_hd = next(iter(merged["histdata"].values()))
            before = float(np.sum(sample_hd.mc_hist))
            _scale_merged_mc(merged, float(mc_pot_fix))
            after = float(np.sum(sample_hd.mc_hist))
            scales = dict(merged_payload.get("exposure_scales") or {})
            if scales.get("scale_mc") is not None:
                scales["scale_mc_before_mcpotfix"] = scales["scale_mc"]
                scales["scale_mc"] = float(scales["scale_mc"]) * float(mc_pot_fix)
            merged_payload = dict(merged_payload)
            merged_payload["merged"] = merged
            merged_payload["exposure_scales"] = scales
            merged_payload["mc_pot_fix_factor"] = float(mc_pot_fix)
            merged_payload["mc_pot_fix_note"] = (
                f"Recorded neutrino-MC POT is high by {MC_POT_FIX}. "
                f"Neutrino-MC histograms and bar MC counts multiplied by "
                f"{MC_POT_FIX}. Dirt and cosmics unchanged."
            )
            with open(local_pkl, "wb") as f:
                pickle.dump(merged_payload, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(
                f"  MC POT fix ×{mc_pot_fix}: example sum(mc) {before:.6g} -> {after:.6g}",
                flush=True,
            )
            print(f"  wrote MC-POT-fixed counts -> {local_pkl}", flush=True)
            if scales.get("scale_mc") is not None:
                print(f"  exposure scales (MC corrected): {scales}", flush=True)

    def _syst_loader(ps, stage_key):
        # Match the exact PlotSpec (incl. not_mu suffix), not the first same-named var.
        from analysis_village.numucc_1p0pi.selection_framework import ChunkRunner

        want = (stage_key, ChunkRunner.plot_key(stage_key, ps))
        hd = merged["histdata"].get(want)
        if hd is None:
            for (sk, _pk), h in merged["histdata"].items():
                if (
                    sk == stage_key
                    and h.var_save_name == ps.var_config.var_save_name
                    and h.breakdown_type == ps.breakdown_type
                ):
                    hd = h
                    break
        if hd is None:
            return None
        return product_a_rate_cov(
            ps.var_config,
            stage_key,
            hd,
            name_suffix=getattr(ps, "name_suffix", None),
            cosmics_mode=cosmics_mode,
            detector_mode=detector_mode,
            syst_disk_root=syst_root,
        )

    print(f"  rendering overlays -> {out}", flush=True)
    agg.render_overlay_plots(
        merged,
        plot_label_map={},
        save_fig_dir=out,
        pot_str=pot_str,
        save_fig=True,
        show_fig=False,
        syst_disk_root=None,  # use custom loader only
        syst_cov_loader=_syst_loader,
        cosmic_estimate="offbeam",
    )
    try:
        agg.render_summary_breakdown_plot(
            merged, out, save_fig=True, show_fig=False, cosmic_estimate="offbeam"
        )
    except Exception as ex:
        print(f"  summary bar skip: {ex}", flush=True)

    # Reportable counts: topology overlays only (signal/bkg well-defined)
    items = []
    for (stage_key, plot_key), hd in merged["histdata"].items():
        if getattr(hd, "breakdown_type", "") != "topology":
            continue
        slug = f"{stage_key}__{hd.var_save_name}"
        items.append((slug, hd))
    export_counts_report(
        out,
        product="A",
        histdata_items=items,
        extra_meta={
            "batches_dir": batches_dir,
            "merged_pkl": cache_pkl if path.isfile(cache_pkl) else local_pkl,
            "syst_disk_root": syst_root,
            "syst_kind": "rate",
            "mcstat": "excluded",
            "detector": f"product_a_detector_{detector_mode}",
            "cosmics": f"product_a_{cosmics_mode}_x_contam_mc_plus_offbeam",
            "chi2": "full_absolute_covariance",
            "data_pot": data_pot,
            "cosmics_mode": cosmics_mode,
            "detector_mode": detector_mode,
            "mc_pot_fix_factor": merged_payload.get("mc_pot_fix_factor") or mc_pot_fix,
        },
    )


# ===========================================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--product",
        choices=("A", "B", "both"),
        default="both",
        help="Which product overlays to build (default: both)",
    )
    p.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Ignore cached overlay / merged counts and refill",
    )
    p.add_argument(
        "--product-a-out",
        default="",
        help="Override Product A overlay output directory (default: production path)",
    )
    p.add_argument(
        "--product-b-out",
        default="",
        help="Override Product B overlay output directory (default: production path)",
    )
    p.add_argument(
        "--product-a-syst-root",
        default="",
        help="Override Product A syst-disk root (default: production path)",
    )
    p.add_argument(
        "--product-b-syst-root",
        default="",
        help="Override Product B syst-disk root (default: production path)",
    )
    p.add_argument(
        "--product-a-cosmics-mode",
        choices=("unisim_rank1", "absdiff_diag"),
        default="unisim_rank1",
        help="Product A cosmics template: coherent offbeam→intime unisim "
        "rank-1 (default / nominal) or uncorrelated |Δ| envelope (test)",
    )
    p.add_argument(
        "--product-a-detector-mode",
        choices=("diag", "unisim_full"),
        default="diag",
        help="Product A detector cov: diagonal of nested WireMod+DENT "
        "(default / nominal) or full on-disk unisim correlations (test)",
    )
    p.add_argument(
        "--product-a-merged-pkl",
        default="",
        help="Explicit Product A merged_histdata.pkl to load for counts "
        "(needed when rebinned counts live outside production productA_sel_all/)",
    )
    p.add_argument(
        "--mc-pot-fix",
        action="store_true",
        help="Multiply neutrino MC by MC_POT_FIX (recorded MC POT is high) "
        "and write overlays under parallel *_mcpotfix directories",
    )
    return p.parse_args()


def main():
    args = parse_args()
    t0 = datetime.now()
    print(f"data_mc_overlay_products start {t0.isoformat()}", flush=True)
    print(f"  OUTPUT_ROOT={OUTPUT_ROOT}", flush=True)
    a_syst = args.product_a_syst_root or PRODUCT_A_SYST_ROOT
    b_syst = args.product_b_syst_root or PRODUCT_B_SYST_ROOT
    print(f"  Product A syst={a_syst}", flush=True)
    print(f"  Product B syst={b_syst}", flush=True)
    makedirs(OUTPUT_ROOT, exist_ok=True)
    mc_pot_fix = MC_POT_FIX if args.mc_pot_fix else None
    if args.mc_pot_fix:
        print(
            f"  MC POT fix: recorded MC POT is high; "
            f"scale neutrino MC by {MC_POT_FIX}",
            flush=True,
        )
    a_out = args.product_a_out or (
        PRODUCT_A_OUT + "_mcpotfix" if args.mc_pot_fix else None
    )
    b_out = args.product_b_out or (
        PRODUCT_B_OUT + "_mcpotfix" if args.mc_pot_fix else None
    )

    if args.product in ("B", "both"):
        run_product_b(
            force_rebuild=args.force_rebuild,
            out_dir=b_out,
            syst_disk_root=args.product_b_syst_root or None,
            mc_pot_fix=mc_pot_fix,
        )
    if args.product in ("A", "both"):
        run_product_a(
            force_rebuild=args.force_rebuild,
            out_dir=a_out,
            cosmics_mode=args.product_a_cosmics_mode,
            detector_mode=args.product_a_detector_mode,
            syst_disk_root=args.product_a_syst_root or None,
            merged_pkl=args.product_a_merged_pkl or None,
            mc_pot_fix=mc_pot_fix,
        )

    print(f"\nDone in {datetime.now() - t0}", flush=True)


if __name__ == "__main__":
    main()
