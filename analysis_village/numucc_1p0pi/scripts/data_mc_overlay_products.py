#!/usr/bin/env python3
"""Data–MC comparison overlays for PRL Product A (sel_all cut-stage) and Product B (sel_mup).

Product B
  - Style / counts cache: same as ``selected_xsec_overlay.py`` / notebook.
  - Data: ``sel_mup`` (not sel_all). After beam-quality cuts, ``n_evt_good`` must be
    **12,804** (used to pin the correct χ²μ / FV campaign among variants).
  - Systematics: PRL ``productB_sel_mup`` CategorySummary, ``syst_kind="rate"``
    (cosmics already contamination-scaled ``SelectedRate``).

Product A
  - Prefer live-PRL batched counts under ``event_selection-batched-live-PRL``;
    aggregate + save ``merged_histdata.pkl`` (May ``…-20260525`` is outdated /
    pre-DQ and is not used).
  - Systematics: PRL ``productA_sel_all`` combined on the fly (no CategorySummary).
    Exclude MCstat; GENIE ``genie_rate``; cosmics raw template × contamination
    from topology layer 0 (Product A NPZ has no ``SelectedRate``).
  - Flat POT / Ntargets: fully correlated (same as CategorySummary).

Outputs under::

    /exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/data_mc_overlays/
      productB_sel_mup/   overlay_histdata.pkl, counts_report.npz, PNGs
      productA_sel_all/   merged_histdata.pkl, counts_report.npz, PNGs

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

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    prl_syst_disk_root,
)
from analysis_village.numucc_1p0pi.categories import topology_labels  # noqa: E402
from analysis_village.numucc_1p0pi.syst_category_summary import (  # noqa: E402
    NTARGETS_FRAC_UNC_PCT,
    POT_FRAC_UNC_PCT,
    _flat_cov_frac,
)
from analysis_village.numucc_1p0pi.syst_cosmics_common import (  # noqa: E402
    flat_uncorrelated_cov_frac,
    scale_cov_frac_by_contamination,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import SYST_DISK_ENV  # noqa: E402
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_category_summary_syst_unc,
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
EXPECTED_DATA_N = 12804

PRODUCT_A_BATCHES_CHI2AVG = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL-chi2avg/batches"
)
PRODUCT_A_BATCHES_LEGACY_I2 = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL/batches"
)


def _resolve_product_a_batches() -> str:
    """Prefer chi2.avg remake batches when complete; else legacy I2 live-PRL."""
    import glob as _glob

    n_avg = len(_glob.glob(path.join(PRODUCT_A_BATCHES_CHI2AVG, "*__batch_*.pkl")))
    if n_avg >= 180:
        return PRODUCT_A_BATCHES_CHI2AVG
    return PRODUCT_A_BATCHES_LEGACY_I2


PRODUCT_A_BATCHES = PRODUCT_A_BATCHES_CHI2AVG  # default target; resolved at run time


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


def _topo_signal_bkg_from_mc_hist(mc_hist: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Split stacked topology MC into signal (last layer) and background (rest)."""
    mh = np.asarray(mc_hist, dtype=float)
    if mh.ndim != 2 or mh.shape[0] == 0:
        z = np.zeros(mh.shape[-1] if mh.ndim else 0, dtype=float)
        return z, z
    sig = mh[_TOPO_SIGNAL_IDX].copy()
    bkg = mh.sum(axis=0) - sig
    return sig, bkg


def export_counts_report(
    out_dir: str,
    *,
    product: str,
    histdata_items: Sequence[Tuple[str, OverlayHistData]],
    extra_meta: Optional[Mapping[str, Any]] = None,
) -> str:
    """Write NPZ + JSON manifest of bin-by-bin data / MC signal / MC background."""
    makedirs(out_dir, exist_ok=True)
    arrays: Dict[str, np.ndarray] = {}
    rows = []
    for slug, hd in histdata_items:
        if not getattr(hd, "has_data", False) and hd.data_hist is None:
            continue
        data = np.asarray(hd.data_hist, dtype=float).ravel()
        bins = np.asarray(hd.bins, dtype=float)
        if getattr(hd, "has_mc", False) and hd.mc_hist is not None:
            sig, bkg = _topo_signal_bkg_from_mc_hist(hd.mc_hist)
            mc_tot = np.asarray(hd.mc_hist, dtype=float).sum(axis=0)
        else:
            sig = bkg = mc_tot = np.zeros_like(data)
        # include intime/dirt/offbeam in "other" prediction if present
        other = np.zeros_like(data)
        for attr in ("intime_hist", "dirt_hist", "offbeam_hist"):
            h = getattr(hd, attr, None)
            if h is not None:
                other = other + np.asarray(h, dtype=float).ravel()
        key = slug.replace("/", "_")
        arrays[f"{key}__bins"] = bins
        arrays[f"{key}__data"] = data
        arrays[f"{key}__mc_signal"] = sig
        arrays[f"{key}__mc_background"] = bkg
        arrays[f"{key}__mc_total"] = mc_tot
        arrays[f"{key}__mc_other"] = other  # intime/dirt/offbeam
        rows.append(
            {
                "slug": slug,
                "var_save_name": getattr(hd, "var_save_name", ""),
                "breakdown_type": getattr(hd, "breakdown_type", ""),
                "n_bins": int(len(data)),
                "data_sum": float(np.nansum(data)),
                "mc_signal_sum": float(np.nansum(sig)),
                "mc_background_sum": float(np.nansum(bkg)),
                "mc_total_sum": float(np.nansum(mc_tot)),
                "mc_other_sum": float(np.nansum(other)),
                "topology_labels": list(topology_labels),
                "signal_layer_index": len(topology_labels) - 1,
                "cosmic_layer_index": 0,
            }
        )

    npz_path = path.join(out_dir, COUNTS_REPORT_NAME)
    np.savez_compressed(npz_path, **arrays)
    manifest = {
        "schema": "data_mc_overlay_counts_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "product": product,
        "npz_path": npz_path,
        "n_variables": len(rows),
        "variables": rows,
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


def _product_b_syst_cov(var_config):
    """CategorySummary total_rate (includes SelectedRate cosmics)."""
    try:
        _unc, cov = get_category_summary_syst_unc(
            var_config,
            syst_kind="rate",
            syst_disk_root=PRODUCT_B_SYST_ROOT,
        )
        if cov is None or not np.any(np.isfinite(cov)):
            return None
        return cov
    except Exception as ex:
        print(f"  [B syst] skip {var_config.var_save_name}: {ex}", flush=True)
        return None


def run_product_b(*, force_rebuild: bool = False) -> None:
    print("\n" + "=" * 72 + "\nProduct B (sel_mup)\n" + "=" * 72, flush=True)
    os.environ[SYST_DISK_ENV] = PRODUCT_B_SYST_ROOT
    makedirs(PRODUCT_B_OUT, exist_ok=True)

    plot_set = {
        "tag": "productB_sel_mup",
        "output_dir": PRODUCT_B_OUT,
        "mc_dir": PRODUCT_B_MC_DIR,
        "mc_filename_str": PRODUCT_B_MC_FN,
        "data_dir": PRODUCT_B_DATA_DIR,
        "data_filename_str": PRODUCT_B_DATA_FN,
    }

    # Align sxo module knobs with this driver
    sxo.APPLY_BEAM_QUALITY = True
    sxo.LOAD_SYST = True
    sxo.SYST_DISK_ROOT = PRODUCT_B_SYST_ROOT
    sxo.FORCE_REBUILD_COUNTS = force_rebuild
    sxo.SAVE_FIG = True
    sxo.PLOT = False
    sxo.APPROVAL = ""
    sxo.TEXTCHI2 = True

    payload = None if force_rebuild else load_overlay_counts(PRODUCT_B_OUT)
    bq_n = None
    if payload is not None:
        print(f"  replot from counts: {histdata_pkl_path(PRODUCT_B_OUT)}", flush=True)
        raw_pot_label = payload.get("pot_label") or "Events / Bin"
        pot_text = format_pot_corner_text(raw_pot_label)
        pot_label = strip_pot_from_ylabel(raw_pot_label) or "Events / Bin"
        histdata_map = payload["histdata"]
        bq_n = (payload.get("plot_set") or {}).get("n_evt_good")
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
        pot_text = format_pot_corner_text(pot_label_raw)
        pot_label = strip_pot_from_ylabel(pot_label_raw) or "Events / Bin"
        histdata_map = build_overlay_histdata_map(
            sxo.VAR_CONFIGS,
            sxo.BREAKDOWN_TYPES,
            mc_df=mc_evt,
            data_df=data_evt,
        )
        plot_set["n_evt_good"] = bq_n
        pkl = save_overlay_counts(
            PRODUCT_B_OUT,
            histdata_map,
            # Keep POT in pickle so corner text survives ylabel stripping on replot.
            pot_label=pot_label_raw,
            plot_set=plot_set,
            var_save_names=[vc.var_save_name for vc in sxo.VAR_CONFIGS],
            breakdown_types=sxo.BREAKDOWN_TYPES,
        )
        print(f"  wrote counts -> {pkl}", flush=True)
        del mc_evt, mc_hdr, data_evt, data_hdr
        gc.collect()

    if bq_n is not None and int(bq_n) != EXPECTED_DATA_N:
        print(
            f"  WARN: cached counts report n_evt_good={bq_n} "
            f"(expected {EXPECTED_DATA_N}); rebuild with --force-rebuild",
            flush=True,
        )

    plot_overlay_counts_map(
        histdata_map,
        sxo.VAR_CONFIGS,
        sxo.BREAKDOWN_TYPES,
        pot_label=pot_label,
        out_dir=PRODUCT_B_OUT,
        get_syst=_product_b_syst_cov,
        ax_ylim_ratio=sxo.AX_YLIM_RATIO,
        ratio=sxo.RATIO,
        textloc=sxo.TEXTLOC,
        approval="",
        pot_text=pot_text,
        save_fig=True,
        plot=False,
        textchi2=True,
    )

    # Prefer topology breakdown for the reportable signal/bkg split
    items = []
    for vc in sxo.VAR_CONFIGS:
        key = (vc.var_save_name, "topology")
        if key in histdata_map:
            items.append((vc.var_save_name, histdata_map[key]))
    export_counts_report(
        PRODUCT_B_OUT,
        product="B",
        histdata_items=items,
        extra_meta={
            "data_dir": PRODUCT_B_DATA_DIR,
            "mc_dir": PRODUCT_B_MC_DIR,
            "syst_disk_root": PRODUCT_B_SYST_ROOT,
            "syst_kind": "rate",
            "n_evt_good": bq_n,
            "expected_n_evt_good": EXPECTED_DATA_N,
        },
    )


# ===========================================================================
# Product A — sel_all cut-stage from live-PRL batches
# ===========================================================================


def _product_a_syst_slug(var_save_name: str, stage_key: str) -> Optional[str]:
    """Map pipeline PlotSpec names onto Product A PRL syst keys."""
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
        if var_save_name == "chi2_mu":
            return f"chi2_avg_mu__at_{stage_key}"
        if var_save_name == "chi2_p":
            return f"chi2_avg_p__at_{stage_key}"
    return None


def _contam_frac_from_histdata(hd: OverlayHistData) -> np.ndarray:
    """Per-bin cosmic contamination from topology MC layer 0."""
    if hd.mc_hist is None:
        n = len(hd.bins) - 1
        return np.zeros(n, dtype=float)
    total = np.asarray(hd.mc_hist, dtype=float).sum(axis=0)
    cosmic = np.asarray(hd.mc_hist[_TOPO_COSMIC_IDX], dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(total > 0, cosmic / total, 0.0)
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


def product_a_rate_cov(
    var_config,
    stage_key: str,
    hd: OverlayHistData,
) -> Optional[np.ndarray]:
    """PRL Product A rate covariance with contamination-scaled cosmics (no MCstat)."""
    slug = _product_a_syst_slug(var_config.var_save_name, stage_key)
    if slug is None:
        return None

    bins = np.asarray(var_config.bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    proxy = _VarProxy(slug, bins, centers)

    # Disk components except mcstat / cosmics / flat (we add flats correlated).
    try:
        _unc, cov = get_syst_unc(
            proxy,
            syst_disk_root=PRODUCT_A_SYST_ROOT,
            syst_components=("flux", "g4", "genie", "detector"),
            genie_cov_frac_key="genie_rate",
            skip_missing_vars=True,
            plot=False,
        )
    except Exception as ex:
        print(f"  [A syst] base skip {slug}: {ex}", flush=True)
        cov = np.zeros((len(centers), len(centers)), dtype=float)

    cov = np.asarray(cov, dtype=float)
    if cov.shape != (len(centers), len(centers)):
        cov = np.zeros((len(centers), len(centers)), dtype=float)
    cov = _sanitize_cov_frac(cov)

    # Cosmics: raw template × contamination (Product A NPZ has no SelectedRate)
    try:
        from analysis_village.numucc_1p0pi.syst_disk_layout import syst_disk_paths

        cpath = syst_disk_paths(PRODUCT_A_SYST_ROOT)["cosmics"]
        blob = dict(np.load(cpath, allow_pickle=True))
        if slug in blob:
            cell = blob[slug].item()
            raw = np.asarray(cell["Cosmics"]["cov_frac"], dtype=float)
            raw = flat_uncorrelated_cov_frac(raw)
            contam = _contam_frac_from_histdata(hd)
            if contam.shape[0] == raw.shape[0]:
                cov = cov + _sanitize_cov_frac(
                    scale_cov_frac_by_contamination(raw, contam)
                )
            else:
                print(
                    f"  [A syst] cosmics bin mismatch for {slug}: "
                    f"contam {contam.shape} vs cov {raw.shape}",
                    flush=True,
                )
    except Exception as ex:
        print(f"  [A syst] cosmics skip {slug}: {ex}", flush=True)

    n = len(centers)
    cov = cov + _flat_cov_frac(n, POT_FRAC_UNC_PCT)
    cov = cov + _flat_cov_frac(n, NTARGETS_FRAC_UNC_PCT)
    cov = _sanitize_cov_frac(cov)
    if not np.any(cov):
        return None
    return cov


def run_product_a(*, force_rebuild: bool = False) -> None:
    print("\n" + "=" * 72 + "\nProduct A (sel_all / cut-stage)\n" + "=" * 72, flush=True)
    os.environ[SYST_DISK_ENV] = PRODUCT_A_SYST_ROOT
    makedirs(PRODUCT_A_OUT, exist_ok=True)

    batches_dir = _resolve_product_a_batches()
    print(f"  batches_dir={batches_dir}", flush=True)

    # Load aggregate helpers from the reduce script
    import importlib.util

    agg_path = path.join(_SCRIPTS, "event_selection_aggregate.py")
    spec = importlib.util.spec_from_file_location("event_selection_aggregate", agg_path)
    assert spec is not None and spec.loader is not None
    agg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(agg)

    merged_pkl = path.join(PRODUCT_A_OUT, "merged_histdata.pkl")
    if (not force_rebuild) and path.isfile(merged_pkl):
        print(f"  loading cached merged counts: {merged_pkl}", flush=True)
        with open(merged_pkl, "rb") as f:
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
        print(f"  exposure scales: {exposure_scales}", flush=True)
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
        }
        with open(merged_pkl, "wb") as f:
            pickle.dump(merged_payload, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"  wrote {merged_pkl}", flush=True)

    def _syst_loader(ps, stage_key):
        hd = None
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
        return product_a_rate_cov(ps.var_config, stage_key, hd)

    print(f"  rendering overlays -> {PRODUCT_A_OUT}", flush=True)
    agg.render_overlay_plots(
        merged,
        plot_label_map={},
        save_fig_dir=PRODUCT_A_OUT,
        pot_str=pot_str,
        save_fig=True,
        show_fig=False,
        syst_disk_root=None,  # use custom loader only
        syst_cov_loader=_syst_loader,
        cosmic_estimate="intime",
    )
    try:
        agg.render_summary_breakdown_plot(
            merged, PRODUCT_A_OUT, save_fig=True, show_fig=False, cosmic_estimate="intime"
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
        PRODUCT_A_OUT,
        product="A",
        histdata_items=items,
        extra_meta={
            "batches_dir": batches_dir,
            "merged_pkl": merged_pkl,
            "syst_disk_root": PRODUCT_A_SYST_ROOT,
            "syst_kind": "rate",
            "mcstat": "excluded",
            "cosmics": "raw_template_x_contamination_from_topology_layer0",
            "data_pot": data_pot,
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
    return p.parse_args()


def main():
    args = parse_args()
    t0 = datetime.now()
    print(f"data_mc_overlay_products start {t0.isoformat()}", flush=True)
    print(f"  OUTPUT_ROOT={OUTPUT_ROOT}", flush=True)
    print(f"  Product A syst={PRODUCT_A_SYST_ROOT}", flush=True)
    print(f"  Product B syst={PRODUCT_B_SYST_ROOT}", flush=True)
    makedirs(OUTPUT_ROOT, exist_ok=True)

    if args.product in ("B", "both"):
        run_product_b(force_rebuild=args.force_rebuild)
    if args.product in ("A", "both"):
        run_product_a(force_rebuild=args.force_rebuild)

    print(f"\nDone in {datetime.now() - t0}", flush=True)


if __name__ == "__main__":
    main()
