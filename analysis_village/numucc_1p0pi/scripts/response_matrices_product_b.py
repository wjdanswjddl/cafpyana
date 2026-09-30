#!/usr/bin/env python3
"""Batch-accumulate Product B response matrices + efficiencies (sel_mup MC).

Streaming over ``*.df`` files avoids loading the full 2k-file MC sample into
RAM. Per file we accumulate:

* ``nevts_allmc`` — truth signal on ``mcnu`` with ``signal_truth_fv='none'``
  (vertex FV + 1μ1p0π topology; **no** μ/p end ``per_tpc`` cut). This is the
  cross-section signal strength; end-FV / selection enter only via efficiency.
* ``nevts_sel_truth`` / ``nevts_sel_reco`` — selected events that are true
  signal under the **selection** definition (``per_tpc``), for migration.
* ``nevts_allsel_reco`` — all selected (signal+bkg) reco
* ``reco_vs_true`` — migration ``histogram2d(truth, reco; wgt_sel_truth)``

Then ``eff = sel_truth / allmc`` and ``R = get_response_matrix(reco_vs_true, eff)``.
Same split as Gen1 ``signal_hists(..., mode='unfold', signal_truth_fv='none')``.

POT scale ``data_tot_pot / mc_tot_pot`` is applied once at the end (flat on MC).
Efficiency and response are scale-invariant; absolute rates are scaled for the
unfold pack.

Overlay data / OffBeam background counts are taken from the Product B
``counts_report.npz`` (already produced by ``data_mc_overlay_products.py``),
not re-filled here.

    Usage::

    /exp/sbnd/app/users/munjung/xsec/freeze/cafpyana/envs/venv_py310_cafpyana/bin/python \\
      analysis_village/numucc_1p0pi/scripts/response_matrices_product_b.py \\
      [--workers 16] [--max-files 0] [--force]

Outputs under::

    /exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL/response_matrices/
"""

from __future__ import annotations

import argparse
import gc
import glob
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from os import makedirs, path
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

_REPO = Path(__file__).resolve().parents[3]
_SCRIPTS = Path(__file__).resolve().parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from pyanalib.split_df_helpers_new import get_n_split, load_dfs  # noqa: E402

from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.categories import (  # noqa: E402
    IsNuInFV_NumuCC_1p0pi,
)
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_clipped_evts,
    get_response_matrix,
    get_topo_category,
    plot_heatmap,
    fig_ext,
    dpi,
)

# Unfold / xsec truth signal (efficiency denominator + Wiener model).
SIGNAL_TRUTH_FV_UNFOLD = "none"
# True-signal label among *selected* events (migration numerator) matches reco FV.
SIGNAL_TRUTH_FV_SELECTED = "per_tpc"
from analysis_village.numucc_1p0pi.event_selection_batch_core import (  # noqa: E402
    prefix_mcnu_columns,
)
from analysis_village.numucc_1p0pi.evt_derived_kinematics import (  # noqa: E402
    ensure_mc_level_phi_mcnu,
)
from pyanalib.variable_calculator import add_mc_cc1p0pi_tki_mcnu  # noqa: E402
import data_mc_overlay_products as dmo  # noqa: E402

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PRL_ROOT = Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL")
RESP_DIR = PRL_ROOT / "response_matrices"
FIG_DIR = RESP_DIR / "plots"
COUNTS_NPZ = PRL_ROOT / "data_mc_overlays/productB_sel_mup/counts_report.npz"

DEFAULT_BATCH_FILES = 20
DEFAULT_WORKERS = 16


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-dir", default=dmo.PRODUCT_B_MC_DIR)
    p.add_argument("--mc-fn", default=dmo.PRODUCT_B_MC_FN)
    p.add_argument("--data-dir", default=dmo.PRODUCT_B_DATA_DIR)
    p.add_argument("--data-fn", default=dmo.PRODUCT_B_DATA_FN)
    p.add_argument(
        "--counts-npz",
        default=str(COUNTS_NPZ),
        help="Product B counts_report.npz for data/bkg ingredients",
    )
    p.add_argument("--out-dir", default=str(RESP_DIR))
    p.add_argument(
        "--batch-files",
        type=int,
        default=DEFAULT_BATCH_FILES,
        help="Unused with --workers>1 (kept for CLI compat); serial flush size when workers=1",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Parallel processes for Pass 1b POT sum + Pass 2 accumulate (1 = serial)",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="If >0, only process the first N MC files (debug)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing response_matrices.npz",
    )
    p.add_argument("--no-plots", action="store_true")
    p.add_argument(
        "--selection-eff-only",
        action="store_true",
        help=(
            "Write selection_efficiency.npz only: generated and selected "
            "1μ1p0π with per_tpc (vertex FV + μ/p containment) in var_nu_col. "
            "Does not rebuild response_matrices.npz."
        ),
    )
    return p.parse_args(argv)


def _list_df_files(sample_dir: str, filename_str: str) -> List[str]:
    pattern = path.join(sample_dir, f"*{filename_str}*.df")
    return sorted(glob.glob(pattern))


def _hdr_pot_one(fp: str) -> float:
    try:
        n_split = int(get_n_split(fp))
        dfs = load_dfs(fp, ["hdr"], n_max_concat=n_split)
        pot = float(dfs["hdr"]["pot"].sum())
        del dfs
        return pot
    except Exception as ex:
        print(f"  WARN: skip POT for {fp}: {ex}", flush=True)
        return 0.0


def _sum_hdr_pot(files: Sequence[str], workers: int = 1) -> float:
    if workers <= 1:
        tot = 0.0
        for fp in tqdm(files, desc="sum MC POT", leave=False):
            tot += _hdr_pot_one(fp)
            gc.collect()
        return tot
    tot = 0.0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(_hdr_pot_one, fp) for fp in files]
        for fut in tqdm(as_completed(futs), total=len(futs), desc="sum MC POT", leave=False):
            tot += float(fut.result())
    return tot


def _data_tot_pot_and_n() -> Tuple[float, int]:
    """Live beam-quality data POT + n_evt_good (same as Product B overlays)."""
    import selected_xsec_overlay as sxo

    sxo.APPLY_BEAM_QUALITY = True
    data_evt, data_hdr = sxo.load_data_sample(dmo.PRODUCT_B_DATA_DIR, dmo.PRODUCT_B_DATA_FN)
    n = int(len(data_evt))
    pot = float(data_hdr["pot"].sum())
    del data_evt, data_hdr
    gc.collect()
    return pot, n


def _empty_acc(nbins: int) -> Dict[str, np.ndarray]:
    return {
        "nevts_allmc": np.zeros(nbins, dtype=np.float64),
        "nevts_sel_truth": np.zeros(nbins, dtype=np.float64),
        "nevts_sel_reco": np.zeros(nbins, dtype=np.float64),
        "nevts_allsel_reco": np.zeros(nbins, dtype=np.float64),
        "reco_vs_true": np.zeros((nbins, nbins), dtype=np.float64),
    }


def _merge_acc(
    dst: Dict[str, Dict[str, np.ndarray]], src: Mapping[str, Mapping[str, np.ndarray]]
) -> None:
    for vsn, pack in src.items():
        for k, arr in pack.items():
            dst[vsn][k] += np.asarray(arr, dtype=np.float64)


def _var_spec_list(var_configs: Sequence[Any]) -> List[Dict[str, Any]]:
    """Pickle-friendly VariableConfig fields for worker processes."""
    out = []
    for vc in var_configs:
        out.append(
            {
                "var_save_name": vc.var_save_name,
                "bins": np.asarray(vc.bins, dtype=float),
                "var_nu_col": vc.var_nu_col,
                "var_evt_truth_col": vc.var_evt_truth_col,
                "var_evt_reco_col": vc.var_evt_reco_col,
            }
        )
    return out


def _prepare_evt_mcnu(
    evt: pd.DataFrame, mcnu: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Prefix/annotate ``mcnu`` like GENIE/syst chunk maps; add topo + unit weights."""
    evt = evt.copy()
    mcnu = mcnu.copy()

    # Raw HDF mcnu has no leading ``mc`` level; categories + var_nu_col expect it.
    prefix_mcnu_columns(mcnu)
    mcnu = ensure_mc_level_phi_mcnu(mcnu)
    mcnu = add_mc_cc1p0pi_tki_mcnu(mcnu)

    if "mc" in evt.columns.get_level_values(0):
        evt.loc[evt.mc.iscc.isna(), ("mc", "iscc")] = 999
    if "mc" in mcnu.columns.get_level_values(0):
        # after prefix, iscc is under mc; pad MultiIndex assign carefully
        try:
            mcnu.loc[mcnu.mc.iscc.isna(), ("mc", "iscc")] = 999
        except Exception:
            iscc_key = None
            for c in mcnu.columns:
                if isinstance(c, tuple) and len(c) >= 2 and c[0] == "mc" and c[1] == "iscc":
                    iscc_key = c
                    break
            if iscc_key is not None:
                mcnu.loc[mcnu[iscc_key].isna(), iscc_key] = 999

    # Selection-matched true signal on evt (per_tpc). Always recompute so a
    # stale topo_categ baked into the DF cannot silently change the response.
    evt.loc[:, "topo_categ"] = get_topo_category(
        evt, signal_truth_fv=SIGNAL_TRUTH_FV_SELECTED
    )

    # unit weights for accumulation; global POT scale applied later
    evt["pot_weight"] = np.ones(len(evt), dtype=np.float64)
    mcnu["pot_weight"] = np.ones(len(mcnu), dtype=np.float64)
    return evt, mcnu


def _accumulate_file_worker(
    fp: str, var_specs: Sequence[Mapping[str, Any]]
) -> Optional[Dict[str, Dict[str, np.ndarray]]]:
    """Process one MC file → per-variable histograms (for process pool)."""
    try:
        n_split = int(get_n_split(fp))
        dfs = load_dfs(fp, ["evt", "mcnu"], n_max_concat=n_split)
    except Exception as ex:
        print(f"  WARN: skip {fp}: {ex}", flush=True)
        return None
    evt, mcnu = _prepare_evt_mcnu(dfs["evt"], dfs["mcnu"])
    del dfs
    evt_sig = evt[evt.topo_categ == 1]
    mcnu_sig = mcnu[
        IsNuInFV_NumuCC_1p0pi(mcnu, signal_truth_fv=SIGNAL_TRUTH_FV_UNFOLD)
    ]

    acc: Dict[str, Dict[str, np.ndarray]] = {
        spec["var_save_name"]: _empty_acc(len(spec["bins"]) - 1) for spec in var_specs
    }
    # Minimal Namespace-like objects for get_clipped_evts
    for spec in var_specs:
        vsn = spec["var_save_name"]
        bins = np.asarray(spec["bins"], dtype=float)
        nb = len(bins) - 1
        a = acc[vsn]

        var_allmc, w_allmc = get_clipped_evts(
            mcnu_sig, spec["var_nu_col"], bins, var_save_name=vsn
        )
        h, _ = np.histogram(var_allmc, bins=bins, weights=w_allmc)
        a["nevts_allmc"] += np.asarray(h, dtype=np.float64)

        var_sel_t, w_sel_t = get_clipped_evts(
            evt_sig, spec["var_evt_truth_col"], bins, var_save_name=vsn
        )
        var_sel_r, w_sel_r = get_clipped_evts(
            evt_sig, spec["var_evt_reco_col"], bins, var_save_name=vsn
        )
        ht, _ = np.histogram(var_sel_t, bins=bins, weights=w_sel_t)
        hr, _ = np.histogram(var_sel_r, bins=bins, weights=w_sel_r)
        a["nevts_sel_truth"] += np.asarray(ht, dtype=np.float64)
        a["nevts_sel_reco"] += np.asarray(hr, dtype=np.float64)

        var_allsel_r, w_allsel_r = get_clipped_evts(
            evt, spec["var_evt_reco_col"], bins, var_save_name=vsn
        )
        ha, _ = np.histogram(var_allsel_r, bins=bins, weights=w_allsel_r)
        a["nevts_allsel_reco"] += np.asarray(ha, dtype=np.float64)

        if nb == 1:
            a["reco_vs_true"][0, 0] += float(np.sum(ht))
        else:
            rvt, _, _ = np.histogram2d(
                np.asarray(var_sel_t, dtype=float),
                np.asarray(var_sel_r, dtype=float),
                weights=np.asarray(w_sel_t, dtype=float),
                bins=[bins, bins],
            )
            a["reco_vs_true"] += np.asarray(rvt, dtype=np.float64)

    del evt, mcnu, evt_sig, mcnu_sig
    gc.collect()
    return acc


def _empty_fvcont(nbins: int) -> Dict[str, np.ndarray]:
    """Selection-efficiency counts: full FV + containment on both sides."""
    return {
        "nevts_all_fvcont": np.zeros(nbins, dtype=np.float64),
        "nevts_sel_fvcont": np.zeros(nbins, dtype=np.float64),
    }


def _accumulate_fvcont_worker(
    fp: str, var_specs: Sequence[Mapping[str, Any]]
) -> Optional[Dict[str, Dict[str, np.ndarray]]]:
    """Generated and selected signal with ``per_tpc`` (Gen-1 + same-TPC containment).

    Both histograms use ``var_nu_col`` so the ratio is the selection efficiency
    in one truth variable. ``per_tpc`` is the same volume as reco
    ``event_contained_per_tpc`` (full ``SBND_Gen1``, including high-YZ).
    """
    try:
        n_split = int(get_n_split(fp))
        dfs = load_dfs(fp, ["evt", "mcnu"], n_max_concat=n_split)
    except Exception as ex:
        print(f"  WARN: skip {fp}: {ex}", flush=True)
        return None
    evt, mcnu = _prepare_evt_mcnu(dfs["evt"], dfs["mcnu"])
    del dfs
    evt_sig = evt[evt.topo_categ == 1]
    mcnu_sig = mcnu[
        IsNuInFV_NumuCC_1p0pi(mcnu, signal_truth_fv=SIGNAL_TRUTH_FV_SELECTED)
    ]
    acc: Dict[str, Dict[str, np.ndarray]] = {
        spec["var_save_name"]: _empty_fvcont(len(spec["bins"]) - 1) for spec in var_specs
    }
    for spec in var_specs:
        vsn = spec["var_save_name"]
        bins = np.asarray(spec["bins"], dtype=float)
        a = acc[vsn]
        var_all, w_all = get_clipped_evts(
            mcnu_sig, spec["var_nu_col"], bins, var_save_name=vsn
        )
        h_all, _ = np.histogram(var_all, bins=bins, weights=w_all)
        a["nevts_all_fvcont"] += np.asarray(h_all, dtype=np.float64)
        var_sel, w_sel = get_clipped_evts(
            evt_sig, spec["var_nu_col"], bins, var_save_name=vsn
        )
        h_sel, _ = np.histogram(var_sel, bins=bins, weights=w_sel)
        a["nevts_sel_fvcont"] += np.asarray(h_sel, dtype=np.float64)
    del evt, mcnu, evt_sig, mcnu_sig
    gc.collect()
    return acc


def _accumulate_fvcont(
    mc_files: Sequence[str],
    var_configs: Sequence[Any],
    workers: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    acc_by_var = {
        vc.var_save_name: _empty_fvcont(len(vc.bins) - 1) for vc in var_configs
    }
    var_specs = _var_spec_list(var_configs)
    if workers <= 1:
        for fp in tqdm(mc_files, desc="selection eff"):
            part = _accumulate_fvcont_worker(fp, var_specs)
            if part is not None:
                _merge_acc(acc_by_var, part)
            gc.collect()
        return acc_by_var
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(_accumulate_fvcont_worker, fp, var_specs): fp for fp in mc_files
        }
        for fut in tqdm(as_completed(futs), total=len(futs), desc="selection eff"):
            part = fut.result()
            if part is not None:
                _merge_acc(acc_by_var, part)
    return acc_by_var


def _write_selection_efficiency(args: argparse.Namespace) -> int:
    """Sidecar histograms for selection efficiency (does not touch the response NPZ)."""
    out_dir = Path(args.out_dir)
    makedirs(out_dir, exist_ok=True)
    out_npz = out_dir / "selection_efficiency.npz"
    resp_npz = out_dir / "response_matrices.npz"
    if not resp_npz.is_file():
        raise FileNotFoundError(
            f"Need {resp_npz} for mc_pot_scale before writing {out_npz.name}"
        )
    resp = np.load(resp_npz, allow_pickle=True)
    resp_meta = json.loads(str(resp["meta_json"][0]))
    mc_scale = float(resp_meta["mc_pot_scale"])

    var_configs = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
    workers = max(1, int(args.workers))
    mc_files = _list_df_files(args.mc_dir, args.mc_fn)
    if args.max_files > 0:
        mc_files = mc_files[: args.max_files]
    if not mc_files:
        raise FileNotFoundError(f"No MC files under {args.mc_dir} matching {args.mc_fn}")
    print(f"Selection efficiency: {len(mc_files)} MC files, workers={workers}", flush=True)
    acc_by_var = _accumulate_fvcont(mc_files, var_configs, workers=workers)

    savez: Dict[str, Any] = {}
    meta = {
        "schema": "prl_productB_selection_efficiency_v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mc_dir": args.mc_dir,
        "mc_fn": args.mc_fn,
        "n_mc_files": len(mc_files),
        "mc_pot_scale": mc_scale,
        "response_npz": str(resp_npz),
        "signal_truth_fv": SIGNAL_TRUTH_FV_SELECTED,
        "truth_column": "var_nu_col",
        "note": (
            "nevts_all_fvcont: generated 1μ1p0π with per_tpc = full SBND_Gen1 "
            "(including high-YZ) on the true vertex and true μ/p start+end, "
            "and all those points in the same TPC. Matches reco event_contained_per_tpc. "
            "nevts_sel_fvcont: selected events with that same truth definition. "
            "Both use var_nu_col. Counts are unit weights times mc_pot_scale."
        ),
    }
    savez["meta_json"] = np.array([json.dumps(meta)], dtype=object)
    for vc in var_configs:
        vsn = vc.var_save_name
        acc = acc_by_var[vsn]
        all_h = acc["nevts_all_fvcont"] * mc_scale
        sel_h = acc["nevts_sel_fvcont"] * mc_scale
        savez[f"{vsn}::nevts_all_fvcont"] = all_h
        savez[f"{vsn}::nevts_sel_fvcont"] = sel_h
        savez[f"{vsn}::bins"] = np.asarray(vc.bins, dtype=float)
        savez[f"{vsn}::bin_centers"] = np.asarray(vc.bin_centers, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.divide(sel_h, all_h, out=np.zeros_like(sel_h), where=all_h > 0)
        print(
            f"  {vsn}: sel_eff_mean={float(np.nanmean(eff)):.4f}  "
            f"all_fvcont={float(np.sum(all_h)):.1f}  sel={float(np.sum(sel_h)):.1f}",
            flush=True,
        )
    np.savez_compressed(out_npz, **savez)
    print("wrote", out_npz, flush=True)
    return 0


def _accumulate_pass2(
    mc_files: Sequence[str],
    var_configs: Sequence[Any],
    workers: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    acc_by_var = {
        vc.var_save_name: _empty_acc(len(vc.bins) - 1) for vc in var_configs
    }
    var_specs = _var_spec_list(var_configs)
    if workers <= 1:
        for fp in tqdm(mc_files, desc="MC files"):
            part = _accumulate_file_worker(fp, var_specs)
            if part is not None:
                _merge_acc(acc_by_var, part)
            gc.collect()
        return acc_by_var

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(_accumulate_file_worker, fp, var_specs): fp for fp in mc_files
        }
        for fut in tqdm(as_completed(futs), total=len(futs), desc="MC files"):
            part = fut.result()
            if part is not None:
                _merge_acc(acc_by_var, part)
    return acc_by_var


def _load_overlay_counts(counts_npz: str) -> Dict[str, Dict[str, np.ndarray]]:
    if not path.isfile(counts_npz):
        raise FileNotFoundError(
            f"Missing {counts_npz}. Run data_mc_overlay_products.py --product B "
            f"--force-rebuild first (with OffBeam)."
        )
    z = np.load(counts_npz)
    out: Dict[str, Dict[str, np.ndarray]] = {}
    slugs = sorted({k.split("__", 1)[0] for k in z.files if "__" in k})
    for slug in slugs:
        pack: Dict[str, np.ndarray] = {}
        for field in (
            "bins",
            "data",
            "mc_signal",
            "mc_background",
            "mc_nu_background",
            "offbeam",
            "mc_total",
            "cosmic_estimate",
        ):
            key = f"{slug}__{field}"
            if key in z.files:
                pack[field] = np.asarray(z[key], dtype=float)
        if "data" in pack and "mc_background" in pack:
            pack["n_sel_data"] = pack["data"] - pack["mc_background"]
        out[slug] = pack
    return out


def _plot_efficiency(vc, eff, nevts_sel, nevts_allmc, save_name: str) -> None:
    bins = vc.bins
    centers = vc.bin_centers
    fig, ax = plt.subplots()
    ax_eff = ax.twinx()
    ax.hist(
        centers,
        bins=bins,
        weights=nevts_allmc,
        histtype="step",
        color="C0",
        label="All signal (mcnu)",
        linewidth=1.5,
    )
    ax.hist(
        centers,
        bins=bins,
        weights=nevts_sel,
        histtype="step",
        color="C1",
        label="Selected signal",
        linewidth=1.5,
    )
    ax_eff.errorbar(centers, eff, fmt="ko-", markersize=4, label="Efficiency")
    ax.set_xlabel(vc.var_labels[0])
    ax.set_ylabel("Events / Bin (POT-scaled)")
    ax_eff.set_ylabel("Efficiency")
    ymax = float(np.nanmax(eff)) if np.any(np.isfinite(eff)) else 0.2
    ax_eff.set_ylim(0, min(1.05, max(0.2, ymax * 1.3)))
    ax.set_xlim(bins[0], bins[-1])
    ax.set_title(f"Signal efficiency — {vc.var_save_name}")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_eff.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=10)
    fig.savefig(save_name + fig_ext, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.selection_eff_only:
        return _write_selection_efficiency(args)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "plots"
    makedirs(out_dir, exist_ok=True)
    makedirs(fig_dir, exist_ok=True)

    out_npz = out_dir / "response_matrices.npz"
    if out_npz.is_file() and not args.force:
        print(f"exists: {out_npz} (pass --force to overwrite)")
        return 0

    var_configs = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
    workers = max(1, int(args.workers))

    mc_files = _list_df_files(args.mc_dir, args.mc_fn)
    if args.max_files > 0:
        mc_files = mc_files[: args.max_files]
    if not mc_files:
        raise FileNotFoundError(f"No MC files under {args.mc_dir} matching {args.mc_fn}")
    print(f"MC files: {len(mc_files)} under {args.mc_dir}", flush=True)
    print(f"workers={workers}", flush=True)

    print("Pass 1: data beam-quality POT + n_evt_good...", flush=True)
    data_tot_pot, n_evt_good = _data_tot_pot_and_n()
    if n_evt_good != dmo.EXPECTED_DATA_N:
        raise RuntimeError(
            f"n_evt_good={n_evt_good} != expected {dmo.EXPECTED_DATA_N}"
        )
    print(f"  data_tot_pot={data_tot_pot:.6e}  n_evt_good={n_evt_good}", flush=True)

    print("Pass 1b: sum MC POT...", flush=True)
    mc_tot_pot = _sum_hdr_pot(mc_files, workers=workers)
    if mc_tot_pot <= 0:
        raise RuntimeError("mc_tot_pot <= 0")
    mc_scale = data_tot_pot / mc_tot_pot
    print(f"  mc_tot_pot={mc_tot_pot:.6e}  mc_scale={mc_scale:.6e}", flush=True)

    print(f"Pass 2: accumulate response ingredients (workers={workers})...", flush=True)
    acc_by_var = _accumulate_pass2(mc_files, var_configs, workers=workers)

    print("Load overlay counts for data/bkg...", flush=True)
    overlay = _load_overlay_counts(args.counts_npz)

    # Apply POT scale to absolute MC rates; build R, eff; merge overlay
    savez: Dict[str, Any] = {}
    meta = {
        "schema": "prl_productB_response_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mc_dir": args.mc_dir,
        "mc_fn": args.mc_fn,
        "n_mc_files": len(mc_files),
        "data_tot_pot": data_tot_pot,
        "mc_tot_pot": mc_tot_pot,
        "mc_pot_scale": mc_scale,
        "n_evt_good": n_evt_good,
        "counts_report": args.counts_npz,
        "cosmic_estimate": "offbeam",
        "batch_files": int(args.batch_files),
        "workers": workers,
        "signal_truth_fv_unfold": SIGNAL_TRUTH_FV_UNFOLD,
        "signal_truth_fv_selected": SIGNAL_TRUTH_FV_SELECTED,
        "note": (
            "nevts_allmc uses signal_truth_fv='none' (xsec signal); "
            "selected-signal migration uses per_tpc (matches reco FV)."
        ),
    }
    savez["meta_json"] = np.array([json.dumps(meta)], dtype=object)

    variables_out: Dict[str, Any] = {}
    for vc in var_configs:
        vsn = vc.var_save_name
        acc = acc_by_var[vsn]
        # scale absolute rates
        for k in (
            "nevts_allmc",
            "nevts_sel_truth",
            "nevts_sel_reco",
            "nevts_allsel_reco",
            "reco_vs_true",
        ):
            acc[k] = acc[k] * mc_scale

        eff = np.divide(
            acc["nevts_sel_truth"],
            acc["nevts_allmc"],
            out=np.zeros_like(acc["nevts_sel_truth"]),
            where=acc["nevts_allmc"] > 0,
        )
        response = get_response_matrix(acc["reco_vs_true"], eff)

        # overlay bookkeeping
        if vsn in overlay:
            n_data = overlay[vsn]["data"]
            n_bkg = overlay[vsn]["mc_background"]
            n_sel_data = overlay[vsn].get("n_sel_data", n_data - n_bkg)
        else:
            # ``integrated`` (and any other single-bin var) is not written by
            # overlay counts_report. Reconstruct by summing a differential
            # variable that covers the same selected sample.
            n_data = np.zeros_like(acc["nevts_sel_reco"])
            n_bkg = np.zeros_like(acc["nevts_sel_reco"])
            n_sel_data = np.zeros_like(acc["nevts_sel_reco"])
            donor = None
            for cand in ("muon-p", "proton-p", "tki-del_alpha", "muon-dir_z"):
                if cand in overlay and "data" in overlay[cand]:
                    donor = cand
                    break
            if donor is not None and n_data.size == 1:
                n_data = np.asarray(
                    [float(np.sum(overlay[donor]["data"]))], dtype=float
                )
                n_bkg = np.asarray(
                    [float(np.sum(overlay[donor]["mc_background"]))], dtype=float
                )
                n_sel_data = overlay[donor].get(
                    "n_sel_data", overlay[donor]["data"] - overlay[donor]["mc_background"]
                )
                n_sel_data = np.asarray([float(np.sum(n_sel_data))], dtype=float)
                print(
                    f"  {vsn}: not in counts_report; "
                    f"data/bkg from sum({donor}) "
                    f"(n_data={n_data[0]:.1f}, n_sel_data={n_sel_data[0]:.1f})",
                    flush=True,
                )
            else:
                print(
                    f"  WARN: {vsn} not in counts_report; data/bkg left zero",
                    flush=True,
                )

        pack = {
            "bins": np.asarray(vc.bins, dtype=float),
            "bin_centers": np.asarray(vc.bin_centers, dtype=float),
            "reco_vs_true": acc["reco_vs_true"],
            "eff": eff,
            "response": response,
            "nevts_allmc": acc["nevts_allmc"],
            "nevts_sel_truth": acc["nevts_sel_truth"],
            "nevts_sel_reco": acc["nevts_sel_reco"],
            "nevts_allsel_reco": acc["nevts_allsel_reco"],
            "n_data": np.asarray(n_data, dtype=float),
            "n_mc_bkg": np.asarray(n_bkg, dtype=float),
            "n_sel_data": np.asarray(n_sel_data, dtype=float),
            "var_labels": list(vc.var_labels),
        }
        if vsn in overlay and "offbeam" in overlay[vsn]:
            pack["offbeam"] = overlay[vsn]["offbeam"]
        if vsn in overlay and "mc_nu_background" in overlay[vsn]:
            pack["mc_nu_background"] = overlay[vsn]["mc_nu_background"]

        for key, arr in pack.items():
            if key == "var_labels":
                savez[f"{vsn}::var_labels"] = np.array(arr, dtype=object)
            else:
                savez[f"{vsn}::{key}"] = np.asarray(arr)

        variables_out[vsn] = {
            "eff_mean": float(np.nanmean(eff)),
            "nevts_allmc_sum": float(np.sum(acc["nevts_allmc"])),
            "n_sel_data_sum": float(np.sum(n_sel_data)),
        }
        print(
            f"  {vsn}: eff_mean={variables_out[vsn]['eff_mean']:.4f}  "
            f"allmc={variables_out[vsn]['nevts_allmc_sum']:.1f}  "
            f"n_sel_data={variables_out[vsn]['n_sel_data_sum']:.1f}",
            flush=True,
        )

        if not args.no_plots:
            _plot_efficiency(
                vc,
                eff,
                acc["nevts_sel_truth"],
                acc["nevts_allmc"],
                str(fig_dir / f"{vsn}-efficiency"),
            )
            if response.shape[0] > 1:
                plot_heatmap(
                    acc["reco_vs_true"],
                    vc.bins,
                    plot_labels=["True", "Reco", f"{vsn} migration"],
                    plot=False,
                    save_fig=True,
                    save_name=str(fig_dir / f"{vsn}-reco_vs_true"),
                )
                plot_heatmap(
                    response,
                    vc.bins,
                    plot_labels=["True", "Reco", f"{vsn} response"],
                    plot=False,
                    save_fig=True,
                    save_name=str(fig_dir / f"{vsn}-response"),
                )

    np.savez_compressed(out_npz, **savez)
    man = {
        "schema": "prl_productB_response_v1",
        "created_utc": meta["created_utc"],
        "npz_path": str(out_npz),
        "variables": sorted(variables_out),
        "meta": meta,
        "summary": variables_out,
        "key_format": "{var_save_name}::{field}",
    }
    man_path = out_dir / "response_matrices_manifest.json"
    with open(man_path, "w") as f:
        json.dump(man, f, indent=2)
    print("wrote", out_npz)
    print("wrote", man_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
