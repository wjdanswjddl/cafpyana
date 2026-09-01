#!/usr/bin/env python
"""Render event-selection overlay plots with slim sel_all GENIE uncertainties.

Uses pre-merged histograms (``merged_histdata.pkl`` or fresh aggregate from
batch pickles) and fractional covariances from the slim GENIE syst disk produced
by ``run_prl_genie_sel_all.sh`` (``PRL_genie_sel_all/syst_disk/GENIE/``).

Cut-stage plots use the combined ``genie_rate`` matrix; final-stage plots use
``genie_rate`` as well (event-count overlays). Chi2 plots map ``chi2_mu`` /
``chi2_p`` (plane I2) to ``chi2_mu_I2__at_<stage>`` / ``chi2_p_I2__at_<stage>``.
"""
from __future__ import annotations

import argparse
import pickle
import sys
from os import makedirs, path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_pipeline
from analysis_village.numucc_1p0pi.selection_framework import ChunkRunner
from analysis_village.numucc_1p0pi.scripts import event_selection_aggregate as agg
from analysis_village.numucc_1p0pi.syst_disk_layout import FILE_GENIE, SUB_GENIE

_DEFAULT_WORK = "/exp/sbnd/data/users/munjung/PRL_data/PRL_genie_sel_all"
_DEFAULT_BATCHES = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/event_selection-batched-20260525/batches"
)
_DEFAULT_MERGED = (
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/event_selection-batched-20260525/plots/merged_histdata.pkl"
)


def resolve_genie_slug(stage_key: str, var_save_name: str, name_suffix: str | None) -> tuple[str, str]:
    """Map a PlotSpec to ``(slug, cov_key)`` in ``cov_mat_dict.pkl``."""
    if name_suffix == "final":
        return var_save_name, "genie_rate"
    if var_save_name == "chi2_mu":
        return f"chi2_mu_I2__at_{stage_key}", "genie_rate"
    if var_save_name == "chi2_p":
        return f"chi2_p_I2__at_{stage_key}", "genie_rate"
    return var_save_name, "genie_rate"


def align_cov_frac_diag(cov: np.ndarray, n_target: int) -> np.ndarray:
    """Resize a fractional covariance to ``n_target`` bins (diagonal-only band display)."""
    cov = np.asarray(cov, dtype=np.float64)
    n_src = cov.shape[0]
    if n_src == n_target:
        return cov
    if n_target <= 0:
        return cov
    diag_frac = np.sqrt(np.maximum(np.diag(cov), 0.0))
    x_src = (np.arange(n_src, dtype=np.float64) + 0.5) / n_src
    x_tgt = (np.arange(n_target, dtype=np.float64) + 0.5) / n_target
    frac_tgt = np.interp(x_tgt, x_src, diag_frac)
    return np.diag(frac_tgt ** 2)


def load_slim_genie_cov_mat(genie_pkl: str) -> dict:
    with open(genie_pkl, "rb") as f:
        return pickle.load(f)


def make_slim_genie_cov_loader(genie_blob: dict, merged: dict):
    """Return ``syst_cov_loader(ps, stage_key)`` for :func:`render_overlay_plots`."""

    def loader(ps, stage_key: str):
        slug, cov_key = resolve_genie_slug(
            stage_key, ps.var_config.var_save_name, ps.name_suffix
        )
        row = genie_blob.get(slug)
        if row is None:
            print(
                f"[slim-genie-plots] no GENIE row for slug={slug!r} "
                f"(stage={stage_key}, var={ps.var_config.var_save_name!r})",
                flush=True,
            )
            return None
        if cov_key not in row:
            print(
                f"[slim-genie-plots] slug={slug!r} missing key {cov_key!r} "
                f"(have {sorted(row.keys())[:6]}…)",
                flush=True,
            )
            return None
        cov = np.asarray(row[cov_key], dtype=np.float64)
        plot_key = (stage_key, ChunkRunner.plot_key(stage_key, ps))
        hd = merged.get("histdata", {}).get(plot_key)
        n_bins = len(hd.bins) - 1 if hd is not None else len(ps.var_config.bin_centers)
        cov = align_cov_frac_diag(cov, n_bins)
        if cov.shape != (n_bins, n_bins):
            print(
                f"[slim-genie-plots] could not align slug={slug!r} to {n_bins} bins",
                flush=True,
            )
            return None
        return cov

    return loader


def load_merged_payload(merged_pkl: str) -> tuple[dict, str, float]:
    with open(merged_pkl, "rb") as f:
        payload = pickle.load(f)
    merged = payload["merged"]
    pot_str = payload.get("pot_str") or agg.get_pot_str(payload.get("data_pot", 1.0))
    data_pot = float(payload.get("data_pot", 1.0))
    return merged, pot_str, data_pot


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out-dir",
        default=path.join(_DEFAULT_WORK, "plots", "event_selection_slim"),
        help="Output directory for PNGs (default: PRL_genie_sel_all/plots/event_selection_slim)",
    )
    p.add_argument(
        "--genie-cov-pkl",
        default=path.join(_DEFAULT_WORK, "syst_disk", SUB_GENIE, FILE_GENIE),
        help="Slim GENIE cov_mat_dict.pkl from syst_genie_aggregate (slim group only)",
    )
    p.add_argument(
        "--merged-histdata",
        default=_DEFAULT_MERGED,
        help="Pre-aggregated merged_histdata.pkl (skips batch aggregate when set)",
    )
    p.add_argument(
        "--in-dir",
        default=None,
        help="Batch pickle directory (mc__*.pkl, …). Used when --merged-histdata is absent.",
    )
    p.add_argument(
        "--cosmic-estimate",
        default="intime",
        choices=("intime", "offbeam"),
    )
    p.add_argument("--show-fig", action="store_true", default=False)
    p.add_argument("--skip-summary", action="store_true", help="Skip event_selection_summary.png")
    p.add_argument("--skip-efficiency", action="store_true", help="Skip efficiency-*.png plots")
    return p.parse_args()


def main():
    args = parse_args()
    makedirs(args.out_dir, exist_ok=True)

    try:
        plt.style.use(path.join(path.dirname(__file__), "presentation.mplstyle"))
    except Exception:
        pass

    if args.merged_histdata and path.isfile(args.merged_histdata):
        print(f"[slim-genie-plots] loading merged histdata from {args.merged_histdata}", flush=True)
        merged, pot_str, _data_pot = load_merged_payload(args.merged_histdata)
    elif args.in_dir:
        print(f"[slim-genie-plots] aggregating batches from {args.in_dir}", flush=True)
        chunk_groups = agg.collect_chunks(args.in_dir)
        samples = {}
        for s, files in chunk_groups.items():
            if files:
                samples[s] = agg.aggregate_chunk_files(files)
        if not samples:
            raise SystemExit(f"[slim-genie-plots] no batch pickles under {args.in_dir}")
        merged = agg.merge_samples(samples)
        agg.sanitize_merged_histdata_finite(merged)
        totals = agg.accumulate_exposure_totals_from_dir(args.in_dir)
        agg.apply_global_exposure_scales(merged, totals, f_offbeam_coincident=0.08)
        data_pot = totals.data_pot if totals.data_pot > 0 else 1.0
        pot_str = agg.get_pot_str(data_pot)
    else:
        raise SystemExit(
            "[slim-genie-plots] need --merged-histdata or --in-dir "
            f"(default merged path missing: {args.merged_histdata})"
        )

    genie_pkl = path.abspath(path.expanduser(args.genie_cov_pkl))
    if not path.isfile(genie_pkl):
        raise SystemExit(f"[slim-genie-plots] GENIE pickle not found: {genie_pkl}")
    print(f"[slim-genie-plots] GENIE cov: {genie_pkl}", flush=True)
    genie_blob = load_slim_genie_cov_mat(genie_pkl)
    cov_loader = make_slim_genie_cov_loader(genie_blob, merged)

    n_plots = sum(len(st.plots) for st in build_pipeline())
    print(f"[slim-genie-plots] rendering {n_plots} overlay plots -> {args.out_dir}", flush=True)

    agg.render_overlay_plots(
        merged,
        plot_label_map={},
        save_fig_dir=args.out_dir,
        pot_str=pot_str,
        save_fig=True,
        show_fig=args.show_fig,
        syst_disk_root=None,
        syst_cov_loader=cov_loader,
    )

    if not args.skip_summary:
        agg.render_summary_breakdown_plot(
            merged, args.out_dir, save_fig=True, show_fig=args.show_fig,
            cosmic_estimate=args.cosmic_estimate,
        )

    if not args.skip_efficiency:
        agg.render_efficiency_plots(
            merged, args.out_dir, pot_str, save_fig=True, show_fig=args.show_fig,
        )

    out_pkl = path.join(args.out_dir, "merged_histdata.pkl")
    with open(out_pkl, "wb") as f:
        pickle.dump({"merged": merged, "pot_str": pot_str, "genie_cov_pkl": genie_pkl}, f)
    print(f"[slim-genie-plots] wrote {out_pkl}", flush=True)
    print(f"[slim-genie-plots] DONE -> {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
