#!/usr/bin/env python
"""Aggregate per-map-shard pickles and produce all event-selection plots.

This is the **reduce** pass of the chunked histogram workflow (one pickle per
*(sample, CAF .df file)* — not an exposure batch in time; see ``exposure_access``).

  1. globs the per-(sample, shard) pickles under ``--in_dir``,
  2. sums histograms across chunks within each sample,
  3. merges sample-level results into one combined histogram dict,
  4. renders each plot through ``overlay_hists_from_histdata`` -- yielding the
     SAME visual output as the original notebook, but starting from the
     pre-binned content,
  5. produces the summary breakdown bar plot and the efficiency curve plots.

**Systematic uncertainty bands (Flux / G4 / GENIE universes)** — entirely compatible with
chunking: each chunk bin-wise histogram for universe ``i`` is a sum of independent event
contributions, so you accumulate ``mc_univ_hist[s][univ, cat, bin]`` exactly like ``mc_hist``,
apply the same global POT scaling, then form a fractional covariance from universe spread
around CV (see ``selection_framework.frac_cov_from_mc_univ_histdata``). Enable map-phase
``event_selection_chunk.py --mc-univ-syst Flux,G4,GENIE``. Aggregation draws syst bands from
those chunks by default; pass ``--no-overlay-syst-from-universes`` to omit them. This complements the notebook workflow that loads
precomputed ``cov_frac`` matrices from disk (``selected_events.ipynb`` / ``utils.get_syst_unc``
via ``--syst-disk-root`` / ``NUMUCC_SYST_DISK_ROOT``).

The plotting step is a thin layer on top of the existing ``overlay_hists`` and
``plot_efficiency`` routines in ``utils.py``; the new precomputed-histogram
entry point is ``overlay_hists_from_histdata``.

Usage
-----
    python event_selection_aggregate.py --in_dir AGG_INPUT_DIR \
                                         --out_dir PLOTS_OUTPUT_DIR

The expected pickle layout in ``--in_dir`` is::

    mc__<chunkname>.pkl
    data__<chunkname>.pkl
    intime__<chunkname>.pkl
    offbeam__<chunkname>.pkl
    dirt__<chunkname>.pkl

i.e. the names produced by ``event_selection_chunk.py``.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from os import path, makedirs
import pickle
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# turn off pandas chatter
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from analysis_village.numucc_1p0pi.event_selection_pipeline_def import (
    build_pipeline, EFFICIENCY_VARS,
)
from analysis_village.numucc_1p0pi.selection_framework import (
    BarBreakdown,
    ExposureTotals, aggregate_chunk_files, merge_samples,
    sanitize_merged_histdata_finite,
    apply_global_exposure_scales,
    frac_cov_from_mc_univ_histdata,
)
from analysis_village.numucc_1p0pi.utils import (
    overlay_hists_from_histdata,
    get_pot_str,
    fig_ext,
    dpi,
    add_approval_text,
    format_singlebin_plot,
    get_syst_unc,
)
from analysis_village.numucc_1p0pi.categories import (
    topology_labels, topology_colors,
    genie_mode_labels, genie_mode_colors,
    pdg_labels,
)
from pyanalib.stat_helpers import return_data_stat_err

# style sheet that the notebook uses
try:
    plt.style.use(path.join(path.dirname(__file__), "presentation.mplstyle"))
except Exception:
    try:
        plt.style.use("presentation.mplstyle")
    except Exception:
        pass


SAMPLES = ("mc", "data", "intime", "offbeam", "dirt")


# ===========================================================================
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in_dir", required=True, help="Dir holding per-chunk pickles")
    p.add_argument("--out_dir", required=True, help="Where to write the plots")
    p.add_argument("--data_pot", type=float, default=None,
                   help="Override data POT for legend (default: sum of chunk_pot from data pickles)")
    p.add_argument("--cosmic_estimate", choices=("intime", "offbeam"), default="intime",
                   help="Which cosmic sample drives the stacked cosmic component")
    p.add_argument("--hide_cosmic_model_unc", action="store_true",
                   help="Do not shade intime-vs-offbeam bin-wise spread")
    p.add_argument("--f_offbeam_frac", type=float, default=0.08,
                   help="Offbeam-coincident-with-BNB fraction used in cosmic gates scaling")
    p.add_argument("--skip_global_exposure", action="store_true",
                   help="Skip POT/gates rescaling (only for legacy chunks already scaled)")
    p.add_argument("--syst_tag", default="", help="Prefix tag for output names")
    p.add_argument("--save_fig", action="store_true", default=True,
                   help="Save figures to disk (default True)")
    p.add_argument("--show_fig", action="store_true", default=False,
                   help="Show figures interactively (off by default)")
    p.add_argument(
        "--no-overlay-syst-from-universes",
        action="store_true",
        help="Skip fractional covariance / hatched syst bands from chunked MC universe histograms "
             "(default is ON when chunks include mc_univ_hist). Chunks need "
             "event_selection_chunk.py --mc-univ-syst Flux,G4,GENIE matching HDF columns.",
    )
    p.add_argument(
        "--syst-disk-root",
        default=None,
        help="Root directory with MCstat/, Flux/, G4/, GENIE/, Cosmics/, Detector/ trees "
             "(see analysis_village.numucc_1p0pi.syst_disk_layout). Env: NUMUCC_SYST_DISK_ROOT.",
    )
    return p.parse_args()


# ===========================================================================
def collect_chunks(in_dir: str) -> Dict[str, List[str]]:
    """Group pickles by sample name (the prefix before ``__``)."""
    out = {s: sorted(glob.glob(path.join(in_dir, f"{s}__*.pkl"))) for s in SAMPLES}
    for s, files in out.items():
        print(f"[aggregate] sample={s} -> {len(files)} chunks")
    return out


# ===========================================================================
def accumulate_exposure_totals_from_dir(in_dir: str) -> ExposureTotals:
    """Sum POT / gates denominators from every per-chunk pickle (map phase metadata)."""
    totals = ExposureTotals()
    pattern = path.join(in_dir, "*__*.pkl")
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f"[aggregate] WARN: no chunk pickles matched {pattern}")
        return totals
    for cf in paths:
        with open(cf, "rb") as f:
            d = pickle.load(f)
        sample = d.get("sample")
        m = d.get("meta", {}) or {}
        if sample == "data":
            totals.data_pot += float(m.get("chunk_pot", 0.0))
            totals.data_gates_bnb += float(m.get("chunk_gates_bnb", 0.0))
        elif sample == "mc":
            totals.mc_pot += float(m.get("chunk_pot", 0.0))
        elif sample == "dirt":
            totals.dirt_pot += float(m.get("chunk_pot", 0.0))
        elif sample == "intime":
            totals.intime_gates += float(m.get("chunk_cosmic_gates_intime", 0.0))
        elif sample == "offbeam":
            totals.offbeam_gates += float(m.get("chunk_cosmic_gates_offbeam", 0.0))
    print(
        f"[aggregate] exposure totals: data_pot={totals.data_pot:.3e} "
        f"bnb_gates={totals.data_gates_bnb:.3e} mc_pot={totals.mc_pot:.3e} "
        f"dirt_pot={totals.dirt_pot:.3e} intime_gates={totals.intime_gates:.3e} "
        f"offbeam_gates={totals.offbeam_gates:.3e}"
    )
    return totals


# ===========================================================================
def render_overlay_plots(
    merged: dict,
    plot_label_map: dict,
    save_fig_dir: str,
    pot_str: str,
    save_fig: bool,
    show_fig: bool,
    cosmic_estimate: str,
    show_cosmic_model_unc: bool,
    overlay_syst_from_universes: bool = True,
    syst_disk_root: str | None = None,
):
    """Render every plot stored in ``merged['histdata']``."""
    # We need the pipeline definition to recover the per-plot kwargs and labels.
    pipeline = build_pipeline()
    spec_lookup = {}  # (stage_key, plot_key) -> PlotSpec (recomputed from build_pipeline)
    from analysis_village.numucc_1p0pi.selection_framework import ChunkRunner
    for stage in pipeline:
        for ps in stage.plots:
            key = (stage.key, ChunkRunner.plot_key(stage.key, ps))
            spec_lookup[key] = ps

    vars_missing_syst = []

    for key, hd in merged["histdata"].items():
        stage_key, plot_key = key
        ps = spec_lookup.get(key)
        if ps is None:
            print(f"[aggregate] WARN: no PlotSpec for {key}; skipping")
            continue

        # build plot labels
        if ps.plot_label_template is not None:
            plot_labels = [
                ps.plot_label_template[0],
                ps.plot_label_template[1].replace("{pot}", pot_str),
                ps.plot_label_template[2].replace("{pot}", pot_str) if len(ps.plot_label_template) >= 3 else "",
            ]
        else:
            plot_labels = [ps.var_config.var_labels[0], f"Events / Bin (POT={pot_str})", ""]

        save_name = path.join(save_fig_dir, f"{stage_key}__{ps.breakdown_type}__{ps.var_config.var_save_name}"
                              + (("_" + ps.name_suffix) if ps.name_suffix else ""))
        kwargs = dict(ps.save_kwargs)
        kwargs.setdefault("ratio", False)
        kwargs.setdefault("save_fig", save_fig)
        kwargs.setdefault("save_name", save_name)
        kwargs.setdefault("plot", show_fig)
        kwargs["plot_labels"] = plot_labels
        kwargs.setdefault("cosmic_estimate", cosmic_estimate)
        kwargs.setdefault("show_cosmic_model_unc", show_cosmic_model_unc)
        kwargs.setdefault("verbose_hist", (ps.name_suffix or "") == "final")
        # Match ``selected_events.ipynb``: combined syst as hatched band (not norm/shape/mixed fill).
        kwargs.setdefault("syst_decomp", False)

        if overlay_syst_from_universes and kwargs.get("syst") is None:
            fc = frac_cov_from_mc_univ_histdata(
                hd, cosmic_estimate=kwargs.get("cosmic_estimate", cosmic_estimate)
            )
            if fc is not None:
                kwargs["syst"] = fc

        # Notebook parity: precomputed fractional covariances on disk (GENIE / flux / …).
        if kwargs.get("syst") is None and syst_disk_root is not None:
            _, cov_disk = get_syst_unc(
                ps.var_config,
                syst_disk_root=syst_disk_root,
                skip_missing_vars=True,
            )
            if np.any(cov_disk):
                kwargs["syst"] = cov_disk

        if kwargs.get("syst") is None:
            vars_missing_syst.append(ps.var_config.var_save_name)

        try:
            overlay_hists_from_histdata(hd, var_config=ps.var_config, **kwargs)
        except Exception as e:
            print(f"[aggregate] WARN: plot {key} failed: {e}")
            plt.close('all')

    if vars_missing_syst:
        uniq = sorted(set(vars_missing_syst))
        print(
            f"[aggregate] overlay plots without syst covariance ({len(vars_missing_syst)} plots, "
            f"{len(uniq)} distinct var_save_name): {', '.join(uniq)}",
            flush=True,
        )


# ===========================================================================
def render_summary_breakdown_plot(merged: dict, save_fig_dir: str,
                                  save_fig: bool, show_fig: bool,
                                  cosmic_estimate: str):
    """Reproduce the cell-79 summary bar plot from the notebook."""
    bar = merged["bar"]
    stage_keys = merged["stage_keys"]
    stage_labels = merged["stage_labels"]

    # only stages that we actually filled the breakdown for
    avail_stages = [k for k in stage_keys if k in bar]
    if not avail_stages:
        print("[aggregate] no breakdown stages -> skipping summary plot")
        return
    stage_label_lookup = dict(zip(stage_keys, stage_labels))

    # Build the percentage matrices (stages x categories)
    def fractions(bb: BarBreakdown) -> np.ndarray:
        # in cuts order; reverse to match the labels/colors used by bar plot
        v = bb.mc_counts.astype(float)
        cosmic_part = bb.offbeam_count if cosmic_estimate == "offbeam" else bb.intime_count
        tot = v.sum() + cosmic_part + bb.dirt_count
        if tot <= 0:
            return np.zeros_like(v)
        return 100.0 * v / tot

    topo_data = np.array([fractions(bar[k]["topology"])[::-1] for k in avail_stages[::-1]])
    genie_data = np.array([fractions(bar[k]["genie"])[::-1] for k in avail_stages[::-1]])

    y = np.arange(len(avail_stages))
    bar_width = 0.3

    def stack_bars(ax, data, yoffset, colors, label):
        left = np.zeros(len(avail_stages))
        for i, color in enumerate(colors[:data.shape[1]]):
            ax.barh(y + yoffset, data[:, i], bar_width, left=left, color=color,
                    label=label if i == 0 else None)
            left += data[:, i]

    fig, ax = plt.subplots(figsize=(10, 10))
    stack_bars(ax, topo_data, -bar_width / 2, topology_colors, "Topology")
    stack_bars(ax, genie_data,  bar_width / 2, genie_mode_colors, "GENIE")

    if not np.any(topo_data > 0) and not np.any(genie_data > 0):
        print(
            "[aggregate] WARN: summary breakdown fractions are all zero "
            "(need MC chunk pickles with ``save_for_breakdown`` filled — "
            "mc_counts sum + cosmic gates term).",
            flush=True,
        )

    ax.set_xlabel("Percentage (%)")
    ax.set_yticks(y)
    ax.set_yticklabels([stage_label_lookup[k] for k in avail_stages[::-1]], fontsize=12)

    common_patches = [Patch(facecolor=c, label=l) for c, l in zip(
        ["gray", "sienna", "crimson", "darkgreen"],
        ["Cosmic", r"Out FV $\nu$", r"In FV other $\nu$", r"In FV $\nu_{\mu}$ NC"]
    )]
    genie_patches = [Patch(facecolor=c, label=l) for c, l in zip(
        ["#BFB17C", "#D88A3B", "#2c7c94", "#390C1E", "#9b5580"],
        [r"In FV $\nu_{\mu}$ CC Other", r"In FV $\nu_{\mu}$ CC SIS/DIS",
         r"In FV $\nu_{\mu}$ CC RES", r"In FV $\nu_{\mu}$ CC MEC", r"In FV $\nu_{\mu}$ CC QE"]
    )]
    topo_patches = [Patch(facecolor=c, label=l) for c, l in zip(
        ["coral", "darkslateblue", "mediumslateblue"],
        [r"In FV $\nu_{\mu}$ CC Other", r"In FV $\nu_{\mu}$ CC Np0$\pi$",
         r"In FV $\nu_{\mu}$ CC 1p0$\pi$"]
    )]
    ax.legend(handles=common_patches, loc='upper left', bbox_to_anchor=(0.01, 1.18),
              ncol=4, fontsize=12, frameon=False)
    for i, handles in enumerate([genie_patches[::-1], topo_patches[::-1]]):
        ax_i = ax.twinx()
        ax_i.legend(handles=handles, loc='upper left',
                    bbox_to_anchor=(0.01, 1.14 - 0.07*i),
                    ncol=3 if i == 0 else 4, fontsize=12, frameon=False)
        ax_i.set_yticks([])

    fig.tight_layout()
    if save_fig:
        plt.savefig(path.join(save_fig_dir, "event_selection_summary.png"),
                    dpi=300, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close()


# ===========================================================================
def render_efficiency_plots(merged: dict, save_fig_dir: str, pot_str: str,
                            save_fig: bool, show_fig: bool):
    """Reproduce ``plot_efficiency`` from the notebook using the eff accumulators."""
    from statsmodels.stats.proportion import proportion_confint
    eff = merged["eff"]
    if not eff:
        print("[aggregate] no efficiency data -> skipping")
        return
    stage_keys = merged["stage_keys"]
    stage_labels = merged["stage_labels"]
    stage_label_lookup = dict(zip(stage_keys, stage_labels))

    # group var-keyed accumulators by var
    vars_seen = set()
    for stage_key, by_v in eff.items():
        vars_seen.update(by_v.keys())

    var_lookup = {vc.var_save_name: vc for vc in EFFICIENCY_VARS}

    eff_dict = {}

    for var_save_name in sorted(vars_seen):
        var_config = var_lookup.get(var_save_name)

        if var_config is None:
            continue

        if var_config.var_save_name == "muon-dir_phi":
            continue

        # Denominator: generated signal on ``mcnu`` × ``var_nu_col`` (filled only on the
        # first efficiency stage, usually ``allreco``). Numerator: ``evt`` × ``var_evt_truth_col``.
        denom_stage = None
        for sk in stage_keys:
            if sk in eff and var_save_name in eff[sk]:
                denom_stage = sk
                break
        if denom_stage is None:
            continue
        denom = eff[denom_stage][var_save_name]
        n_tot_pot = np.asarray(
            getattr(denom, "n_truth_nu_pot", denom.n_signal_pot), dtype=float
        )
        n_tot_raw = np.asarray(
            getattr(denom, "n_truth_nu_raw", denom.n_signal_raw), dtype=float
        )
        # Legacy chunk pickles (pre mcnu denominator): fall back to evt-only first stage.
        if float(np.sum(n_tot_raw)) <= 0.0 and float(np.sum(n_tot_pot)) <= 0.0:
            n_tot_pot = np.asarray(denom.n_signal_pot, dtype=float)
            n_tot_raw = np.asarray(denom.n_signal_raw, dtype=float)
        denom_int_pot = float(np.sum(n_tot_pot))

        bins = var_config.bins
        bin_centers = 0.5 * (bins[:-1] + bins[1:])

        fig, ax = plt.subplots()
        ax_eff = ax.twinx()

        eff_list, eff_err_list = [], []
        ymax_hist = 0.0
        plot_idx = 0
        for stage_key in stage_keys:
            if stage_key not in eff or var_save_name not in eff[stage_key]:
                continue
            ea = eff[stage_key][var_save_name]
            n_pot = np.asarray(ea.n_signal_pot, dtype=float)
            n_int = ea.n_total_signal_int
            ymax_hist = max(ymax_hist, float(np.max(n_pot)) if n_pot.size else 0.0)
            with np.errstate(divide='ignore', invalid='ignore'):
                this_eff = np.where(n_tot_pot > 0, n_pot / n_tot_pot, 0.0)
            # raw counts for Wilson interval (avoid div-zero)
            n_succ = ea.n_signal_raw
            err_low, err_high = [], []
            for i in range(len(this_eff)):
                if n_tot_raw[i] > 0:
                    interval = proportion_confint(int(n_succ[i]), int(n_tot_raw[i]), method='wilson')
                    err_low.append(abs(this_eff[i] - interval[0]))
                    err_high.append(abs(this_eff[i] - interval[1]))
                else:
                    err_low.append(0.0)
                    err_high.append(0.0)
            this_eff_err = [err_low, err_high]
            eff_int_pct = (n_int / denom_int_pot * 100) if denom_int_pot > 0 else 0.0
            label = stage_label_lookup.get(stage_key, stage_key) + f" ({eff_int_pct:.2f}%)"
            color = plt.cm.tab10(plot_idx % 10)
            # Pre-binned spectra: draw as a step outline, matching the curve color.
            # (Closer to notebook ``histtype='step'``; avoid filled patches obscuring cuts.)
            ax.stairs(n_pot, bins, fill=False, alpha=0.5, color=color, linewidth=1.5, zorder=1)
            ax_eff.errorbar(bin_centers, this_eff, yerr=this_eff_err,
                            fmt="o-", color=color, label=label, markersize=4, zorder=4)
            eff_list.append(this_eff)
            eff_err_list.append(this_eff_err)
            plot_idx += 1

        ax.set_xlabel(var_config.var_labels[0])
        ax.set_ylabel("Events")
        ax_eff.set_ylabel("Efficiency")
        ax_eff.set_ylim(0, 1.05)
        ax.set_xlim(bins[0], bins[-1])
        if ymax_hist > 0:
            ax.set_ylim(0.0, ymax_hist * 1.15)
        ax_eff.set_zorder(3)
        ax_eff.patch.set_visible(False)

        fig.subplots_adjust(top=0.88)

        ax_eff.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=2,
            fontsize=8,
            frameon=False,
        )

        add_approval_text("internal", 0.98, 0.97, "right", fontsize=14)
        if var_config.var_save_name == "integrated":
            format_singlebin_plot()
        if save_fig:
            plt.savefig(path.join(save_fig_dir, f"efficiency-{var_save_name}{fig_ext}"),
                        bbox_inches="tight", dpi=dpi)
        if show_fig:
            plt.show()
        else:
            plt.close()

        eff_dict[var_save_name] = {"eff_list": eff_list, "eff_err_list": eff_err_list}

    # Integrated purity is identical for every efficiency variable (same event counts).
    for vs in sorted(vars_seen):
        if vs not in var_lookup:
            continue
        last_stage = next((s for s in reversed(stage_keys) if s in eff and vs in eff[s]), None)
        if last_stage is None:
            continue
        ea_last = eff[last_stage][vs]
        if ea_last.n_at_stage_int <= 0:
            continue
        purity = ea_last.n_total_signal_int / ea_last.n_at_stage_int * 100.0
        print(f"[aggregate] final selection purity: {purity:.2f}%", flush=True)
        break

    # Save the eff dict for downstream tools
    out_pkl = path.join(save_fig_dir, "eff_dict.pkl")
    with open(out_pkl, "wb") as f:
        pickle.dump(eff_dict, f)
    print(f"[aggregate] wrote {out_pkl}")


# ===========================================================================
def main():
    args = parse_args()
    makedirs(args.out_dir, exist_ok=True)
    save_fig_dir = args.out_dir

    # ---- discover chunk pickles by sample
    chunk_groups = collect_chunks(args.in_dir)

    # ---- aggregate per sample
    samples = {}
    for s, files in chunk_groups.items():
        if not files:
            continue
        print(f"[aggregate] aggregating {len(files)} chunks for sample={s}")
        samples[s] = aggregate_chunk_files(files)

    if not samples:
        print("[aggregate] no chunks found; nothing to do")
        return

    # ---- merge across samples
    merged = merge_samples(samples)
    n_hd_fixed, n_bar_fixed = sanitize_merged_histdata_finite(merged)
    if n_hd_fixed or n_bar_fixed:
        print(
            f"[aggregate] sanitized NaN/inf histogram bins (legacy weights): "
            f"{n_hd_fixed} OverlayHistData keys, {n_bar_fixed} bar breakdown rows",
            flush=True,
        )

    # ---- exposure denominators from chunk metadata (every pickle), then global scales
    totals = accumulate_exposure_totals_from_dir(args.in_dir)
    exposure_scales = None
    if not args.skip_global_exposure:
        exposure_scales = apply_global_exposure_scales(
            merged, totals, f_offbeam_coincident=args.f_offbeam_frac
        )
        print(f"[aggregate] applied global scales: {exposure_scales}")
    else:
        print("[aggregate] --skip_global_exposure: histograms left as in pickles")

    # ---- POT string for axis labels / legend
    if args.data_pot is not None:
        data_pot = args.data_pot
    else:
        data_pot = totals.data_pot
        if data_pot <= 0:
            data_pot = 1.0
    pot_str = get_pot_str(data_pot)
    print(f"[aggregate] data_pot (legend)={data_pot:.3e} -> POT label={pot_str}")

    show_cosmic_unc = not args.hide_cosmic_model_unc

    # ---- render
    overlay_syst = not args.no_overlay_syst_from_universes
    if overlay_syst:
        print(
            "[aggregate] overlay MC syst bands from chunked universe histograms when available "
            "(flux/G4/GENIE weights on MC slices; cosmic/dirt held fixed per univ). "
            "Use --no-overlay-syst-from-universes to disable.",
            flush=True,
        )
    render_overlay_plots(
        merged, plot_label_map={}, save_fig_dir=save_fig_dir,
        pot_str=pot_str, save_fig=args.save_fig, show_fig=args.show_fig,
        cosmic_estimate=args.cosmic_estimate,
        show_cosmic_model_unc=show_cosmic_unc,
        overlay_syst_from_universes=overlay_syst,
        syst_disk_root=args.syst_disk_root or os.environ.get("NUMUCC_SYST_DISK_ROOT"),
    )
    render_summary_breakdown_plot(
        merged, save_fig_dir,
        save_fig=args.save_fig, show_fig=args.show_fig,
        cosmic_estimate=args.cosmic_estimate,
    )
    render_efficiency_plots(merged, save_fig_dir, pot_str,
                            save_fig=args.save_fig, show_fig=args.show_fig)

    # save the merged dict so downstream stuff can pull histdata directly
    out_pkl = path.join(save_fig_dir, "merged_histdata.pkl")
    with open(out_pkl, "wb") as f:
        pickle.dump({
            "merged": merged,
            "data_pot": data_pot,
            "pot_str": pot_str,
            "exposure_totals": totals,
            "exposure_scales": exposure_scales,
            "cosmic_estimate": args.cosmic_estimate,
            "f_offbeam_frac": args.f_offbeam_frac,
        }, f)
    print(f"[aggregate] wrote {out_pkl}")


if __name__ == "__main__":
    main()
