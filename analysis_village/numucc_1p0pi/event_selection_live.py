"""Live / interactive batched event selection for notebooks.

Processes ``.df`` files in chunks of ``files_per_job`` (default = plot-update
cadence), concatenating each chunk like legacy ``dfs_from_dir`` / ``n_max_concat``,
then refreshes a large multi-panel figure as statistics accumulate.
"""
from __future__ import annotations

import gc
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from analysis_village.numucc_1p0pi.categories import (
    genie_mode_colors,
    genie_mode_labels,
    pdg_labels,
    topology_colors,
    topology_labels,
)
from analysis_village.numucc_1p0pi.event_selection_batch_core import run_batch_selection
from analysis_village.numucc_1p0pi.event_selection_batched import (
    BatchJob,
    EventSelectionBatchedConfig,
    FileRecord,
    SAMPLES,
    default_event_selection_batched_work_root,
    print_survey_summary,
    survey_files,
    survey_files_from_dirs,
    write_manifest,
)
from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_pipeline
from analysis_village.numucc_1p0pi.selection_framework import ChunkRunner
from analysis_village.numucc_1p0pi.utils import get_pot_str, pdg_colors

# Aggregate helpers live in the scripts/ module (same as run_aggregate).
import importlib.util

_AGG_PATH = Path(__file__).resolve().parent / "scripts" / "event_selection_aggregate.py"


def _load_aggregate_module():
    # Importing the aggregate script sets MPLBACKEND=Agg; keep notebook backend intact.
    prev_backend = plt.get_backend()
    spec = importlib.util.spec_from_file_location("event_selection_aggregate_live", _AGG_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {_AGG_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    try:
        if plt.get_backend() != prev_backend:
            plt.switch_backend(prev_backend)
    except Exception:
        pass
    return mod


@dataclass
class LiveAccumulateConfig:
    """Notebook-friendly knobs for the live accumulating run."""

    max_files_per_sample: int | None = 2
    """Cap ``.df`` files per sample. ``None`` = all files."""

    update_every_round: bool = True
    """If False, only draw the panel once at the end."""

    update_every_n_files: int = 20
    """Refresh the live panel every N processed files (always also on the last file)."""

    files_per_job: int | None = None
    """How many ``.df`` files to concat+process per map job.

    ``None`` → use ``update_every_n_files`` (legacy-style ``n_max_concat`` batches).
    Set to ``1`` for one-file-at-a-time loading.
    """

    one_file_per_job: bool = False
    """Deprecated alias: if True, forces ``files_per_job=1``."""

    concat_load: bool = True
    """Concat files in a job before one pipeline pass (legacy ``dfs_from_dir`` style)."""

    # Notebook file layout (same as event_selection.ipynb config cell).
    base_dir: Path | str | None = None
    sample_dirs: Dict[str, str] | None = None
    filename_str: str = "sel_all"

    f_offbeam_frac: float = 0.08
    cosmic_estimate: str = "intime"
    use_mc_genweight: bool = False
    mc_univ_syst: Sequence[str] = ()
    skip_existing: bool = False
    work_base: Path | str | None = None
    batches_dir: Path | str | None = None
    plots_dir: Path | str | None = None
    samples: Sequence[str] = field(default_factory=lambda: tuple(SAMPLES))
    figsize: Tuple[float, float] = (28, 40)
    ncols: int = 5
    include_summary_bar: bool = True
    include_efficiency: bool = True
    show_in_notebook: bool = True
    """Inline-update the accumulation panel in the notebook cell."""

    save_merged_payload: bool = True
    """Always write ``merged_histdata.pkl`` (reloadable counts for re-plotting)."""

    save_final_panel: bool = False
    """If True, also write ``live_accumulation_panel.png`` under ``plots_dir``."""

    panel_dpi: int = 120
    legend_fontsize: float = 6.5
    trace: bool = False


@dataclass
class LiveAccumulateResult:
    work_base: Path
    batches_dir: Path
    plots_dir: Path
    merged_payload: dict | None
    pot_str: str
    data_pot: float
    n_files_done: Dict[str, int]
    failed: List[Tuple[str, str]]
    panel_path: Path | None
    payload_path: Path | None = None


def _jobs_chunked(
    records: Sequence[FileRecord],
    *,
    files_per_job: int,
) -> List[BatchJob]:
    """Pack sorted files into jobs of at most ``files_per_job`` per sample."""
    n = max(1, int(files_per_job))
    by_sample: Dict[str, List[FileRecord]] = {}
    for rec in records:
        by_sample.setdefault(rec.sample, []).append(rec)
    jobs: List[BatchJob] = []
    for sample in sorted(by_sample):
        files = sorted(by_sample[sample], key=lambda r: r.path)
        job_id = 0
        for i in range(0, len(files), n):
            chunk = files[i : i + n]
            jobs.append(
                BatchJob(
                    sample=sample,
                    job_id=job_id,
                    files=[r.path for r in chunk],
                    total_bytes=sum(r.size_bytes for r in chunk),
                )
            )
            job_id += 1
    return jobs


def _jobs_one_file_each(records: Sequence[FileRecord]) -> List[BatchJob]:
    return _jobs_chunked(records, files_per_job=1)


def _resolve_files_per_job(cfg: LiveAccumulateConfig) -> int:
    if cfg.one_file_per_job:
        return 1
    if cfg.files_per_job is not None:
        return max(1, int(cfg.files_per_job))
    return max(1, int(cfg.update_every_n_files or 1))


def _batch_out_path(batches_dir: Path, sample: str, job: BatchJob) -> Path:
    return batches_dir / f"{sample}__{job.tag}.pkl"


def process_job_inprocess(
    job: BatchJob,
    batches_dir: Path,
    *,
    cfg: LiveAccumulateConfig,
) -> Path:
    """Run selection on one job in-process (no subprocess)."""
    batches_dir.mkdir(parents=True, exist_ok=True)
    out_pkl = _batch_out_path(batches_dir, job.sample, job)
    if cfg.skip_existing and out_pkl.is_file():
        return out_pkl

    mc_univ = tuple(cfg.mc_univ_syst) if job.sample == "mc" else ()
    pipeline_trace = (lambda msg: print(msg, flush=True)) if cfg.trace else None
    run_batch_selection(
        job.sample,
        job.files,
        str(out_pkl),
        job_id=job.tag,
        use_mc_genweight=cfg.use_mc_genweight,
        mc_univ_syst_tags=mc_univ,
        concat_load=cfg.concat_load,
        pipeline_trace=pipeline_trace,
    )
    return out_pkl


def aggregate_batches_so_far(
    batches_dir: Path,
    *,
    f_offbeam_frac: float = 0.08,
    quiet: bool = True,
) -> Tuple[dict, float, str, dict, object]:
    """Merge all pickles currently under ``batches_dir`` and apply exposure scales.

    Returns ``(merged, data_pot, pot_str, scales, totals)``. When ``quiet=True``,
    suppresses aggregate chunk-count spam (progress bar already covers that).
    """
    import contextlib
    import io

    agg = _load_aggregate_module()
    sink = io.StringIO() if quiet else None
    ctx = contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()
    with ctx:
        chunk_groups = agg.collect_chunks(str(batches_dir))
        samples = {}
        for s, files in chunk_groups.items():
            if not files:
                continue
            samples[s] = agg.aggregate_chunk_files(files)
        if not samples:
            raise RuntimeError(f"No batch pickles under {batches_dir}")

        merged = agg.merge_samples(samples)
        agg.sanitize_merged_histdata_finite(merged)
        totals = agg.accumulate_exposure_totals(chunk_groups)
        scales = agg.apply_global_exposure_scales(
            merged, totals, f_offbeam_coincident=f_offbeam_frac
        )
    data_pot = totals.data_pot if totals.data_pot > 0 else 1.0
    pot_str = get_pot_str(data_pot)
    return merged, data_pot, pot_str, scales, totals


def _breakdown_style(breakdown_type: str) -> Tuple[List[str], List[str]]:
    """Signal-first labels/colors (same starting point as ``overlay_hists_from_histdata``)."""
    if breakdown_type == "pdg":
        return list(pdg_labels), list(pdg_colors)
    if breakdown_type == "topology":
        return list(topology_labels), list(topology_colors)
    if breakdown_type == "genie":
        return list(genie_mode_labels), list(genie_mode_colors)
    return [f"cat{i}" for i in range(8)], [f"C{i}" for i in range(8)]


def _draw_compact_overlay(
    ax,
    hd,
    *,
    xlabel: str,
    title: str,
    pot_str: str,
    vlines=None,
    show_legend: bool = True,
    legend_fontsize: float = 6.5,
) -> None:
    """Stack MC (+ cosmics folded into cat-0, dirt) like ``overlay_hists_from_histdata``."""
    from matplotlib.patches import Patch

    bins = np.asarray(hd.bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    widths = np.diff(bins)
    labels, colors = _breakdown_style(hd.breakdown_type)

    weights: List[np.ndarray] = []
    if hd.mc_hist is not None and (
        hd.has_mc
        or float(np.sum(hd.mc_hist)) != 0.0
        or hd.has_intime
        or getattr(hd, "has_offbeam", False)
        or hd.has_dirt
    ):
        weights = [np.asarray(hd.mc_hist[i], dtype=float).copy() for i in range(hd.mc_hist.shape[0])]

    cosmic = None
    if hd.has_intime and hd.intime_hist is not None:
        cosmic = np.asarray(hd.intime_hist, dtype=float)
    elif getattr(hd, "has_offbeam", False) and hd.offbeam_hist is not None:
        cosmic = np.asarray(hd.offbeam_hist, dtype=float)
    if cosmic is not None and weights:
        weights[0] = weights[0] + cosmic

    plot_labels = list(labels)
    plot_colors = list(colors)
    if hd.has_dirt and hd.dirt_hist is not None and weights:
        dirt = np.asarray(hd.dirt_hist, dtype=float)
        weights = [dirt] + weights
        plot_colors = plot_colors + ["black"]
        plot_labels = plot_labels + ["Dirt"]

    plot_colors = plot_colors[::-1]
    plot_labels = plot_labels[::-1]

    layer_integrals = [float(np.sum(np.asarray(w, dtype=float))) for w in weights]
    tot_int = float(sum(layer_integrals))
    fracs = [li / tot_int if tot_int > 0 else 0.0 for li in layer_integrals]

    bottom = np.zeros(len(centers), dtype=float)
    legend_handles = []
    for w, col, lab, frac in zip(weights, plot_colors, plot_labels, fracs):
        vals = np.asarray(w, dtype=float)
        if float(np.sum(vals)) == 0.0:
            bottom = bottom + vals
            continue
        ax.bar(
            centers,
            vals,
            width=widths,
            bottom=bottom,
            color=col,
            align="center",
            linewidth=0,
            zorder=2,
        )
        legend_handles.append(
            Patch(facecolor=col, edgecolor="none", label=f"{lab} ({frac * 100:.1f}%)")
        )
        bottom = bottom + vals

    n_data = 0.0
    if hd.has_data and hd.data_hist is not None:
        data = np.asarray(hd.data_hist, dtype=float)
        n_data = float(np.sum(data))
        err = np.sqrt(np.maximum(np.asarray(hd.data_err2, dtype=float), 0.0))
        ax.errorbar(
            centers,
            data,
            yerr=err,
            fmt="o",
            color="black",
            markersize=2.5,
            elinewidth=0.8,
            capsize=0,
            zorder=5,
            label=f"Data (N={n_data:.0f})",
        )
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="black",
                linestyle="None",
                markersize=3,
                label=f"Data (N={n_data:.0f})",
            )
        )

    if vlines:
        for item in vlines:
            try:
                x = float(item[0]) if isinstance(item, (list, tuple)) else float(item)
            except Exception:
                continue
            ax.axvline(x, color="crimson", ls="--", lw=0.9, alpha=0.85)

    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_title(title, fontsize=8, pad=2)
    ax.tick_params(labelsize=6)
    ymax = float(np.nanmax(bottom)) if bottom.size else 0.0
    if hd.has_data and hd.data_hist is not None and len(hd.data_hist):
        ymax = max(ymax, float(np.nanmax(hd.data_hist)))
    if ymax > 0:
        ax.set_ylim(0, ymax * 1.25)
    if show_legend and legend_handles:
        ax.legend(
            handles=legend_handles,
            fontsize=legend_fontsize,
            loc="upper right",
            frameon=False,
            borderpad=0.2,
            handlelength=0.9,
            labelspacing=0.2,
        )


def _draw_summary_bar(
    ax,
    merged: dict,
    cosmic_estimate: str = "intime",
    *,
    legend_fontsize: float = 6.5,
) -> None:
    """Per-stage topology composition including cosmics + dirt (bars sum to 100%)."""

    bar = merged.get("bar") or {}
    stage_keys = list(merged.get("stage_keys") or [])
    stage_labels = list(merged.get("stage_labels") or [])
    label_lookup = dict(zip(stage_keys, stage_labels))
    avail = [k for k in stage_keys if k in bar]
    if not avail:
        ax.set_axis_off()
        ax.set_title("summary (no stages yet)", fontsize=8)
        return

    n_topo = len(topology_colors)
    rows = []
    final_fracs = None
    for k in avail[::-1]:
        bb = bar[k]["topology"]
        v = np.asarray(bb.mc_counts, dtype=float).copy()
        cosmic = float(
            bb.offbeam_count if cosmic_estimate == "offbeam" else bb.intime_count
        )
        dirt = float(bb.dirt_count)
        if len(v) > 0:
            v[0] += cosmic
        tot = float(v.sum() + dirt)
        if tot <= 0:
            rows.append(np.zeros(n_topo + 1, dtype=float))
            continue
        mc_disp = (100.0 * v / tot)[::-1]
        if len(mc_disp) < n_topo:
            mc_disp = np.pad(mc_disp, (0, n_topo - len(mc_disp)))
        elif len(mc_disp) > n_topo:
            mc_disp = mc_disp[:n_topo]
        dirt_frac = 100.0 * dirt / tot
        row = np.concatenate([mc_disp, [dirt_frac]])
        rows.append(row)
        # First plotted row is the latest stage (avail[::-1][0] == avail[-1]).
        if final_fracs is None:
            final_fracs = row

    data = np.asarray(rows, dtype=float)
    y = np.arange(len(avail))
    left = np.zeros(len(avail))
    stack_colors = list(topology_colors) + ["sienna"]
    stack_labels = list(topology_labels) + ["Dirt"]
    for i, (color, lab) in enumerate(zip(stack_colors, stack_labels)):
        if i >= data.shape[1]:
            break
        vals = data[:, i]
        if float(np.sum(vals)) == 0.0:
            left = left + vals
            continue
        pct = float(final_fracs[i]) if final_fracs is not None else float(np.mean(vals))
        ax.barh(
            y,
            vals,
            left=left,
            color=color,
            height=0.7,
            label=f"{lab} ({pct:.1f}%)",
        )
        left = left + vals

    # Data count at final stage (if recorded).
    final_key = avail[-1]
    n_data_final = float(getattr(bar[final_key]["topology"], "data_count", 0.0) or 0.0)
    title = "Stage topology fractions"
    if n_data_final > 0:
        title += f"  |  final data N={n_data_final:.0f}"

    ax.set_yticks(y)
    ax.set_yticklabels([label_lookup.get(k, k) for k in avail[::-1]], fontsize=6)
    ax.set_xlabel("Percentage (%)", fontsize=7)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=6)
    ax.set_xlim(0, 100)
    ax.legend(
        fontsize=legend_fontsize,
        loc="lower right",
        frameon=False,
        borderpad=0.2,
        handlelength=0.9,
        labelspacing=0.2,
    )


def ordered_plot_keys(merged: dict) -> List[Tuple[Tuple[str, str], object]]:
    """Pipeline order for histdata keys present in ``merged``."""
    pipeline = build_pipeline()
    histdata = merged.get("histdata") or {}
    ordered: List[Tuple[Tuple[str, str], object]] = []
    seen = set()
    for stage in pipeline:
        for ps in stage.plots:
            key = (stage.key, ChunkRunner.plot_key(stage.key, ps))
            if key in histdata and key not in seen:
                ordered.append((key, ps))
                seen.add(key)
    for key, hd in histdata.items():
        if key not in seen:
            ordered.append((key, None))
            seen.add(key)
    return ordered


def ordered_efficiency_vars(merged: dict) -> List[object]:
    """Efficiency VariableConfigs present in ``merged['eff']`` (pipeline order)."""
    from analysis_village.numucc_1p0pi.event_selection_pipeline_def import EFFICIENCY_VARS

    eff = merged.get("eff") or {}
    vars_seen = set()
    for by_v in eff.values():
        vars_seen.update(by_v.keys())
    out = []
    for vc in EFFICIENCY_VARS:
        if vc.var_save_name == "muon-dir_phi":
            continue
        if vc.var_save_name in vars_seen:
            out.append(vc)
    return out


def _draw_efficiency_curve(
    ax,
    merged: dict,
    var_config,
    *,
    legend_fontsize: float = 6.5,
) -> None:
    """Compact per-stage efficiency curves for one variable (integrated % in legend)."""
    try:
        from statsmodels.stats.proportion import proportion_confint
    except ImportError:
        proportion_confint = None

    eff = merged.get("eff") or {}
    stage_keys = list(merged.get("stage_keys") or [])
    stage_labels = list(merged.get("stage_labels") or [])
    label_lookup = dict(zip(stage_keys, stage_labels))
    var_save_name = var_config.var_save_name

    denom_stage = None
    for sk in stage_keys:
        if sk in eff and var_save_name in eff[sk]:
            denom_stage = sk
            break
    if denom_stage is None:
        ax.set_axis_off()
        ax.set_title(f"eff: {var_save_name} (no data)", fontsize=8)
        return

    denom = eff[denom_stage][var_save_name]
    n_tot_pot = np.asarray(getattr(denom, "n_truth_nu_pot", denom.n_signal_pot), dtype=float)
    n_tot_raw = np.asarray(getattr(denom, "n_truth_nu_raw", denom.n_signal_raw), dtype=float)
    if float(np.sum(n_tot_raw)) <= 0.0 and float(np.sum(n_tot_pot)) <= 0.0:
        n_tot_pot = np.asarray(denom.n_signal_pot, dtype=float)
        n_tot_raw = np.asarray(denom.n_signal_raw, dtype=float)
    denom_int_pot = float(np.sum(n_tot_pot))

    bins = np.asarray(var_config.bins, dtype=float)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    ax_eff = ax.twinx()

    plot_idx = 0
    final_eff_pct = None
    for stage_key in stage_keys:
        if stage_key not in eff or var_save_name not in eff[stage_key]:
            continue
        ea = eff[stage_key][var_save_name]
        n_pot = np.asarray(ea.n_signal_pot, dtype=float)
        n_int = float(ea.n_total_signal_int)
        with np.errstate(divide="ignore", invalid="ignore"):
            this_eff = np.where(n_tot_pot > 0, n_pot / n_tot_pot, 0.0)
        err_low, err_high = [], []
        n_succ = np.asarray(ea.n_signal_raw, dtype=float)
        for i in range(len(this_eff)):
            if proportion_confint is not None and n_tot_raw[i] > 0:
                interval = proportion_confint(int(n_succ[i]), int(n_tot_raw[i]), method="wilson")
                err_low.append(abs(this_eff[i] - interval[0]))
                err_high.append(abs(this_eff[i] - interval[1]))
            else:
                err_low.append(0.0)
                err_high.append(0.0)
        eff_int_pct = (n_int / denom_int_pot * 100.0) if denom_int_pot > 0 else 0.0
        final_eff_pct = eff_int_pct
        label = f"{label_lookup.get(stage_key, stage_key)} ({eff_int_pct:.2f}%)"
        color = plt.cm.tab10(plot_idx % 10)
        ax.stairs(n_pot, bins, fill=False, alpha=0.45, color=color, linewidth=1.2, zorder=1)
        ax_eff.errorbar(
            bin_centers,
            this_eff,
            yerr=[err_low, err_high],
            fmt="o-",
            color=color,
            label=label,
            markersize=3,
            linewidth=1.0,
            zorder=4,
        )
        plot_idx += 1

    xlabel = var_config.var_labels[0] if getattr(var_config, "var_labels", None) else var_save_name
    title = f"eff: {var_save_name}"
    if final_eff_pct is not None:
        title += f" ({final_eff_pct:.2f}%)"
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel("Events", fontsize=7)
    ax_eff.set_ylabel("Efficiency", fontsize=7)
    ax_eff.set_ylim(0, 1.05)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_title(title, fontsize=8, pad=2)
    ax.tick_params(labelsize=6)
    ax_eff.tick_params(labelsize=6)
    ax_eff.set_zorder(3)
    ax_eff.patch.set_visible(False)
    ax_eff.legend(
        fontsize=legend_fontsize,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=2,
        frameon=False,
        borderpad=0.15,
        handlelength=1.0,
        labelspacing=0.15,
        columnspacing=0.8,
    )


def make_accumulation_panel(
    merged: dict,
    pot_str: str,
    *,
    status: str = "",
    figsize: Tuple[float, float] = (28, 40),
    ncols: int = 5,
    include_summary_bar: bool = True,
    include_efficiency: bool = True,
    cosmic_estimate: str = "intime",
    legend_fontsize: float = 6.5,
) -> plt.Figure:
    """Build one large figure with overlays + optional efficiency curves."""
    items = ordered_plot_keys(merged)
    eff_vars = ordered_efficiency_vars(merged) if include_efficiency else []
    n_hist = len(items)
    n_eff = len(eff_vars)
    n_extra = 1 if include_summary_bar else 0
    n_panels = max(n_hist + n_extra + n_eff, 1)
    ncols = max(1, int(ncols))
    nrows = int(math.ceil(n_panels / ncols))
    fig_w = float(figsize[0])
    fig_h = max(float(figsize[1]), 5.2 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), squeeze=False)
    axes_flat = list(axes.ravel())

    panel_i = 0
    if include_summary_bar:
        _draw_summary_bar(
            axes_flat[panel_i],
            merged,
            cosmic_estimate=cosmic_estimate,
            legend_fontsize=legend_fontsize,
        )
        panel_i += 1

    for key, ps in items:
        ax = axes_flat[panel_i]
        panel_i += 1
        hd = merged["histdata"][key]
        stage_key, plot_key = key
        if ps is not None and ps.plot_label_template is not None:
            xlabel = ps.plot_label_template[0]
            title = stage_key
            if ps.name_suffix:
                title = f"{stage_key} [{ps.name_suffix}]"
        elif ps is not None:
            xlabel = ps.var_config.var_labels[0]
            title = stage_key
        else:
            xlabel = getattr(hd, "var_save_name", plot_key)
            title = stage_key
        vlines = (ps.save_kwargs or {}).get("vline") if ps is not None else None
        _draw_compact_overlay(
            ax,
            hd,
            xlabel=xlabel,
            title=title,
            pot_str=pot_str,
            vlines=vlines,
            legend_fontsize=legend_fontsize,
        )

    for vc in eff_vars:
        ax = axes_flat[panel_i]
        panel_i += 1
        _draw_efficiency_curve(
            ax, merged, vc, legend_fontsize=max(5.5, legend_fontsize - 0.5)
        )

    for j in range(panel_i, len(axes_flat)):
        axes_flat[j].set_axis_off()

    fig.suptitle(
        status or f"Event selection accumulation (POT={pot_str})",
        fontsize=12,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    return fig


def build_live_payload(
    merged: dict,
    *,
    data_pot: float,
    pot_str: str,
    scales: dict,
    totals=None,
    cosmic_estimate: str = "intime",
    f_offbeam_frac: float = 0.08,
    n_files_done: Dict[str, int] | None = None,
    extra: dict | None = None,
) -> dict:
    """Serializable payload with everything needed to re-plot the live panel."""
    payload = {
        "merged": merged,
        "data_pot": data_pot,
        "pot_str": pot_str,
        "exposure_scales": scales,
        "exposure_totals": totals,
        "cosmic_estimate": cosmic_estimate,
        "f_offbeam_frac": f_offbeam_frac,
        "n_files_done": dict(n_files_done or {}),
        "format": "numucc_1p0pi.event_selection_live.v1",
    }
    if extra:
        payload.update(extra)
    return payload


def save_live_payload(path: Path | str, payload: dict) -> Path:
    import pickle

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def load_live_payload(path: Path | str) -> dict:
    import pickle

    with open(path, "rb") as f:
        return pickle.load(f)


def replot_from_payload(
    payload: dict | Path | str,
    *,
    figsize: Tuple[float, float] = (28, 40),
    ncols: int = 5,
    include_summary_bar: bool = True,
    include_efficiency: bool = True,
    legend_fontsize: float = 6.5,
    status: str | None = None,
    cosmic_estimate: str | None = None,
) -> plt.Figure:
    """Rebuild the accumulation panel from a saved ``merged_histdata.pkl`` payload."""
    if not isinstance(payload, dict):
        payload = load_live_payload(payload)
    merged = payload["merged"]
    pot_str = payload.get("pot_str") or get_pot_str(float(payload.get("data_pot") or 1.0))
    return make_accumulation_panel(
        merged,
        pot_str,
        status=status
        or f"Replot  POT={pot_str}  files={payload.get('n_files_done')}",
        figsize=figsize,
        ncols=ncols,
        include_summary_bar=include_summary_bar,
        include_efficiency=include_efficiency,
        cosmic_estimate=cosmic_estimate or payload.get("cosmic_estimate", "intime"),
        legend_fontsize=legend_fontsize,
    )


def resolve_merged_payload(
    *,
    merged_payload: dict | None = None,
    live_result=None,
    plots_dir: Path | str | None = None,
    merged_histdata_pkl: Path | str | None = None,
) -> tuple[dict, Path]:
    """Return ``(payload, payload_path)`` from memory or ``merged_histdata.pkl``.

    Resolution order: explicit ``merged_histdata_pkl`` → in-memory
    ``merged_payload`` / ``live_result.merged_payload`` →
    ``{plots_dir}/merged_histdata.pkl``.
    """
    if merged_histdata_pkl is not None:
        pkl = Path(merged_histdata_pkl)
        return load_live_payload(pkl), pkl

    if merged_payload is not None:
        pkl = None
        if live_result is not None and getattr(live_result, "payload_path", None):
            pkl = Path(live_result.payload_path)
        elif plots_dir is not None:
            pkl = Path(plots_dir) / "merged_histdata.pkl"
        else:
            pkl = Path(".")
        return merged_payload, pkl

    if live_result is not None and getattr(live_result, "merged_payload", None) is not None:
        pkl = Path(live_result.payload_path) if live_result.payload_path else (
            Path(live_result.plots_dir) / "merged_histdata.pkl"
        )
        return live_result.merged_payload, pkl

    if plots_dir is not None:
        pkl = Path(plots_dir) / "merged_histdata.pkl"
        if pkl.is_file():
            return load_live_payload(pkl), pkl

    raise FileNotFoundError(
        "No merged payload available. Re-run the live batched cell or set "
        "MERGED_HISTDATA_PKL to a merged_histdata.pkl path."
    )


def render_followup_from_payload(
    payload: dict | Path | str,
    *,
    save_fig_dir: Path | str | None = None,
    cosmic_estimate: str = "intime",
    save_fig: bool = True,
    show_fig: bool = True,
    render_overlays: bool = True,
    render_summary: bool = True,
    render_efficiency: bool = True,
    syst_disk_root: str | None = None,
) -> dict:
    """Render summary / overlay / efficiency plots from a live/aggregate payload.

    Parameters
    ----------
    cosmic_estimate
        ``"intime"`` or ``"offbeam"`` — which cosmic sample to use in the
        summary bars and distribution overlays.
    """
    if not isinstance(payload, dict):
        payload = load_live_payload(payload)

    merged = payload["merged"]
    pot_str = payload.get("pot_str") or get_pot_str(float(payload.get("data_pot") or 1.0))
    if save_fig_dir is None:
        save_fig_dir = Path(".")
    save_fig_dir = Path(save_fig_dir)
    save_fig_dir.mkdir(parents=True, exist_ok=True)

    cosmic_estimate = str(cosmic_estimate).lower()
    if cosmic_estimate not in ("intime", "offbeam"):
        raise ValueError("cosmic_estimate must be 'intime' or 'offbeam'")

    agg = _load_aggregate_module()
    if render_overlays:
        agg.render_overlay_plots(
            merged,
            plot_label_map={},
            save_fig_dir=str(save_fig_dir),
            pot_str=pot_str,
            save_fig=save_fig,
            show_fig=show_fig,
            syst_disk_root=syst_disk_root,
            cosmic_estimate=cosmic_estimate,
        )
    if render_summary:
        agg.render_summary_breakdown_plot(
            merged,
            str(save_fig_dir),
            save_fig=save_fig,
            show_fig=show_fig,
            cosmic_estimate=cosmic_estimate,
        )
    if render_efficiency:
        agg.render_efficiency_plots(
            merged,
            str(save_fig_dir),
            pot_str,
            save_fig=save_fig,
            show_fig=show_fig,
        )

    return {
        "merged": merged,
        "pot_str": pot_str,
        "save_fig_dir": str(save_fig_dir),
        "cosmic_estimate": cosmic_estimate,
        "data_pot": payload.get("data_pot"),
        "n_files_done": payload.get("n_files_done"),
    }


def _fig_to_png_bytes(fig: plt.Figure, dpi: int = 120) -> bytes:
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


class LiveNotebookDisplay:
    """Update one inline image in the notebook cell (no file write)."""

    def __init__(self, dpi: int = 120):
        self.dpi = dpi
        self._display_id = None

    def __call__(self, fig: plt.Figure, status: str, _merged: dict | None = None) -> None:
        from IPython.display import Image, display, update_display

        png = _fig_to_png_bytes(fig, dpi=self.dpi)
        plt.close(fig)
        img = Image(data=png)
        if self._display_id is None:
            handle = display(img, display_id=True)
            self._display_id = handle.display_id
        else:
            update_display(img, display_id=self._display_id)


def _display_panel(fig: plt.Figure, status: str) -> None:
    """Fallback one-shot inline display (clears prior cell output)."""
    try:
        from IPython.display import Image, clear_output, display

        png = _fig_to_png_bytes(fig)
        clear_output(wait=True)
        print(status, flush=True)
        display(Image(data=png))
    except Exception:
        print(status, flush=True)
        try:
            plt.show()
        except Exception:
            pass
    finally:
        plt.close(fig)


def run_live_accumulate(
    live_cfg: LiveAccumulateConfig | None = None,
    *,
    on_update: Optional[Callable[[plt.Figure, str, dict], None]] = None,
) -> LiveAccumulateResult:
    """Run map jobs file-by-file and refresh the accumulation panel each round.

    Parameters
    ----------
    on_update
        Optional callback ``(fig, status, merged)``. Default displays inline in
        IPython via ``clear_output`` + ``display``.
    """
    live_cfg = live_cfg or LiveAccumulateConfig()
    today = __import__("datetime").datetime.now().strftime("%Y%m%d")
    work = Path(
        live_cfg.work_base or default_event_selection_batched_work_root(f"live-{today}")
    ).expanduser()
    batches_dir = Path(live_cfg.batches_dir or work / "batches").expanduser()
    plots_dir = Path(live_cfg.plots_dir or work / "plots").expanduser()
    work.mkdir(parents=True, exist_ok=True)
    batches_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    if live_cfg.base_dir is not None and live_cfg.sample_dirs:
        records = survey_files_from_dirs(
            live_cfg.base_dir,
            live_cfg.sample_dirs,
            samples=live_cfg.samples,
            max_files_per_sample=live_cfg.max_files_per_sample,
            filename_str=live_cfg.filename_str,
        )
        if not records:
            raise RuntimeError(
                f"No input .df files under base_dir={live_cfg.base_dir!r} "
                f"sample_dirs={live_cfg.sample_dirs!r} "
                f"(filename_str={live_cfg.filename_str!r})"
            )
    else:
        records = survey_files(live_cfg.samples, live_cfg.max_files_per_sample)
        if not records:
            raise RuntimeError("No input .df files matched EVENT_SELECTION_GLOBS")

    files_per_job = _resolve_files_per_job(live_cfg)
    jobs = _jobs_chunked(records, files_per_job=files_per_job)
    write_manifest(work, records, jobs, max_job_bytes=1 << 30)
    print(
        f"[live] files_per_job={files_per_job}  concat_load={live_cfg.concat_load}  "
        f"update_every_n_files={live_cfg.update_every_n_files}  n_jobs={len(jobs)}",
        flush=True,
    )
    print_survey_summary(records, jobs)

    by_sample: Dict[str, List[BatchJob]] = {s: [] for s in live_cfg.samples}
    for job in jobs:
        by_sample.setdefault(job.sample, []).append(job)
    for s in by_sample:
        by_sample[s].sort(key=lambda j: j.job_id)

    max_rounds = max((len(v) for v in by_sample.values()), default=0)
    n_done: Dict[str, int] = {s: 0 for s in live_cfg.samples}
    n_files_total = 0
    n_files_planned = sum(len(j.files) for j in jobs)
    n_jobs_planned = len(jobs)
    failed: List[Tuple[str, str]] = []
    merged_payload = None
    pot_str = ""
    data_pot = 0.0
    panel_path = None
    payload_path = None

    display_fn = on_update
    if display_fn is None and live_cfg.show_in_notebook:
        display_fn = LiveNotebookDisplay(dpi=live_cfg.panel_dpi)
    elif display_fn is None:
        display_fn = lambda fig, _status, _merged: plt.close(fig)

    update_every = max(1, int(live_cfg.update_every_n_files or 1))

    try:
        from tqdm.auto import tqdm
    except ImportError:
        from tqdm import tqdm  # type: ignore

    def _should_update_panel(*, is_last: bool) -> bool:
        if not live_cfg.update_every_round:
            return is_last
        if is_last:
            return True
        return n_files_total > 0 and (n_files_total % update_every == 0)

    def _sample_counts_str() -> str:
        return " ".join(f"{s}:{n_done.get(s, 0)}" for s in live_cfg.samples)

    def _refresh_panel(round_i: int) -> None:
        nonlocal merged_payload, data_pot, pot_str, payload_path
        try:
            merged, data_pot, pot_str, scales, totals = aggregate_batches_so_far(
                batches_dir, f_offbeam_frac=live_cfg.f_offbeam_frac, quiet=True
            )
        except RuntimeError as exc:
            tqdm.write(f"[live] skip panel ({exc})")
            return

        tqdm.write(
            f"[live] exposure  files={n_files_total}/{n_files_planned}  "
            f"data_pot={totals.data_pot:.3e}  mc_pot={totals.mc_pot:.3e}  "
            f"dirt_pot={totals.dirt_pot:.3e}  "
            f"intime_gates={totals.intime_gates:.3e}  "
            f"scales={{mc:{scales.get('scale_mc', 0):.4g}, "
            f"dirt:{scales.get('scale_dirt', 0):.4g}, "
            f"intime:{scales.get('scale_intime', 0):.4g}}}"
        )

        status = (
            f"Live accumulation  files={n_files_total}/{n_files_planned}  "
            f"{_sample_counts_str()}  POT={pot_str}"
        )
        fig = make_accumulation_panel(
            merged,
            pot_str,
            status=status,
            figsize=live_cfg.figsize,
            ncols=live_cfg.ncols,
            include_summary_bar=live_cfg.include_summary_bar,
            include_efficiency=live_cfg.include_efficiency,
            cosmic_estimate=live_cfg.cosmic_estimate,
            legend_fontsize=live_cfg.legend_fontsize,
        )
        display_fn(fig, status, merged)
        merged_payload = build_live_payload(
            merged,
            data_pot=data_pot,
            pot_str=pot_str,
            scales=scales,
            totals=totals,
            cosmic_estimate=live_cfg.cosmic_estimate,
            f_offbeam_frac=live_cfg.f_offbeam_frac,
            n_files_done=n_done,
        )
        if live_cfg.save_merged_payload:
            payload_path = save_live_payload(
                plots_dir / "merged_histdata.pkl", merged_payload
            )

    pbar = tqdm(
        total=n_files_planned,
        desc="live selection",
        unit="file",
        leave=True,
        dynamic_ncols=True,
    )
    try:
        for round_i in range(max_rounds):
            for sample in live_cfg.samples:
                sample_jobs = by_sample.get(sample) or []
                if round_i >= len(sample_jobs):
                    continue
                job = sample_jobs[round_i]
                tag = f"{sample}/{job.tag}"
                n_job_files = len(job.files)
                pbar.set_postfix_str(
                    f"{_sample_counts_str()} | {sample}×{n_job_files}",
                    refresh=False,
                )
                try:
                    process_job_inprocess(job, batches_dir, cfg=live_cfg)
                    n_done[sample] = n_done.get(sample, 0) + n_job_files
                except Exception as exc:
                    msg = f"{type(exc).__name__}: {exc}"
                    tqdm.write(f"[live] FAILED {tag}: {msg}")
                    failed.append((tag, msg))
                n_files_total += n_job_files
                pbar.update(n_job_files)
                gc.collect()

                is_last = n_files_total >= n_files_planned
                if _should_update_panel(is_last=is_last):
                    _refresh_panel(round_i)
                    pbar.set_postfix_str(_sample_counts_str(), refresh=True)
    finally:
        pbar.close()

    if merged_payload is not None:
        if live_cfg.save_merged_payload:
            payload_path = save_live_payload(
                plots_dir / "merged_histdata.pkl", merged_payload
            )
            print(f"[live] wrote reloadable counts → {payload_path}", flush=True)
        if live_cfg.save_final_panel:
            panel_path = plots_dir / "live_accumulation_panel.png"
            fig = make_accumulation_panel(
                merged_payload["merged"],
                merged_payload["pot_str"],
                status=(
                    f"Final live panel  files/sample={dict(n_done)}  "
                    f"POT={merged_payload['pot_str']}"
                ),
                figsize=live_cfg.figsize,
                ncols=live_cfg.ncols,
                include_summary_bar=live_cfg.include_summary_bar,
                include_efficiency=live_cfg.include_efficiency,
                cosmic_estimate=live_cfg.cosmic_estimate,
                legend_fontsize=live_cfg.legend_fontsize,
            )
            fig.savefig(panel_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"[live] wrote {panel_path}", flush=True)

    return LiveAccumulateResult(
        work_base=work,
        batches_dir=batches_dir,
        plots_dir=plots_dir,
        merged_payload=merged_payload,
        pot_str=pot_str,
        data_pot=data_pot,
        n_files_done=n_done,
        failed=failed,
        panel_path=panel_path,
        payload_path=payload_path,
    )


def batched_cfg_from_live(live_cfg: LiveAccumulateConfig) -> EventSelectionBatchedConfig:
    """Map live knobs onto the non-interactive ``EventSelectionBatchedConfig``."""
    files_per_job = _resolve_files_per_job(live_cfg)
    return EventSelectionBatchedConfig(
        work_base=live_cfg.work_base,
        batches_dir=live_cfg.batches_dir,
        plots_dir=live_cfg.plots_dir,
        max_files_per_sample=live_cfg.max_files_per_sample,
        max_job_bytes=1 if files_per_job <= 1 else (1 << 30),
        mc_univ_syst=live_cfg.mc_univ_syst,
        use_mc_genweight=live_cfg.use_mc_genweight,
        skip_existing_batches=live_cfg.skip_existing,
        f_offbeam_frac=live_cfg.f_offbeam_frac,
        cosmic_estimate=live_cfg.cosmic_estimate,
        samples=live_cfg.samples,
        base_dir=live_cfg.base_dir,
        sample_dirs=live_cfg.sample_dirs,
        filename_str=live_cfg.filename_str,
        trace=live_cfg.trace,
    )


# keep deepcopy available for callers who mutate merged
__all__ = [
    "LiveAccumulateConfig",
    "LiveAccumulateResult",
    "LiveNotebookDisplay",
    "run_live_accumulate",
    "make_accumulation_panel",
    "aggregate_batches_so_far",
    "build_live_payload",
    "save_live_payload",
    "load_live_payload",
    "replot_from_payload",
    "resolve_merged_payload",
    "render_followup_from_payload",
    "batched_cfg_from_live",
    "process_job_inprocess",
]
