"""Notebook-friendly map/reduce driver for numuCC 1p0pi event selection.

Processes one ``.df`` file at a time (see ``scripts/event_selection_chunk.py``),
writes per-shard histogram pickles, then merges and plots via
``scripts/event_selection_aggregate.py``. Peak memory is one input file plus
accumulated histograms — not the full sample concatenated in RAM.

Typical notebook usage::

    from analysis_village.numucc_1p0pi.event_selection_chunked import (
        EventSelectionChunkedConfig,
        run_full,
        show_saved_plots,
    )

    cfg = EventSelectionChunkedConfig()
    result = run_full(cfg)
    show_saved_plots(result.plots_dir)

Shell equivalent: ``scripts/run_event_selection_chunked.sh``.
"""
from __future__ import annotations

import importlib.util
import os
import pickle
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from analysis_village.numucc_1p0pi.dataset_locations import (
    EVENT_SELECTION_GLOBS,
    PLOTS_BASE,
    default_event_selection_work_root,
    default_syst_disk_root,
    iter_event_selection_df_paths,
)

SAMPLES: Tuple[str, ...] = ("mc", "data", "intime", "offbeam", "dirt")
_SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
_CHUNK_SCRIPT = _SCRIPTS_DIR / "event_selection_chunk.py"
_AGG_SCRIPT = _SCRIPTS_DIR / "event_selection_aggregate.py"


@dataclass
class EventSelectionChunkedConfig:
    """Configuration for the chunked event-selection workflow."""

    work_base: Path | str | None = None
    chunks_dir: Path | str | None = None
    plots_dir: Path | str | None = None
    mc_univ_syst: Sequence[str] = ("Flux", "G4", "GENIE")
    load_mode: str = "splits"
    max_splits: int = 0
    use_mc_genweight: bool = False
    skip_existing_chunks: bool = True
    aggregate_only: bool = False
    skip_aggregate: bool = False
    cosmic_estimate: str = "intime"
    hide_cosmic_model_unc: bool = False
    f_offbeam_frac: float = 0.08
    save_fig: bool = True
    show_fig: bool = False
    overlay_syst_from_universes: bool = True
    syst_disk_root: Path | str | None = None
    syst_tag: str = ""
    samples: Sequence[str] = SAMPLES
    max_files_per_sample: int | None = None
    trace: bool = False
    python_executable: str = field(default_factory=lambda: sys.executable)

    def resolve_paths(self) -> Tuple[Path, Path, Path]:
        today = datetime.now().strftime("%Y%m%d")
        work = Path(self.work_base or default_event_selection_work_root(today)).expanduser()
        chunks = Path(self.chunks_dir or work / "chunks").expanduser()
        plots = Path(self.plots_dir or work / "plots").expanduser()
        return work, chunks, plots


@dataclass
class EventSelectionChunkedResult:
    work_base: Path
    chunks_dir: Path
    plots_dir: Path
    chunk_paths: Dict[str, List[str]]
    failed: List[Tuple[str, str, str]]
    merged_payload: dict | None
    pot_str: str
    data_pot: float


def discover_jobs(
    samples: Sequence[str] = SAMPLES,
    max_files_per_sample: int | None = None,
) -> List[Tuple[str, str]]:
    """Return ``(sample, df_path)`` jobs from :data:`EVENT_SELECTION_GLOBS`."""
    jobs: List[Tuple[str, str]] = []
    for sample in samples:
        if sample not in EVENT_SELECTION_GLOBS:
            raise KeyError(
                f"unknown sample {sample!r}; choose from {tuple(EVENT_SELECTION_GLOBS)}"
            )
        paths = list(iter_event_selection_df_paths(sample))
        if max_files_per_sample is not None:
            paths = paths[: max(0, int(max_files_per_sample))]
        if not paths:
            print(
                f"[chunked] sample={sample} no files matched {EVENT_SELECTION_GLOBS[sample]}",
                flush=True,
            )
            continue
        print(
            f"[chunked] sample={sample}  {len(paths)} file(s)  ({EVENT_SELECTION_GLOBS[sample]})",
            flush=True,
        )
        for p in paths:
            jobs.append((sample, p))
    return jobs


def _chunk_out_path(chunks_dir: Path, sample: str, df_file: str) -> Path:
    stem = Path(df_file).stem
    return chunks_dir / f"{sample}__{stem}.pkl"


def process_one_file(
    df_file: str,
    sample: str,
    out_dir: Path | str,
    *,
    cfg: EventSelectionChunkedConfig | None = None,
    skip_existing: bool = True,
) -> Optional[str]:
    """Run ``event_selection_chunk.py`` on one ``.df`` file. Returns output pickle path."""
    cfg = cfg or EventSelectionChunkedConfig()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pkl = _chunk_out_path(out_dir, sample, df_file)
    if skip_existing and out_pkl.is_file():
        print(f"[chunked] {out_pkl} exists; skipping", flush=True)
        return str(out_pkl)

    cmd = [
        cfg.python_executable,
        str(_CHUNK_SCRIPT),
        "--df_file",
        df_file,
        "--sample",
        sample,
        "--out_dir",
        str(out_dir),
        "--load_mode",
        cfg.load_mode,
    ]
    if cfg.max_splits > 0:
        cmd.extend(["--max_splits", str(cfg.max_splits)])
    if cfg.use_mc_genweight:
        cmd.append("--use-mc-genweight")
    if cfg.trace:
        cmd.append("--trace")
    if sample == "mc" and cfg.mc_univ_syst:
        cmd.extend(["--mc-univ-syst", ",".join(cfg.mc_univ_syst)])

    print(f"[chunked] sample={sample} file={df_file}", flush=True)
    subprocess.run(cmd, check=True)
    return str(out_pkl)


def run_map(
    cfg: EventSelectionChunkedConfig,
    jobs: Sequence[Tuple[str, str]] | None = None,
) -> Tuple[Path, Dict[str, List[str]], List[Tuple[str, str, str]]]:
    """Map phase: one pickle per *(sample, .df file)*."""
    work, chunks_dir, _ = cfg.resolve_paths()
    chunks_dir.mkdir(parents=True, exist_ok=True)
    jobs = list(jobs if jobs is not None else discover_jobs(cfg.samples, cfg.max_files_per_sample))
    if not jobs:
        raise RuntimeError("No input .df files matched EVENT_SELECTION_GLOBS")

    chunk_paths: Dict[str, List[str]] = {s: [] for s in cfg.samples}
    failed: List[Tuple[str, str, str]] = []

    print(f"[chunked] WORK_BASE={work}", flush=True)
    print(f"[chunked] CHUNKS_DIR={chunks_dir}", flush=True)
    print(f"[chunked] MC_UNIV_SYST={','.join(cfg.mc_univ_syst) if cfg.mc_univ_syst else '<disabled>'}", flush=True)

    for sample, df_path in jobs:
        if not os.path.isfile(df_path):
            msg = f"missing file: {df_path}"
            print(f"[chunked] sample={sample} {msg}", flush=True)
            failed.append((sample, df_path, msg))
            continue
        try:
            out = process_one_file(
                df_path,
                sample,
                chunks_dir,
                cfg=cfg,
                skip_existing=cfg.skip_existing_chunks,
            )
            if out:
                chunk_paths.setdefault(sample, []).append(out)
        except subprocess.CalledProcessError as exc:
            msg = f"chunk failed exit={exc.returncode}"
            print(f"[chunked] FAILED sample={sample} file={df_path} ({msg})", flush=True)
            failed.append((sample, df_path, msg))

    for sample, paths in chunk_paths.items():
        print(f"[chunked] sample={sample} wrote {len(paths)} chunk pickle(s)", flush=True)
    return chunks_dir, chunk_paths, failed


def _load_aggregate_module():
    spec = importlib.util.spec_from_file_location(
        "event_selection_aggregate", _AGG_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {_AGG_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_aggregate(cfg: EventSelectionChunkedConfig, chunks_dir: Path | str | None = None) -> EventSelectionChunkedResult:
    """Reduce phase: merge chunk pickles and render all event-selection plots."""
    work, default_chunks, plots_dir = cfg.resolve_paths()
    chunks_dir = Path(chunks_dir or default_chunks)
    plots_dir.mkdir(parents=True, exist_ok=True)

    syst_root = cfg.syst_disk_root
    if syst_root is None:
        syst_root = os.environ.get("NUMUCC_SYST_DISK_ROOT")
    if syst_root is None:
        syst_root = default_syst_disk_root()
    syst_root = Path(syst_root).expanduser()

    agg = _load_aggregate_module()
    chunk_groups = agg.collect_chunks(str(chunks_dir))

    samples = {}
    for s, files in chunk_groups.items():
        if not files:
            continue
        print(f"[chunked] aggregating {len(files)} chunks for sample={s}", flush=True)
        samples[s] = agg.aggregate_chunk_files(files)

    if not samples:
        raise RuntimeError(f"No chunk pickles found under {chunks_dir}")

    merged = agg.merge_samples(samples)
    n_hd_fixed, n_bar_fixed = agg.sanitize_merged_histdata_finite(merged)
    if n_hd_fixed or n_bar_fixed:
        print(
            f"[chunked] sanitized NaN/inf bins: {n_hd_fixed} hist keys, {n_bar_fixed} bar rows",
            flush=True,
        )

    totals = agg.accumulate_exposure_totals_from_dir(str(chunks_dir))
    exposure_scales = agg.apply_global_exposure_scales(
        merged, totals, f_offbeam_coincident=cfg.f_offbeam_frac
    )
    print(f"[chunked] applied global scales: {exposure_scales}", flush=True)

    data_pot = totals.data_pot if totals.data_pot > 0 else 1.0
    pot_str = agg.get_pot_str(data_pot)
    print(f"[chunked] data_pot={data_pot:.3e} -> POT label={pot_str}", flush=True)

    overlay_syst = cfg.overlay_syst_from_universes
    syst_disk_arg = str(syst_root) if syst_root.is_dir() else None
    if syst_disk_arg is None:
        print(
            f"[chunked] WARN: syst disk not found ({syst_root}); "
            "overlay bands from MC universes in chunks only",
            flush=True,
        )

    agg.render_overlay_plots(
        merged,
        plot_label_map={},
        save_fig_dir=str(plots_dir),
        pot_str=pot_str,
        save_fig=cfg.save_fig,
        show_fig=cfg.show_fig,
        cosmic_estimate=cfg.cosmic_estimate,
        show_cosmic_model_unc=not cfg.hide_cosmic_model_unc,
        overlay_syst_from_universes=overlay_syst,
        syst_disk_root=syst_disk_arg,
    )
    agg.render_summary_breakdown_plot(
        merged,
        str(plots_dir),
        save_fig=cfg.save_fig,
        show_fig=cfg.show_fig,
        cosmic_estimate=cfg.cosmic_estimate,
    )
    agg.render_efficiency_plots(
        merged,
        str(plots_dir),
        pot_str,
        save_fig=cfg.save_fig,
        show_fig=cfg.show_fig,
    )

    merged_payload = {
        "merged": merged,
        "data_pot": data_pot,
        "pot_str": pot_str,
        "exposure_totals": totals,
        "exposure_scales": exposure_scales,
        "cosmic_estimate": cfg.cosmic_estimate,
        "f_offbeam_frac": cfg.f_offbeam_frac,
    }
    out_pkl = plots_dir / "merged_histdata.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(merged_payload, f)
    print(f"[chunked] wrote {out_pkl}", flush=True)

    chunk_paths = {s: sorted(str(p) for p in chunks_dir.glob(f"{s}__*.pkl")) for s in SAMPLES}

    return EventSelectionChunkedResult(
        work_base=work,
        chunks_dir=chunks_dir,
        plots_dir=plots_dir,
        chunk_paths=chunk_paths,
        failed=[],
        merged_payload=merged_payload,
        pot_str=pot_str,
        data_pot=data_pot,
    )


def load_merged_results(plots_dir: Path | str) -> dict:
    """Load ``merged_histdata.pkl`` written by :func:`run_aggregate`."""
    pkl = Path(plots_dir) / "merged_histdata.pkl"
    with open(pkl, "rb") as f:
        return pickle.load(f)


def run_full(cfg: EventSelectionChunkedConfig | None = None) -> EventSelectionChunkedResult:
    """Run map (unless ``aggregate_only``) then reduce (unless ``skip_aggregate``)."""
    cfg = cfg or EventSelectionChunkedConfig()
    work, chunks_dir, plots_dir = cfg.resolve_paths()

    if not cfg.aggregate_only:
        _, _, failed = run_map(cfg)
        if failed:
            print(f"[chunked] WARN: {len(failed)} chunk job(s) failed", flush=True)
    else:
        print("[chunked] aggregate_only=1 — skipping map phase", flush=True)

    if cfg.skip_aggregate:
        print("[chunked] skip_aggregate=1 — map pickles only", flush=True)
        chunk_paths = {s: sorted(str(p) for p in chunks_dir.glob(f"{s}__*.pkl")) for s in SAMPLES}
        return EventSelectionChunkedResult(
            work_base=work,
            chunks_dir=chunks_dir,
            plots_dir=plots_dir,
            chunk_paths=chunk_paths,
            failed=[],
            merged_payload=None,
            pot_str="",
            data_pot=0.0,
        )

    result = run_aggregate(cfg, chunks_dir=chunks_dir)
    print(f"[chunked] DONE plots -> {result.plots_dir}", flush=True)
    return result


def show_saved_plots(
    plots_dir: Path | str,
    *,
    patterns: Sequence[str] = ("*.png",),
    max_images: int | None = None,
) -> None:
    """Display PNGs from the aggregate step inside a Jupyter notebook."""
    try:
        from IPython.display import Image, display
    except ImportError as exc:
        raise ImportError(
            "show_saved_plots requires IPython (run inside a Jupyter notebook)"
        ) from exc

    plots_dir = Path(plots_dir)
    paths: List[Path] = []
    for pat in patterns:
        paths.extend(sorted(plots_dir.glob(pat)))
    paths = sorted(set(paths))
    if max_images is not None:
        paths = paths[: max(0, int(max_images))]

    if not paths:
        print(f"[chunked] no plots matched under {plots_dir}", flush=True)
        return

    print(f"[chunked] displaying {len(paths)} plot(s) from {plots_dir}", flush=True)
    for p in paths:
        display(Image(filename=str(p)))
