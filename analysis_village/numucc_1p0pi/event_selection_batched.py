"""Batched event selection — notebook pipeline over ≤1 GB file groups.

Survey input ``.df`` files from :data:`dataset_locations.EVENT_SELECTION_GLOBS`,
pack them into jobs under a size budget, run the same selection pipeline as
``notebooks/event_selection.ipynb`` (via ``build_runner`` / ``ChunkRunner``),
save histogram + breakdown pickles per job, then aggregate with
``scripts/event_selection_aggregate.py``.

Typical notebook usage::

    from analysis_village.numucc_1p0pi.event_selection_batched import (
        EventSelectionBatchedConfig,
        run_full,
    )

    result = run_full(EventSelectionBatchedConfig())
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from analysis_village.numucc_1p0pi.dataset_locations import (
    EVENT_SELECTION_GLOBS,
    PLOTS_BASE,
    default_syst_disk_root,
    iter_event_selection_df_paths,
)

SAMPLES: Tuple[str, ...] = ("mc", "data", "intime", "offbeam", "dirt")
_SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
_MAP_SCRIPT = _SCRIPTS_DIR / "event_selection_batch_map.py"
_AGG_SCRIPT = _SCRIPTS_DIR / "event_selection_aggregate.py"

DEFAULT_MAX_JOB_BYTES = 1 << 30  # 1 GiB


@dataclass
class FileRecord:
    sample: str
    path: str
    size_bytes: int

    @property
    def size_gb(self) -> float:
        return self.size_bytes / (1024.0 ** 3)


@dataclass
class BatchJob:
    sample: str
    job_id: int
    files: List[str]
    total_bytes: int

    @property
    def tag(self) -> str:
        return f"batch_{self.job_id:04d}"


@dataclass
class EventSelectionBatchedConfig:
    work_base: Path | str | None = None
    batches_dir: Path | str | None = None
    plots_dir: Path | str | None = None
    max_job_bytes: int = DEFAULT_MAX_JOB_BYTES
    # Optional map-phase bookkeeping only (not used to form overlay bands).
    mc_univ_syst: Sequence[str] = ()
    use_mc_genweight: bool = False
    skip_existing_batches: bool = True
    aggregate_only: bool = False
    skip_aggregate: bool = False
    cosmic_estimate: str = "intime"
    f_offbeam_frac: float = 0.08
    save_fig: bool = True
    show_fig: bool = False
    syst_disk_root: Path | str | None = None
    syst_tag: str = ""
    samples: Sequence[str] = SAMPLES
    max_files_per_sample: int | None = None
    # Notebook-style inputs (preferred when set). Else fall back to EVENT_SELECTION_GLOBS.
    base_dir: Path | str | None = None
    sample_dirs: Dict[str, str] | None = None
    filename_str: str = "sel_all"
    trace: bool = False
    python_executable: str = field(default_factory=lambda: sys.executable)

    def resolve_paths(self) -> Tuple[Path, Path, Path]:
        today = datetime.now().strftime("%Y%m%d")
        work = Path(self.work_base or default_event_selection_batched_work_root(today)).expanduser()
        batches = Path(self.batches_dir or work / "batches").expanduser()
        plots = Path(self.plots_dir or work / "plots").expanduser()
        return work, batches, plots


@dataclass
class EventSelectionBatchedResult:
    work_base: Path
    batches_dir: Path
    plots_dir: Path
    batch_paths: Dict[str, List[str]]
    failed: List[Tuple[str, str, str]]
    merged_payload: dict | None
    pot_str: str
    data_pot: float
    manifest_path: Path


def default_event_selection_batched_work_root(tag: str | None = None) -> Path:
    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_EVENT_SELECTION_WORK_BASE")
    if base:
        return Path(base)
    user = os.environ.get("USER", "user")
    return Path(f"/exp/sbnd/data/users/{user}/xsec/numucc_1p0pi/event_selection-batched-{t}")


def survey_files(
    samples: Sequence[str] = SAMPLES,
    max_files_per_sample: int | None = None,
) -> List[FileRecord]:
    """List input ``.df`` files with on-disk sizes from ``EVENT_SELECTION_GLOBS``."""
    records: List[FileRecord] = []
    for sample in samples:
        if sample not in EVENT_SELECTION_GLOBS:
            raise KeyError(
                f"unknown sample {sample!r}; choose from {tuple(EVENT_SELECTION_GLOBS)}"
            )
        paths = list(iter_event_selection_df_paths(sample))
        if max_files_per_sample is not None:
            paths = paths[: max(0, int(max_files_per_sample))]
        for p in paths:
            if not os.path.isfile(p):
                continue
            records.append(
                FileRecord(sample=sample, path=p, size_bytes=os.path.getsize(p))
            )
    return records


def survey_files_from_dirs(
    base_dir: Path | str,
    sample_dirs: Dict[str, str],
    *,
    samples: Sequence[str] | None = None,
    max_files_per_sample: int | None = None,
    filename_str: str = "sel_all",
) -> List[FileRecord]:
    """List ``.df`` files under notebook-style ``base_dir/<subdir>/*<filename_str>*.df``.

    ``sample_dirs`` maps sample name → subdirectory under ``base_dir`` (same layout
    as the event-selection notebook config cell).
    """
    import glob

    base_dir = Path(base_dir)
    sample_order = list(samples) if samples is not None else list(sample_dirs.keys())
    records: List[FileRecord] = []
    for sample in sample_order:
        if sample not in sample_dirs:
            raise KeyError(
                f"sample {sample!r} missing from sample_dirs; have {tuple(sample_dirs)}"
            )
        search_dir = base_dir / sample_dirs[sample]
        pattern = str(search_dir / f"*{filename_str}*.df")
        paths = sorted(glob.glob(pattern))
        if max_files_per_sample is not None:
            paths = paths[: max(0, int(max_files_per_sample))]
        for p in paths:
            if not os.path.isfile(p):
                continue
            records.append(
                FileRecord(sample=sample, path=p, size_bytes=os.path.getsize(p))
            )
    return records


def group_files_into_jobs(
    records: Sequence[FileRecord],
    max_bytes: int = DEFAULT_MAX_JOB_BYTES,
) -> List[BatchJob]:
    """Greedy bin-packing: group files per sample without exceeding ``max_bytes``."""
    by_sample: Dict[str, List[FileRecord]] = {}
    for rec in records:
        by_sample.setdefault(rec.sample, []).append(rec)

    jobs: List[BatchJob] = []
    for sample in sorted(by_sample):
        files = sorted(by_sample[sample], key=lambda r: r.path)
        batch_files: List[str] = []
        batch_bytes = 0
        job_id = 0
        for rec in files:
            if rec.size_bytes > max_bytes:
                if batch_files:
                    jobs.append(
                        BatchJob(sample=sample, job_id=job_id, files=batch_files, total_bytes=batch_bytes)
                    )
                    job_id += 1
                    batch_files, batch_bytes = [], 0
                jobs.append(
                    BatchJob(
                        sample=sample,
                        job_id=job_id,
                        files=[rec.path],
                        total_bytes=rec.size_bytes,
                    )
                )
                job_id += 1
                continue
            if batch_files and batch_bytes + rec.size_bytes > max_bytes:
                jobs.append(
                    BatchJob(sample=sample, job_id=job_id, files=batch_files, total_bytes=batch_bytes)
                )
                job_id += 1
                batch_files, batch_bytes = [], 0
            batch_files.append(rec.path)
            batch_bytes += rec.size_bytes
        if batch_files:
            jobs.append(
                BatchJob(sample=sample, job_id=job_id, files=batch_files, total_bytes=batch_bytes)
            )
    return jobs


def write_manifest(
    work_dir: Path | str,
    records: Sequence[FileRecord],
    jobs: Sequence[BatchJob],
    *,
    max_job_bytes: int = DEFAULT_MAX_JOB_BYTES,
) -> Path:
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created": datetime.now().isoformat(),
        "max_job_bytes": max_job_bytes,
        "file_survey": [asdict(r) for r in records],
        "jobs": [
            {
                "sample": j.sample,
                "job_id": j.job_id,
                "tag": j.tag,
                "files": j.files,
                "total_bytes": j.total_bytes,
                "total_gb": j.total_bytes / (1024.0 ** 3),
            }
            for j in jobs
        ],
    }
    out = work_dir / "manifest.json"
    with open(out, "w") as f:
        json.dump(manifest, f, indent=2)
    return out


def print_survey_summary(records: Sequence[FileRecord], jobs: Sequence[BatchJob]) -> None:
    by_sample: Dict[str, List[FileRecord]] = {}
    for r in records:
        by_sample.setdefault(r.sample, []).append(r)
    for sample, recs in sorted(by_sample.items()):
        tot = sum(r.size_bytes for r in recs)
        loc = EVENT_SELECTION_GLOBS.get(sample) if sample in EVENT_SELECTION_GLOBS else None
        if loc is None and recs:
            loc = str(Path(recs[0].path).parent)
        print(
            f"[batched] sample={sample}  files={len(recs)}  total={tot / (1024**3):.2f} GiB  "
            f"loc={loc}",
            flush=True,
        )
    print(f"[batched] {len(jobs)} job(s) under size budget", flush=True)
    for j in jobs[:12]:
        print(
            f"  {j.sample} {j.tag}: {len(j.files)} file(s), {j.total_bytes / (1024**3):.3f} GiB",
            flush=True,
        )
    if len(jobs) > 12:
        print(f"  ... and {len(jobs) - 12} more", flush=True)


def survey_files_for_config(cfg: EventSelectionBatchedConfig) -> List[FileRecord]:
    """Resolve input files from notebook dirs if set, else ``EVENT_SELECTION_GLOBS``."""
    if cfg.base_dir is not None and cfg.sample_dirs:
        return survey_files_from_dirs(
            cfg.base_dir,
            cfg.sample_dirs,
            samples=cfg.samples,
            max_files_per_sample=cfg.max_files_per_sample,
            filename_str=cfg.filename_str,
        )
    return survey_files(cfg.samples, cfg.max_files_per_sample)


def discover_jobs(cfg: EventSelectionBatchedConfig) -> Tuple[List[FileRecord], List[BatchJob], Path]:
    work, _, _ = cfg.resolve_paths()
    records = survey_files_for_config(cfg)
    if not records:
        if cfg.base_dir is not None and cfg.sample_dirs:
            raise RuntimeError(
                f"No input .df files under base_dir={cfg.base_dir!r} "
                f"sample_dirs={cfg.sample_dirs!r} (filename_str={cfg.filename_str!r})"
            )
        raise RuntimeError("No input .df files matched EVENT_SELECTION_GLOBS")
    jobs = group_files_into_jobs(records, max_bytes=cfg.max_job_bytes)
    manifest_path = write_manifest(work, records, jobs, max_job_bytes=cfg.max_job_bytes)
    print_survey_summary(records, jobs)
    return records, jobs, manifest_path


def _batch_out_path(batches_dir: Path, sample: str, job: BatchJob) -> Path:
    return batches_dir / f"{sample}__{job.tag}.pkl"


def process_one_job(
    job: BatchJob,
    out_dir: Path | str,
    *,
    cfg: EventSelectionBatchedConfig | None = None,
    skip_existing: bool = True,
) -> Optional[str]:
    cfg = cfg or EventSelectionBatchedConfig()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pkl = _batch_out_path(out_dir, job.sample, job)
    if skip_existing and out_pkl.is_file():
        print(f"[batched] {out_pkl} exists; skipping", flush=True)
        return str(out_pkl)

    cmd = [
        cfg.python_executable,
        str(_MAP_SCRIPT),
        "--sample",
        job.sample,
        "--job_id",
        str(job.job_id),
        "--out_dir",
        str(out_dir),
    ]
    for f in job.files:
        cmd.extend(["--df_file", f])
    if cfg.use_mc_genweight:
        cmd.append("--use-mc-genweight")
    if cfg.trace:
        cmd.append("--trace")
    if job.sample == "mc" and cfg.mc_univ_syst:
        cmd.extend(["--mc-univ-syst", ",".join(cfg.mc_univ_syst)])

    print(
        f"[batched] sample={job.sample} {job.tag}  files={len(job.files)}  "
        f"size={job.total_bytes / (1024**3):.3f} GiB",
        flush=True,
    )
    subprocess.run(cmd, check=True)
    return str(out_pkl)


def run_map(
    cfg: EventSelectionBatchedConfig,
    jobs: Sequence[BatchJob] | None = None,
) -> Tuple[Path, Dict[str, List[str]], List[Tuple[str, str, str]], Path]:
    work, batches_dir, _ = cfg.resolve_paths()
    batches_dir.mkdir(parents=True, exist_ok=True)
    if jobs is None:
        _, jobs, manifest_path = discover_jobs(cfg)
    else:
        records = survey_files_for_config(cfg)
        manifest_path = write_manifest(work, records, jobs, max_job_bytes=cfg.max_job_bytes)

    batch_paths: Dict[str, List[str]] = {s: [] for s in cfg.samples}
    failed: List[Tuple[str, str, str]] = []

    print(f"[batched] WORK_BASE={work}", flush=True)
    print(f"[batched] BATCHES_DIR={batches_dir}", flush=True)

    for job in jobs:
        tag = f"{job.sample}/{job.tag}"
        try:
            out = process_one_job(job, batches_dir, cfg=cfg, skip_existing=cfg.skip_existing_batches)
            if out:
                batch_paths.setdefault(job.sample, []).append(out)
        except subprocess.CalledProcessError as exc:
            msg = f"batch job failed exit={exc.returncode}"
            print(f"[batched] FAILED {tag} ({msg})", flush=True)
            failed.append((job.sample, tag, msg))

    for sample, paths in batch_paths.items():
        print(f"[batched] sample={sample} wrote {len(paths)} batch pickle(s)", flush=True)
    return batches_dir, batch_paths, failed, manifest_path


def _load_aggregate_module():
    spec = importlib.util.spec_from_file_location("event_selection_aggregate", _AGG_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {_AGG_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_aggregate(cfg: EventSelectionBatchedConfig, batches_dir: Path | str | None = None) -> EventSelectionBatchedResult:
    work, default_batches, plots_dir = cfg.resolve_paths()
    batches_dir = Path(batches_dir or default_batches)
    plots_dir.mkdir(parents=True, exist_ok=True)

    syst_root = cfg.syst_disk_root
    if syst_root is None:
        syst_root = os.environ.get("NUMUCC_SYST_DISK_ROOT")
    if syst_root is None:
        syst_root = default_syst_disk_root()
    syst_root = Path(syst_root).expanduser()

    agg = _load_aggregate_module()
    chunk_groups = agg.collect_chunks(str(batches_dir))

    samples = {}
    for s, files in chunk_groups.items():
        if not files:
            continue
        print(f"[batched] aggregating {len(files)} batches for sample={s}", flush=True)
        samples[s] = agg.aggregate_chunk_files(files)

    if not samples:
        raise RuntimeError(f"No batch pickles found under {batches_dir}")

    merged = agg.merge_samples(samples)
    n_hd_fixed, n_bar_fixed = agg.sanitize_merged_histdata_finite(merged)
    if n_hd_fixed or n_bar_fixed:
        print(
            f"[batched] sanitized NaN/inf bins: {n_hd_fixed} hist keys, {n_bar_fixed} bar rows",
            flush=True,
        )

    totals = agg.accumulate_exposure_totals(chunk_groups)
    exposure_scales = agg.apply_global_exposure_scales(
        merged, totals, f_offbeam_coincident=cfg.f_offbeam_frac
    )
    print(f"[batched] applied global scales: {exposure_scales}", flush=True)

    data_pot = totals.data_pot if totals.data_pot > 0 else 1.0
    pot_str = agg.get_pot_str(data_pot)
    print(f"[batched] data_pot={data_pot:.3e} -> POT label={pot_str}", flush=True)

    syst_disk_arg = str(syst_root) if syst_root.is_dir() else None
    print(f"[batched] systematics disk root: {syst_disk_arg}", flush=True)

    agg.render_overlay_plots(
        merged,
        plot_label_map={},
        save_fig_dir=str(plots_dir),
        pot_str=pot_str,
        save_fig=cfg.save_fig,
        show_fig=cfg.show_fig,
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
    import pickle

    out_pkl = plots_dir / "merged_histdata.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(merged_payload, f)
    print(f"[batched] wrote {out_pkl}", flush=True)

    batch_paths = {s: sorted(str(p) for p in batches_dir.glob(f"{s}__*.pkl")) for s in SAMPLES}
    manifest_path = work / "manifest.json"
    return EventSelectionBatchedResult(
        work_base=work,
        batches_dir=batches_dir,
        plots_dir=plots_dir,
        batch_paths=batch_paths,
        failed=[],
        merged_payload=merged_payload,
        pot_str=pot_str,
        data_pot=data_pot,
        manifest_path=manifest_path,
    )


def run_full(cfg: EventSelectionBatchedConfig | None = None) -> EventSelectionBatchedResult:
    cfg = cfg or EventSelectionBatchedConfig()
    work, batches_dir, plots_dir = cfg.resolve_paths()
    manifest_path = work / "manifest.json"

    failed: List[Tuple[str, str, str]] = []
    if not cfg.aggregate_only:
        _, _, failed, manifest_path = run_map(cfg)
        if failed:
            print(f"[batched] WARN: {len(failed)} batch job(s) failed", flush=True)
    else:
        print("[batched] aggregate_only=1 — skipping map phase", flush=True)

    if cfg.skip_aggregate:
        print("[batched] skip_aggregate=1 — map pickles only", flush=True)
        return EventSelectionBatchedResult(
            work_base=work,
            batches_dir=batches_dir,
            plots_dir=plots_dir,
            batch_paths={},
            failed=failed,
            merged_payload=None,
            pot_str="",
            data_pot=0.0,
            manifest_path=manifest_path,
        )

    result = run_aggregate(cfg, batches_dir)
    result.failed = failed
    print(f"[batched] DONE plots -> {result.plots_dir}", flush=True)
    return result


def show_saved_plots(plots_dir: Path | str, max_images: int = 20) -> None:
    from IPython.display import Image, display

    plots_dir = Path(plots_dir)
    paths = sorted(plots_dir.glob("*.png"))[:max_images]
    if not paths:
        print(f"[batched] no plots matched under {plots_dir}", flush=True)
        return
    print(f"[batched] displaying {len(paths)} plot(s) from {plots_dir}", flush=True)
    for p in paths:
        display(Image(filename=str(p)))
