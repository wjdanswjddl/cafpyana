#!/usr/bin/env python
"""CLI for the batched event-selection workflow (notebook equivalent).

Surveys ``EVENT_SELECTION_GLOBS``, packs ``.df`` files into ≤1 GiB jobs, runs the
notebook selection pipeline per job, aggregates histograms, and renders plots.

This replaces the orchestration formerly only in
``notebooks/event_selection.ipynb`` and ``run_event_selection_batched.sh``.

Examples
--------
    # Full sample (default 1 GiB jobs)
    python run_event_selection_batched.py --plots-dir /path/to/plots

    # Smoke test: 2 files per sample
    python run_event_selection_batched.py --max-files-per-sample 2

    # Aggregate only (map pickles already present)
    python run_event_selection_batched.py --aggregate-only
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.event_selection_batched import (  # noqa: E402
    EventSelectionBatchedConfig,
    discover_jobs,
    run_full,
    default_event_selection_batched_work_root,
)
from analysis_village.numucc_1p0pi.dataset_locations import PLOTS_BASE  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Batched event selection: map → aggregate → plots "
        "(same workflow as notebooks/event_selection.ipynb batched cells)."
    )
    p.add_argument(
        "--work-base",
        default=None,
        help="Working directory for manifest + batches "
        "(default: /exp/sbnd/data/users/$USER/.../event_selection-batched-$TODAY)",
    )
    p.add_argument(
        "--batches-dir",
        default=None,
        help="Directory for per-job pickles (default: <work-base>/batches)",
    )
    p.add_argument(
        "--plots-dir",
        default=None,
        help=f"Output plots directory (default: <work-base>/plots or under {PLOTS_BASE})",
    )
    p.add_argument(
        "--max-job-gb",
        type=float,
        default=1.0,
        help="Max input size per map job in GiB (default: 1.0)",
    )
    p.add_argument(
        "--max-files-per-sample",
        type=int,
        default=None,
        help="Cap number of .df files per sample (smoke / partial runs)",
    )
    p.add_argument(
        "--mc-univ-syst",
        default="",
        help="Comma-separated MC universe tags for map bookkeeping (e.g. Flux,G4,GENIE)",
    )
    p.add_argument("--use-mc-genweight", action="store_true")
    p.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="Re-run map jobs even if output pickles exist",
    )
    p.add_argument("--aggregate-only", action="store_true")
    p.add_argument("--skip-aggregate", action="store_true")
    p.add_argument("--save-fig", action="store_true", default=True)
    p.add_argument("--no-save-fig", action="store_true")
    p.add_argument("--show-fig", action="store_true")
    p.add_argument("--syst-disk-root", default=None)
    p.add_argument("--syst-tag", default="")
    p.add_argument("--f-offbeam-frac", type=float, default=0.08)
    p.add_argument("--cosmic-estimate", choices=("intime", "offbeam"), default="intime")
    p.add_argument("--trace", action="store_true")
    p.add_argument(
        "--discover-only",
        action="store_true",
        help="Survey files, write manifest, print jobs, then exit",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    today = __import__("datetime").datetime.now().strftime("%Y%m%d")
    work_base = args.work_base or str(default_event_selection_batched_work_root(today))
    plots_dir = args.plots_dir
    if plots_dir is None:
        plots_dir = str(Path(work_base) / "plots")

    mc_univ = tuple(x.strip() for x in (args.mc_univ_syst or "").split(",") if x.strip())
    cfg = EventSelectionBatchedConfig(
        work_base=work_base,
        batches_dir=args.batches_dir,
        plots_dir=plots_dir,
        max_job_bytes=int(float(args.max_job_gb) * 1024**3),
        max_files_per_sample=args.max_files_per_sample,
        mc_univ_syst=mc_univ,
        use_mc_genweight=args.use_mc_genweight,
        skip_existing_batches=not args.no_skip_existing,
        aggregate_only=args.aggregate_only,
        skip_aggregate=args.skip_aggregate,
        save_fig=(False if args.no_save_fig else args.save_fig),
        show_fig=args.show_fig,
        syst_disk_root=args.syst_disk_root,
        syst_tag=args.syst_tag,
        f_offbeam_frac=args.f_offbeam_frac,
        cosmic_estimate=args.cosmic_estimate,
        trace=args.trace,
    )

    records, jobs, manifest_path = discover_jobs(cfg)
    print(f"[run_event_selection_batched] {len(jobs)} job(s)  manifest={manifest_path}")
    print(f"[run_event_selection_batched] files surveyed: {len(records)}")
    if args.discover_only:
        return 0

    result = run_full(cfg)
    print("[run_event_selection_batched] DONE")
    print("  batches :", result.batches_dir)
    print("  plots   :", result.plots_dir)
    print("  POT     :", result.pot_str)
    if result.failed:
        print(f"  FAILED  : {len(result.failed)} job(s)")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
