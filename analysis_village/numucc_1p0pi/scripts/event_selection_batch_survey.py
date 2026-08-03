#!/usr/bin/env python
"""Survey event-selection input files and write a batch job manifest."""
from __future__ import annotations

import argparse
import sys
from os import path

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from analysis_village.numucc_1p0pi.event_selection_batched import (  # noqa: E402
    DEFAULT_MAX_JOB_BYTES,
    EventSelectionBatchedConfig,
    discover_jobs,
    group_files_into_jobs,
    print_survey_summary,
    survey_files,
    write_manifest,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--work_dir",
        required=True,
        help="Directory for manifest.json (and default batch output root parent)",
    )
    p.add_argument(
        "--max_job_gb",
        type=float,
        default=DEFAULT_MAX_JOB_BYTES / (1024.0 ** 3),
        help="Maximum total input size per job (GiB)",
    )
    p.add_argument("--max_files_per_sample", type=int, default=0, help="0 = all files")
    return p.parse_args()


def main():
    args = parse_args()
    max_bytes = int(args.max_job_gb * (1024.0 ** 3))
    max_files = None if args.max_files_per_sample <= 0 else args.max_files_per_sample

    records = survey_files(max_files_per_sample=max_files)
    jobs = group_files_into_jobs(records, max_bytes=max_bytes)
    manifest = write_manifest(args.work_dir, records, jobs)
    print_survey_summary(records, jobs)
    print(f"[survey] wrote {manifest}", flush=True)


if __name__ == "__main__":
    main()
