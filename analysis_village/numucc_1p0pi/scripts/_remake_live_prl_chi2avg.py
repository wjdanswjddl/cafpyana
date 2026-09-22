#!/usr/bin/env python3
"""Remake live-PRL event-selection batches with VariableConfig chi2.avg (not I2).

Reads the existing live-PRL manifest (same input .df files) and re-runs map jobs
into a new batches directory. Parallel workers speed up the wall clock.

After completion, point Product A overlays at this batches dir and re-aggregate.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.event_selection_batched import (  # noqa: E402
    BatchJob,
    EventSelectionBatchedConfig,
    process_one_job,
)

SRC_MANIFEST = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL/manifest.json"
)
DST_WORK = Path(
    "/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/"
    "event_selection-batched-live-PRL-chi2avg"
)


def _run_one(job_dict: dict, out_dir: str, skip_existing: bool) -> tuple[str, str, float, str]:
    t0 = time.time()
    job = BatchJob(
        sample=job_dict["sample"],
        job_id=int(job_dict["job_id"]),
        files=list(job_dict["files"]),
        total_bytes=int(job_dict.get("total_bytes") or 0),
    )
    cfg = EventSelectionBatchedConfig(skip_existing_batches=skip_existing)
    try:
        out = process_one_job(
            job, Path(out_dir), cfg=cfg, skip_existing=skip_existing
        )
        return job.sample, job.tag, time.time() - t0, ("ok:" + str(out))
    except Exception as ex:  # noqa: BLE001
        return job.sample, job.tag, time.time() - t0, f"FAIL:{ex}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--skip-existing", action="store_true", default=True)
    p.add_argument("--force", action="store_true", help="Re-run even if pickle exists")
    p.add_argument("--max-jobs", type=int, default=None)
    p.add_argument(
        "--samples",
        default="",
        help="Comma-separated sample filter (default: all)",
    )
    args = p.parse_args()
    skip_existing = False if args.force else True

    dst_batches = DST_WORK / "batches"
    dst_batches.mkdir(parents=True, exist_ok=True)
    (DST_WORK / "manifest.json").write_text(SRC_MANIFEST.read_text())

    jobs = json.loads(SRC_MANIFEST.read_text())["jobs"]
    if args.samples:
        keep = {s.strip() for s in args.samples.split(",") if s.strip()}
        jobs = [j for j in jobs if j["sample"] in keep]
    if args.max_jobs is not None:
        jobs = jobs[: max(0, args.max_jobs)]

    print(
        f"[chi2avg remake] n_jobs={len(jobs)} workers={args.workers} "
        f"out={dst_batches} skip_existing={skip_existing}",
        flush=True,
    )
    t0 = time.time()
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = [
            ex.submit(_run_one, j, str(dst_batches), skip_existing) for j in jobs
        ]
        for i, fut in enumerate(as_completed(futs), 1):
            sample, tag, dt, status = fut.result()
            if status.startswith("ok"):
                ok += 1
            else:
                fail += 1
                print(f"  FAIL {sample}/{tag}: {status}", flush=True)
            if i % 10 == 0 or i == len(futs):
                print(
                    f"  progress {i}/{len(futs)} ok={ok} fail={fail} "
                    f"last={sample}/{tag} {dt:.0f}s elapsed={time.time()-t0:.0f}s",
                    flush=True,
                )
    print(
        f"[chi2avg remake] done ok={ok} fail={fail} elapsed={time.time()-t0:.0f}s",
        flush=True,
    )
    if fail:
        sys.exit(1)


if __name__ == "__main__":
    # Map subprocesses need the repo on PYTHONPATH
    os.environ.setdefault("PYTHONPATH", str(_REPO))
    main()
