#!/usr/bin/env bash
# Batched event selection: group input .df files into ≤1 GiB jobs, run notebook
# pipeline per job, aggregate histograms, render final plots.
#
# Map:  event_selection_batch_map.py  (one job = multiple files, sequential load)
# Reduce: event_selection_aggregate.py (unchanged pickle format)
#
# Override paths:
#   WORK_BASE=/path/to/out bash run_event_selection_batched.sh
#   MAX_JOB_GB=0.5 bash run_event_selection_batched.sh
#   SKIP_AGGREGATE=1 bash run_event_selection_batched.sh
#   AGGREGATE_ONLY=1 bash run_event_selection_batched.sh

set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
TODAY="${TODAY:-$(date +%Y%m%d)}"
WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-batched-${TODAY}}"
BATCHES_DIR="${BATCHES_DIR:-${WORK_BASE}/batches}"
PLOTS_DIR="${PLOTS_DIR:-${WORK_BASE}/plots}"
MAX_JOB_GB="${MAX_JOB_GB:-1.0}"
PYTHON="${PYTHON:-${REPO_ROOT}/envs/venv_py310_cafpyana/bin/python}"
WORKERS="${WORKERS:-1}"
SKIP_AGGREGATE="${SKIP_AGGREGATE:-0}"
AGGREGATE_ONLY="${AGGREGATE_ONLY:-0}"

mkdir -p "$WORK_BASE" "$BATCHES_DIR" "$PLOTS_DIR"

echo "[batched] WORK_BASE=$WORK_BASE"
echo "[batched] BATCHES_DIR=$BATCHES_DIR"
echo "[batched] PLOTS_DIR=$PLOTS_DIR"
echo "[batched] MAX_JOB_GB=$MAX_JOB_GB"

if [[ "$AGGREGATE_ONLY" != "1" ]]; then
  "$PYTHON" "$THIS_DIR/event_selection_batch_survey.py" \
    --work_dir "$WORK_BASE" \
    --max_job_gb "$MAX_JOB_GB"

  "$PYTHON" - <<PY
from analysis_village.numucc_1p0pi.event_selection_batched import (
    EventSelectionBatchedConfig,
    run_map,
)
import json
from pathlib import Path

work = Path("${WORK_BASE}")
with open(work / "manifest.json") as f:
    manifest = json.load(f)

from analysis_village.numucc_1p0pi.event_selection_batched import BatchJob

jobs = [
    BatchJob(
        sample=j["sample"],
        job_id=j["job_id"],
        files=j["files"],
        total_bytes=j["total_bytes"],
    )
    for j in manifest["jobs"]
]

cfg = EventSelectionBatchedConfig(
    work_base=work,
    batches_dir=Path("${BATCHES_DIR}"),
    max_job_bytes=int(float("${MAX_JOB_GB}") * 1024**3),
    skip_existing_batches=True,
)
run_map(cfg, jobs=jobs)
PY
fi

if [[ "$SKIP_AGGREGATE" == "1" ]]; then
  echo "[batched] SKIP_AGGREGATE=1 — map pickles only"
  exit 0
fi

"$PYTHON" "$THIS_DIR/event_selection_aggregate.py" \
  --in_dir "$BATCHES_DIR" \
  --out_dir "$PLOTS_DIR" \
  --save_fig

echo "[batched] DONE plots -> $PLOTS_DIR"
