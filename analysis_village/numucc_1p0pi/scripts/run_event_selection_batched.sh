#!/usr/bin/env bash
# Batched event selection (map → aggregate → plots).
#
# Thin wrapper around run_event_selection_batched.py (same workflow as the
# notebook's non-interactive batched path).
#
# Override paths / knobs:
#   WORK_BASE=/path/to/out bash run_event_selection_batched.sh
#   MAX_JOB_GB=0.5 bash run_event_selection_batched.sh
#   MAX_FILES_PER_SAMPLE=2 bash run_event_selection_batched.sh
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

mkdir -p "$WORK_BASE" "$BATCHES_DIR" "$PLOTS_DIR"

ARGS=(
  --work-base "$WORK_BASE"
  --batches-dir "$BATCHES_DIR"
  --plots-dir "$PLOTS_DIR"
  --max-job-gb "$MAX_JOB_GB"
  --save-fig
)

if [[ -n "${MAX_FILES_PER_SAMPLE:-}" ]]; then
  ARGS+=(--max-files-per-sample "$MAX_FILES_PER_SAMPLE")
fi
if [[ "${AGGREGATE_ONLY:-0}" == "1" ]]; then
  ARGS+=(--aggregate-only)
fi
if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
  ARGS+=(--skip-aggregate)
fi
if [[ -n "${MC_UNIV_SYST:-}" ]]; then
  ARGS+=(--mc-univ-syst "$MC_UNIV_SYST")
fi
if [[ "${USE_MC_GENWEIGHT:-0}" == "1" ]]; then
  ARGS+=(--use-mc-genweight)
fi
if [[ "${TRACE:-0}" == "1" ]]; then
  ARGS+=(--trace)
fi

echo "[batched] WORK_BASE=$WORK_BASE"
echo "[batched] BATCHES_DIR=$BATCHES_DIR"
echo "[batched] PLOTS_DIR=$PLOTS_DIR"
echo "[batched] MAX_JOB_GB=$MAX_JOB_GB"

exec "$PYTHON" "$THIS_DIR/run_event_selection_batched.py" "${ARGS[@]}"
