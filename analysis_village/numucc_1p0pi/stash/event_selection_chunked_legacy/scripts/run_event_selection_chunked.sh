#!/usr/bin/env bash
#
# Driver script for the map/reduce event-selection framework.
# -----------------------------------------------------------------------------
# Phase 1 (map):   processes each CAF .df shard with event_selection_chunk.py,
#                  writing one pickle per (sample, file stem). HDF5 splits inside a
#                  file are read sequentially (not fully concatenated in RAM).
#                  MC chunks pass --mc-univ-syst (default Flux,G4,GENIE) for hatched
#                  syst bands; see event_selection_chunked_lib.sh.
# Phase 2 (reduce): aggregates pickles and applies exposure normalization from
#                   summed shard metadata, then runs event_selection_aggregate.py.
#                   Uses NUMUCC_SYST_DISK_ROOT (default: dataset_locations default)
#                   as NPZ fallback via --syst-disk-root.
#
# Sample globs: ``analysis_village.numucc_1p0pi.dataset_locations.EVENT_SELECTION_GLOBS``.
# Override ``SPRING_GEN1_ROOT`` / ``NUMUCC_SPRING_GEN1_ROOT`` there or via env before import.
# Discovers every ``*.df`` under each directory in EVENT_SELECTION_GLOBS
# (see dataset_locations.py). For a hand-picked subset, use run_event_selection_paths.sh.
# For grid jobs, submit one chunk job per file.
# -----------------------------------------------------------------------------
#
# Environment overrides: WORK_BASE, PLOTS_DIR, MC_UNIV_SYST, NUMUCC_SYST_DISK_ROOT,
#   AGG_EXTRA_ARGS, AGGREGATE_ONLY=1, SKIP_AGGREGATE=1
#
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-"/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-chunked-$TODAY"}
CHUNKS_DIR="${CHUNKS_DIR:-$WORK_BASE/chunks}"
#PLOTS_DIR="$WORK_BASE/plots-nominal"
#PLOTS_DIR="$WORK_BASE/plots-EField_R00"
PLOTS_DIR="${PLOTS_DIR:-$WORK_BASE/plots-EField_R30_Short}"

source "$THIS_DIR/event_selection_chunked_lib.sh"

declare -a JOB_ENTRIES=()
while IFS= read -r line; do
    JOB_ENTRIES+=("$line")
done < <(event_selection_discover_jobs_from_dataset_locations)

if [[ ${#JOB_ENTRIES[@]} -eq 0 ]]; then
    echo "[run] no input .df files matched EVENT_SELECTION_GLOBS — nothing to do" >&2
    exit 1
fi

if [[ "${AGGREGATE_ONLY:-0}" != "1" ]]; then
    event_selection_run_map "${JOB_ENTRIES[@]}"
fi

if [[ "${SKIP_AGGREGATE:-0}" != "1" ]]; then
    event_selection_run_aggregate
fi
