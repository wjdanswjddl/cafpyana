#!/usr/bin/env bash
#
# Chunked event selection on an explicit subset of .df files.
#
# For the full sample (all .df under dataset_locations.py directories), use:
#   bash run_event_selection_chunked.sh
#
# This script is for manifests or CLI sample/path pairs only.
#
# Usage
# -----
#   bash run_event_selection_paths.sh /path/to/manifest.txt
#   bash run_event_selection_paths.sh mc /path/a.df mc /path/b.df data /path/c.df
#
# Manifest format (one per line):  sample<TAB or space>path/to/file.df
#
#   AGGREGATE_ONLY=1   skip map (chunks already in CHUNKS_DIR)
#   SKIP_AGGREGATE=1   chunks only
#
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -eq 0 ]]; then
    echo "This script needs a manifest file or sample/path pairs." >&2
    echo "For all files under dataset_locations.EVENT_SELECTION_GLOBS, run:" >&2
    echo "  bash run_event_selection_chunked.sh" >&2
    exit 1
fi

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-"/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-paths-$TODAY"}
CHUNKS_DIR="${CHUNKS_DIR:-$WORK_BASE/chunks}"
PLOTS_DIR="${PLOTS_DIR:-$WORK_BASE/plots}"

source "$THIS_DIR/event_selection_chunked_lib.sh"

declare -a JOB_ENTRIES=()

_parse_manifest_lines() {
    local line sample f
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line#"${line%%[![:space:]]*}"}"
        [[ -z "$line" ]] && continue
        if [[ "$line" == *$'\t'* ]]; then
            sample="${line%%$'\t'*}"
            f="${line#*$'\t'}"
        else
            sample="${line%% *}"
            f="${line#* }"
        fi
        sample="${sample#"${sample%%[![:space:]]*}"}"
        sample="${sample%"${sample##*[![:space:]]}"}"
        f="${f#"${f%%[![:space:]]*}"}"
        f="${f%"${f##*[![:space:]]}"}"
        [[ -z "$sample" || -z "$f" ]] && continue
        JOB_ENTRIES+=("${sample}|${f}")
    done
}

if [[ $# -eq 1 && -f "$1" ]]; then
    _parse_manifest_lines <"$1"
elif [[ $# -ge 2 && $(( $# % 2 )) -eq 0 ]]; then
    while [[ $# -ge 2 ]]; do
        JOB_ENTRIES+=("$1|$2")
        shift 2
    done
else
    echo "Usage: $0 manifest.txt" >&2
    echo "   or: $0 sample path.df [sample path.df ...]" >&2
    exit 1
fi

if [[ ${#JOB_ENTRIES[@]} -eq 0 ]]; then
    echo "[run] no jobs in manifest" >&2
    exit 1
fi

if [[ "${AGGREGATE_ONLY:-0}" != "1" ]]; then
    event_selection_run_map "${JOB_ENTRIES[@]}"
fi

if [[ "${SKIP_AGGREGATE:-0}" != "1" ]]; then
    event_selection_run_aggregate
fi
