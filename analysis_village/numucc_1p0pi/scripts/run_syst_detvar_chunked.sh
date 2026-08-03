#!/usr/bin/env bash
#
# Driver for the chunked detector (calorimetry) unisim systematic.
# -----------------------------------------------------------------------------
# Phase 1 (map):   for each WireMod-tag input directory of variation .df files,
#                  run syst_detvar_chunk.py once per file. Each pickle is
#                  prefixed with the WireMod tag so the aggregator can group
#                  WireMod models without re-parsing paths.
# Phase 2 (reduce): syst_detvar_aggregate.py builds the envelope plots,
#                   the unisim covariance per (WireMod, calo), the per-WireMod
#                   quadrature combination, and the global detector covariance,
#                   and writes ``detector_syst_dict.npz`` -- the same format
#                   loaded by analysis_village.numucc_1p0pi.utils.get_syst_unc,
#                   which is what scripts/event_selection_aggregate.py picks up
#                   at plot time. The NPZ also includes ``detector_by_wiremod``
#                   (``var -> {wiremod_tag -> pack}``), analogous to Flux/G4
#                   ``*_by_knob`` breakdowns. Drop the npz into the path baked into
#                   ``utils.get_syst_unc`` and re-run the main driver to see
#                   detector error bars.
# -----------------------------------------------------------------------------
#
# Edit WIREMOD_DIRS for your inputs. Each entry is "tag|glob" where files
# matching the glob are processed under that WireMod tag.
#
# Override output root (chunk pickles); aggregate NPZ root defaults separately:
#   WORK_BASE=/path/to/out bash run_syst_detvar_chunked.sh
#
# ``NUMUCC_SYST_DISK_ROOT`` overrides where ``Detector/detector_syst_dict.npz`` is written
# (default: ``dataset_locations.default_syst_disk_root()``).
#
# Skip aggregation:
#   SKIP_AGGREGATE=1 bash run_syst_detvar_chunked.sh
#
# Progress: file-level lines match other ``run_syst_*`` drivers (overall k/N, BEGIN/END):
#   [detvar-run] progress map overall 42/3000 (1%) tag=wiremod_yz  BEGIN … stem.df
# Enable tqdm bars inside Python (splits / merge / plots):
#   PYTHON_PROGRESS_BARS=1 bash run_syst_detvar_chunked.sh
#
# Cap how many .df files are processed (**globally**, after MAX_FILES_PER_TAG):
#   MAX_FILES=100 bash run_syst_detvar_chunked.sh
#   bash run_syst_detvar_chunked.sh --max-files 100   # same (-n 100)
# Jobs stay in WIREMOD_DIRS order; within each glob paths are sorted. Phase 2
# still merges every pickle under CHUNKS_DIR — use a fresh WORK_BASE for a clean
# partial-sample aggregation.
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

_usage() {
    cat <<'EOF'
run_syst_detvar_chunked.sh — chunked detector (calo) unisim driver.

Environment:
  WORK_BASE              Chunk pickle root (default: ``default_detvar_syst_work_root``)
  NUMUCC_SYST_DISK_ROOT  Aggregate / syst_disk_layout root (default: ``default_syst_disk_root``)
  MAX_FILES_PER_TAG      Cap files per WireMod glob (0 = all)
  MAX_FILES              Global cap on queued .df jobs after per-tag cap (0 = all)
  MAX_FILES also accepts:  --max-files N   or   -n N   (CLI wins over env)
  SKIP_AGGREGATE=1       Map phase only (skip aggregation / plots)
  PYTHON_PROGRESS_BARS=1 Enable tqdm inside Python (default: off)

Examples:
  MAX_FILES=50 bash run_syst_detvar_chunked.sh
  bash run_syst_detvar_chunked.sh --max-files 50
EOF
}

_CLI_MAX_FILES=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-files|-n)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "[detvar-run] ERROR: $1 requires a non-negative integer" >&2
                exit 1
            fi
            _CLI_MAX_FILES="$2"
            shift 2
            ;;
        --help|-h)
            _usage
            exit 0
            ;;
        *)
            echo "[detvar-run] Unknown option: $1 (use --help)" >&2
            exit 1
            ;;
    esac
done

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_detvar_syst_work_root
print(default_detvar_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR="$WORK_BASE/chunks"
FAILED_LOG="$WORK_BASE/failed_df_files.log"
SYST_DISK_ROOT="$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_root
print(default_syst_disk_root())
")"

# tag|glob (one entry per WireMod model).
declare -a WIREMOD_DIRS=(
    "wiremod_yz|/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModYZ/*.df"
    "wiremod_xtxw|/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_11_103733__sel_2prong-mc-BNB_cosmics-WireModXTXW/*df"
    # add more WireMod variants here, e.g.:
    # "wiremod_xtxw|/path/to/WireModXThetaXW/*.df"
)

SKIP_AGGREGATE=${SKIP_AGGREGATE:-0}
MAX_FILES_PER_TAG=${MAX_FILES_PER_TAG:-0}   # 0 = all (per WireMod glob)
# Global cap after building the job queue (--max-files / -n overrides MAX_FILES env).
MAX_FILES=${_CLI_MAX_FILES:-${MAX_FILES:-0}} # 0 = no extra global limit
PYTHON_PROGRESS_BARS=${PYTHON_PROGRESS_BARS:-0}

PYTHON_PROGRESS_ARGS=()
if [[ "$PYTHON_PROGRESS_BARS" != "1" ]]; then
    PYTHON_PROGRESS_ARGS=(--no-progress)
fi

mkdir -p "$CHUNKS_DIR" "$SYST_DISK_ROOT"
echo "[detvar-run] WORK_BASE=$WORK_BASE  CHUNKS_DIR=$CHUNKS_DIR"
echo "[detvar-run] SYST_DISK_ROOT=$SYST_DISK_ROOT  (Detector/ aggregate target)"
echo "[detvar-run] logging chunk failures to $FAILED_LOG"

# Collect all (tag, path) jobs so we can print global File k/N progress.
declare -a JOB_TAGS=()
declare -a JOB_FILES=()
for entry in "${WIREMOD_DIRS[@]}"; do
    tag="${entry%%|*}"
    glob_pat="${entry#*|}"
    shopt -s nullglob
    files=( $glob_pat )
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "[detvar-run] tag=$tag no files matched ($glob_pat) -- skipping"
        continue
    fi

    mapfile -t sorted_files < <(printf '%s\n' "${files[@]}" | sort)
    files=( "${sorted_files[@]}" )

    if [[ "$MAX_FILES_PER_TAG" != "0" ]]; then
        files=( "${files[@]:0:$MAX_FILES_PER_TAG}" )
    fi

    for f in "${files[@]}"; do
        JOB_TAGS+=("$tag")
        JOB_FILES+=("$f")
    done
done

N_TOTAL=${#JOB_FILES[@]}
if [[ "$N_TOTAL" -eq 0 ]]; then
    echo "[detvar-run] no input .df files matched any WIREMOD_DIRS entry; exiting"
    exit 1
fi

if [[ "$MAX_FILES" =~ ^[0-9]+$ ]] && [[ "$MAX_FILES" -gt 0 ]] && [[ "$MAX_FILES" -lt "$N_TOTAL" ]]; then
    echo "[detvar-run] MAX_FILES=${MAX_FILES}: trimming queue from ${N_TOTAL} file(s)"
    JOB_TAGS=( "${JOB_TAGS[@]:0:MAX_FILES}" )
    JOB_FILES=( "${JOB_FILES[@]:0:MAX_FILES}" )
    N_TOTAL=${#JOB_FILES[@]}
fi

echo "[detvar-run] chunk-map queue: ${N_TOTAL} .df file(s) (Phase 1 map)"

# ---- Phase 1: map ---------------------------------------------------------
for ((i = 0; i < N_TOTAL; i++)); do
    tag="${JOB_TAGS[i]}"
    f="${JOB_FILES[i]}"
    k=$((i + 1))
    pct=$((100 * k / N_TOTAL))
    [[ "$pct" -gt 100 ]] && pct=100
    out_pkl_base="$(basename "$f" .df)"
    out_pkl="$CHUNKS_DIR/${tag}__${out_pkl_base}.pkl"

    if [[ -f "$out_pkl" ]]; then
        echo "[detvar-run] progress map overall ${k}/${N_TOTAL} (${pct}%) tag=${tag}  (skip existing) ${out_pkl_base}.df"
        continue
    fi

    echo "[detvar-run] progress map overall ${k}/${N_TOTAL} (${pct}%) tag=${tag}  BEGIN $(date -Is) ${out_pkl_base}.df"
    if ! python "$THIS_DIR/syst_detvar_chunk.py" \
        --df_file "$f" \
        --wiremod_tag "$tag" \
        --out_dir "$CHUNKS_DIR" \
        "${PYTHON_PROGRESS_ARGS[@]}"; then
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[detvar-run] progress map overall ${k}/${N_TOTAL} (${pct}%) tag=${tag}  FAILED $(date -Is) $f" >&2
        echo "[detvar-run] FAILED chunk (see $FAILED_LOG): $f" >&2
        printf '%s\t%s\t%s\n' "$ts" "$tag" "$f" >> "$FAILED_LOG"
        continue
    fi
    echo "[detvar-run] progress map overall ${k}/${N_TOTAL} (${pct}%) tag=${tag}  END $(date -Is) ${out_pkl_base}.df"
done

echo "[detvar-run] Phase 1 complete (${N_TOTAL} file job(s) considered)"

# ---- Phase 2: reduce -------------------------------------------------------
if [[ "$SKIP_AGGREGATE" != "1" ]]; then
    echo "[detvar-run] Phase 2 (aggregate): merging pickles under $CHUNKS_DIR"
    echo "[detvar-run] progress aggregate 1/1  BEGIN $(date -Is) syst_detvar_aggregate.py"
    python "$THIS_DIR/syst_detvar_aggregate.py" \
        --in_dir "$CHUNKS_DIR" \
        --syst-disk-root "$SYST_DISK_ROOT" \
        "${PYTHON_PROGRESS_ARGS[@]}"
    echo "[detvar-run] progress aggregate 1/1  END $(date -Is)"
    echo "[detvar-run] DONE detector NPZ -> $SYST_DISK_ROOT/Detector/detector_syst_dict.npz"

    NPZ_SRC="$SYST_DISK_ROOT/Detector/detector_syst_dict.npz"
    if [[ -f "$NPZ_SRC" ]]; then
        echo "[detvar-run] Merge with multisim/cosmics/GENIE outputs under one NUMUCC_SYST_DISK_ROOT for utils.get_syst_unc."
    fi
else
    echo "[detvar-run] SKIP_AGGREGATE=1 -> chunks only in $CHUNKS_DIR"
fi
