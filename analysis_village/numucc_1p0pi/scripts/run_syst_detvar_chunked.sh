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
#                   at plot time. Drop the npz into the path baked into
#                   ``utils.get_syst_unc`` and re-run the main driver to see
#                   detector error bars.
# -----------------------------------------------------------------------------
#
# Edit WIREMOD_DIRS for your inputs. Each entry is "tag|glob" where files
# matching the glob are processed under that WireMod tag.
#
# Override output root:
#   WORK_BASE=/path/to/out bash run_syst_detvar_chunked.sh
#
# Skip aggregation:
#   SKIP_AGGREGATE=1 bash run_syst_detvar_chunked.sh
#
# Progress: by default this script prints only **file-level** progress:
#   [detvar-run] File 42/3000 (2%) tag=wiremod_yz  stem.df
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

_usage() {
    cat <<'EOF'
run_syst_detvar_chunked.sh — chunked detector (calo) unisim driver.

Environment:
  WORK_BASE              Output root (default: ~/xsec/.../detvar_systematics-<date>)
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
WORK_BASE=${WORK_BASE:-"/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/detvar_systematics-$TODAY"}
CHUNKS_DIR="$WORK_BASE/chunks"
FAILED_LOG="$WORK_BASE/failed_df_files.log"

# tag|glob (one entry per WireMod model).
declare -a WIREMOD_DIRS=(
    "wiremod_yz|/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModYZ/*.df"
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

mkdir -p "$CHUNKS_DIR"
echo "[detvar-run] WORK_BASE=$WORK_BASE"
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

echo "[detvar-run] Phase 1 (map): ${N_TOTAL} .df file(s) queued"

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
        echo "[detvar-run] File ${k}/${N_TOTAL} (${pct}%) tag=${tag}  SKIP (exists) ${out_pkl_base}.df"
        continue
    fi

    echo "[detvar-run] File ${k}/${N_TOTAL} (${pct}%) tag=${tag}  RUN ${out_pkl_base}.df"
    if ! python "$THIS_DIR/syst_detvar_chunk.py" \
        --df_file "$f" \
        --wiremod_tag "$tag" \
        --out_dir "$CHUNKS_DIR" \
        "${PYTHON_PROGRESS_ARGS[@]}"; then
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[detvar-run] FAILED chunk (see $FAILED_LOG): $f" >&2
        printf '%s\t%s\t%s\n' "$ts" "$tag" "$f" >> "$FAILED_LOG"
        continue
    fi
done

echo "[detvar-run] Phase 1 complete (${N_TOTAL} file job(s) considered)"

# ---- Phase 2: reduce -------------------------------------------------------
if [[ "$SKIP_AGGREGATE" != "1" ]]; then
    echo "[detvar-run] Phase 2 (aggregate): merging pickles under $CHUNKS_DIR"
    python "$THIS_DIR/syst_detvar_aggregate.py" \
        --in_dir "$CHUNKS_DIR" \
        --syst-disk-root "$WORK_BASE" \
        "${PYTHON_PROGRESS_ARGS[@]}"
    echo "[detvar-run] DONE detector NPZ -> $WORK_BASE/Detector/detector_syst_dict.npz"

    NPZ_SRC="$WORK_BASE/Detector/detector_syst_dict.npz"
    if [[ -f "$NPZ_SRC" ]]; then
        echo "[detvar-run] Merge with multisim/cosmics/GENIE outputs under one NUMUCC_SYST_DISK_ROOT for utils.get_syst_unc."
    fi
else
    echo "[detvar-run] SKIP_AGGREGATE=1 -> chunks only in $CHUNKS_DIR"
fi
