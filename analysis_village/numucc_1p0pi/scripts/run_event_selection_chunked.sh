#!/usr/bin/env bash
#
# Driver script for the chunked event-selection framework.
# -----------------------------------------------------------------------------
# Phase 1 (map):   processes each .df file with event_selection_chunk.py,
#                  writing one pickle per (sample, chunk). HDF5 splits inside a
#                  file are read sequentially (not fully concatenated in RAM).
#                  Each pickle stores intrinsic weights plus POT/gates metadata.
# Phase 2 (reduce): aggregates pickles and applies exposure normalization from
#                   summed chunk metadata, then runs event_selection_aggregate.py.
#
# Edit SAMPLE_DIRS for your paths. For grid jobs, submit one chunk job per file.
# -----------------------------------------------------------------------------

set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-"/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-chunked-$TODAY"}
CHUNKS_DIR="$WORK_BASE/chunks"
PLOTS_DIR="$WORK_BASE/plots"
# Failed chunk inputs (unreadable/missing HDF, etc.) are appended here; the driver keeps going.
FAILED_LOG="$WORK_BASE/failed_df_files.log"

declare -a SAMPLE_DIRS=(
    "mc|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/BNB_cosmics/*-sel_all-wgts.df"
    "data|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/_Fixed_all.df"
    "intime|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/intime/*_all.df"
    "offbeam|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/OffBeam/*_all.df"
    "dirt|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/lowE/*_all.df"
)

mkdir -p "$CHUNKS_DIR" "$PLOTS_DIR"

# ---- Phase 1: map ---------------------------------------------------------
echo "[run] logging chunk failures to $FAILED_LOG"
for entry in "${SAMPLE_DIRS[@]}"; do
    sample="${entry%%|*}"
    glob_pat="${entry#*|}"
    files=( $glob_pat )
    if [[ ${#files[@]} -eq 0 || ! -e "${files[0]}" ]]; then
        echo "[run] sample=$sample no files matched ($glob_pat) -- skipping"
        continue
    fi
    echo "[run] sample=$sample  ${#files[@]} files"
    for f in "${files[@]}"; do
        out_pkl_base="$(basename "$f" .df)"
        out_pkl="$CHUNKS_DIR/${sample}__${out_pkl_base}.pkl"
        if [[ -f "$out_pkl" ]]; then
            echo "[run] $out_pkl exists; skipping"
            continue
        fi
        if ! python "$THIS_DIR/event_selection_chunk.py" \
            --df_file "$f" \
            --sample "$sample" \
            --out_dir "$CHUNKS_DIR"; then
            ts="$(date '+%Y-%m-%d %H:%M:%S')"
            echo "[run] FAILED chunk (see $FAILED_LOG): $f" >&2
            printf '%s\t%s\t%s\n' "$ts" "$sample" "$f" >> "$FAILED_LOG"
            continue
        fi
    done
done

# ---- Phase 2: reduce + plot -----------------------------------------------
# Optional: --cosmic_estimate offbeam  --hide_cosmic_model_unc  etc.
echo "[run] aggregating + plotting"
python "$THIS_DIR/event_selection_aggregate.py" \
    --in_dir "$CHUNKS_DIR" \
    --out_dir "$PLOTS_DIR"

echo "[run] DONE -> $PLOTS_DIR"
