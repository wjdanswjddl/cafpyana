#!/usr/bin/env bash
#
# Smoke-test driver: run event_selection_chunk.py on exactly ONE .df file per
# sample type (mc, data, intime, offbeam, dirt), then aggregate + plot.
#
# Same globs as run_event_selection_chunked.sh; the first match after sorting
# paths lexicographically is used so the choice is stable.
#
# Override output root:
#   WORK_BASE=/path/to/out bash run_event_selection_chunked_example.sh
#
# Extra aggregator flags:
#   AGG_EXTRA_ARGS="--cosmic_estimate offbeam" bash run_event_selection_chunked_example.sh
#
# Skip aggregation (only write chunk pickles):
#   SKIP_AGGREGATE=1 bash run_event_selection_chunked_example.sh
#
# Chunk logs RSS/VmHWM by default (--no_mem_diag to silence). Notebook-like load:
#   .../event_selection_chunk.py ... --load_mode concat
# Cap splits: ... --max_splits 3
# Mystery kills / no traceback: add  --trace  (writes chunk_trace__*.log under chunks dir).
# -----------------------------------------------------------------------------

set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-"/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-example-$TODAY"}
CHUNKS_DIR="$WORK_BASE/chunks"
PLOTS_DIR="$WORK_BASE/plots_example"

declare -a SAMPLE_DIRS=(
    "mc|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/BNB_cosmics/*-sel_all-wgts.df"
    "data|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/_Fixed_all.df"
    "intime|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/intime/*_all.df"
    "offbeam|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/OffBeam/*_all.df"
    "dirt|/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/lowE/*_all.df"
)

SKIP_AGGREGATE=${SKIP_AGGREGATE:-0}
AGG_EXTRA_ARGS=${AGG_EXTRA_ARGS:-}

mkdir -p "$CHUNKS_DIR" "$PLOTS_DIR"

echo "[example] WORK_BASE=$WORK_BASE"

# ---- Phase 1: one chunk per sample ----------------------------------------
for entry in "${SAMPLE_DIRS[@]}"; do
    sample="${entry%%|*}"
    glob_pat="${entry#*|}"
    shopt -s nullglob
    matches=( $glob_pat )
    shopt -u nullglob

    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "[example] sample=$sample no files matched ($glob_pat) -- skipping"
        continue
    fi

    mapfile -t sorted_matches < <(printf '%s\n' "${matches[@]}" | sort)
    f="${sorted_matches[0]}"
    echo "[example] sample=$sample using ONE file: $f"

    out_pkl_base="$(basename "$f" .df)"
    out_pkl="$CHUNKS_DIR/${sample}__${out_pkl_base}.pkl"
    if [[ -f "$out_pkl" ]]; then
        echo "[example] $out_pkl exists; skipping chunk (delete to re-run)"
        continue
    fi

    python "$THIS_DIR/event_selection_chunk.py" \
        --df_file "$f" \
        --sample "$sample" \
        --out_dir "$CHUNKS_DIR"
done

# ---- Phase 2: reduce + plot -----------------------------------------------
if [[ "$SKIP_AGGREGATE" != "1" ]]; then
    echo "[example] aggregating + plotting -> $PLOTS_DIR"
    # shellcheck disable=SC2086
    python "$THIS_DIR/event_selection_aggregate.py" \
        --in_dir "$CHUNKS_DIR" \
        --out_dir "$PLOTS_DIR" \
        $AGG_EXTRA_ARGS
    echo "[example] DONE plots + merged_histdata.pkl -> $PLOTS_DIR"
else
    echo "[example] SKIP_AGGREGATE=1 -> chunks only in $CHUNKS_DIR"
fi
