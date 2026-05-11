#!/usr/bin/env bash
#
# GENIE systematic driver (map → merge per knob group).
# -----------------------------------------------------------------------------
# Phase 1: for each ``(GENIE_GROUP, .df)`` from ``dataset_locations.iter_genie_chunk_map_tasks``,
#          run ``get_systematics_genie.py chunk-map`` with matching ``--genie-group``.
# Phase 2: for each GENIE group, ``chunk-merge`` → ``<MERGE_ROOT>/<group>/`` (NPZ optional).
#
# Paths: ``analysis_village.numucc_1p0pi.dataset_locations`` (``GENIE_GROUP_GLOBS``).
#
# Environment:
#   WORK_BASE     Output root (default: ``dataset_locations.default_genie_syst_work_root``)
#   CHUNKS_DIR    Chunk pickle dir (default ``$WORK_BASE/chunks``)
#   MERGE_ROOT    Per-group merge output (default ``$WORK_BASE/merged``)
#   SKIP_MERGE    1 → map only
#   XSEC_UNIT     Passed to chunk-merge (default 1.0)
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_genie_syst_work_root
print(default_genie_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR=${CHUNKS_DIR:-"$WORK_BASE/chunks"}
MERGE_ROOT=${MERGE_ROOT:-"$WORK_BASE/merged"}
XSEC_UNIT=${XSEC_UNIT:-1.0}

mkdir -p "$CHUNKS_DIR" "$MERGE_ROOT"

genie_py="$THIS_DIR/get_systematics_genie.py"
FAILED_LOG="$WORK_BASE/failed_genie_df_files.log"
touch "$FAILED_LOG"

echo "[genie-run] WORK_BASE=$WORK_BASE  CHUNKS_DIR=$CHUNKS_DIR  MERGE_ROOT=$MERGE_ROOT"

while IFS=$'\t' read -r grp f; do
    [[ -z "${grp:-}" ]] && continue
    out_stem="$(basename "$f" .df)"
    out_pkl="$CHUNKS_DIR/genie__${grp}__${out_stem}.pkl"
    if [[ -f "$out_pkl" ]]; then
        echo "[genie-run] skip existing $out_pkl"
        continue
    fi
    if python "$genie_py" chunk-map \
        --df-file "$f" \
        --out-dir "$CHUNKS_DIR" \
        --genie-group "$grp"; then
        :
    else
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[genie-run] FAILED: group=$grp file=$f" >&2
        printf '%s\t%s\t%s\n' "$ts" "$grp" "$f" >> "$FAILED_LOG"
    fi
done < <(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_genie_chunk_map_tasks
for grp, p in iter_genie_chunk_map_tasks():
    print('%s\t%s' % (grp, p))
")

if [[ "${SKIP_MERGE:-0}" == "1" ]]; then
    echo "[genie-run] SKIP_MERGE=1 — chunks only → $CHUNKS_DIR"
    exit 0
fi

for grp in CCQE MEC DIS Other Ar23p; do
    out_sub="$MERGE_ROOT/$grp"
    mkdir -p "$out_sub"
    echo "[genie-run] chunk-merge group=$grp -> $out_sub"
    python "$genie_py" chunk-merge \
        --chunks-dir "$CHUNKS_DIR" \
        --genie-group "$grp" \
        --out-dir "$out_sub" \
        --xsec-unit "$XSEC_UNIT"
done

echo "[genie-run] DONE map=$CHUNKS_DIR merge=$MERGE_ROOT"
