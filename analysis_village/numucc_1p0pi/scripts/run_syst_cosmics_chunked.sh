#!/usr/bin/env bash
#
# Cosmics systematic driver (offbeam + intime map → aggregate).
# -----------------------------------------------------------------------------
# Phase 1: ``syst_cosmics_chunk.py`` for every ``.df`` under
# ``dataset_locations.EVENT_SELECTION_GLOBS`` offbeam / intime.
# Phase 2: ``syst_cosmics_aggregate.py`` → ``<SYST_DISK_ROOT>/Cosmics/cosmics_syst_dict.npz``.
#
# Environment:
#   WORK_BASE          Default chunk dir root (``default_cosmics_syst_work_root``)
#   CHUNKS_DIR         (default ``$WORK_BASE/chunks``)
#   NUMUCC_SYST_DISK_ROOT  Required for aggregate unless passed as first extra arg
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_cosmics_syst_work_root
print(default_cosmics_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR=${CHUNKS_DIR:-"$WORK_BASE/chunks"}
SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-${1:-}}"

mkdir -p "$CHUNKS_DIR"

chunk_py="$THIS_DIR/syst_cosmics_chunk.py"
agg_py="$THIS_DIR/syst_cosmics_aggregate.py"

echo "[cosmics-run] CHUNKS_DIR=$CHUNKS_DIR"

for sample in offbeam intime; do
    while IFS= read -r f; do
        [[ -z "$f" ]] && continue
        out_stem="$(basename "$f" .df)"
        out_pkl="$CHUNKS_DIR/cosmics__${sample}__${out_stem}.pkl"
        if [[ -f "$out_pkl" ]]; then
            echo "[cosmics-run] skip existing $out_pkl"
            continue
        fi
        python "$chunk_py" --sample "$sample" --df_file "$f" --out_dir "$CHUNKS_DIR"
    done < <(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_cosmics_chunk_df_paths
for p in iter_cosmics_chunk_df_paths('${sample}'):
    print(p)
")
done

if [[ -z "$SYST_DISK_ROOT" ]]; then
    echo "Set NUMUCC_SYST_DISK_ROOT or pass syst disk root as \$1 for aggregate." >&2
    exit 1
fi

python "$agg_py" --chunks_dir "$CHUNKS_DIR" --syst-disk-root "$SYST_DISK_ROOT"
echo "[cosmics-run] DONE -> ${SYST_DISK_ROOT}/Cosmics/"
