#!/usr/bin/env bash
#
# GENIE systematic driver (map → merge per knob group).
# -----------------------------------------------------------------------------
# Phase 1: for each ``(GENIE_GROUP, .df)`` from ``dataset_locations.iter_genie_chunk_map_tasks``,
#          run ``get_systematics_genie.py chunk-map`` with matching ``--genie-group`` and
#          ``--input-stage`` (driven by ``MC_DF_STAGE``).
# Phase 2: for each group in ``GENIE_GROUP_ORDER`` that has chunk pickles under ``CHUNKS_DIR``,
#          ``chunk-merge`` → ``<MERGE_ROOT>/<group>/`` (NPZ optional).
#
# Paths: ``analysis_village.numucc_1p0pi.dataset_locations``:
#   ``MC_DF_STAGE=final|sel_all`` → ``GENIE_GROUP_GLOBS`` vs ``GENIE_GROUP_GLOBS_SEL_ALL``.
#
# ``MC_DF_STAGE`` also drives ``--input-stage`` on the Python driver:
#   ``final``   → chunk reads ``evt``+``mcnu`` only (already selected); full rate+xsec
#                 stack for configured variables (legacy GENIE xsec path unchanged).
#   ``sel_all`` → chunk reads ``evt``+``trk``+``hdr``+``mcnu``, re-runs the selection
#                 pipeline; **cut-stage** observables get **rate** systematics only (no
#                 xsec tensors). Final variables still use the existing GENIE xsec recipe.
#
# Environment:
#   WORK_BASE     Output root (default: ``dataset_locations.default_genie_syst_work_root``)
#   CHUNKS_DIR    Chunk pickle dir (default ``$WORK_BASE/chunks``)
#   MERGE_ROOT    Per-group merge output (default ``$WORK_BASE/merged``)
#   MC_DF_STAGE   ``final`` (default) or ``sel_all``
#   SKIP_MERGE    1 → map only (default: run per-group ``chunk-merge`` after map)
#   XSEC_UNIT     Passed to chunk-merge (default 1.0)
#
# Logs ``progress map overall`` / ``progress merge`` lines (k/total, BEGIN/END timestamps).
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

MC_DF_STAGE=${MC_DF_STAGE:-final}
case "$MC_DF_STAGE" in
    final|sel_all) ;;
    *)
        echo "[genie-run] ERROR: MC_DF_STAGE must be 'final' or 'sel_all', got '$MC_DF_STAGE'" >&2
        exit 2
        ;;
esac
export MC_DF_STAGE
INPUT_STAGE="$MC_DF_STAGE"

mkdir -p "$CHUNKS_DIR" "$MERGE_ROOT"

genie_py="$THIS_DIR/get_systematics_genie.py"
FAILED_LOG="$WORK_BASE/failed_genie_df_files.log"
touch "$FAILED_LOG"

echo "[genie-run] WORK_BASE=$WORK_BASE  CHUNKS_DIR=$CHUNKS_DIR  MERGE_ROOT=$MERGE_ROOT"
echo "[genie-run] MC_DF_STAGE=$MC_DF_STAGE  INPUT_STAGE=$INPUT_STAGE  XSEC_UNIT=$XSEC_UNIT"

_genie_map_total="$(python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_genie_chunk_map_tasks
stage = os.environ.get('MC_DF_STAGE', 'final')
print(sum(1 for _ in iter_genie_chunk_map_tasks(mc_df_stage=stage)))
")"
echo "[genie-run] chunk-map queue: ${_genie_map_total} (.df, group) job(s) for MC_DF_STAGE=$MC_DF_STAGE"

mapfile -t _genie_map_jobs < <(python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_genie_chunk_map_tasks
stage = os.environ.get('MC_DF_STAGE', 'final')
for grp, p in iter_genie_chunk_map_tasks(mc_df_stage=stage):
    print('%s\t%s' % (grp, p))
")
_genie_map_n="${#_genie_map_jobs[@]}"
_genie_map_done=0

for ((_genie_mi = 0; _genie_mi < _genie_map_n; _genie_mi++)); do
    line="${_genie_map_jobs[_genie_mi]}"
    [[ -z "$line" ]] && continue
    IFS=$'\t' read -r grp f <<<"$line"
    if [[ -z "${grp:-}" ]] || [[ -z "${f:-}" ]]; then
        continue
    fi
    _genie_cur=$((_genie_mi + 1))
    out_stem="$(basename "$f" .df)"
    out_pkl="$CHUNKS_DIR/genie__${grp}__${out_stem}.pkl"
    if [[ -f "$out_pkl" ]]; then
        ((_genie_map_done++)) || true
        echo "[genie-run] progress map overall ${_genie_map_done}/${_genie_map_total}  job ${_genie_cur}/${_genie_map_n}  group=$grp  (skip existing) $out_pkl"
        continue
    fi
    echo "[genie-run] progress map overall $((_genie_map_done + 1))/${_genie_map_total}  job ${_genie_cur}/${_genie_map_n}  group=$grp  BEGIN $(date -Is) df_file=$f"
    if python "$genie_py" chunk-map \
        --df-file "$f" \
        --out-dir "$CHUNKS_DIR" \
        --genie-group "$grp" \
        --input-stage "$INPUT_STAGE"; then
        ((_genie_map_done++)) || true
        echo "[genie-run] progress map overall ${_genie_map_done}/${_genie_map_total}  job ${_genie_cur}/${_genie_map_n}  group=$grp  END $(date -Is) df_file=$f"
    else
        ((_genie_map_done++)) || true
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[genie-run] progress map overall ${_genie_map_done}/${_genie_map_total}  job ${_genie_cur}/${_genie_map_n}  group=$grp  FAILED $(date -Is) df_file=$f" >&2
        printf '%s\t%s\t%s\t%s\n' "$ts" "$MC_DF_STAGE" "$grp" "$f" >> "$FAILED_LOG"
    fi
done

if [[ "${SKIP_MERGE:-0}" == "1" ]]; then
    echo "[genie-run] SKIP_MERGE=1 — chunks only → $CHUNKS_DIR"
    exit 0
fi

GENIE_GROUPS=$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_ORDER
print(' '.join(GENIE_GROUP_ORDER))
")

_genie_merge_targets=()
for grp in $GENIE_GROUPS; do
    if compgen -G "$CHUNKS_DIR/genie__${grp}__"*.pkl > /dev/null; then
        _genie_merge_targets+=("$grp")
    else
        echo "[genie-run] skip chunk-merge (no genie__${grp}__*.pkl under $CHUNKS_DIR)"
    fi
done
_genie_merge_n="${#_genie_merge_targets[@]}"
echo "[genie-run] chunk-merge queue: ${_genie_merge_n} knob group(s) with chunk pickles"

for ((_genie_mj = 0; _genie_mj < _genie_merge_n; _genie_mj++)); do
    grp="${_genie_merge_targets[_genie_mj]}"
    _genie_mc=$((_genie_mj + 1))
    out_sub="$MERGE_ROOT/$grp"
    mkdir -p "$out_sub"
    echo "[genie-run] progress merge ${_genie_mc}/${_genie_merge_n}  group=$grp  BEGIN $(date -Is) -> $out_sub"
    python "$genie_py" chunk-merge \
        --chunks-dir "$CHUNKS_DIR" \
        --genie-group "$grp" \
        --input-stage "$INPUT_STAGE" \
        --out-dir "$out_sub" \
        --xsec-unit "$XSEC_UNIT"
    echo "[genie-run] progress merge ${_genie_mc}/${_genie_merge_n}  group=$grp  END $(date -Is)"
done

echo "[genie-run] DONE map=$CHUNKS_DIR merge=$MERGE_ROOT"
