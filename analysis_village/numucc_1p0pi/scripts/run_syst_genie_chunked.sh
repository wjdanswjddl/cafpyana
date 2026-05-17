#!/usr/bin/env bash
#
# GENIE systematic driver (map → merge per knob group).
# -----------------------------------------------------------------------------
# Phase 1: for each ``(GENIE_GROUP, .df)`` from ``dataset_locations.iter_genie_chunk_map_tasks``,
#          run ``get_systematics_genie.py chunk-map`` with matching ``--genie-group`` and
#          ``--input-stage`` (driven by ``MC_DF_STAGE``).
# Phase 2: for each group in the active glob map (optionally filtered by ``GENIE_RUN_GROUPS``)
#          that has chunk pickles under ``CHUNKS_DIR``, ``chunk-merge`` → ``<MERGE_ROOT>/<group>/``.
#          Merge output dict / NPZ ``syst`` includes per-knob entries plus
#          ``__GENIE_group_combined__`` → per-variable ``{rate, xsec?}`` with cov / cov_frac / corr
#          (independent-sum across knobs in that group, same recipe as multisim Flux/G4 aggregate).
#
# Phase 3: ``syst_genie_aggregate.py`` folds chunk pickles into the unified syst-disk tree:
#          ``<SYST_DISK_ROOT>/GENIE/cov_mat_dict.pkl`` (same role as ``syst_multisim_aggregate.py``
#          for Flux/G4). Always runs after phases 1–2.
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
#   MAX_FILES     Max ``.df`` paths per GENIE knob group for chunk-map (0 = all in each group).
#                 Overridden by ``--max-files N`` / ``-n N``.
#   GENIE_RUN_GROUPS  Comma-separated subset of ``GENIE_GROUP_GLOBS`` keys (e.g. ``CCQE,MEC``).
#                     Overridden by ``--genie-groups LIST`` / ``-g LIST``. Empty = all groups
#                     that appear in the active glob map for ``MC_DF_STAGE``.
#   WORKERS           Map-phase parallel worker count for ``syst_genie_parallel.py``
#                     (default: ``min(nproc, 8)``). Overridden by ``--workers N`` / ``-j N``.
#                     Each worker calls ``get_systematics_genie.run_chunk_map`` directly inside
#                     a forked process — no Python re-launch per file — so heavy imports are
#                     paid ``WORKERS`` times instead of once per ``.df``.
#   MERGE_WORKERS     Per-group merge worker count (default: ``min(N_groups, WORKERS)``).
#                     0 = use the default.
#   NUMUCC_SYST_DISK_ROOT     Syst disk root for phase 3 (default: ``dataset_locations.default_syst_disk_root()``).
#
# Logs ``progress map overall`` / ``progress merge`` lines (k/total, BEGIN/END timestamps).
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

_usage() {
    cat <<'EOF'
run_syst_genie_chunked.sh [OPTIONS]

Options:
  --max-files|-n N       Run at most N chunk-map jobs per GENIE knob group (each .df × group). 0 = no limit.
  --genie-groups|-g LIST Comma-separated GENIE groups (keys of ``GENIE_GROUP_GLOBS`` for ``MC_DF_STAGE=final``).
                         Default: all groups that have entries in the active glob map.
  --workers|-j N         Number of parallel worker processes for the chunk-map phase
                         (default: min(nproc, 8)). Imports are shared via fork — primary speed-up.
  --merge-workers M      Number of parallel workers for the per-group chunk-merge phase
                         (default: min(N_groups, --workers); 0 = default).
  -h, --help             Show this message

Environment:
  MAX_FILES            Same as --max-files: per GENIE knob group, not total (CLI wins).
  GENIE_RUN_GROUPS     Same as --genie-groups (CLI wins).
  WORKERS              Same as --workers (CLI wins).
  MERGE_WORKERS        Same as --merge-workers (CLI wins).
  NUMUCC_SYST_DISK_ROOT     Syst disk root for phase 3 (default: default_syst_disk_root()).
  See script header for other variables.
EOF
}

_CLI_MAX_FILES=""
_CLI_GENIE_GROUPS=""
_CLI_WORKERS=""
_CLI_MERGE_WORKERS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-files|-n)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "[genie-run] ERROR: $1 requires a non-negative integer" >&2
                exit 2
            fi
            _CLI_MAX_FILES="$2"
            shift 2
            ;;
        --genie-groups|-g)
            if [[ -z "${2:-}" ]]; then
                echo "[genie-run] ERROR: $1 requires a comma-separated list of GENIE group tags" >&2
                exit 2
            fi
            _CLI_GENIE_GROUPS="$2"
            shift 2
            ;;
        --workers|-j)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]]; then
                echo "[genie-run] ERROR: $1 requires a positive integer" >&2
                exit 2
            fi
            _CLI_WORKERS="$2"
            shift 2
            ;;
        --merge-workers)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "[genie-run] ERROR: $1 requires a non-negative integer" >&2
                exit 2
            fi
            _CLI_MERGE_WORKERS="$2"
            shift 2
            ;;
        -h|--help)
            _usage
            exit 0
            ;;
        *)
            echo "[genie-run] ERROR: unknown option: $1 (use --help)" >&2
            exit 2
            ;;
    esac
done

MAX_FILES="${_CLI_MAX_FILES:-${MAX_FILES:-0}}"
GENIE_RUN_GROUPS="${_CLI_GENIE_GROUPS:-${GENIE_RUN_GROUPS:-}}"
_DEFAULT_WORKERS=$(python3 -c "import os; print(min(os.cpu_count() or 8, 8))")
WORKERS="${_CLI_WORKERS:-${WORKERS:-$_DEFAULT_WORKERS}}"
MERGE_WORKERS="${_CLI_MERGE_WORKERS:-${MERGE_WORKERS:-0}}"
export MAX_FILES GENIE_RUN_GROUPS WORKERS MERGE_WORKERS

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_genie_syst_work_root
print(default_genie_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR=${CHUNKS_DIR:-"$WORK_BASE/chunks"}
MERGE_ROOT=${MERGE_ROOT:-"$WORK_BASE/merged_perTPC"}
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

parallel_py="$THIS_DIR/syst_genie_parallel.py"
FAILED_LOG="$WORK_BASE/failed_genie_df_files.log"
touch "$FAILED_LOG"

echo "[genie-run] WORK_BASE=$WORK_BASE  CHUNKS_DIR=$CHUNKS_DIR  MERGE_ROOT=$MERGE_ROOT"
echo "[genie-run] MC_DF_STAGE=$MC_DF_STAGE  INPUT_STAGE=$INPUT_STAGE  XSEC_UNIT=$XSEC_UNIT"
echo "[genie-run] MAX_FILES=${MAX_FILES} (0 = no per-group cap; else max .df files per GENIE knob group)"
if [[ -n "${GENIE_RUN_GROUPS// /}" ]]; then
    echo "[genie-run] GENIE_RUN_GROUPS=${GENIE_RUN_GROUPS} (subset mode)"
else
    echo "[genie-run] GENIE_RUN_GROUPS=(empty) — all groups in active GENIE glob map"
fi
echo "[genie-run] WORKERS=${WORKERS} parallel chunk-map worker processes (MERGE_WORKERS=${MERGE_WORKERS})"

# ``syst_genie_parallel.py`` imports the GENIE chunk-map code once and forks
# workers. Per-file pickle outputs (``genie__<GROUP>__<stem>.pkl``) match the
# legacy serial path bit-for-bit, so ``chunk-merge`` semantics are unchanged.
# Atomic ``.tmp`` rename inside ``run_chunk_map`` makes resume safe.
parallel_args=(
    --mc-df-stage "$MC_DF_STAGE"
    --chunks-dir "$CHUNKS_DIR"
    --merge-root "$MERGE_ROOT"
    --failed-log "$FAILED_LOG"
    --workers "$WORKERS"
    --merge-workers "$MERGE_WORKERS"
    --max-files "$MAX_FILES"
    --xsec-unit "$XSEC_UNIT"
)
if [[ -n "${GENIE_RUN_GROUPS// /}" ]]; then
    parallel_args+=(--genie-groups "$GENIE_RUN_GROUPS")
fi
if [[ "${SKIP_MERGE:-0}" == "1" ]]; then
    parallel_args+=(--skip-merge)
fi

echo "[genie-run] progress map BEGIN $(date -Is) (workers=${WORKERS})"
python3 "$parallel_py" "${parallel_args[@]}"
echo "[genie-run] DONE map=$CHUNKS_DIR merge=$MERGE_ROOT"

agg_py="$THIS_DIR/syst_genie_aggregate.py"
SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_root
print(default_syst_disk_root())
")}"
mkdir -p "$SYST_DISK_ROOT"
echo "[genie-run] SYST_DISK_ROOT=$SYST_DISK_ROOT  (GENIE/ aggregate target)"
agg_cmd=(
    python3 "$agg_py"
    --chunks-dir "$CHUNKS_DIR"
    --syst-disk-root "$SYST_DISK_ROOT"
    --mc-df-stage "$MC_DF_STAGE"
    --xsec-unit "$XSEC_UNIT"
)
if [[ -n "${GENIE_RUN_GROUPS// /}" ]]; then
    agg_cmd+=(--genie-groups "$GENIE_RUN_GROUPS")
fi
echo "[genie-run] progress syst-disk 1/1  BEGIN $(date -Is)"
"${agg_cmd[@]}"
echo "[genie-run] progress syst-disk 1/1  END $(date -Is)"
echo "[genie-run] $(date -Is) DONE GENIE -> ${SYST_DISK_ROOT}/GENIE/cov_mat_dict.pkl"
