#!/usr/bin/env bash
#
# Joint-bin GENIE **rate** map → aggregate into PRL Product B ``JointCC/JointGenie/``.
# -----------------------------------------------------------------------------
# Phase 1: same (GENIE group, ``.df``) queue as ``run_syst_genie_chunked.sh``, but each job
# runs ``syst_cc_joint_genie_chunk.py`` → ``nu__joint_cc_genie_stack__<GROUP>__<stem>.pkl``
# with one inclusive stacked universe vector (muon p, muon cosθ, proton p, proton cosθ;
# ``bkgd_subtract=False``). ``JOINT_PAIRS`` switches to pairwise ``nu__joint_cc_genie__*`` shards.
#
# Phase 2: ``syst_cc_joint_genie_aggregate.py`` merges chunks and writes
# ``JointGenie/joint_genie_combined.npz`` with ``JointGenie`` (sum over groups) plus
# ``JointGenie_by_knob`` (per reweight knob matrices; used by :mod:`cc_joint_cov` with joint multisim).
#
# This path uses GENIE **rate** reweights (``get_univ_rates(..., cov_type=rate)``), not the
# response-matrix **xsec** tensors in marginal ``cov_mat_dict.pkl`` ``["genie"]``.
#
# Environment (mirrors marginal GENIE driver where applicable):
#   JOINT_GENIE_WORK_BASE   Chunk root (default: ``default_joint_genie_cc_work_root``);
#                           or ``NUMUCC_JOINT_GENIE_CC_WORK_BASE``.
#   MC_DF_STAGE             ``final`` only for this pipeline.
#   GENIE_RUN_GROUPS        Default ``FSI_compare`` (``GENIE_slim_v3`` only). Do not add
#                           VecFF/Ar23p/CCQE/MEC — those knobs are already inside slim_v3.
#                           After aggregate, ``align_joint_cc_genie.py`` sets combined
#                           JointGenie to that slim_v3 joint (full off-diagonals).
#                           Pass ``--mec-splice`` on the align script to add May−Sep interpolators.
#   JOINT_CC_MODE           ``stack`` (default, inclusive 4-var vector) or ``pairs``.
#   JOINT_PAIRS             If set, switches to pairwise legacy mode.
#   MAX_FILES, WORKERS      Per-group cap and parallel pool size.
#   SYST_DISK_CC_ROOT       Output tree (default: PRL Product B ``…/JointCC``).
#   JOINT_CC_CHUNKS_SUBDIR  Optional segment under ``.../chunks/`` (default ``CC_joint``). Set to
#                           an empty string before invoking the script to write directly under
#                           ``chunks/`` (pickle prefix ``nu__joint_cc_genie__`` still busts
#                           ``skip_existing`` vs legacy ``nu__joint_genie__*``).
#   CHUNKS_DIR              Override full chunk output directory (default:
#                           ``$JOINT_GENIE_WORK_BASE/chunks/$JOINT_CC_CHUNKS_SUBDIR``).
#   JOINT_PAIRS             Optional CSV for preset pair slugs (passed as ``--pairs``). For a
#                           multi-Y constraint of proton_p by (muon_p, muon_costheta) you need:
#                           ``muon_p__proton_p,muon_costheta__proton_p,muon_p__muon_costheta``.
#   NUMUCC_JOINT_GENIE_SHAPES  1 (default) | all | 0 — log stacked u_j/cv_j shapes from chunk.py
#
# Exports ``NUMUCC_SYST_DISK_CC_ROOT`` when unset so downstream matches aggregate target.
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
VENV="${REPO_ROOT}/envs/venv_py310_cafpyana/bin/activate"
if [[ -f "$VENV" ]]; then
    # shellcheck disable=SC1090
    source "$VENV"
fi
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg

_CLI_MAX_FILES=""
_CLI_GENIE_GROUPS=""
_CLI_WORKERS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-files|-n)
            [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] && { echo "[cc-joint-genie-run] ERROR: $1 requires a non-negative integer" >&2; exit 2; }
            _CLI_MAX_FILES="$2"
            shift 2
            ;;
        --genie-groups|-g)
            [[ -z "${2:-}" ]] && { echo "[cc-joint-genie-run] ERROR: $1 requires a value" >&2; exit 2; }
            _CLI_GENIE_GROUPS="$2"
            shift 2
            ;;
        --workers|-j)
            [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]] && { echo "[cc-joint-genie-run] ERROR: $1 requires a positive integer" >&2; exit 2; }
            _CLI_WORKERS="$2"
            shift 2
            ;;
        -h|--help)
            head -n 28 "$0" | tail -n +2
            exit 0
            ;;
        *)
            echo "[cc-joint-genie-run] ERROR: unknown option: $1 (use --help)" >&2
            exit 2
            ;;
    esac
done

MAX_FILES="${_CLI_MAX_FILES:-${MAX_FILES:-0}}"
_DEFAULT_WORKERS=$(python3 -c "import os; print(min(os.cpu_count() or 8, 8))")
WORKERS="${_CLI_WORKERS:-${WORKERS:-$_DEFAULT_WORKERS}}"
GENIE_RUN_GROUPS="${_CLI_GENIE_GROUPS:-${GENIE_RUN_GROUPS:-FSI_compare}}"

TODAY=$(date +%Y%m%d)
JOINT_GENIE_WORK_BASE="${JOINT_GENIE_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_joint_genie_cc_work_root
print(default_joint_genie_cc_work_root('${TODAY}'))
")}"

JOINT_CC_CHUNKS_SUBDIR="${JOINT_CC_CHUNKS_SUBDIR-CC_joint}"
_cc_g_chunks_suffix=""
if [[ -n "${JOINT_CC_CHUNKS_SUBDIR}" ]]; then
    _cc_g_chunks_suffix="/${JOINT_CC_CHUNKS_SUBDIR}"
fi
CHUNKS_DIR="${CHUNKS_DIR:-$JOINT_GENIE_WORK_BASE/chunks${_cc_g_chunks_suffix}}"

if [[ -z "${SYST_DISK_CC_ROOT:-}" ]]; then
    SYST_DISK_CC_ROOT="$(python3 -c "
import sys
from pathlib import Path
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import prl_syst_disk_root
print((prl_syst_disk_root('B') / 'JointCC').resolve())
")"
fi
export NUMUCC_SYST_DISK_CC_ROOT="$SYST_DISK_CC_ROOT"

MC_DF_STAGE="${MC_DF_STAGE:-final}"
if [[ "$MC_DF_STAGE" != "final" ]]; then
    echo "[cc-joint-genie-run] ERROR: MC_DF_STAGE must be 'final' for joint GENIE (got '$MC_DF_STAGE')" >&2
    exit 2
fi

FAILED_LOG="$JOINT_GENIE_WORK_BASE/failed_cc_joint_genie_df_files.log"
mkdir -p "$CHUNKS_DIR"
touch "$FAILED_LOG"

parallel_py="$THIS_DIR/syst_cc_joint_genie_parallel.py"
agg_py="$THIS_DIR/syst_cc_joint_genie_aggregate.py"

_optional_pairs=()
_optional_mode=(--mode "${JOINT_CC_MODE:-stack}")
if [[ -n "${JOINT_PAIRS:-}" ]]; then
    _optional_pairs=(--pairs "${JOINT_PAIRS}")
    _optional_mode=(--mode pairs)
fi

_optional_groups=()
if [[ -n "${GENIE_RUN_GROUPS// /}" ]]; then
    _optional_groups=(--genie-groups "$GENIE_RUN_GROUPS")
fi

echo "[cc-joint-genie-run] JOINT_GENIE_WORK_BASE=$JOINT_GENIE_WORK_BASE"
echo "[cc-joint-genie-run] CHUNKS_DIR=$CHUNKS_DIR  JOINT_CC_CHUNKS_SUBDIR=${JOINT_CC_CHUNKS_SUBDIR:-"(empty)"}  SYST_DISK_CC_ROOT=$SYST_DISK_CC_ROOT"
echo "[cc-joint-genie-run] MAX_FILES=${MAX_FILES} WORKERS=${WORKERS}"

MAX_SPLITS="${MAX_SPLITS:-0}"

parallel_args=(
    --mc-df-stage "$MC_DF_STAGE"
    --chunks-dir "$CHUNKS_DIR"
    --failed-log "$FAILED_LOG"
    --max-files "$MAX_FILES"
    --max-splits "$MAX_SPLITS"
    --workers "$WORKERS"
)
if ((${#_optional_groups[@]})); then
    parallel_args+=("${_optional_groups[@]}")
fi
if ((${#_optional_pairs[@]})); then
    parallel_args+=("${_optional_pairs[@]}")
fi
if ((${#_optional_mode[@]})); then
    parallel_args+=("${_optional_mode[@]}")
fi

echo "[cc-joint-genie-run] progress map BEGIN $(date -Is)"
set +e
python3 "$parallel_py" "${parallel_args[@]}"
_rc=$?
set -e
echo "[cc-joint-genie-run] progress map END $(date -Is) rc=${_rc}"
if [[ "$_rc" -ne 0 ]]; then
    echo "[cc-joint-genie-run] map phase failed entirely (rc=${_rc}); aborting before aggregate" >&2
    exit "$_rc"
fi

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[cc-joint-genie-run] SKIP_AGGREGATE=1 — chunks: $CHUNKS_DIR"
    exit 0
fi

agg_cmd=(
    python3 "$agg_py"
    --chunks_dir "$CHUNKS_DIR"
    --syst-disk-cc-root "$SYST_DISK_CC_ROOT"
    --mc-df-stage "$MC_DF_STAGE"
)
echo "[cc-joint-genie-run] aggregate BEGIN $(date -Is) ${agg_cmd[*]}"
"${agg_cmd[@]}"
echo "[cc-joint-genie-run] aggregate END $(date -Is)"
if [[ "${SKIP_ALIGN_GENIE:-0}" != "1" ]]; then
    echo "[cc-joint-genie-run] align BEGIN $(date -Is)"
    python3 "$THIS_DIR/align_joint_cc_genie.py" --cc-root "$SYST_DISK_CC_ROOT"
    echo "[cc-joint-genie-run] align END $(date -Is)"
fi
echo "[cc-joint-genie-run] DONE → $SYST_DISK_CC_ROOT (JointGenie/joint_genie_combined.npz)"
