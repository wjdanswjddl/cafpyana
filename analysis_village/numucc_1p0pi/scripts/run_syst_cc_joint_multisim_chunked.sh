#!/usr/bin/env bash
#
# Joint-bin (cross-variable) multisim map → aggregate into ``syst_disk_CC``.
# -----------------------------------------------------------------------------
# Phase 1: same (syst, ``.df``) queue as ``run_syst_multisim_chunked.sh``, but each job
# runs ``syst_cc_joint_multisim_chunk.py`` → stacked ``(n_X+n_Y)`` universes per kinematic
# pair (X bins first, then Y). Pickles: ``nu__joint_cc__<tag>__<stem>.pkl`` (``joint_cc`` prefix
# so parallel ``skip_existing`` does not reuse legacy ``nu__joint__*`` shards after pair-list updates).
# Phase 2: ``syst_cc_joint_multisim_aggregate.py`` → ``JointMCstat/``, ``JointFlux/``, ``JointG4/``
# (one ``joint_*_combined.npz`` per category; ``cc_joint_cov`` sums them for the total multisim term).
#
# Pair coverage: by default the chunk script iterates **every** preset pair returned by
# ``syst_cc_joint_multisim_common.default_kinematic_joint_pairs`` — including the same-side
# pairs (``muon_p__muon_costheta``, ``proton_p__proton_costheta``) needed for the
# **multi-variable Y** conditional constraint (``Y_i × Y_j`` and ``X_i × X_j`` cross blocks of
# the joint covariance ``Σ``). Restrict via ``JOINT_PAIRS`` CSV when re-running a subset.
#
# Environment (mirrors marginal multisim where applicable):
#   JOINT_CC_WORK_BASE   Chunk work root (default: ``default_joint_multisim_cc_work_root``).
#                        Or set ``NUMUCC_JOINT_MULTISIM_CC_WORK_BASE``.
#   JOINT_CC_CHUNKS_SUBDIR  Optional path segment under each chunk root's ``chunks/`` directory
#                           (default ``CC_joint``), e.g. ``.../g4_syst-.../chunks/CC_joint/`` so CC
#                           joint shards do not sit next to unrelated pickles. Set
#                           ``JOINT_CC_CHUNKS_SUBDIR=`` (empty, exported before the script) to
#                           use the legacy flat ``chunks/`` layout.
#   MCSTAT_WORK_BASE, G4_WORK_BASE, FLUX_WORK_BASE — same defaults as marginal script.
#   MULTISIM_SYST_TYPES  all (Flux,G4) | full | comma subset (MCstat,Flux,G4).
#   SYST_DISK_CC_ROOT    Output tree for CC joint NPZ (default: ``default_syst_disk_cc_root``).
#   MC_DF_STAGE          Must be ``final`` (joint chunk does not implement sel_all).
#   SKIP_AGGREGATE, MAX_FILES, WORKERS, G4_MODE, FLUX_MODE, FLUX_KNOB_GROUPS, NO_PLOTS (unused)
#   JOINT_PAIRS          Optional CSV passed as ``--pairs`` (pair slugs). Examples:
#                        ``muon_p__proton_p`` (single cross pair),
#                        ``muon_p__muon_costheta`` (Y-Y same-side pair for multi-Y),
#                        ``muon_p__proton_p,muon_p__muon_costheta,muon_costheta__proton_p``
#                        (all pairs needed for X=proton_p constrained by both muon vars).
#   MAX_SPLITS           HDF splits cap per ``.df`` (0 = all); passed as ``--max-splits``.
#   N_UNIVERSE           Universe count for rate histograms (default 100); ``--n-universe``.
#   NUMUCC_JOINT_MULTISIM_SHAPES  1 (default) | all | 0 — log stacked u_j/cv_j shapes in chunk.py
#
# Exports ``NUMUCC_SYST_DISK_CC_ROOT`` to the resolved aggregate target so downstream
# loaders pick up the same tree without extra env setup.
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

_CLI_SYST_TYPES=""
_CLI_MAX_FILES=""
_CLI_WORKERS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --syst-types)
            [[ -z "${2:-}" ]] && { echo "[cc-joint-multisim-run] ERROR: --syst-types requires a value" >&2; exit 2; }
            _CLI_SYST_TYPES="$2"
            shift 2
            ;;
        --max-files|-n)
            [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] && { echo "[cc-joint-multisim-run] ERROR: $1 requires a non-negative integer" >&2; exit 2; }
            _CLI_MAX_FILES="$2"
            shift 2
            ;;
        --workers|-j)
            [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]] && { echo "[cc-joint-multisim-run] ERROR: $1 requires a positive integer" >&2; exit 2; }
            _CLI_WORKERS="$2"
            shift 2
            ;;
        -h|--help)
            head -n 35 "$0" | tail -n +2
            exit 0
            ;;
        *)
            echo "[cc-joint-multisim-run] ERROR: unknown option: $1 (use --help)" >&2
            exit 2
            ;;
    esac
done

MAX_FILES="${_CLI_MAX_FILES:-${MAX_FILES:-0}}"
_DEFAULT_WORKERS=$(python3 -c "import os; print(min(os.cpu_count() or 8, 8))")
WORKERS="${_CLI_WORKERS:-${WORKERS:-$_DEFAULT_WORKERS}}"

MULTISIM_SYST_TYPES="${_CLI_SYST_TYPES:-${MULTISIM_SYST_TYPES:-all}}"
MULTISIM_SYST_TYPES="$(echo "${MULTISIM_SYST_TYPES}" | tr '[:upper:]' '[:lower:]')"
_NEUTRINO_ORDER=('MCstat' 'Flux' 'G4')

if [[ "$MULTISIM_SYST_TYPES" == "all" || -z "$MULTISIM_SYST_TYPES" ]]; then
    ONLY_SYSTS=('Flux' 'G4')
elif [[ "${MULTISIM_SYST_TYPES}" == "full" ]]; then
    ONLY_SYSTS=("${_NEUTRINO_ORDER[@]}")
else
    IFS=',' read -ra _raw_systs <<< "${MULTISIM_SYST_TYPES}"
    _raw_norm=()
    for _s in "${_raw_systs[@]}"; do
        _sn="$(echo "${_s}" | tr '[:lower:]' '[:upper:]')"
        case "${_sn}" in
            MCSTAT|MC_STAT|MC) _raw_norm+=('MCstat') ;;
            FLUX) _raw_norm+=('Flux') ;;
            G4) _raw_norm+=('G4') ;;
            *) echo "[cc-joint-multisim-run] ERROR: unknown syst type '${_s}'" >&2; exit 2 ;;
        esac
    done
    _uniq=()
    for _sn in "${_raw_norm[@]}"; do
        if [[ " ${_uniq[*]} " != *" ${_sn} "* ]]; then _uniq+=("${_sn}"); fi
    done
    ONLY_SYSTS=()
    for _sn in "${_NEUTRINO_ORDER[@]}"; do
        if [[ " ${_uniq[*]} " == *" ${_sn} "* ]]; then ONLY_SYSTS+=("${_sn}"); fi
    done
fi

ONLY_SYSTS_CSV="$(IFS=','; echo "${ONLY_SYSTS[*]}")"

TODAY=$(date +%Y%m%d)
JOINT_CC_WORK_BASE="${JOINT_CC_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_joint_multisim_cc_work_root
print(default_joint_multisim_cc_work_root('${TODAY}'))
")}"

G4_WORK_BASE=${G4_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_g4_syst_work_root
print(default_g4_syst_work_root('${TODAY}'))
")}
FLUX_WORK_BASE=${FLUX_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_flux_syst_work_root
print(default_flux_syst_work_root('${TODAY}'))
")}
MCSTAT_WORK_BASE=${MCSTAT_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_mcstat_syst_work_root
print(default_mcstat_syst_work_root('${TODAY}'))
")}

JOINT_CC_CHUNKS_SUBDIR="${JOINT_CC_CHUNKS_SUBDIR-CC_joint}"
_cc_ms_chunks_suffix=""
if [[ -n "${JOINT_CC_CHUNKS_SUBDIR}" ]]; then
    _cc_ms_chunks_suffix="/${JOINT_CC_CHUNKS_SUBDIR}"
fi

CHUNKS_MULTISIM="$JOINT_CC_WORK_BASE/chunks${_cc_ms_chunks_suffix}"
CHUNKS_MCSTAT="$MCSTAT_WORK_BASE/chunks${_cc_ms_chunks_suffix}"
CHUNKS_G4="$G4_WORK_BASE/chunks${_cc_ms_chunks_suffix}"
CHUNKS_FLUX="$FLUX_WORK_BASE/chunks${_cc_ms_chunks_suffix}"

SYST_DISK_CC_ROOT="$(python3 -c "
import sys
from pathlib import Path
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_cc_root
print(Path(default_syst_disk_cc_root()).resolve())
")"
export NUMUCC_SYST_DISK_CC_ROOT="${NUMUCC_SYST_DISK_CC_ROOT:-$SYST_DISK_CC_ROOT}"

FAILED_LOG="$JOINT_CC_WORK_BASE/failed_cc_joint_multisim_df_files.log"
MC_DF_STAGE="${MC_DF_STAGE:-final}"
MAX_SPLITS="${MAX_SPLITS:-0}"
N_UNIVERSE="${N_UNIVERSE:-100}"
if [[ "$MC_DF_STAGE" != "final" ]]; then
    echo "[cc-joint-multisim-run] ERROR: MC_DF_STAGE must be 'final' for joint multisim (got '$MC_DF_STAGE')" >&2
    exit 2
fi

G4_MODE="${G4_MODE:-knobs}"
FLUX_MODE="${FLUX_MODE:-knobs}"
FLUX_KNOB_GROUPS="${FLUX_KNOB_GROUPS:-all}"

parallel_py="$THIS_DIR/syst_cc_joint_multisim_parallel.py"
agg_py="$THIS_DIR/syst_cc_joint_multisim_aggregate.py"

mkdir -p "${CHUNKS_MULTISIM}/Combined" "$CHUNKS_MCSTAT" "$CHUNKS_G4" "$CHUNKS_FLUX" "$SYST_DISK_CC_ROOT"
touch "$FAILED_LOG"

_optional_pairs=()
if [[ -n "${JOINT_PAIRS:-}" ]]; then
    _optional_pairs=(--pairs "${JOINT_PAIRS}")
fi

echo "[cc-joint-multisim-run] JOINT_CC_WORK_BASE=$JOINT_CC_WORK_BASE"
echo "[cc-joint-multisim-run] CHUNKS_MULTISIM=$CHUNKS_MULTISIM  JOINT_CC_CHUNKS_SUBDIR=${JOINT_CC_CHUNKS_SUBDIR:-"(empty)"}  SYST_DISK_CC_ROOT=$SYST_DISK_CC_ROOT"
echo "[cc-joint-multisim-run] ONLY_SYSTS_CSV=$ONLY_SYSTS_CSV  MC_DF_STAGE=$MC_DF_STAGE"
echo "[cc-joint-multisim-run] MAX_FILES=${MAX_FILES} WORKERS=${WORKERS} MAX_SPLITS=${MAX_SPLITS} N_UNIVERSE=${N_UNIVERSE}"

parallel_args=(
    --mc-df-stage "$MC_DF_STAGE"
    --syst-types "$MULTISIM_SYST_TYPES"
    --max-files "$MAX_FILES"
    --workers "$WORKERS"
    --max-splits "$MAX_SPLITS"
    --n-universe "$N_UNIVERSE"
    --chunks-multisim "$CHUNKS_MULTISIM"
    --chunks-mcstat "$CHUNKS_MCSTAT"
    --chunks-flux "$CHUNKS_FLUX"
    --chunks-g4 "$CHUNKS_G4"
    --failed-log "$FAILED_LOG"
    --g4-mode "$G4_MODE"
    --flux-mode "$FLUX_MODE"
    --flux-knob-groups "$FLUX_KNOB_GROUPS"
)
if ((${#_optional_pairs[@]})); then
    parallel_args+=("${_optional_pairs[@]}")
fi

echo "[cc-joint-multisim-run] progress map BEGIN $(date -Is)"
set +e
python3 "$parallel_py" "${parallel_args[@]}"
_ms_map_rc=$?
set -e
echo "[cc-joint-multisim-run] progress map END $(date -Is) rc=${_ms_map_rc}"
if [[ "$_ms_map_rc" -ne 0 ]]; then
    echo "[cc-joint-multisim-run] map phase failed entirely (rc=${_ms_map_rc}); aborting before aggregate" >&2
    exit "$_ms_map_rc"
fi

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[cc-joint-multisim-run] SKIP_AGGREGATE=1 — chunks under $CHUNKS_MULTISIM"
    exit 0
fi

_agg_chunks=(
    "$CHUNKS_MULTISIM"
    "$CHUNKS_MCSTAT"
    "$CHUNKS_G4"
    "$CHUNKS_FLUX"
)
[[ -d "${CHUNKS_MULTISIM}/G4" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/G4")
[[ -d "${CHUNKS_MULTISIM}/Flux" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/Flux")
[[ -d "${CHUNKS_MULTISIM}/MCstat" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/MCstat")

agg_cmd=(
    python3 "$agg_py"
    --chunks_dir "${_agg_chunks[@]}"
    --syst-disk-cc-root "$SYST_DISK_CC_ROOT"
    --syst-types "$ONLY_SYSTS_CSV"
)

echo "[cc-joint-multisim-run] aggregate BEGIN $(date -Is) ${agg_cmd[*]}"
"${agg_cmd[@]}"
echo "[cc-joint-multisim-run] aggregate END $(date -Is)"
echo "[cc-joint-multisim-run] DONE → $SYST_DISK_CC_ROOT (export NUMUCC_SYST_DISK_CC_ROOT=$SYST_DISK_CC_ROOT)"
