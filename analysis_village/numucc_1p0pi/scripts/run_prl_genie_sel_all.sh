#!/usr/bin/env bash
#
# GENIE systematics on sel_all dfs (event-selection cut variables + final vars).
# -----------------------------------------------------------------------------
# Writes under WORK_BASE = PRL_data/PRL_genie_sel_all (sibling of the sel_mup
# PRL_genie_per_knob_Ar23p tree) so existing sel_mup chunks/merged/per_knob are
# never overwritten.
#
# Process order: slim first, then per-mode knob groups.
#
# Usage:
#   bash /exp/sbnd/data/users/munjung/PRL_data/run_prl_genie_sel_all.sh
#   MAX_FILES=2 bash .../run_prl_genie_sel_all.sh          # smoke
#   PHASE=slim bash .../run_prl_genie_sel_all.sh           # slim only
#   PHASE=modes bash .../run_prl_genie_sel_all.sh          # modes only
#   GENIE_RUN_GROUPS=CCQE,MEC bash .../run_prl_genie_sel_all.sh
# -----------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana"
PRL_DATA="/exp/sbnd/data/users/munjung/PRL_data"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

VENV="${REPO_ROOT}/envs/venv_py310_cafpyana/bin/activate"
if [[ -f "$VENV" ]]; then
    # shellcheck disable=SC1090
    source "$VENV"
fi

DRIVER="${REPO_ROOT}/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh"
if [[ ! -f "$DRIVER" ]]; then
    echo "[sel_all-GENIE] ERROR: missing $DRIVER" >&2
    exit 1
fi

WORK_BASE="${WORK_BASE:-${PRL_DATA}/PRL_genie_sel_all}"
export WORK_BASE
export CHUNKS_DIR="${CHUNKS_DIR:-${WORK_BASE}/chunks}"
export MERGE_ROOT="${MERGE_ROOT:-${WORK_BASE}/merged}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-${WORK_BASE}/syst_disk}"
export MC_DF_STAGE=sel_all
export XSEC_UNIT="${XSEC_UNIT:-1.0}"
export MAX_FILES="${MAX_FILES:-0}"
export WORKERS="${WORKERS:-16}"
export MERGE_WORKERS="${MERGE_WORKERS:-4}"

PHASE="${PHASE:-all}"   # all | slim | modes
MODE_GROUPS="${GENIE_RUN_GROUPS:-CCQE,MEC,RES,nonRES,Other}"

mkdir -p "$WORK_BASE" "$CHUNKS_DIR" "$MERGE_ROOT" "$NUMUCC_SYST_DISK_ROOT"
LOG="$WORK_BASE/sel_all_genie_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "========================================================================"
echo "[sel_all-GENIE] BEGIN $(date -Is)"
echo "[sel_all-GENIE] WORK_BASE=$WORK_BASE"
echo "[sel_all-GENIE] CHUNKS_DIR=$CHUNKS_DIR  MERGE_ROOT=$MERGE_ROOT"
echo "[sel_all-GENIE] SYST_DISK=$NUMUCC_SYST_DISK_ROOT"
echo "[sel_all-GENIE] PHASE=$PHASE  MAX_FILES=$MAX_FILES  WORKERS=$WORKERS"
echo "========================================================================"

# Refuse to write into the sel_mup PRL tree.
case "$WORK_BASE" in
    */PRL_genie_per_knob_Ar23p|*/PRL_genie_per_knob_Ar23p/*)
        echo "[sel_all-GENIE] ERROR: refusing WORK_BASE under sel_mup tree: $WORK_BASE" >&2
        exit 1
        ;;
esac

run_groups() {
    local groups="$1"
    echo "[sel_all-GENIE] === GENIE_RUN_GROUPS=$groups  $(date -Is) ==="
    GENIE_RUN_GROUPS="$groups" bash "$DRIVER"
}

case "$PHASE" in
    slim)
        run_groups "slim"
        ;;
    modes)
        run_groups "$MODE_GROUPS"
        ;;
    all)
        run_groups "slim"
        run_groups "$MODE_GROUPS"
        ;;
    *)
        echo "[sel_all-GENIE] ERROR: PHASE must be all|slim|modes, got $PHASE" >&2
        exit 2
        ;;
esac

echo "[sel_all-GENIE] DONE $(date -Is)"
echo "[sel_all-GENIE] merged → $MERGE_ROOT"
echo "[sel_all-GENIE] syst_disk → $NUMUCC_SYST_DISK_ROOT/GENIE"
