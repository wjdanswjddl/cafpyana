#!/usr/bin/env bash
#
# Re-process sel_all GENIE dfs with the current event-selection thresholds
# (``make_pandora_evtdf`` / ``MU_CHI2MU_TH`` etc. from ``makedf.selections``)
# and write covariances + per-knob uncertainties under PRL_corrected.
#
# Does **not** overwrite PRL_genie_sel_all.
#
# Usage:
#   bash analysis_village/numucc_1p0pi/scripts/run_prl_corrected.sh
#   MAX_FILES=2 bash .../run_prl_corrected.sh          # smoke
#
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

export WORK_BASE="${WORK_BASE:-${PRL_DATA}/PRL_corrected}"
export CHUNKS_DIR="${CHUNKS_DIR:-${WORK_BASE}/chunks}"
export MERGE_ROOT="${MERGE_ROOT:-${WORK_BASE}/merged}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-${WORK_BASE}/syst_disk}"
export MC_DF_STAGE=sel_all
export XSEC_UNIT="${XSEC_UNIT:-1.0}"
export MAX_FILES="${MAX_FILES:-0}"
export WORKERS="${WORKERS:-20}"
export MERGE_WORKERS="${MERGE_WORKERS:-4}"
export PHASE="${PHASE:-all}"

mkdir -p "$WORK_BASE" "$CHUNKS_DIR" "$MERGE_ROOT" "$NUMUCC_SYST_DISK_ROOT"
LOG="$WORK_BASE/prl_corrected_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "========================================================================"
echo "[PRL_corrected] BEGIN $(date -Is)"
echo "[PRL_corrected] WORK_BASE=$WORK_BASE"
echo "[PRL_corrected] PHASE=$PHASE  MAX_FILES=$MAX_FILES  WORKERS=$WORKERS"
python3 - <<'PY'
from analysis_village.numucc_1p0pi.makedf.selections import MU_CHI2MU_TH, P_CHI2P_TH, NU_SCORE_TH, VTXDIST_TH
print("[PRL_corrected] selection thresholds: NU_SCORE_TH=%s  VTXDIST_TH=%s  MU_CHI2MU_TH=%s  P_CHI2P_TH=%s"
      % (NU_SCORE_TH, VTXDIST_TH, MU_CHI2MU_TH, P_CHI2P_TH))
PY
echo "========================================================================"

case "$WORK_BASE" in
    */PRL_genie_sel_all|*/PRL_genie_sel_all/*|*/PRL_genie_per_knob_Ar23p|*/PRL_genie_per_knob_Ar23p/*)
        echo "[PRL_corrected] ERROR: refusing to write into existing PRL tree: $WORK_BASE" >&2
        exit 1
        ;;
esac

bash "${PRL_DATA}/run_prl_genie_sel_all.sh"

echo "[PRL_corrected] map/merge/aggregate done; writing per-knob NPZ + combined cov $(date -Is)"
python3 "${REPO_ROOT}/analysis_village/numucc_1p0pi/scripts/prl_genie_sel_all_finalize.py" \
    --work-dir "$WORK_BASE" \
    --groups "${GENIE_FINAL_GROUPS:-slim,CCQE,MEC,RES,nonRES,Other}" \
    --no-wait \
    --force-merge

echo "[PRL_corrected] DONE $(date -Is)"
echo "[PRL_corrected] combined cov → ${WORK_BASE}/syst_disk/GENIE/cov_mat_dict.pkl"
echo "[PRL_corrected] per-knob     → ${WORK_BASE}/per_knob/per_knob_all_vars.json"
