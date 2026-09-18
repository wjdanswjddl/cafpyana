#!/usr/bin/env bash
# Final-stage GENIE syst for the VecFF group (single knob: SBN_v1_multisigma_VecFFCCQEshape).
# Uses GENIE_GROUP_GLOBS["VecFF"] — the sel_mup-wgts_genie_VecFF jobs from
# submit_mc_jobs_GENIE_vecff.sh. Work dirs are parallel to the Ar23p ones, so nothing
# from the earlier GENIE campaigns is touched.
#
# Usage:
#   nohup bash analysis_village/numucc_1p0pi/scripts/launch_sel_mup_VecFF.sh > run_VecFF.nohup.out 2>&1 &
#
# Optional env: WORKERS (default 16), MERGE_WORKERS, REPO, WORK_BASE, NUMUCC_SYST_DISK_ROOT
#
set -euo pipefail

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

export WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_mup_VecFF}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_mup_VecFF}"
export MC_DF_STAGE=final
export GENIE_RUN_GROUPS=VecFF
export GENIE_VAR_SAVE_NAMES="${GENIE_VAR_SAVE_NAMES:-integrated,tki-del_Tp,tki-del_alpha,tki-del_phi}"
export GENIE_MULTISIGMA_DIVIDE_BY_CV="${GENIE_MULTISIGMA_DIVIDE_BY_CV:-1}"

WORKERS="${WORKERS:-16}"
MERGE_WORKERS="${MERGE_WORKERS:-1}"

mkdir -p "$WORK_BASE/chunks" "$WORK_BASE/merged" "$NUMUCC_SYST_DISK_ROOT"
cd "$WORK_BASE"
LOG="$WORK_BASE/run_VecFF.log"

echo "===== LAUNCH $(date -Is) final VecFF vars=${GENIE_VAR_SAVE_NAMES} workers=${WORKERS} =====" | tee -a "$LOG"
echo "[info] WORK_BASE=$WORK_BASE" | tee -a "$LOG"
echo "[info] NUMUCC_SYST_DISK_ROOT=$NUMUCC_SYST_DISK_ROOT" | tee -a "$LOG"
echo "[info] input from dataset_locations.GENIE_GROUP_GLOBS['VecFF']" | tee -a "$LOG"

bash "$REPO/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh" \
  -g VecFF -j "$WORKERS" --merge-workers "$MERGE_WORKERS" >>"$LOG" 2>&1

echo "===== DONE $(date -Is) -> $NUMUCC_SYST_DISK_ROOT/GENIE/cov_mat_dict.pkl =====" | tee -a "$LOG"
