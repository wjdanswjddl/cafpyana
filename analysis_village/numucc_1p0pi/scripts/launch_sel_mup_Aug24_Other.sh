#!/usr/bin/env bash
# Final-stage GENIE syst on Aug24 sel_mup Other:
#   2026_08_24_130336__sel_mup-wgts_genie_Other  (GENIE_GROUP_GLOBS["Other"])
#
# Usage:
#   nohup bash launch_sel_mup_Aug24_Other.sh > run_Other.nohup.out 2>&1 &
#
# Optional env: WORKERS (default 16), MERGE_WORKERS, REPO, WORK_BASE, NUMUCC_SYST_DISK_ROOT
#
set -euo pipefail

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

export WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_mup_Aug24_Other}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_mup_Aug24_Other}"
export MC_DF_STAGE=final
export GENIE_RUN_GROUPS=Other
export GENIE_VAR_SAVE_NAMES="${GENIE_VAR_SAVE_NAMES:-integrated,tki-del_Tp,tki-del_alpha,tki-del_phi}"
export GENIE_MULTISIGMA_DIVIDE_BY_CV="${GENIE_MULTISIGMA_DIVIDE_BY_CV:-1}"

WORKERS="${WORKERS:-16}"
MERGE_WORKERS="${MERGE_WORKERS:-1}"

mkdir -p "$WORK_BASE/chunks" "$WORK_BASE/merged" "$NUMUCC_SYST_DISK_ROOT"
cd "$WORK_BASE"
LOG="$WORK_BASE/run_Other.log"

echo "===== LAUNCH $(date -Is) final Other Aug24 vars=${GENIE_VAR_SAVE_NAMES} workers=${WORKERS} =====" | tee -a "$LOG"
echo "[info] WORK_BASE=$WORK_BASE" | tee -a "$LOG"
echo "[info] NUMUCC_SYST_DISK_ROOT=$NUMUCC_SYST_DISK_ROOT" | tee -a "$LOG"
echo "[info] input glob: .../2026_08_24_130336__sel_mup-wgts_genie_Other/*.df" | tee -a "$LOG"

bash "$REPO/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh" \
  -g Other -j "$WORKERS" --merge-workers "$MERGE_WORKERS" >>"$LOG" 2>&1

echo "===== DONE $(date -Is) -> $NUMUCC_SYST_DISK_ROOT/GENIE/cov_mat_dict.pkl =====" | tee -a "$LOG"
