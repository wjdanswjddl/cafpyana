#!/usr/bin/env bash
# Final-stage GENIE syst on Aug24 sel_mup MEC + Other (same recipe as CCQE/Ar23p).
#
# Intended for a *second* machine while Ar23p runs elsewhere. Uses shared
# /exp/sbnd data paths by default; override REPO / WORK_BASE / SYST_DISK if needed.
#
# Usage:
#   nohup bash launch_sel_mup_Aug24_MEC_Other.sh > run_MEC_Other.nohup.out 2>&1 &
#
# Optional env:
#   WORKERS=16          map workers (default 16)
#   MERGE_WORKERS=1     merge workers (default 1)
#   REPO=...            cafpyana checkout
#   WORK_BASE=...       chunk/merge work root
#   NUMUCC_SYST_DISK_ROOT=...   cov_mat_dict output root
#   GENIE_VAR_SAVE_NAMES=...    default: integrated + 3 TKI vars
#   GENIE_MULTISIGMA_DIVIDE_BY_CV=1   MvA-style /cv (harmless for MEC/Other)
#
set -euo pipefail

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

export WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_mup_Aug24_MEC_Other}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_mup_Aug24_MEC_Other}"
export MC_DF_STAGE=final
export GENIE_RUN_GROUPS=MEC,Other
export GENIE_VAR_SAVE_NAMES="${GENIE_VAR_SAVE_NAMES:-integrated,tki-del_Tp,tki-del_alpha,tki-del_phi}"
export GENIE_MULTISIGMA_DIVIDE_BY_CV="${GENIE_MULTISIGMA_DIVIDE_BY_CV:-1}"

WORKERS="${WORKERS:-16}"
MERGE_WORKERS="${MERGE_WORKERS:-1}"

mkdir -p "$WORK_BASE/chunks" "$WORK_BASE/merged" "$NUMUCC_SYST_DISK_ROOT"
cd "$WORK_BASE"
LOG="$WORK_BASE/run_MEC_Other.log"

echo "===== LAUNCH $(date -Is) final MEC,Other Aug24 vars=${GENIE_VAR_SAVE_NAMES} workers=${WORKERS} =====" | tee -a "$LOG"
echo "[info] WORK_BASE=$WORK_BASE" | tee -a "$LOG"
echo "[info] NUMUCC_SYST_DISK_ROOT=$NUMUCC_SYST_DISK_ROOT" | tee -a "$LOG"
echo "[info] globs from dataset_locations.GENIE_GROUP_GLOBS (Aug24 sel_mup MEC + Other)" | tee -a "$LOG"

bash "$REPO/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh" \
  -g MEC,Other -j "$WORKERS" --merge-workers "$MERGE_WORKERS" >>"$LOG" 2>&1

echo "===== DONE $(date -Is) -> $NUMUCC_SYST_DISK_ROOT/GENIE/cov_mat_dict.pkl =====" | tee -a "$LOG"
