#!/usr/bin/env bash
# FSI compare df campaigns (Product A + B): BASE ∪ retired v1 _N ∪ v3 nucleon FSI,
# with GENIE_slim_v1 / GENIE_slim_v3 / GENIE_slim_both (morph/multisigma throws, recipe B).
#
# Must use the *_xrootd.list (root://…), NOT flatcaf_sbnd.list (/pnfs/…).
# Requires jobsub_submit on PATH and a valid analysis token (e.g. ``. ~/get_token.sh``).
set -euo pipefail

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_DIR="$REPO"
export CAFPYANA_WD="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"
mkdir -p "$CAFPYANA_GRID_OUT_DIR"

ar23_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
ar23_list="${ar23_inputdir}/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list"

NGRID="${NGRID:-2000}"
export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-6h}"
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-10GB}"
export JOBSUB_DISK="${JOBSUB_DISK:-20GB}"

echo "[info] CAFPYANA_WD=$CAFPYANA_WD"
echo "[info] CAFPYANA_GRID_OUT_DIR=$CAFPYANA_GRID_OUT_DIR"
echo "[info] python=$(which python) lifetime=$JOBSUB_LIFETIME mem=$JOBSUB_MEMORY"

echo "===== FSI_compare sel_mup (Product B) ngrid=${NGRID} $(date -Is) ====="
python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_mup-geniewgts-fsi_compare.py \
  -l "$ar23_list" \
  -o "sel_mup-wgts_genie_FSI_compare" \
  -ngrid "$NGRID"

echo "===== FSI_compare sel_all (Product A) ngrid=${NGRID} $(date -Is) ====="
python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_all-geniewgts-fsi_compare.py \
  -l "$ar23_list" \
  -o "sel_all-wgts_genie_FSI_compare" \
  -ngrid "$NGRID"

echo "===== SUBMIT DONE $(date -Is) ====="
