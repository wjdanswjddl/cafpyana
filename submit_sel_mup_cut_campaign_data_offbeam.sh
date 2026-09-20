#!/usr/bin/env bash
# Resubmit only data + offbeam for one cut tag (higher ngrid to avoid ARG_MAX).
# Usage: CUT_TAG=nu_score0 bash submit_sel_mup_cut_campaign_data_offbeam.sh
set -euo pipefail

CUT_TAG="${CUT_TAG:?set CUT_TAG}"

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export CAFPYANA_DIR="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"
export CAFPYANA_GRID_SUBDIR="$CUT_TAG"

# Higher than the failed 40/30 — keep files/job closer to historical working counts.
NGRID_OFFBEAM="${NGRID_OFFBEAM:-80}"
NGRID_DATA="${NGRID_DATA:-150}"

export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-8h}"
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-6GB}"
export JOBSUB_DISK="${JOBSUB_DISK:-20GB}"
export JOBSUB_CPU="${JOBSUB_CPU:-7}"

OFFBEAM_LIST=/exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list
DATA_LIST=/exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list

echo "[data/offbeam resubmit] CUT_TAG=$CUT_TAG ngrid offbeam=$NGRID_OFFBEAM data=$NGRID_DATA"

# shellcheck disable=SC1090
source ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l "$OFFBEAM_LIST" \
  -o sel_mup-data-OffBeamLight -ngrid "$NGRID_OFFBEAM" -gsubdir "$CUT_TAG"

source ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l "$DATA_LIST" \
  -o sel_mup-data-1e20 -ngrid "$NGRID_DATA" -gsubdir "$CUT_TAG"

echo "[data/offbeam resubmit] done CUT_TAG=$CUT_TAG"
