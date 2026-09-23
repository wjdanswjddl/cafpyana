#!/usr/bin/env bash
# Submit sel_mup dfs for all overlay samples under one cut-campaign tag.
#
# Usage:
#   CUT_TAG=nu_score0 bash submit_sel_mup_cut_campaign.sh
#
# Outputs:
#   $CAFPYANA_GRID_OUT_DIR/dfs/$CUT_TAG/<timestamp>__sel_mup-<sample>/
#
# ngrid is deliberately lower than the historical 2000/500/... shard counts so
# each job does more work and queue overhead stays manageable across many cuts.
set -euo pipefail

CUT_TAG="${CUT_TAG:?set CUT_TAG (e.g. nu_score0, chi2mu15, vz_exclude_200_300, fv_z_lt_200)}"

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export CAFPYANA_DIR="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"
export CAFPYANA_GRID_SUBDIR="$CUT_TAG"

# Fewer jobs than prior sel_mup campaigns (~4x fewer for BNB MC).
NGRID_MC="${NGRID_MC:-500}"
NGRID_DIRT="${NGRID_DIRT:-150}"
NGRID_INTIME="${NGRID_INTIME:-200}"
NGRID_OFFBEAM="${NGRID_OFFBEAM:-30}"
NGRID_DATA="${NGRID_DATA:-40}"

export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-8h}"
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-6GB}"
export JOBSUB_DISK="${JOBSUB_DISK:-20GB}"
export JOBSUB_CPU="${JOBSUB_CPU:-7}"

MC_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
DIRT_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/lowE/mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
INTIME_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/intime/mc_MCP2025B_1e20__v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
OFFBEAM_LIST=/exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list
DATA_LIST=/exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list

echo "[cut campaign] CUT_TAG=$CUT_TAG"
echo "[cut campaign] CAFPYANA_GRID_SUBDIR=$CAFPYANA_GRID_SUBDIR"
echo "[cut campaign] ngrid mc=$NGRID_MC dirt=$NGRID_DIRT intime=$NGRID_INTIME offbeam=$NGRID_OFFBEAM data=$NGRID_DATA"
echo "[cut campaign] lifetime=$JOBSUB_LIFETIME mem=$JOBSUB_MEMORY"

# shellcheck disable=SC1090
. ~/get_token.sh

submit_one() {
  local cfg="$1" list="$2" out="$3" ngrid="$4"
  echo "---- submitting $out (ngrid=$ngrid) ----"
  python run_df_maker.py -c "$cfg" -l "$list" -o "$out" -ngrid "$ngrid" -gsubdir "$CUT_TAG"
}

submit_one configs/numucc_1p0pi/sel_mup.py      "$MC_LIST"      sel_mup-mc-BNB_cosmics   "$NGRID_MC"
. ~/get_token.sh
submit_one configs/numucc_1p0pi/sel_mup.py      "$DIRT_LIST"    sel_mup-mc-dirt          "$NGRID_DIRT"
. ~/get_token.sh
submit_one configs/numucc_1p0pi/sel_mup.py      "$INTIME_LIST"  sel_mup-mc-Intime        "$NGRID_INTIME"
. ~/get_token.sh
submit_one configs/numucc_1p0pi/sel_mup-data.py "$OFFBEAM_LIST" sel_mup-data-OffBeamLight "$NGRID_OFFBEAM"
. ~/get_token.sh
submit_one configs/numucc_1p0pi/sel_mup-data.py "$DATA_LIST"    sel_mup-data-1e20        "$NGRID_DATA"

echo "[cut campaign] done submitting CUT_TAG=$CUT_TAG"
echo "[cut campaign] outputs under: $CAFPYANA_GRID_OUT_DIR/dfs/$CUT_TAG/"
