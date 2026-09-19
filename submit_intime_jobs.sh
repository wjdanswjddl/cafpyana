#!/bin/bash
# Nest sample outputs under dfs/<CUT_TAG>/ when set, e.g.:
#   CUT_TAG=fvfix bash submit_intime_jobs.sh
if [ -n "${CUT_TAG:-}" ]; then
  export CAFPYANA_GRID_SUBDIR="$CUT_TAG"
  echo "[cut campaign] CAFPYANA_GRID_SUBDIR=$CAFPYANA_GRID_SUBDIR"
fi

. ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_all-mc.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/intime/mc_MCP2025B_1e20__v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_xrootd.list -o sel_all-mc-Intime -ngrid 500
. ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/intime/mc_MCP2025B_1e20__v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_xrootd.list -o sel_mup-mc-Intime -ngrid 500
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-mc.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/intime/mc_MCP2025B_1e20__v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_xrootd.list -o sel_2prong-mc-Intime -ngrid 200
