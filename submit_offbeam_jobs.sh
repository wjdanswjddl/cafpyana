#!/bin/bash
# Nest sample outputs under dfs/<CUT_TAG>/ when set, e.g.:
#   CUT_TAG=fvfix bash submit_offbeam_jobs.sh
if [ -n "${CUT_TAG:-}" ]; then
  export CAFPYANA_GRID_SUBDIR="$CUT_TAG"
  echo "[cut campaign] CAFPYANA_GRID_SUBDIR=$CAFPYANA_GRID_SUBDIR"
fi

. ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-OffBeamLight -ngrid 500
. ~/get_token.sh
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-OffBeamLight -ngrid 50
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_2prong-data-OffBeamLight -ngrid 200
