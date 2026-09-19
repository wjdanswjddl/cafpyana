#!/bin/bash
# Nest sample outputs under dfs/<CUT_TAG>/ when set, e.g.:
#   CUT_TAG=fvfix bash submit_data_jobs.sh
if [ -n "${CUT_TAG:-}" ]; then
  export CAFPYANA_GRID_SUBDIR="$CUT_TAG"
  echo "[cut campaign] CAFPYANA_GRID_SUBDIR=$CAFPYANA_GRID_SUBDIR"
fi

#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Fixed -ngrid 10
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Fixed -ngrid 10
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Rolling -ngrid 10
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Rolling -ngrid 10

#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-1e20 -ngrid 500
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_2prong-data-1e20 -ngrid 100
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-1e20-fvfix-chi2fix-real -ngrid 100
#list="/exp/sbnd/app/users/munjung/misc/filelists/data/sanity_data_xrootd.list"
list="/exp/sbnd/app/users/munjung/misc/filelists/data/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list"
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l $list -o sel_mup-data-FixedDev -ngrid 10

