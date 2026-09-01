#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Fixed -ngrid 10
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Fixed -ngrid 10
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Rolling -ngrid 10
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Rolling -ngrid 10

#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-1e20 -ngrid 500
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_2prong-data-1e20 -ngrid 100
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-1e20-fvfix-chi2fix -ngrid 100

