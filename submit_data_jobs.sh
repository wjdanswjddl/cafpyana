python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Fixed -ngrid 10
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Fixed -ngrid 10

python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Rolling -ngrid 10
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Rolling -ngrid 10

#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-BNB_cosmics -ngrid 500
#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-Gen1 -ngrid 100
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-Gen1 -ngrid 50
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong_wcandidates-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_2prong_wcandidates-data-BNB_cosmics -ngrid 50
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_FullData1e20_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_2prong-data-BNB_cosmics -ngrid 50

#python run_df_maker.py -c configs/numucc_1p0pi/sel_all-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_Dev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_all-data-BNB_cosmics-Dev -ngrid 10
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-data.py -l /exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/BNB/data_MCP2025C_Spring25_reprocess_Dev_bnblight_v10_06_00_09_flatcaf_sbnd_xrootd.list -o sel_mup-data-BNB_cosmics-Dev -ngrid 10

