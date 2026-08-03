#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list  -o sel_2prong-mc-BNB_cosmics-WireModYZ -ngrid 3000
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/WireMod_XTXW_testprod_xrootd.list  -o sel_2prong-mc-BNB_cosmics-WireModXTXW -ngrid 3000
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list  -o sel_2prong-mc-BNB_cosmics-WireModXTXW -ngrid 3000
#
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd_xrootd.list -o sel_2prong-mc-BNB_cosmics-CV -ngrid 1000

#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_0xSCE_caf_flat_caf_sbnd_xrootd.list -o sel_mup-mc-BNB_cosmics-0xSCE -ngrid 1000
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_2xSCE_caf_flat_caf_sbn_xrootd.list -o sel_mup-mc-BNB_cosmics-2xSCE -ngrid 1000
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updatecalo.py -l /exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd_xrootd.list -o sel_mup-mc-BNB_cosmics-CV -ngrid 1000

# filextxw="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list"
# python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updateefield.py -l $filextxw -o sel_mup-mc-WireModXTXW-efield -ngrid 2000
# #
# fileyz="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list"
# python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updateefield.py -l $fileyz -o sel_mup-mc-WireModYZ-efield -ngrid 2000

# filextxw="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list"
# python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filextxw -o sel_2prong-mc-WireModXTXW-efield -ngrid 2000
# #
# fileyz="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list"
# python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $fileyz -o sel_2prong-mc-WireModYZ-efield -ngrid 2000

filecv="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filecv -o sel_2prong-mc-BNB_cosmics-efieldvar -ngrid 2000
#python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-updatecalo.py -l $filecv -o sel_mup-mc-BNB_cosmics-calovar -ngrid 3000
#
#filecv="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_test_100"
#python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filecv -o sel_2prong-mc-BNB_cosmics-efieldvar 
