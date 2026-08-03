filecv="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_ab"
python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filecv -o sel_2prong-mc-BNB_cosmics-efieldvar_ab

filecv="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_ac"
python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filecv -o sel_2prong-mc-BNB_cosmics-efieldvar_ac

filecv="/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_ad"
python run_df_maker.py -c configs/numucc_1p0pi/sel_2prong-updateefield.py -l $filecv -o sel_2prong-mc-BNB_cosmics-efieldvar_ad
