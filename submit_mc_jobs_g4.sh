inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
list="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"

# All three categories in one .df
python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup-g4wgts.py \
    -l "$list" \
    -o sel_mup-wgts_g4 \
    -ngrid 2000

# One category (same as old thin configs)
#FLUX_GROUP=beam python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py \
#  -l "$mc_list" -o sel_mup-wgts_flux_beam -ngrid 200
