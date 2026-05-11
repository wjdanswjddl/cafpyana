inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
list="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"

# All BNB flux knobs in one .df (HDF keys evt, hdr)
python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py \
    -l "$list" \
    -o sel_mup-wgts_flux \
    -ngrid 2000

# Flux subset only: sel_mup-wgts_flux_beam.py (etc.)
