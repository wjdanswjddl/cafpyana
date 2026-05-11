inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
list="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd_test.list"

# GENIE_KNOB_GROUP is exported into each grid worker by run_df_maker.py (-ngrid); without that,
# workers would load all knob groups (evt_CCQE, evt_MEC, ...). One group -> HDF keys evt, mcnu, hdr.

#for g in CCQE MEC RES nonRES DIS Other Ar23p ZExp; do
for g in MEC; do
  GENIE_KNOB_GROUP=$g python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
    -l "$list" \
    -o "sel_mup-wgts_genie_${g}" \
    -ngrid 1
done
