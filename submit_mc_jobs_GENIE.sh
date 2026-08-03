#inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#list="${inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#for g in MEC RES nonRES DIS Other; do
#  GENIE_KNOB_GROUP=$g python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
#    -l "$list" \
#    -o "sel_mup-wgts_genie_${g}" \
#    -ngrid 1000
#done

inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/Ar23+
list="${inputdir}/ar23p_respin-xrootd.list"

for g in Ar23p; do
  GENIE_KNOB_GROUP=$g python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
    -l "$list" \
    -o "sel_mup-wgts_genie_${g}" \
    -ngrid 4000
done
