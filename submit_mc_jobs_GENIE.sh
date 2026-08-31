# GENIE knob-group dfs via configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py
# Knob lists: makedf.geniesyst.GENIE_KNOB_GROUPS
# Valid GENIE_KNOB_GROUP keys: CCQE ZExp MEC RES nonRES DIS Other Ar23p
#
# Spring CV sample for standard GENIE groups; Ar23+ respin for Ar23p-only knobs.

#spring_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#spring_list="${spring_inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#for g in CCQE MEC RES nonRES DIS Other; do
#  GENIE_KNOB_GROUP=$g python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
#    -l "$spring_list" \
#    -o "sel_mup-wgts_genie_${g}" \
#    -ngrid 1000
#done

ar23_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
ar23_list="${ar23_inputdir}/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list"

for g in Ar23p; do
  GENIE_KNOB_GROUP=$g python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
    -l "$ar23_list" \
    -o "sel_mup-wgts_genie_${g}" \
    -ngrid 4000
done
