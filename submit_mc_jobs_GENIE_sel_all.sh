# GENIE weights at loose sel_all via configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py
# Set GENIE_KNOB_GROUP=slim for bundled mc.GENIE.univ_* (see makedf/getsyst.py),
# or CCQE|MEC|RES|nonRES|DIS|Other|Ar23p for per-group tables.
#
# Output tables per split: evt, trk, mcnu, hdr (+ histpot/histgenevt from run_df_maker)
# Use with get_systematics_genie.py chunk-map --input-stage sel_all

#spring_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#spring_list="${spring_inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#GENIE_KNOB_GROUP=slim python run_df_maker.py \
#  -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
#  -l "$spring_list" \
#  -o "sel_all-wgts_genie_slim" \
#  -ngrid 3000

#spring_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
#spring_list="${spring_inputdir}/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list"
#
#for g in CCQE MEC RES nonRES DIS Other; do
#  GENIE_KNOB_GROUP=$g python run_df_maker.py \
#    -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
#    -l "$spring_list" \
#    -o "sel_all-wgts_genie_${g}" \
#    -ngrid 4000
#done

# Must use the *_xrootd.list (root://…), NOT flatcaf_sbnd.list (/pnfs/…).
# 2026_08_31_142802__sel_all-wgts_genie_Ar23p used the pnfs list → all 6000 jobs
# failed at xrdcp ("no such file or directory") and wrote empty 6.9K stub .df files.
ar23_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
ar23_list="${ar23_inputdir}/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list"

for g in Ar23p; do
  GENIE_KNOB_GROUP=$g python run_df_maker.py \
    -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
    -l "$ar23_list" \
    -o "sel_all-wgts_genie_${g}" \
    -ngrid 2000
done
