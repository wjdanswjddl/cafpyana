# VecFFCCQEshape GENIE knob on the Ar23+ sample at loose sel_all
# (knob list: makedf.geniesyst.GENIE_KNOB_GROUPS["VecFF"]).
#
# Same chain as the Ar23p sel_all knob-group jobs: sel_all-geniewgts-knobgroups.py
# with GENIE_KNOB_GROUP set, so worker jobs build HDF keys evt / trk / mcnu / hdr.
#
# Reuses the xrootd list from the Ar23p run (worker scripts xrdcp each line verbatim,
# so the dcache:-prefixed list cannot be used here).
#
# One knob instead of Ar23p's 47: 1000 jobs x ~377 files (vs Ar23p sel_all's 2000
# at ~188 files/job). sel_all is heavier than sel_mup (trk + unselected events),
# so more jobs than the 500-job sel_mup VecFF pass and an 8h lifetime.
#
# Run from the repo root after ``source setup.sh`` and ``. ~/get_token.sh``.

ar23_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
ar23_list="${ar23_inputdir}/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list"

JOBSUB_LIFETIME=8h GENIE_KNOB_GROUP=VecFF python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
  -l "$ar23_list" \
  -o "sel_all-wgts_genie_VecFF" \
  -ngrid 1000
