# VecFFCCQEshape GENIE knob on the Ar23+ sample (knob list: makedf.geniesyst.GENIE_KNOB_GROUPS["VecFF"]).
#
# Same chain as the other Ar23p knob-group jobs: sel_mup-geniewgts-knobgroups.py with
# GENIE_KNOB_GROUP set, so worker jobs build HDF keys evt / mcnu / hdr.
#
# Reuses the xrootd list from the Ar23p run (worker scripts xrdcp each line verbatim,
# so the dcache:-prefixed list cannot be used here).
#
# One knob instead of Ar23p's 47, so fewer/fatter jobs than that run's -ngrid 2000:
# 500 jobs x ~753 files. The 2000-job Ar23p pass took ~65 min worst case at 188
# files/job, so the longer lifetime keeps 4x-bigger jobs clear of the 3h default.
#
# Run from the repo root after ``source setup.sh``.

ar23_inputdir=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics
ar23_list="${ar23_inputdir}/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list"

JOBSUB_LIFETIME=6h GENIE_KNOB_GROUP=VecFF python run_df_maker.py \
  -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py \
  -l "$ar23_list" \
  -o "sel_mup-wgts_genie_VecFF" \
  -ngrid 500
