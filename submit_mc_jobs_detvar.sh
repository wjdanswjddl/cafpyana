#!/bin/bash
# Detector-variation dataframe jobs (WireMod / calo+efield / DENT / SCE).
# Always use sel_all configs — match and walk selection downstream.
#
# sel_all-updatecalo.py writes: hdr, evt_cv/trk_cv, eight calo ± universes,
# and evt_efield/trk_efield (E-field redo) in one job.
set -euo pipefail

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
cd "$REPO"
# shellcheck disable=SC1091
source envs/venv_py310_cafpyana/bin/activate
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"

# Sized from updatecalo test jobs (~0.5 GB RSS, ~2 GB disk, few min/CAF).
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-8GB}"
export JOBSUB_DISK="${JOBSUB_DISK:-20GB}"
export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-12h}"
export JOBSUB_CPU="${JOBSUB_CPU:-7}"

CFG=configs/numucc_1p0pi/sel_all-updatecalo.py
YZ_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list
XTXW_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list

echo "jobsub resources: disk=$JOBSUB_DISK mem=$JOBSUB_MEMORY life=$JOBSUB_LIFETIME cpu=$JOBSUB_CPU"
wc -l "$YZ_LIST" "$XTXW_LIST"

# shellcheck disable=SC1090
source ~/get_token.sh
python run_df_maker.py -c "$CFG" -l "$YZ_LIST" -o sel_all-mc-BNB_cosmics-WireModYZ -ngrid 3000

# shellcheck disable=SC1090
source ~/get_token.sh
python run_df_maker.py -c "$CFG" -l "$XTXW_LIST" -o sel_all-mc-BNB_cosmics-WireModXTXW -ngrid 3000

echo "Submitted WireMod YZ + XTXW. Monitor: jobsub_q -G sbnd --user munjung"
