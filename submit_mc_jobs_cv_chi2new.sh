#!/bin/bash
# Submit CV-only updatecalo (chi2_*_new recalculation) on the Sep-4 MCP2025C
# CV sample. No calo ± universes, no efield.
set -euo pipefail

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
cd "$REPO"
# shellcheck disable=SC1091
source envs/venv_py310_cafpyana/bin/activate
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"

# Lighter than full updatecalo (1 remake vs 10); keep same headroom as detvar.
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-8GB}"
export JOBSUB_DISK="${JOBSUB_DISK:-20GB}"
export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-12h}"
export JOBSUB_CPU="${JOBSUB_CPU:-7}"

CFG=configs/numucc_1p0pi/sel_all-updatecalo-cvonly.py
# Same CAFs as Sep-4 sel_all-mc-CV (2026_09_04_172912__sel_all-mc-CV)
CV_LIST=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list

echo "jobsub resources: disk=$JOBSUB_DISK mem=$JOBSUB_MEMORY life=$JOBSUB_LIFETIME cpu=$JOBSUB_CPU"
wc -l "$CV_LIST"

# shellcheck disable=SC1090
source ~/get_token.sh
python run_df_maker.py -c "$CFG" -l "$CV_LIST" -o sel_all-mc-CV-updatecalo-cvonly -ngrid 2000

echo "Submitted CV chi2_new-only. Monitor: jobsub_q -G sbnd --user munjung"
