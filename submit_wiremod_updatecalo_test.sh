#!/usr/bin/env bash
# Submit a small WireMod sel_all-updatecalo test (YZ / XTXW / CV).
# Uses a few CAFs per sample so we can verify HDF products + resource usage
# before a full production.
set -euo pipefail

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
cd "$REPO"

# shellcheck disable=SC1091
source envs/venv_py310_cafpyana/bin/activate
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export CAFPYANA_GRID_OUT_DIR="${CAFPYANA_GRID_OUT_DIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out}"

# updatecalo + efield is heavy (10 remakes of calo/PID per CAF).
export JOBSUB_DISK="${JOBSUB_DISK:-40GB}"
export JOBSUB_MEMORY="${JOBSUB_MEMORY:-16GB}"
export JOBSUB_LIFETIME="${JOBSUB_LIFETIME:-12h}"
export JOBSUB_CPU="${JOBSUB_CPU:-7}"

NFILES="${NFILES:-3}"   # CAFs per sample
NGRID="${NGRID:-3}"     # jobs per sample (1 CAF/job when NFILES==NGRID)

LIST_DIR="$REPO/tmp/wiremod_updatecalo_test_lists"
mkdir -p "$LIST_DIR"

YZ_FULL=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list
XTXW_FULL=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list
CV_FULL=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd_xrootd.list

head -n "$NFILES" "$YZ_FULL"   > "$LIST_DIR/wiremod_yz_${NFILES}.list"
head -n "$NFILES" "$XTXW_FULL" > "$LIST_DIR/wiremod_xtxw_${NFILES}.list"
head -n "$NFILES" "$CV_FULL"   > "$LIST_DIR/wiremod_cv_${NFILES}.list"

echo "Lists under $LIST_DIR (NFILES=$NFILES NGRID=$NGRID)"
wc -l "$LIST_DIR"/*.list
echo "jobsub resources: disk=$JOBSUB_DISK mem=$JOBSUB_MEMORY life=$JOBSUB_LIFETIME cpu=$JOBSUB_CPU"

# Refresh token for xrdcp / jobsub
# shellcheck disable=SC1090
source ~/get_token.sh

CFG=configs/numucc_1p0pi/sel_all-updatecalo.py

submit_one() {
  local name="$1" list="$2"
  echo "==== submit $name ===="
  python run_df_maker.py -c "$CFG" -l "$list" -o "$name" -ngrid "$NGRID"
}

submit_one "sel_all-mc-BNB_cosmics-WireModYZ-updatecalo-test"   "$LIST_DIR/wiremod_yz_${NFILES}.list"
submit_one "sel_all-mc-BNB_cosmics-WireModXTXW-updatecalo-test" "$LIST_DIR/wiremod_xtxw_${NFILES}.list"
submit_one "sel_all-mc-BNB_cosmics-calovar-updatecalo-test"     "$LIST_DIR/wiremod_cv_${NFILES}.list"

echo "Submitted. Monitor with: jobsub_q -G sbnd --user munjung"
echo "Outputs under: $CAFPYANA_GRID_OUT_DIR/dfs/*updatecalo-test"
