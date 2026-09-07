#!/usr/bin/env bash
# Full-statistics syst_histcounts grid campaign (DRY-RUN by default).
#
# Creates a fresh campaign directory under pnfs scratch, prints the planned
# jobsub invocations, and only submits if you pass --submit.
#
# Physics layout (do not split GENIE by CCQE/MEC/… — except Ar23p):
#   * hist_mc_genie       — non-Ar23p GENIE knobs + slim_multisim + slim (Spring CV CAFs)
#   * hist_mc_genie_Ar23p — Ar23p knobs only (AR23plus knobs CAF list; separate files)
#   * hist_mc_flux / g4   — Flux / G4 + slim products
#   * dirt: hist_dirt_genie / _flux / _g4 (no Ar23p dirt sample)
#   * unisim nowgt: WireMod, DENT, intime, offbeam + mc CV for WireMod pairing
#
# IMPORTANT
# ---------
# * Grid workers ``git clone`` ``release/numucc_1p0pi`` from GitHub
#   (``bin/grid_executable.sh``). Push the histcounts / VariableConfig-freeze
#   commits before --submit, or workers will run stale code.
# * Every campaign uses a NEW ``CAFPYANA_GRID_OUT_DIR`` so outputs never mix.
#
# Usage:
#   bash submit_syst_histcounts_full.sh              # print plan only
#   bash submit_syst_histcounts_full.sh --submit     # actually jobsub (after OK)
#
# Tunables (env):
#   NGRID_GENIE=1500 NGRID_FLUX=2500 NGRID_UNISIM=1000
#   JOBSUB_DISK JOBSUB_MEMORY JOBSUB_LIFETIME JOBSUB_CPU
#     (defaults: 10GB / 10GB / 6h / 7)
set -euo pipefail

DO_SUBMIT=0
for arg in "$@"; do
  case "$arg" in
    --submit) DO_SUBMIT=1 ;;
    -h|--help)
      sed -n '1,45p' "$0"
      exit 0
      ;;
  esac
done

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
cd "$REPO"

LIST_MC=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
LIST_MC_AR23P=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/BNB_cosmics/mc_SBND2026A_AR23plus_knobs_BNBLight_CV_v1_00_01_flatcaf_sbnd_xrootd.list
LIST_DIRT=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/lowE/mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
LIST_WM_SV=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_SV_v10_06_00_10_flatcaf_sbnd_xrootd.list
LIST_WM_XTXW=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/WireMod/mc_SBND2026A_prodgenie_corsika_proton_rockbox_sbnd_wiremod_X-ThetaXW_v10_06_00_10_flatcaf_sbnd_xrootd.list
LIST_DENT_CV=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/CV_highstats_xrootd.list
LIST_DENT_VAR=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/DENT_highstats_xrootd.list
LIST_INTIME=/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/2025Spring_v10_06_00_09/intime/mc_MCP2025B_1e20__v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_xrootd.list
LIST_OFFBEAM=/exp/sbnd/app/users/munjung/misc/filelists/data/2025Spring_v10_06_00_09/OffBeamLight/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd_xrootd.list

nfiles() { grep -cvE '^\s*(#|$)' "$1"; }

N_MC=$(nfiles "$LIST_MC")
N_MC_AR23P=$(nfiles "$LIST_MC_AR23P")
N_DIRT=$(nfiles "$LIST_DIRT")
N_WM_SV=$(nfiles "$LIST_WM_SV")
N_WM_XTXW=$(nfiles "$LIST_WM_XTXW")
N_DENT_CV=$(nfiles "$LIST_DENT_CV")
N_DENT_VAR=$(nfiles "$LIST_DENT_VAR")
N_INTIME=$(nfiles "$LIST_INTIME")
N_OFFBEAM=$(nfiles "$LIST_OFFBEAM")

# Job counts (capped at nfiles by run_df_maker).
# Tuned for 10GB / 6h: fewer files per job is fine (more -N is OK).
NGRID_GENIE=${NGRID_GENIE:-${NGRID_MC:-1500}}
NGRID_GENIE_AR23P=${NGRID_GENIE_AR23P:-${NGRID_GENIE}}
NGRID_FLUX=${NGRID_FLUX:-2500}
NGRID_G4=${NGRID_G4:-2500}
NGRID_DIRT=${NGRID_DIRT:-800}
NGRID_UNISIM=${NGRID_UNISIM:-1000}
NGRID_COSMICS=${NGRID_COSMICS:-600}

# Full universe counts. Unset smoke leftovers from the same shell.
if [[ "${SYST_HIST_FORCE_PROD_NUNIV:-1}" == "1" ]]; then
  export SYST_HIST_GENIE_NUNIV=100
  export SYST_HIST_FLUX_NUNIV=1000
  export SYST_HIST_G4_NUNIV=1000
else
  export SYST_HIST_GENIE_NUNIV=${SYST_HIST_GENIE_NUNIV:-100}
  export SYST_HIST_FLUX_NUNIV=${SYST_HIST_FLUX_NUNIV:-1000}
  export SYST_HIST_G4_NUNIV=${SYST_HIST_G4_NUNIV:-1000}
fi
# Ar23p knobs live only on the AR23plus CAF sample — never on Spring CV.
# Main hist_*_genie jobs exclude Ar23p; WAVE=ar23p / hist_mc_genie_Ar23p is separate.
export SYST_HIST_EXCLUDE_AR23P=1

# Match typical cafpyana weight-job resources; more -N instead of long walltime.
export JOBSUB_DISK=${JOBSUB_DISK:-10GB}
export JOBSUB_MEMORY=${JOBSUB_MEMORY:-10GB}
export JOBSUB_LIFETIME=${JOBSUB_LIFETIME:-6h}
export JOBSUB_CPU=${JOBSUB_CPU:-7}

STAMP=$(date +%Y_%m_%d_%H%M%S)
CAMPAIGN_ROOT_DEFAULT=/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp
if [[ -n "${CAMPAIGN:-}" ]]; then
  # Reuse an existing campaign dir (append more waves without mixing stamps).
  :
elif [[ "$DO_SUBMIT" -eq 1 ]]; then
  if ! mkdir -p "$CAMPAIGN_ROOT_DEFAULT" 2>/dev/null; then
    echo "ERROR: cannot create $CAMPAIGN_ROOT_DEFAULT (needed for grid outputs)" >&2
    exit 1
  fi
  CAMPAIGN=$CAMPAIGN_ROOT_DEFAULT/syst_histcounts_full_${STAMP}
else
  # Dry-run: keep plan artifacts in-repo so pnfs write is not required.
  CAMPAIGN=$REPO/tmp/syst_histcounts_full_plan_${STAMP}
fi
# Isolate this campaign's dfs/logs from other cafpyana_out traffic.
export CAFPYANA_GRID_OUT_DIR="$CAMPAIGN"
export CAFPYANA_WD="$REPO"
export CAFPYANA_DIR="$REPO"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export BEARER_TOKEN_FILE="${BEARER_TOKEN_FILE:-/tmp/bt_u$(id -u)}"

mkdir -p "$CAMPAIGN"
# shellcheck disable=SC1091
source envs/venv_py310_cafpyana/bin/activate

python - <<PY || true
from analysis_village.numucc_1p0pi.syst_histcounts import write_var_config_snapshot_json
try:
    write_var_config_snapshot_json("$CAMPAIGN/variable_configs.json")
    print("wrote $CAMPAIGN/variable_configs.json")
except Exception as ex:
    print("WARNING: could not refresh variable_configs.json (%s); continuing" % ex)
PY

CFG=configs/numucc_1p0pi/syst_histcounts.py

files_per() {
  local nfiles=$1 ngrid=$2
  if (( ngrid < 1 )); then ngrid=1; fi
  if (( ngrid > nfiles )); then ngrid=$nfiles; fi
  echo $(( (nfiles + ngrid - 1) / ngrid ))
}

n_eff() {
  local nfiles=$1 ngrid=$2
  if (( ngrid > nfiles )); then echo "$nfiles"; else echo "$ngrid"; fi
}

echo "============================================================"
echo " syst_histcounts FULL campaign plan"
echo " CAMPAIGN=$CAMPAIGN"
echo " DO_SUBMIT=$DO_SUBMIT"
echo " NUNIV GENIE=$SYST_HIST_GENIE_NUNIV FLUX=$SYST_HIST_FLUX_NUNIV G4=$SYST_HIST_G4_NUNIV"
echo " Ar23p: EXCLUDED from hist_*_genie (SYST_HIST_EXCLUDE_AR23P=$SYST_HIST_EXCLUDE_AR23P);"
echo "        separate hist_mc_genie_Ar23p on AR23plus CAF list (N=$N_MC_AR23P)"
echo " Resources: disk=$JOBSUB_DISK mem=$JOBSUB_MEMORY life=$JOBSUB_LIFETIME cpu=$JOBSUB_CPU"
echo "============================================================"
printf '%-28s %8s %6s %10s\n' "job" "nfiles" "ngrid" "files/job"
printf '%-28s %8s %6s %10s\n' "----------------------------" "--------" "------" "----------"

plan_line() {
  local name=$1 nfiles=$2 ngrid=$3
  local fp ng
  fp=$(files_per "$nfiles" "$ngrid")
  ng=$(n_eff "$nfiles" "$ngrid")
  printf '%-28s %8d %6d %10d\n' "$name" "$nfiles" "$ng" "$fp"
}

TOTAL_JOBS=0
add_jobs() {
  local nfiles=$1 ngrid=$2
  local ng
  ng=$(n_eff "$nfiles" "$ngrid")
  TOTAL_JOBS=$((TOTAL_JOBS + ng))
}

# Family split only — never split GENIE by mode (slim needs all non-Ar23p knobs).
# Ar23p is a separate sample/job (WAVE=ar23p).
plan_line "hist_mc_genie" "$N_MC" "$NGRID_GENIE"
add_jobs "$N_MC" "$NGRID_GENIE"
plan_line "hist_mc_genie_Ar23p" "$N_MC_AR23P" "$NGRID_GENIE_AR23P"
add_jobs "$N_MC_AR23P" "$NGRID_GENIE_AR23P"
plan_line "hist_mc_flux" "$N_MC" "$NGRID_FLUX"
add_jobs "$N_MC" "$NGRID_FLUX"
plan_line "hist_mc_g4" "$N_MC" "$NGRID_G4"
add_jobs "$N_MC" "$NGRID_G4"
plan_line "hist_mc_cv_nowgt" "$N_MC" "$NGRID_UNISIM"
add_jobs "$N_MC" "$NGRID_UNISIM"

plan_line "hist_dirt_genie" "$N_DIRT" "$NGRID_DIRT"
add_jobs "$N_DIRT" "$NGRID_DIRT"
plan_line "hist_dirt_flux" "$N_DIRT" "$NGRID_DIRT"
add_jobs "$N_DIRT" "$NGRID_DIRT"
plan_line "hist_dirt_g4" "$N_DIRT" "$NGRID_DIRT"
add_jobs "$N_DIRT" "$NGRID_DIRT"

plan_line "hist_wiremod_sv" "$N_WM_SV" "$NGRID_UNISIM"
add_jobs "$N_WM_SV" "$NGRID_UNISIM"
plan_line "hist_wiremod_xtxw" "$N_WM_XTXW" "$NGRID_UNISIM"
add_jobs "$N_WM_XTXW" "$NGRID_UNISIM"
plan_line "hist_dent_cv" "$N_DENT_CV" "$NGRID_UNISIM"
add_jobs "$N_DENT_CV" "$NGRID_UNISIM"
plan_line "hist_dent_var" "$N_DENT_VAR" "$NGRID_UNISIM"
add_jobs "$N_DENT_VAR" "$NGRID_UNISIM"
plan_line "hist_intime" "$N_INTIME" "$NGRID_COSMICS"
add_jobs "$N_INTIME" "$NGRID_COSMICS"
plan_line "hist_offbeam" "$N_OFFBEAM" "$NGRID_COSMICS"
add_jobs "$N_OFFBEAM" "$NGRID_COSMICS"

echo
echo "Approx total jobsub processes if everything is submitted at once: ~$TOTAL_JOBS"
echo "Optional waves (family-level only):"
echo "  WAVE=1      hist_mc_genie + flux + g4 + mc_cv_nowgt  (Ar23p excluded)"
echo "  WAVE=2      hist_dirt_genie + flux + g4"
echo "  WAVE=3      WireMod / DENT / cosmics"
echo "  WAVE=genie   hist_mc_genie + hist_dirt_genie (non-Ar23p only)"
echo "  WAVE=ar23p   hist_mc_genie_Ar23p only (AR23plus CAF list)"
echo "  WAVE=flux_g4 hist_mc/dirt flux + g4 only"
echo "Reuse campaign: CAMPAIGN=/path/to/syst_histcounts_full_... WAVE=flux_g4 bash $0 --submit"
echo "Override e.g. NGRID_GENIE=2000 NGRID_FLUX=3000 bash $0"
echo

WAVE=${WAVE:-all}

submit_one() {
  local name=$1 mode=$2 sample=$3 list=$4 ngrid=$5
  local extra="${6:-}"
  echo "---- $name  mode=$mode sample=$sample ngrid=$ngrid $extra"
  if [[ "$DO_SUBMIT" -ne 1 ]]; then
    return 0
  fi
  # Default: exclude Ar23p (Spring CV). Override via $extra for Ar23p jobs.
  # shellcheck disable=SC2086
  env SYST_HIST_MODE="$mode" SYST_HIST_SAMPLE="$sample" \
    SYST_HIST_EXCLUDE_AR23P="${SYST_HIST_EXCLUDE_AR23P:-1}" $extra \
    python run_df_maker.py -c "$CFG" -l "$list" -o "$name" -ngrid "$ngrid"
}

want_wave() {
  local w=$1
  [[ "$WAVE" == "all" || "$WAVE" == "$w" ]]
}

if [[ "$DO_SUBMIT" -ne 1 ]]; then
  echo "DRY RUN only. Re-run with --submit after you approve the plan."
  echo "Example wave-1 submit:"
  echo "  WAVE=1 bash submit_syst_histcounts_full.sh --submit"
  echo "Campaign dir (snapshot already written): $CAMPAIGN"
  exit 0
fi

echo "Submitting WAVE=$WAVE (token + jobsub required)..."
# shellcheck disable=SC1090
source ~/get_token.sh
which jobsub_submit >/dev/null 2>&1 || {
  # shellcheck disable=SC1091
  source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
  setup jobsub_client
}

# Unset any leftover group filter from the submit shell (except Ar23p wave).
unset GENIE_KNOB_GROUP || true

if want_wave 1; then
  submit_one hist_mc_genie genie mc "$LIST_MC" "$NGRID_GENIE"
  submit_one hist_mc_flux flux mc "$LIST_MC" "$NGRID_FLUX"
  submit_one hist_mc_g4 g4 mc "$LIST_MC" "$NGRID_G4"
  submit_one hist_mc_cv_nowgt nowgt mc "$LIST_MC" "$NGRID_UNISIM"
fi

if want_wave 2; then
  submit_one hist_dirt_genie genie dirt "$LIST_DIRT" "$NGRID_DIRT"
  submit_one hist_dirt_flux flux dirt "$LIST_DIRT" "$NGRID_DIRT"
  submit_one hist_dirt_g4 g4 dirt "$LIST_DIRT" "$NGRID_DIRT"
fi

# GENIE non-Ar23p only (full family + slim). Not part of WAVE=all.
if [[ "$WAVE" == "genie" ]]; then
  submit_one hist_mc_genie genie mc "$LIST_MC" "$NGRID_GENIE"
  submit_one hist_dirt_genie genie dirt "$LIST_DIRT" "$NGRID_DIRT"
fi

# Ar23p knobs on AR23plus CAF sample only (not Spring CV).
if [[ "$WAVE" == "all" || "$WAVE" == "ar23p" ]]; then
  submit_one hist_mc_genie_Ar23p genie mc "$LIST_MC_AR23P" "$NGRID_GENIE_AR23P" \
    "GENIE_KNOB_GROUP=Ar23p SYST_HIST_EXCLUDE_AR23P=0"
fi

# Flux + G4 only (MC + dirt). Not part of WAVE=all.
if [[ "$WAVE" == "flux_g4" ]]; then
  submit_one hist_mc_flux flux mc "$LIST_MC" "$NGRID_FLUX"
  submit_one hist_mc_g4 g4 mc "$LIST_MC" "$NGRID_G4"
  submit_one hist_dirt_flux flux dirt "$LIST_DIRT" "$NGRID_DIRT"
  submit_one hist_dirt_g4 g4 dirt "$LIST_DIRT" "$NGRID_DIRT"
fi

if want_wave 3; then
  # CV for WireMod pairing (nowgt) — remaining with unisim samples.
  submit_one hist_mc_cv_nowgt nowgt mc "$LIST_MC" "$NGRID_UNISIM"
  submit_one hist_wiremod_sv nowgt mc "$LIST_WM_SV" "$NGRID_UNISIM"
  submit_one hist_wiremod_xtxw nowgt mc "$LIST_WM_XTXW" "$NGRID_UNISIM"
  submit_one hist_dent_cv nowgt mc "$LIST_DENT_CV" "$NGRID_UNISIM"
  submit_one hist_dent_var nowgt mc "$LIST_DENT_VAR" "$NGRID_UNISIM"
  submit_one hist_intime nowgt intime "$LIST_INTIME" "$NGRID_COSMICS"
  submit_one hist_offbeam nowgt offbeam "$LIST_OFFBEAM" "$NGRID_COSMICS"
fi

echo
echo "Submitted WAVE=$WAVE under CAFPYANA_GRID_OUT_DIR=$CAFPYANA_GRID_OUT_DIR"
echo "Notebook: export NUMUCC_SYST_HIST_ROOT=$CAMPAIGN/dfs"
echo "(grid outputs land in dfs/<timestamp>__<tag>/ )"
