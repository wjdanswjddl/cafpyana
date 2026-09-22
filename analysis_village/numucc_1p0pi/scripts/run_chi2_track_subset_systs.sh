#!/usr/bin/env bash
# Focused Product A campaign: chi2_avg len>50 + not_mu @ 2prong-vtxdist only.
#
# Does NOT recompute existing all-track packs. Set NUMUCC_CUT_STAGE_SLUGS +
# NUMUCC_SKIP_FINAL_STAGE so walks only fill the four new slugs.
#
# Usage:
#   bash run_chi2_track_subset_systs.sh              # launch all sources
#   FAMILY=Flux bash run_chi2_track_subset_systs.sh  # one family
#   STAGE=merge bash run_chi2_track_subset_systs.sh # merge when ready
#
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="${REPO}${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONWARNINGS=ignore

SLUGS="chi2_avg_mu_len50__at_2prong-vtxdist,chi2_avg_p_len50__at_2prong-vtxdist,chi2_avg_mu_not_mu__at_2prong-vtxdist,chi2_avg_p_not_mu__at_2prong-vtxdist"
export NUMUCC_CUT_STAGE_SLUGS="$SLUGS"
export NUMUCC_SKIP_FINAL_STAGE=1
export GENIE_VAR_SAVE_NAMES="$SLUGS"

STAMP="${STAMP:-$(date +%Y_%m_%d_%H%M%S)}"
OUT_ROOT="${OUT_ROOT:-/pnfs/sbnd/scratch/users/${USER}/cafpyana_tmp/chi2_track_subset_${STAMP}}"
SYST_DISK="${SYST_DISK:-/exp/sbnd/data/users/${USER}/plots/numucc1p0pi/systematics-chi2subset-${STAMP}}"
WORKERS="${WORKERS:-8}"
N_UNIVERSE="${N_UNIVERSE:-1000}"
MAX_FILES="${MAX_FILES:-0}"
FAMILY="${FAMILY:-all}"   # Flux|G4|Cosmics|GENIE|Detector|all|merge
STAGE="${STAGE:-run}"     # run|merge

FLUX_GLOB="${FLUX_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_172234__sel_all-wgts_flux-corrected_updated/*.df}"
G4_GLOB="${G4_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_175940__sel_all-wgts_g4-corrected_updated/*.df}"

mkdir -p "$OUT_ROOT/logs" "$SYST_DISK"
echo "$OUT_ROOT" > /tmp/chi2_track_subset_out.txt
echo "$SYST_DISK" > /tmp/chi2_track_subset_syst_disk.txt
echo "===== chi2 track subset  $(date -Is)  OUT=$OUT_ROOT ====="

run_flux_g4() {
  local fam="$1"
  local glob="$2"
  local tag
  tag="$(echo "$fam" | tr '[:upper:]' '[:lower:]')"
  local out_dir="${OUT_ROOT}/dfs/hist_mc_${tag}"
  mkdir -p "$out_dir"
  echo "===== $fam histcounts $(date -Is) ====="
  python analysis_village/numucc_1p0pi/scripts/syst_histcounts_from_df_parallel.py \
    --input-glob "$glob" \
    --out-dir "$out_dir" \
    --family "$fam" \
    --workers "$WORKERS" \
    --n-universe "$N_UNIVERSE" \
    --sample mc \
    --max-files "$MAX_FILES" \
    --failed-log "${OUT_ROOT}/logs/failed_${fam}.txt" \
    --include-slim

  local summed="${OUT_ROOT}/summed/hist_mc_${tag}__fullstat.df"
  mkdir -p "${OUT_ROOT}/summed"
  python analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum.py \
    --input-glob "${out_dir}/*.df" \
    --out-df "$summed" \
    --family "$fam" \
    --mode dense-rate \
    --build-covs \
    --skip-hist-df \
    --checkpoint-every 25 \
    --slim-skip "slim,slim_multisim,Flux_slim,Flux_slim_multisim,G4_slim,G4_slim_multisim,Flux,G4"

  python analysis_village/numucc_1p0pi/scripts/histcounts_covs_to_syst_npzs.py \
    --cov-pkl "${summed}.rate_covs.pkl" \
    --syst-disk-root "$SYST_DISK" \
    --no-plots
}

run_cosmics() {
  echo "===== Cosmics $(date -Is) ====="
  COSMICS_WORK="${OUT_ROOT}/cosmics" \
  COSMICS_WORKERS="$WORKERS" \
  python analysis_village/numucc_1p0pi/scripts/_launch_cosmics_sel_all_updated.py \
    2>&1 | tee "${OUT_ROOT}/logs/cosmics.log"
  # launcher writes under its SYST_DISK; copy note in log
}

run_genie() {
  echo "===== GENIE FSI_compare slim_v3 rate $(date -Is) ====="
  export MC_DF_STAGE=sel_all
  export GENIE_RUN_GROUPS=FSI_compare
  export WORK_BASE="${OUT_ROOT}/genie"
  export WORKERS
  export MAX_FILES
  export SKIP_MERGE=0
  # Point syst disk away from PRL; we merge manually
  export NUMUCC_SYST_DISK_ROOT="${SYST_DISK}/genie_work"
  mkdir -p "$NUMUCC_SYST_DISK_ROOT"
  bash analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh \
    --genie-groups FSI_compare \
    --workers "$WORKERS" \
    ${MAX_FILES:+--max-files "$MAX_FILES"} \
    2>&1 | tee "${OUT_ROOT}/logs/genie.log"
}

run_detector() {
  echo "===== Detector WireMod+DENT $(date -Is) ====="
  python analysis_village/numucc_1p0pi/scripts/run_chi2_subset_detector.py \
    --out-dir "${OUT_ROOT}/detector" \
    --workers "$WORKERS" \
    --max-files "$MAX_FILES" \
    2>&1 | tee "${OUT_ROOT}/logs/detector.log"
}

do_merge() {
  echo "===== merge into Product A $(date -Is) ====="
  local flux_npz g4_npz
  flux_npz="${SYST_DISK}/productA/Flux/flux_syst_dict.npz"
  g4_npz="${SYST_DISK}/productA/G4/g4_syst_dict.npz"
  args=()
  [[ -f "$flux_npz" ]] && args+=(--flux-npz "$flux_npz")
  [[ -f "$g4_npz" ]] && args+=(--g4-npz "$g4_npz")
  # cosmics / detector / genie paths filled by their runners
  [[ -f "${OUT_ROOT}/cosmics_syst_dict.npz" ]] && args+=(--cosmics-npz "${OUT_ROOT}/cosmics_syst_dict.npz")
  [[ -f "${OUT_ROOT}/detector/detector_sel_syst_dict.npz" ]] && args+=(--detector-npz "${OUT_ROOT}/detector/detector_sel_syst_dict.npz")
  [[ -f "${SYST_DISK}/genie_work/GENIE/cov_mat_dict.pkl" ]] && args+=(--genie-cov-pkl "${SYST_DISK}/genie_work/GENIE/cov_mat_dict.pkl" --genie-slim-npz "${SYST_DISK}/genie_work/GENIE/genie_slim_v3.npz")
  python analysis_village/numucc_1p0pi/scripts/merge_chi2_subset_into_product_a.py "${args[@]}"
}

if [[ "$STAGE" == "merge" ]]; then
  do_merge
  exit 0
fi

case "$(echo "$FAMILY" | tr '[:upper:]' '[:lower:]')" in
  flux) run_flux_g4 Flux "$FLUX_GLOB" ;;
  g4) run_flux_g4 G4 "$G4_GLOB" ;;
  cosmics) run_cosmics ;;
  genie) run_genie ;;
  detector) run_detector ;;
  all)
    # Flux + G4 sequential (IO heavy); others can be started separately in parallel shells
    run_flux_g4 Flux "$FLUX_GLOB"
    run_flux_g4 G4 "$G4_GLOB"
    run_cosmics
    run_genie
    run_detector
    do_merge
    ;;
  *)
    echo "Unknown FAMILY=$FAMILY" >&2
    exit 1
    ;;
esac

echo "===== DONE $(date -Is) ====="
echo "OUT_ROOT=$OUT_ROOT"
echo "SYST_DISK=$SYST_DISK"
echo "Merge with: STAGE=merge OUT_ROOT=$OUT_ROOT SYST_DISK=$SYST_DISK bash $0"
