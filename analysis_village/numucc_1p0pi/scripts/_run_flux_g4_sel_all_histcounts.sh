#!/usr/bin/env bash
# Product A+B Flux/G4 from sel_all weight DFs: batched histcounts walk → stream-sum → NPZs.
#
# One sel_all pipeline walk fills cut-stage (Product A) and final (Product B) rate
# histcounts. Flux/G4 are rate-only; the same fractional cov is used for both.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="${REPO}${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONWARNINGS=ignore

STAMP="${STAMP:-$(date +%Y_%m_%d_%H%M%S)}"
OUT_ROOT="${OUT_ROOT:-/pnfs/sbnd/scratch/users/${USER}/cafpyana_tmp/syst_histcounts_from_df_${STAMP}}"
SYST_DISK="${SYST_DISK:-/exp/sbnd/data/users/${USER}/plots/numucc1p0pi/systematics-histcounts-fluxg4-${STAMP}}"
WORKERS="${WORKERS:-8}"
N_UNIVERSE="${N_UNIVERSE:-1000}"
MAX_FILES="${MAX_FILES:-0}"

FLUX_GLOB="${FLUX_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_172234__sel_all-wgts_flux-corrected/*.df}"
G4_GLOB="${G4_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_175940__sel_all-wgts_g4-corrected/*.df}"

SLIM_SKIP="slim,slim_multisim,Flux_slim,Flux_slim_multisim,G4_slim,G4_slim_multisim,Flux,G4"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "$LOG_DIR" "$SYST_DISK" "${OUT_ROOT}/summed"
echo "$OUT_ROOT" > /tmp/syst_histcounts_from_df_campaign.txt
echo "$SYST_DISK" > /tmp/syst_histcounts_from_df_syst_disk.txt

echo "===== histcounts from-df  $(date -Is) ====="
echo "OUT_ROOT=$OUT_ROOT"
echo "SYST_DISK=$SYST_DISK"
echo "WORKERS=$WORKERS N_UNIVERSE=$N_UNIVERSE MAX_FILES=$MAX_FILES"

FAMILY=both FLUX_GLOB="$FLUX_GLOB" G4_GLOB="$G4_GLOB" \
  OUT_ROOT="$OUT_ROOT" WORKERS="$WORKERS" N_UNIVERSE="$N_UNIVERSE" \
  MAX_FILES="$MAX_FILES" INCLUDE_SLIM=1 SAMPLE=mc STAMP="$STAMP" \
  bash "$REPO/analysis_village/numucc_1p0pi/scripts/run_syst_histcounts_from_df.sh"

sum_family() {
  local fam="$1"
  local tag
  tag="$(echo "$fam" | tr '[:upper:]' '[:lower:]')"
  local glob="${OUT_ROOT}/dfs/hist_mc_${tag}/*.df"
  local out="${OUT_ROOT}/summed/hist_mc_${tag}__fullstat.df"
  echo "===== stream-sum $fam $(date -Is) ====="
  python "$REPO/analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum.py" \
    --input-glob "$glob" \
    --out-df "$out" \
    --family "$fam" \
    --mode dense-rate \
    --build-covs \
    --skip-hist-df \
    --checkpoint-every 25 \
    --slim-skip "$SLIM_SKIP"
}

sum_family Flux
sum_family G4

echo "===== NPZs + plots $(date -Is) ====="
python "$REPO/analysis_village/numucc_1p0pi/scripts/histcounts_covs_to_syst_npzs.py" \
  --cov-pkl "${OUT_ROOT}/summed/hist_mc_flux__fullstat.df.rate_covs.pkl" \
  --cov-pkl "${OUT_ROOT}/summed/hist_mc_g4__fullstat.df.rate_covs.pkl" \
  --syst-disk-root "$SYST_DISK"

echo "===== DONE $(date -Is) ====="
echo "histcounts: $OUT_ROOT"
echo "Product A/B NPZs: $SYST_DISK/productA  $SYST_DISK/productB"
