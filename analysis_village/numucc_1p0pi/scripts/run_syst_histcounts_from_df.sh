#!/usr/bin/env bash
# Run Flux / G4 histcounts from cached sel_all weight DFs (no CAF).
#
# Defaults point at the corrected Spring sel_all caches. Override with env vars.
#
# For notebook / kernel-host runs prefer the Python twin:
#   analysis_village/numucc_1p0pi/scripts/run_syst_histcounts_from_df.py
#   (import run_histcounts_from_df — see systematics-histcounts.ipynb).
#
# Examples:
#   # Smoke (2 files, 2 workers)
#   FAMILY=Flux MAX_FILES=2 WORKERS=2 bash run_syst_histcounts_from_df.sh
#
#   # Full flux then g4
#   FAMILY=Flux WORKERS=4 bash run_syst_histcounts_from_df.sh
#   FAMILY=G4   WORKERS=4 bash run_syst_histcounts_from_df.sh
#
#   # Both
#   FAMILY=both WORKERS=4 bash run_syst_histcounts_from_df.sh
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
chmod +x analysis_village/numucc_1p0pi/scripts/run_syst_histcounts_from_df.sh 2>/dev/null || true

# shellcheck disable=SC1091
if [[ -f envs/venv_py310_cafpyana/bin/activate ]]; then
  # Prefer project venv when present
  # shellcheck disable=SC1091
  source envs/venv_py310_cafpyana/bin/activate
fi

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MPLBACKEND=Agg

STAMP="${STAMP:-$(date +%Y_%m_%d_%H%M%S)}"
OUT_ROOT="${OUT_ROOT:-/pnfs/sbnd/scratch/users/${USER}/cafpyana_tmp/syst_histcounts_from_df_${STAMP}}"
FAMILY="${FAMILY:-both}"   # Flux | G4 | both
WORKERS="${WORKERS:-2}"
N_UNIVERSE="${N_UNIVERSE:-1000}"
MAX_FILES="${MAX_FILES:-0}"
INCLUDE_SLIM="${INCLUDE_SLIM:-1}"
SAMPLE="${SAMPLE:-mc}"

FLUX_GLOB="${FLUX_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_172234__sel_all-wgts_flux-corrected_updated/*.df}"
G4_GLOB="${G4_GLOB:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_04_175940__sel_all-wgts_g4-corrected_updated/*.df}"

SLIM_FLAG=(--include-slim)
if [[ "$INCLUDE_SLIM" == "0" || "$INCLUDE_SLIM" == "false" ]]; then
  SLIM_FLAG=(--no-include-slim)
fi

run_one() {
  local fam="$1"
  local glob="$2"
  local out_dir="${OUT_ROOT}/dfs/hist_mc_$(echo "$fam" | tr '[:upper:]' '[:lower:]')"
  mkdir -p "$out_dir" "${OUT_ROOT}/logs"
  local failed_log="${OUT_ROOT}/logs/failed_${fam}.txt"
  echo "===== FAMILY=$fam  OUT=$out_dir  WORKERS=$WORKERS  MAX_FILES=$MAX_FILES ====="
  python analysis_village/numucc_1p0pi/scripts/syst_histcounts_from_df_parallel.py \
    --input-glob "$glob" \
    --out-dir "$out_dir" \
    --family "$fam" \
    --workers "$WORKERS" \
    --n-universe "$N_UNIVERSE" \
    --sample "$SAMPLE" \
    --max-files "$MAX_FILES" \
    --failed-log "$failed_log" \
    "${SLIM_FLAG[@]}"
}

echo "OUT_ROOT=$OUT_ROOT"
echo "$OUT_ROOT" > /tmp/syst_histcounts_from_df_campaign.txt
mkdir -p "$OUT_ROOT"

case "$(echo "$FAMILY" | tr '[:upper:]' '[:lower:]')" in
  flux)
    run_one Flux "$FLUX_GLOB"
    ;;
  g4)
    run_one G4 "$G4_GLOB"
    ;;
  both|all)
    run_one Flux "$FLUX_GLOB"
    run_one G4 "$G4_GLOB"
    ;;
  *)
    echo "Unknown FAMILY=$FAMILY (use Flux, G4, or both)" >&2
    exit 1
    ;;
esac

echo "===== DONE  campaign=$OUT_ROOT ====="
echo "Point notebook NUMUCC_SYST_HIST_EXTRA_ROOTS at this campaign (or merge under HIST_ROOT)."
