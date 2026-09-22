#!/usr/bin/env bash
# Memory-capped parallel stream-sum: Flux then G4 (never concurrent).
# Keep total host used memory under ~80 GiB (Available stays > ~50 GiB on 131 GiB box).
# Default WORKERS=5 → peak ~20–25 GiB for our jobs.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="${REPO}${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONWARNINGS=ignore

OUT_ROOT="${OUT_ROOT:-/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp/syst_histcounts_from_df_2026_09_18_031124}"
SYST_DISK="${SYST_DISK:-/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-histcounts-fluxg4-2026_09_18_031124}"
WORKERS="${WORKERS:-5}"
# Abort if host used memory would exceed this (GiB). Used = MemTotal - MemAvailable.
MAX_USED_GB="${MAX_USED_GB:-80}"
SLIM_SKIP="slim,slim_multisim,Flux_slim,Flux_slim_multisim,G4_slim,G4_slim_multisim,Flux,G4"

FLUX_OUT="${OUT_ROOT}/summed/hist_mc_flux__fullstat.df"
G4_OUT="${OUT_ROOT}/summed/hist_mc_g4__fullstat.df"
FLUX_SEED_DENSE="${FLUX_OUT}.stream_state.json.dense.pkl"
FLUX_SEED_JSON="${FLUX_OUT}.stream_state.json"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "$LOG_DIR" "${OUT_ROOT}/summed" "$SYST_DISK"

LOG="${LOG_DIR}/parallel_reduce_memcap_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

used_gb() {
  awk '/MemTotal:/ {t=$2} /MemAvailable:/ {a=$2} END {printf "%.1f", (t-a)/1024/1024}' /proc/meminfo
}

watch_mem() {
  # background watchdog: kill children if used memory exceeds cap
  local cap="$1"
  while true; do
    u=$(used_gb)
    # compare as integers via awk
    over=$(awk -v u="$u" -v c="$cap" 'BEGIN {print (u+0 > c+0) ? 1 : 0}')
    if [[ "$over" == "1" ]]; then
      echo "[watchdog] USED=${u}GiB exceeds cap ${cap}GiB — killing stream_sum_parallel" >&2
      pkill -f 'syst_histcounts_stream_sum_parallel.py' 2>/dev/null || true
      exit 99
    fi
    sleep 15
  done
}

echo "===== memcap parallel reduce $(date -Is) ====="
echo "OUT_ROOT=$OUT_ROOT WORKERS=$WORKERS MAX_USED_GB=$MAX_USED_GB"
echo "host used now: $(used_gb) GiB"
echo "LOG=$LOG"

# Clear any leftover worker pools only (never pkill this bash script by name).
pkill -f 'syst_histcounts_stream_sum_parallel.py' 2>/dev/null || true
sleep 2

# Clear incomplete partials from the aborted 12+12 run
rm -rf "${FLUX_OUT}.partials" "${G4_OUT}.partials"
mkdir -p "${FLUX_OUT}.partials" "${G4_OUT}.partials"

watch_mem "$MAX_USED_GB" &
WATCH_PID=$!
trap 'kill $WATCH_PID 2>/dev/null || true' EXIT

echo "===== Flux (seeded, workers=$WORKERS) used=$(used_gb)GiB ====="
python "$REPO/analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum_parallel.py" \
  --input-glob "${OUT_ROOT}/dfs/hist_mc_flux/*.df" \
  --out-df "$FLUX_OUT" \
  --family Flux \
  --workers "$WORKERS" \
  --seed-dense "$FLUX_SEED_DENSE" \
  --seed-done-json "$FLUX_SEED_JSON" \
  --build-covs \
  --slim-skip "$SLIM_SKIP" \
  --progress-every 10 \
  --gc-every 10
echo "Flux DONE $(date -Is) used=$(used_gb)GiB"

echo "===== G4 (full, workers=$WORKERS) used=$(used_gb)GiB ====="
python "$REPO/analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum_parallel.py" \
  --input-glob "${OUT_ROOT}/dfs/hist_mc_g4/*.df" \
  --out-df "$G4_OUT" \
  --family G4 \
  --workers "$WORKERS" \
  --build-covs \
  --slim-skip "$SLIM_SKIP" \
  --progress-every 10 \
  --gc-every 10
echo "G4 DONE $(date -Is) used=$(used_gb)GiB"

echo "===== NPZs + plots $(date -Is) ====="
python "$REPO/analysis_village/numucc_1p0pi/scripts/histcounts_covs_to_syst_npzs.py" \
  --cov-pkl "${FLUX_OUT}.rate_covs.pkl" \
  --cov-pkl "${G4_OUT}.rate_covs.pkl" \
  --syst-disk-root "$SYST_DISK"

kill "$WATCH_PID" 2>/dev/null || true
echo "===== ALL DONE $(date -Is) used=$(used_gb)GiB ====="
echo "Product A/B: $SYST_DISK/productA  $SYST_DISK/productB"
