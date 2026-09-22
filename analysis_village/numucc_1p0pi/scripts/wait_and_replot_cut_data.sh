#!/usr/bin/env bash
# Wait for cut-campaign data-1e20 resubmits to finish, then remake Product B + cut-var overlays.
set -euo pipefail

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
OUT=/exp/sbnd/data/users/munjung/FixedDev/selected_xsec_overlay_cuts
LOG="$OUT/wait_and_replot_data_fix.log"
EXPECT_N=150
# Require nearly all shards present and almost none empty.
MIN_GOOD=140
MAX_TINY_FRAC=0.05
POLL_SEC=300

DATA_DIRS=(
  /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/nu_score0/2026_09_20_022614__sel_mup-data-1e20
  /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/chi2mu15/2026_09_20_022713__sel_mup-data-1e20
  /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/chi2mu45/2026_09_20_022838__sel_mup-data-1e20
  /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/mcs_range_diff1p0/2026_09_20_023002__sel_mup-data-1e20
  /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/vz_exclude_200_300/2026_09_20_023058__sel_mup-data-1e20
)

cd "$REPO"
# shellcheck disable=SC1091
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

echo "[$(date -Is)] wait_and_replot start" | tee -a "$LOG"

ready=0
while [ "$ready" -eq 0 ]; do
  ready=1
  echo "[$(date -Is)] poll:" | tee -a "$LOG"
  for d in "${DATA_DIRS[@]}"; do
    tag=$(basename "$(dirname "$d")")
    n=$(ls "$d"/*.df 2>/dev/null | wc -l || true)
    tiny=$(find "$d" -maxdepth 1 -name '*.df' -size -2k 2>/dev/null | wc -l || true)
    n=${n// /}; tiny=${tiny// /}
    n=${n:-0}; tiny=${tiny:-0}
    good=$((n - tiny))
    assert_n=0
    if compgen -G "$d"/log_*.log >/dev/null 2>&1; then
      assert_n=$(rg -l 'Length of new_levels' "$d"/log_*.log 2>/dev/null | wc -l || true)
      assert_n=${assert_n// /}
      assert_n=${assert_n:-0}
    fi
    if [ "$n" -gt 0 ]; then
      frac_tiny=$(python3 -c "print($tiny/float($n))")
    else
      frac_tiny=1.0
    fi
    echo "  $tag n=$n good=$good tiny=$tiny tiny_frac=$frac_tiny assert_logs=$assert_n" | tee -a "$LOG"
    if ! python3 -c "import sys; sys.exit(0 if ($good >= $MIN_GOOD and float('$frac_tiny') <= $MAX_TINY_FRAC and $assert_n == 0) else 1)"; then
      ready=0
    fi
  done
  if [ "$ready" -eq 0 ]; then
    sleep "$POLL_SEC"
  fi
done

echo "[$(date -Is)] data ready — removing old overlay counts and replotting" | tee -a "$LOG"
for tag in nu_score0 chi2mu15 chi2mu45 mcs_range_diff1p0 vz_exclude_200_300; do
  rm -f "$OUT/$tag/overlay_histdata.pkl" "$OUT/$tag/cutvar_overlay_histdata.pkl"
done

echo "[$(date -Is)] Product B overlays..." | tee -a "$LOG"
python -u analysis_village/numucc_1p0pi/scripts/selected_xsec_overlay.py 2>&1 | tee -a "$LOG"

echo "[$(date -Is)] Cut-var overlays..." | tee -a "$LOG"
python -u analysis_village/numucc_1p0pi/scripts/selected_xsec_overlay_cut_vars.py 2>&1 | tee -a "$LOG"

# Restore FORCE flags to False so later accidental runs don't rebuild
python3 - <<'PY'
from pathlib import Path
for rel in (
    "analysis_village/numucc_1p0pi/scripts/selected_xsec_overlay.py",
    "analysis_village/numucc_1p0pi/scripts/selected_xsec_overlay_cut_vars.py",
):
    p = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana") / rel
    t = p.read_text()
    t2 = t.replace("FORCE_REBUILD_COUNTS = True", "FORCE_REBUILD_COUNTS = False", 1)
    if t2 != t:
        p.write_text(t2)
        print("reset", rel)
PY

echo "[$(date -Is)] DONE" | tee -a "$LOG"
