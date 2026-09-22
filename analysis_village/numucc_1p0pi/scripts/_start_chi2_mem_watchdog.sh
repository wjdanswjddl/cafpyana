#!/usr/bin/env bash
# Start 80 GiB absolute RSS watchdog for the chi2 track-subset campaign.
set -euo pipefail
REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
LOGDIR=/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-chi2subset-watchdog
mkdir -p "$LOGDIR"
LOG="$LOGDIR/mem_watchdog.log"
NOHUP="$LOGDIR/mem_watchdog.nohup.out"
PY="$REPO/envs/venv_py310_cafpyana/bin/python"
MATCH='syst_histcounts_from_df|syst_histcounts_stream_sum|syst_cosmics_(?:chunk|aggregate)|get_systematics_genie|syst_genie_parallel|run_syst_genie_chunked|run_chi2_subset_detector|run_chi2_track_subset|chi2_track_subset|_launch_cosmics'

pkill -f "$REPO/analysis_village/numucc_1p0pi/scripts/mem_watchdog.py" 2>/dev/null || true
sleep 1

nohup "$PY" "$REPO/analysis_village/numucc_1p0pi/scripts/mem_watchdog.py" \
  --threshold-gb 80 \
  --interval 10 \
  --match "$MATCH" \
  --log "$LOG" \
  > "$NOHUP" 2>&1 &
echo $! | tee "$LOGDIR/mem_watchdog.pid"
echo "started pid=$(cat "$LOGDIR/mem_watchdog.pid") log=$LOG"
sleep 12
head -n 20 "$LOG"
ps -p "$(cat "$LOGDIR/mem_watchdog.pid")" -o pid,etime,rss,cmd
