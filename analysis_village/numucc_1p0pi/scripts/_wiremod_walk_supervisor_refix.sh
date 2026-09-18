#!/bin/bash
set -e
CACHE="/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod/cache"
OUT="/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod"
PY="/exp/sbnd/app/users/munjung/env/bin/python"
N_SHARDS=6
BATCH=25
DROP_MAP=$CACHE/wiremod_xtxw_drop_map.pkl
SUP=$CACHE/wiremod_walk_supervisor.log
cd /exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
export BEARER_TOKEN_FILE=${BEARER_TOKEN_FILE:-/tmp/bt_u$(id -u)}
echo "supervisor_refix start $(date -Is)" > "$SUP"

wait_pids() {
  local pf=$1
  while true; do
    alive=0
    for pid in $(cat "$pf" 2>/dev/null); do
      kill -0 "$pid" 2>/dev/null && alive=1
    done
    [ "$alive" -eq 0 ] && return 0
    sleep 60
  done
}

echo "waiting YZ $(date -Is)" >> "$SUP"
wait_pids "$CACHE/wiremod_walk_yz_pids.txt"
echo "YZ done $(date -Is)" >> "$SUP"
for sid in $(seq 0 $((N_SHARDS-1))); do
  grep -q 'done pot=' "$CACHE/wiremod_walk_yz_s${sid}.log" || { echo "YZ $sid incomplete" >> "$SUP"; exit 1; }
done

echo "waiting XTXW drop-map $(date -Is)" >> "$SUP"
DPID=$(cat "$CACHE/wiremod_dedupe_xtxw.pid" 2>/dev/null || true)
if [ -n "$DPID" ]; then
  while kill -0 "$DPID" 2>/dev/null; do sleep 30; done
fi
# accept either unique_claimed= (scan done) or drop map file
for _i in $(seq 1 120); do
  if [ -f "$DROP_MAP" ] && grep -q 'unique_claimed=' "$CACHE/wiremod_dedupe_xtxw.log" 2>/dev/null; then
    break
  fi
  sleep 30
done
tail -8 "$CACHE/wiremod_dedupe_xtxw.log" >> "$SUP"
[ -f "$DROP_MAP" ] || { echo "drop map missing: $DROP_MAP" >> "$SUP"; exit 1; }
grep -q 'unique_claimed=' "$CACHE/wiremod_dedupe_xtxw.log" || { echo "dedupe incomplete" >> "$SUP"; exit 1; }
echo "dedupe/drop-map done $(date -Is)" >> "$SUP"

# refresh file list without importing syst_detvar_common (avoids matplotlib hang)
$PY - <<'PY'
import pickle
from pathlib import Path
out = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
xtxw = sorted(
    str(p)
    for p in (out / "matched" / "xtxw").glob("*_matched.df")
    if "sel_all" in p.name
)
with open(out / "cache" / "wiremod_xtxw_all_files.pkl", "wb") as fh:
    pickle.dump(xtxw, fh, protocol=pickle.HIGHEST_PROTOCOL)
print("xtxw files", len(xtxw))
PY

XTXW_PIDS=""
for sid in $(seq 0 $((N_SHARDS-1))); do
  LOG=$CACHE/wiremod_walk_xtxw_s${sid}.log
  CK=$CACHE/wiremod_walk_xtxw_s${sid}.pkl
  rm -f "$CK"
  : > "$LOG"
  nohup $PY -u analysis_village/numucc_1p0pi/scripts/wiremod_walk_shard.py \
    --files-pkl "$CACHE/wiremod_xtxw_all_files.pkl" \
    --checkpoint "$CK" --shard-id "$sid" --n-shards "$N_SHARDS" \
    --batch-size "$BATCH" --rss-limit-gb 20 --label XTXW \
    --drop-map-pkl "$DROP_MAP" >> "$LOG" 2>&1 &
  XTXW_PIDS="$XTXW_PIDS $!"
  echo "XTXW shard $sid pid=$!" >> "$SUP"
  disown $! 2>/dev/null || true
done
echo "$XTXW_PIDS" > "$CACHE/wiremod_walk_xtxw_pids.txt"
: > "$CACHE/wiremod_mem_monitor_walk_xtxw.log"
nohup bash "$CACHE/_mem_monitor_walk_shards.sh" "$CACHE/wiremod_walk_xtxw_pids.txt" "$CACHE/wiremod_mem_monitor_walk_xtxw.log" >/dev/null 2>&1 &

echo "waiting XTXW $(date -Is)" >> "$SUP"
wait_pids "$CACHE/wiremod_walk_xtxw_pids.txt"
echo "XTXW done $(date -Is)" >> "$SUP"
for sid in $(seq 0 $((N_SHARDS-1))); do
  grep -q 'done pot=' "$CACHE/wiremod_walk_xtxw_s${sid}.log" || { echo "XTXW $sid incomplete" >> "$SUP"; exit 1; }
done

YZ_CKPTS=""; XTXW_CKPTS=""; CV_CKPTS=""
for sid in $(seq 0 $((N_SHARDS-1))); do
  YZ_CKPTS="$YZ_CKPTS $CACHE/wiremod_walk_yz_s${sid}.pkl"
  XTXW_CKPTS="$XTXW_CKPTS $CACHE/wiremod_walk_xtxw_s${sid}.pkl"
  CV_CKPTS="$CV_CKPTS $CACHE/wiremod_walk_cv_s${sid}.pkl"
done
echo "merge $(date -Is)" >> "$SUP"
: > "$CACHE/wiremod_merge_walk.log"
# shellcheck disable=SC2086
$PY -u analysis_village/numucc_1p0pi/scripts/wiremod_merge_walk_shards.py \
  --out-base "$OUT" \
  --yz-shard-ckpts $YZ_CKPTS \
  --xtxw-shard-ckpts $XTXW_CKPTS \
  --cv-shard-ckpts $CV_CKPTS \
  --cv-campaign 2026_09_04_172912__sel_all-mc-CV >> "$CACHE/wiremod_merge_walk.log" 2>&1
echo "merge_rc=$? " >> "$SUP"
tail -8 "$CACHE/wiremod_merge_walk.log" >> "$SUP"

echo "plots $(date -Is)" >> "$SUP"
: > "$CACHE/wiremod_make_plots.log"
$PY -u analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py \
  --out-base "$OUT" --skip-inspect >> "$CACHE/wiremod_make_plots.log" 2>&1
echo "plots_rc=$? $(date -Is)" >> "$SUP"
tail -5 "$CACHE/wiremod_make_plots.log" >> "$SUP"
