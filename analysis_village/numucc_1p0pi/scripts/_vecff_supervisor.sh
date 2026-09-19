#!/bin/bash
# Wait for the VecFF grid jobs (cluster 72035188), report stragglers, then run the
# GENIE syst chain via launch_sel_mup_VecFF.sh.
#
# Completion is decided from the .df count alone, because the bearer token that
# jobsub_q needs expires in a few hours and this waits longer than that; jobsub_q
# is only used for log lines when it happens to still work.
#
# Usage:
#   setsid nohup bash analysis_village/numucc_1p0pi/scripts/_vecff_supervisor.sh &
#
set -u
REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
CLUSTER=${CLUSTER:-72035188}
DFDIR=${DFDIR:-/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_18_184018__sel_mup-wgts_genie_VecFF}
NJOBS=${NJOBS:-500}
WORK=/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_mup_VecFF
SUP=$WORK/vecff_supervisor.log
# A run that loses a few jobs is still worth aggregating; bail out only if it lost a lot.
MIN_FRAC=${MIN_FRAC:-90}
POLL=${POLL:-300}
# Treat the run as finished once no new .df has appeared for this long.
STALL_POLLS=${STALL_POLLS:-9}
MAX_POLLS=${MAX_POLLS:-288}

export PATH=/opt/jobsub_lite/bin:$PATH
mkdir -p "$WORK"
cd "$REPO"
echo "vecff supervisor start $(date -Is) cluster=$CLUSTER pid=$$" >> "$SUP"

n_df() { ls "$DFDIR"/*.df 2>/dev/null | wc -l; }
# Queued count, or empty when jobsub_q cannot be read (expired token, timeout).
n_queued() {
    local out
    out=$(timeout 120 jobsub_q --group sbnd --user munjung 2>/dev/null) || return 1
    echo "$out" | grep -c "^${CLUSTER}\."
}

last=-1
stall=0
for poll in $(seq 1 "$MAX_POLLS"); do
    d=$(n_df)
    q=$(n_queued || echo "?")
    if [ "$d" != "$last" ]; then
        echo "$(date -Is) dfs=$d/$NJOBS queued=${q:-?}" >> "$SUP"
        last=$d
        stall=0
    else
        stall=$((stall + 1))
    fi

    if [ "$d" -ge "$NJOBS" ]; then
        echo "$(date -Is) all $NJOBS dfs present" >> "$SUP"
        break
    fi
    if [ "$stall" -ge "$STALL_POLLS" ]; then
        echo "$(date -Is) no new dfs for $((stall * POLL / 60))min; treating run as finished (dfs=$d/$NJOBS)" >> "$SUP"
        break
    fi
    sleep "$POLL"
done

# Straggler report: job N is missing when its .df never landed.
d=$(n_df)
pct=$((100 * d / NJOBS))
echo "$(date -Is) final dfs=$d/$NJOBS (${pct}%)" >> "$SUP"
missing=""
for i in $(seq 0 $((NJOBS - 1))); do
    [ -f "$DFDIR/sel_mup-wgts_genie_VecFF_${i}.df" ] || missing="$missing $i"
done
[ -n "$missing" ] && echo "$(date -Is) missing job ids:$missing" >> "$SUP"

if [ "$pct" -lt "$MIN_FRAC" ]; then
    echo "$(date -Is) ABORT: only ${pct}% of dfs present (< ${MIN_FRAC}%); not aggregating" >> "$SUP"
    exit 1
fi

# Phases 1-3: chunk-map, per-group merge, syst-disk aggregate.
echo "$(date -Is) launching syst chain" >> "$SUP"
bash "$REPO/analysis_village/numucc_1p0pi/scripts/launch_sel_mup_VecFF.sh" >> "$SUP" 2>&1
rc=$?
echo "$(date -Is) syst chain rc=$rc" >> "$SUP"
ls -la /exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_mup_VecFF/GENIE/ >> "$SUP" 2>&1
echo "$(date -Is) supervisor done" >> "$SUP"
