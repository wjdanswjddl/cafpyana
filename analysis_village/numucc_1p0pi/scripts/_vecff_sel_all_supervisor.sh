#!/bin/bash
# Wait for the VecFF sel_all grid jobs, report stragglers, then run the GENIE
# syst chain via launch_sel_all_VecFF.sh (Product A).
#
# Completion is decided from the .df count alone (bearer token for jobsub_q
# expires before this wait finishes). jobsub_q is only for optional log lines.
#
# Usage (after submit; set CLUSTER / DFDIR to the new stamp):
#   CLUSTER=... DFDIR=... NJOBS=1000 setsid nohup \
#     bash analysis_village/numucc_1p0pi/scripts/_vecff_sel_all_supervisor.sh &
#
set -u
REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
CLUSTER=${CLUSTER:?set CLUSTER to the jobsub cluster id}
DFDIR=${DFDIR:?set DFDIR to .../dfs/<stamp>__sel_all-wgts_genie_VecFF}
NJOBS=${NJOBS:-1000}
WORK=/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_all_VecFF
SUP=$WORK/vecff_sel_all_supervisor.log
MIN_FRAC=${MIN_FRAC:-90}
POLL=${POLL:-300}
STALL_POLLS=${STALL_POLLS:-9}
MAX_POLLS=${MAX_POLLS:-288}

export PATH=/opt/jobsub_lite/bin:$PATH
mkdir -p "$WORK"
cd "$REPO"
echo "vecff sel_all supervisor start $(date -Is) cluster=$CLUSTER pid=$$" >> "$SUP"

n_df() { ls "$DFDIR"/*.df 2>/dev/null | wc -l; }
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

d=$(n_df)
pct=$((100 * d / NJOBS))
echo "$(date -Is) final dfs=$d/$NJOBS (${pct}%)" >> "$SUP"
missing=""
for i in $(seq 0 $((NJOBS - 1))); do
    [ -f "$DFDIR/sel_all-wgts_genie_VecFF_${i}.df" ] || missing="$missing $i"
done
[ -n "$missing" ] && echo "$(date -Is) missing job ids:$missing" >> "$SUP"

if [ "$pct" -lt "$MIN_FRAC" ]; then
    echo "$(date -Is) ABORT: only ${pct}% of dfs present (< ${MIN_FRAC}%); not aggregating" >> "$SUP"
    exit 1
fi

echo "$(date -Is) launching syst chain" >> "$SUP"
bash "$REPO/analysis_village/numucc_1p0pi/scripts/launch_sel_all_VecFF.sh" >> "$SUP" 2>&1
rc=$?
echo "$(date -Is) syst chain rc=$rc" >> "$SUP"
ls -la /exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_all_VecFF/GENIE/ >> "$SUP" 2>&1
echo "$(date -Is) supervisor done" >> "$SUP"
