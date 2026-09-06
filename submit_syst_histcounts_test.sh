#!/usr/bin/env bash
# Fast local-pool histcounts campaign for notebook testing (new sentinel + var_configs).
#
# Creates a fresh timestamped output directory under pnfs scratch, runs 1-file
# jobs with reduced universes. Does NOT use jobsub (grid clones GitHub and would
# miss uncommitted local fixes).
#
# Usage:
#   bash submit_syst_histcounts_test.sh
#   NUMUCC_SYST_HIST_ROOT=<campaign>/out  # point the notebook here when done
set -euo pipefail

REPO=/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
cd "$REPO"
# shellcheck disable=SC1091
source envs/venv_py310_cafpyana/bin/activate
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export CAFPYANA_WD="$REPO"
export BEARER_TOKEN_FILE="${BEARER_TOKEN_FILE:-/tmp/bt_u$(id -u)}"

STAMP=$(date +%Y_%m_%d_%H%M%S)
# Prefer pnfs scratch; fall back to repo tmp if pnfs is not writable from this host.
CAMPAIGN_ROOT_DEFAULT=/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp
if mkdir -p "$CAMPAIGN_ROOT_DEFAULT/.probe_write" 2>/dev/null; then
  rmdir "$CAMPAIGN_ROOT_DEFAULT/.probe_write" 2>/dev/null || true
  CAMPAIGN_ROOT=$CAMPAIGN_ROOT_DEFAULT
else
  CAMPAIGN_ROOT="$REPO/tmp"
  mkdir -p "$CAMPAIGN_ROOT"
  echo "[warn] pnfs scratch not writable; using $CAMPAIGN_ROOT"
fi
CAMPAIGN=$CAMPAIGN_ROOT/syst_histcounts_test_${STAMP}
LISTS=$CAMPAIGN/lists
OUT=$CAMPAIGN/out
LOG=$CAMPAIGN/logs
mkdir -p "$LISTS" "$OUT" "$LOG"

# Reuse the known-good 1-file smoke lists when present; else copy from prior smoke.
PRIOR_LISTS=/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp/syst_histcounts_smoke/lists
if [[ -d "$PRIOR_LISTS" ]]; then
  cp -a "$PRIOR_LISTS"/. "$LISTS"/
else
  echo "ERROR: missing prior smoke lists at $PRIOR_LISTS" >&2
  exit 1
fi
# Prefer pnfs-local Spring CAF (exp/data path breaks xrootd glob on some nodes).
SPRING_PNFS=/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp/syst_histcounts_smoke/spring_mc_one.flat.caf.root
if [[ -f "$SPRING_PNFS" ]]; then
  printf '%s\n' "$SPRING_PNFS" > "$LISTS/mc_bnb.list"
fi

# Freeze VariableConfig into the campaign root (also written into each .df).
python - <<PY
from analysis_village.numucc_1p0pi.syst_histcounts import write_var_config_snapshot_json
write_var_config_snapshot_json("$CAMPAIGN/variable_configs.json")
write_var_config_snapshot_json("$OUT/variable_configs.json")
print("wrote variable_configs.json")
PY

CFG=configs/numucc_1p0pi/syst_histcounts.py
export SYST_HIST_GENIE_NUNIV=10
export SYST_HIST_FLUX_NUNIV=20
export SYST_HIST_G4_NUNIV=20
# Ar23p is a normal GENIE group (include in genie / all modes).
export SYST_HIST_EXCLUDE_AR23P=0

echo "CAMPAIGN=$CAMPAIGN"
echo "OUT=$OUT"
echo "reduced NUNIV: GENIE=$SYST_HIST_GENIE_NUNIV FLUX=$SYST_HIST_FLUX_NUNIV G4=$SYST_HIST_G4_NUNIV"

run_one() {
  local name="$1" mode="$2" sample="$3" list="$4"
  local extra="${5:-}"
  echo "[$(date +%H:%M:%S)] START $name mode=$mode sample=$sample"
  # shellcheck disable=SC2086
  env SYST_HIST_MODE="$mode" SYST_HIST_SAMPLE="$sample" $extra \
    python run_df_maker.py -c "$CFG" -l "$list" -o "$OUT/$name" -ncpu 1 -nfile 1 \
    >"$LOG/${name}.log" 2>&1
  local ec=$?
  echo "[$(date +%H:%M:%S)] DONE  $name exit=$ec"
  return $ec
}

# Unisim / CV (fast)
run_one hist_mc_cv_nowgt nowgt mc "$LISTS/mc_bnb.list" &
run_one hist_wiremod_sv nowgt mc "$LISTS/wiremod_sv.list" &
run_one hist_wiremod_xtxw nowgt mc "$LISTS/wiremod_xtxw.list" &
run_one hist_dent_cv nowgt mc "$LISTS/dent_cv.list" &
run_one hist_dent_var nowgt mc "$LISTS/dent_var.list" &
run_one hist_intime nowgt intime "$LISTS/intime.list" &
run_one hist_offbeam nowgt offbeam "$LISTS/offbeam.list" &
wait

# GENIE: all knobs incl. Ar23p + slim in one pass (do not split by mode).
unset GENIE_KNOB_GROUP || true
run_one hist_mc_genie genie mc "$LISTS/mc_bnb.list" &
run_one hist_mc_flux flux mc "$LISTS/mc_bnb.list" &
run_one hist_mc_g4 g4 mc "$LISTS/mc_bnb.list" &

# Dirt: full GENIE+Flux+G4 on 1 file is OK with reduced univ
run_one hist_dirt_all all dirt "$LISTS/dirt.list" &
wait

echo
echo "=== campaign complete ==="
echo "CAMPAIGN=$CAMPAIGN"
echo "Point notebook:  export NUMUCC_SYST_HIST_ROOT=$OUT"
ls -lh "$OUT"/*.df 2>/dev/null || true
# quick integrity: var_configs key present?
python - <<PY
import glob, pandas as pd
out="$OUT"
for p in sorted(glob.glob(out+"/*.df")):
    with pd.HDFStore(p,"r") as s:
        keys=[k.lstrip("/") for k in s.keys()]
        has_vc=any(k.startswith("var_configs") for k in keys)
        has_sh=any(k.startswith("syst_hists") for k in keys)
        print(f"{p.split('/')[-1]:28s} syst_hists={has_sh} var_configs={has_vc}")
PY
