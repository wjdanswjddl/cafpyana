#!/usr/bin/env bash
# sel_all GENIE syst for the VecFF group (single knob: SBN_v1_multisigma_VecFFCCQEshape).
# Uses GENIE_GROUP_GLOBS_SEL_ALL["VecFF"] — the sel_all-wgts_genie_VecFF jobs from
# submit_mc_jobs_GENIE_vecff_sel_all.sh. Work dirs are parallel to the Ar23p sel_all
# ones and to the sel_mup VecFF Product B run, so nothing existing is touched.
#
# Fills cut-stage (Product A, rate) + final-stage vars. After aggregate, Product B
# xsec is pruned (same as syst_disk_sel_all_20260913_Ar23p) so the canonical pkl
# cannot be confused with the sel_mup Product B numbers.
#
# Usage:
#   nohup bash analysis_village/numucc_1p0pi/scripts/launch_sel_all_VecFF.sh > run_VecFF_sel_all.nohup.out 2>&1 &
#
# Optional env: WORKERS (default 16), MERGE_WORKERS, REPO, WORK_BASE, NUMUCC_SYST_DISK_ROOT,
#               GENIE_VAR_SAVE_NAMES (default: all vars; empty = no filter)
#
set -euo pipefail

REPO="${REPO:-/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana}"
source "$REPO/envs/venv_py310_cafpyana/bin/activate"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

export WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/genie_syst-chunked-sel_all_VecFF}"
export NUMUCC_SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/syst_disk_sel_all_VecFF}"
export MC_DF_STAGE=sel_all
export GENIE_RUN_GROUPS=VecFF
export GENIE_VAR_SAVE_NAMES="${GENIE_VAR_SAVE_NAMES:-}"
export GENIE_MULTISIGMA_DIVIDE_BY_CV="${GENIE_MULTISIGMA_DIVIDE_BY_CV:-1}"

WORKERS="${WORKERS:-16}"
MERGE_WORKERS="${MERGE_WORKERS:-1}"

mkdir -p "$WORK_BASE/chunks" "$WORK_BASE/merged" "$NUMUCC_SYST_DISK_ROOT"
cd "$WORK_BASE"
LOG="$WORK_BASE/run_VecFF_sel_all.log"

echo "===== LAUNCH $(date -Is) SEL_ALL VecFF vars=${GENIE_VAR_SAVE_NAMES:-<all>} workers=${WORKERS} =====" | tee -a "$LOG"
echo "[info] WORK_BASE=$WORK_BASE" | tee -a "$LOG"
echo "[info] NUMUCC_SYST_DISK_ROOT=$NUMUCC_SYST_DISK_ROOT" | tee -a "$LOG"
echo "[info] input from dataset_locations.GENIE_GROUP_GLOBS_SEL_ALL['VecFF']" | tee -a "$LOG"

bash "$REPO/analysis_village/numucc_1p0pi/scripts/run_syst_genie_chunked.sh" \
  -g VecFF -j "$WORKERS" --merge-workers "$MERGE_WORKERS" >>"$LOG" 2>&1

GDIR="$NUMUCC_SYST_DISK_ROOT/GENIE"
PKL="$GDIR/cov_mat_dict.pkl"
if [ -f "$PKL" ]; then
  echo "===== prune Product B xsec from sel_all product $(date -Is) =====" | tee -a "$LOG"
  python3 "$REPO/analysis_village/numucc_1p0pi/scripts/prune_sel_all_genie_xsec.py" \
    "$PKL" "$GDIR/cov_mat_dict_rate_only.pkl" \
    --manifest "$GDIR/genie_covariance_manifest.json" >>"$LOG" 2>&1
  mv -n "$PKL" "$GDIR/cov_mat_dict.pkl.orig_with_productB_xsec"
  mv -n "$GDIR/genie_covariance_manifest.json" \
    "$GDIR/genie_covariance_manifest.json.orig_with_productB_xsec"
  mv -n "$GDIR/cov_mat_dict_rate_only.pkl" "$PKL"
  mv -n "$GDIR/cov_mat_dict_rate_only_manifest.json" "$GDIR/genie_covariance_manifest.json"
  python3 -c "
import json, pathlib
p = pathlib.Path('$GDIR/genie_covariance_manifest.json')
m = json.loads(p.read_text())
m['output_pkl'] = '$PKL'
m['productB_xsec_original_pkl'] = '$GDIR/cov_mat_dict.pkl.orig_with_productB_xsec'
p.write_text(json.dumps(m, indent=2) + '\n')
"
fi

echo "===== DONE $(date -Is) -> $PKL =====" | tee -a "$LOG"
