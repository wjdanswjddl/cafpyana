#!/bin/bash
# Prefer notebooks/systematics-detector-match.ipynb (this script remains a CLI backend).
# SCE matching must use sel_all productions (same as WireMod / DENT).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DFS="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

python dent_match_common_events.py \
  --format sel_all \
  --variation 0xSCE "${DFS}/${SCE_0X_SEL_ALL:-SET_ME__sel_all-mc-BNB_cosmics-0xSCE}" \
  --variation 2xSCE "${DFS}/${SCE_2X_SEL_ALL:-SET_ME__sel_all-mc-BNB_cosmics-2xSCE}" \
  --variation cv    "${DFS}/${SCE_CV_SEL_ALL:-SET_ME__sel_all-mc-BNB_cosmics-CV}" \
  --filename-str sel_all \
  --summary-csv "${DFS}/sce_matched_summary-sel_all.csv"
