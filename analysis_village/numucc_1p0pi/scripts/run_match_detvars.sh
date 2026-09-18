#!/bin/bash
# Prefer notebooks/systematics-detector-match.ipynb (sel_all match before selection).
#
# Legacy: this used to match WireMod at sel_mup. Matching must be done at sel_all.
# Point --variation dirs at sel_all WireMod / calovar productions when available.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DFS="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

# YZ/XTXW: updatecalo. CV: plain nominal (match + POT reference only).
python dent_match_common_events.py \
  --format sel_all \
  --variation yz   "${DFS}/${WIREMOD_YZ_SEL_ALL:-2026_09_14_025216__sel_all-mc-BNB_cosmics-WireModYZ}" \
  --variation xtxw "${DFS}/${WIREMOD_XTXW_SEL_ALL:-2026_09_14_024629__sel_all-mc-BNB_cosmics-WireModXTXW}" \
  --variation cv   "${DFS}/${WIREMOD_CV_SEL_ALL:-2026_09_04_172912__sel_all-mc-CV}" \
  --filename-str sel_all \
  --summary-csv "${DFS}/wiremod_matched_summary-sel_all.csv"
