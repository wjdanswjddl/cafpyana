#!/bin/bash
# Prefer notebooks/systematics-detector-match.ipynb (sel_all match before selection).
#
# Legacy: this used to match WireMod at sel_mup. Matching must be done at sel_all.
# Point --variation dirs at sel_all WireMod / calovar productions when available.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DFS="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

# Example — replace SET_ME paths with real sel_all WireMod campaign dirs:
python dent_match_common_events.py \
  --format sel_all \
  --variation yz   "${DFS}/${WIREMOD_YZ_SEL_ALL:-SET_ME__sel_all-mc-WireModYZ}" \
  --variation xtxw "${DFS}/${WIREMOD_XTXW_SEL_ALL:-SET_ME__sel_all-mc-WireModXTXW}" \
  --variation cv   "${DFS}/${WIREMOD_CV_SEL_ALL:-SET_ME__sel_all-mc-calovar}" \
  --filename-str sel_all \
  --summary-csv "${DFS}/wiremod_matched_summary-sel_all.csv"
