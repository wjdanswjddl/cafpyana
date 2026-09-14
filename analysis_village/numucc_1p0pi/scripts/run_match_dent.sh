#!/bin/bash
# Prefer notebooks/systematics-detector-match.ipynb (this script remains a CLI backend).
# Match events common to CV and DENT at sel_all ONLY (before any event selection).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DFS="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

python dent_match_common_events.py \
  --format sel_all \
  --variation cv   "${DFS}/2026_08_19_031254__sel_all-mc-CV" \
  --variation dent "${DFS}/2026_08_18_120607__sel_all-mc-DENT" \
  --filename-str sel_all \
  --summary-csv "${DFS}/dent_matched_summary-sel_all.csv"
