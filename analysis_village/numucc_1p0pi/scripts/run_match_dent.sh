#!/bin/bash
# Match events common to CV and DENT; write *_matched.df in each input directory.
#
# Primary comparison starts at sel_all (early selection impact).  A sel_mup pass
# is included for final-selected variable comparisons.
#
# Note: production DENT dirs may use the same timestamp suffix as CV
# (2026_08_18_125707 / 130158).  Override with --variation if needed.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DFS="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

# sel_all — compare from the beginning of the selection chain
python dent_match_common_events.py \
  --format sel_all \
  --variation cv   "${DFS}/2026_08_19_031254__sel_all-mc-CV" \
  --variation dent "${DFS}/2026_08_18_120607__sel_all-mc-DENT" \
  --filename-str sel_all \
  --summary-csv "${DFS}/dent_matched_summary-sel_all.csv"

# sel_mup — final-selected variables (optional second pass)
python dent_match_common_events.py \
  --format sel_mup \
  --variation cv   "${DFS}/2026_08_19_031423__sel_mup-mc-CV" \
  --variation dent "${DFS}/2026_08_18_120753__sel_mup-mc-DENT" \
  --filename-str sel_mup \
  --summary-csv "${DFS}/dent_matched_summary-sel_mup.csv"
