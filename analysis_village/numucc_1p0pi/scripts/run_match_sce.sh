#!/bin/bash
# Match events common to 0xSCE and 2xSCE; write *_matched.df under each merged_perTPC dir.
#
python sce_match_common_events.py \
  --variation 0xSCE /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_12_113644__sel_2prong-mc-BNB_cosmics-0xSCE/merged_perTPC \
  --variation 2xSCE /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_12_114000__sel_2prong-mc-BNB_cosmics-2xSCE/merged_perTPC \
  --variation cv /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_17_232429__sel_2prong-mc-BNB_cosmics-CV/merged_perTPC  \
  --filename-str sel_2prong \
  --summary-csv /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/sce_matched_summary.csv

python sce_match_common_events.py \
  --variation 0xSCE /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_17_172521__sel_mup-mc-BNB_cosmics-0xSCE/merged_perTPC \
  --variation 2xSCE /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_17_172736__sel_mup-mc-BNB_cosmics-2xSCE/merged_perTPC \
  --variation cv /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_17_184445__sel_mup-mc-BNB_cosmics-CV/merged_perTPC  \
  --filename-str sel_mup \
  --summary-csv /pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/sce_matched_summary-mup.csv
