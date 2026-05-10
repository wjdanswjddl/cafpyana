OUT_CHUNKS=/exp/sbnd/data/users/munjung/xsec/event_selection_chunks   

# create if not exists: mkdir -p "$OUT_CHUNKS"
mkdir -p "$OUT_CHUNKS"
echo "Saving chunks to $OUT_CHUNKS"

echo "Running event selection chunk for mc"
python analysis_village/numucc_1p0pi/scripts/event_selection_chunk.py \
  --df_file /exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/BNB_cosmics/ab-sel_all-wgts.df \
  --sample mc \
  --out_dir "$OUT_CHUNKS"

echo "Running event selection chunk for data"
python analysis_village/numucc_1p0pi/scripts/event_selection_chunk.py \
  --df_file /exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/_Fixed_all.df \
  --sample data \
  --out_dir "$OUT_CHUNKS"

# repeat with --sample intime | offbeam | dirt and matching .df paths
