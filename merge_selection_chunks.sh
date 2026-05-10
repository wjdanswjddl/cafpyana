OUT_CHUNKS=/path/to/your/chunks_dir
OUT_PLOTS=/path/to/your/plots_dir    # mkdir -p "$OUT_PLOTS"

python analysis_village/numucc_1p0pi/scripts/event_selection_aggregate.py \
  --in_dir "$OUT_CHUNKS" \
  --out_dir "$OUT_PLOTS"

#  --cosmic_estimate intime          # or offbeam
#  --hide_cosmic_model_unc           # turn off grey cosmic band
#  --f_offbeam_frac 0.08             # default
#  --data_pot 1.42e20                # optional legend override only
#  --skip_global_exposure            # only for old per-chunk-scaled pickles
