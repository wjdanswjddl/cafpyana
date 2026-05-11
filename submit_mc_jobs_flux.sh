# All three categories in one .df
python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py \
  -l "$mc_list" -o sel_mup-wgts_flux_all -ngrid 200

# One category (same as old thin configs)
FLUX_GROUP=beam python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py \
  -l "$mc_list" -o sel_mup-wgts_flux_beam -ngrid 200
