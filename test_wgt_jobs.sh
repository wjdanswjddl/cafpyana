filename="/pnfs/sbnd/scratch/users/gputnam/Ar23+_iterB/SBNDSpringMC/27118490_9/out20.flat.caf.root"
python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --list-cases
python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --dry-run
# Pool-mode smoke on one CAF: includes mup-mcstatwgts (Poisson MC stat), mup-g4wgts, single flux/GENIE groups, …
python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --caf "$filename"
python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --caf "$filename" --full-multisim
