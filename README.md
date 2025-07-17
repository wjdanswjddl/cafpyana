# cafpyana

A set of scripts for analyzing SBN CAF files with the python.
Goal is to provide an easy starting point for SBN analysis to everyone.
For more details and instructions, please check [wiki](https://github.com/sungbinoh/cafpyana/wiki).

## Version compatibility

Current `main` branch is based on python v3.9.15.
It is for running the repository without any issue at gpvm servers with the `spack`.
There is no need to open an SL7 image.
If python version is updated in gpvm servers, compatibility issue should be revisited.

# for the dune plotting style, follow the instructions here
https://github.com/DUNE/dune_plot_style.git

# to run the df maker
python run_df_maker.py -c configs/pandora.py -l data/sample_lists/sbnd/prod_2025B/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/test_data_beamon_2025B -nfile 100