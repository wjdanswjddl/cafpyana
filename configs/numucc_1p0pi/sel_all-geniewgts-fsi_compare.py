# Product A (sel_all): one CAF pass with BASE ∪ FSI_v1_N ∪ FSI_v3_N.
# Writes GENIE_base, FSI_v1_N, FSI_v3_N, GENIE_slim_{v1,v3,both} + atomic FSI ±σ.
#
#   python run_df_maker.py -c configs/numucc_1p0pi/sel_all-geniewgts-fsi_compare.py \
#     -l /path/to/ar23_xrootd.list -o sel_all-wgts_genie_FSI_compare -ngrid 2000
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS, ARGS, NAMES = build_genie_fsi_compare_config_sel_all()
