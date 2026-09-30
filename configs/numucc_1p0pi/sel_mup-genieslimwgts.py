# Product B GENIE weights. The systematic is GENIE_slim_v3 = base × FSI v3.
#
#   python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-genieslimwgts.py \
#     -l /path/to/ar23_xrootd.list -o sel_mup-wgts_genie -ngrid 2000
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS, ARGS, NAMES = build_genie_fsi_compare_config_sel_mup()
