# BNB flux multisim by category (bundles in makedf.bnbsyst.BNB_FLUX_GROUPS).
#
# Default: all categories in one output (evt_beam, evt_hadron, evt_xsec, hdr).
# Single category (lighter jobs / keys evt, hdr):
#   FLUX_GROUP=beam python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py ...
# Valid FLUX_GROUP keys: beam, hadron, xsec
import os

from analysis_village.numucc_1p0pi.makedf.makedf import *

_sel = os.environ.get("FLUX_GROUP", "").strip()
DFS, ARGS, NAMES = build_flux_knobgroup_config(group_filter=_sel if _sel else None)
