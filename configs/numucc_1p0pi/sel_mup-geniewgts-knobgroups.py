# GENIE multisim weights by knob group (lists in makedf.geniesyst.GENIE_KNOB_GROUPS).
#
# Default: all groups in one output (evt_<Group>, mcnu_<Group>, ..., hdr).
# Single group (lighter jobs / old HDF key names evt, mcnu, hdr):
#   GENIE_KNOB_GROUP=Ar23p python run_df_maker.py -c configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py ...
import os

from analysis_village.numucc_1p0pi.makedf.makedf import *

_sel = os.environ.get("GENIE_KNOB_GROUP", "").strip()
DFS, ARGS, NAMES = build_genie_knobgroup_config(group_filter=_sel if _sel else None)
