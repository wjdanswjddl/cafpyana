# GENIE weights at loose ``sel_all`` (evt / trk / mcnu / hdr).
#
# --- Job 1: slim GENIE (recommended first pass) ---
# Fold all multisim knobs into one ``mc.GENIE.univ_*`` column set (universe-wise
# product via ``makedf/getsyst.py``). Multisigma / morph knobs stay per-knob.
# Knob list: ``regen_systematics`` (Spring regen CAFs).
#
#   GENIE_KNOB_GROUP=slim python run_df_maker.py \
#     -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
#     -l /path/to/mc.list -o sel_all-wgts_genie_slim -ngrid 1000
#
# --- Job 2+: per-knob-group (for cut-stage syst breakdown) ---
# Knob lists: ``makedf.geniesyst.GENIE_KNOB_GROUPS``.
#
#   GENIE_KNOB_GROUP=CCQE python run_df_maker.py \
#     -c configs/numucc_1p0pi/sel_all-geniewgts-knobgroups.py \
#     -l /path/to/mc.list -o sel_all-wgts_genie_CCQE -ngrid 1000
import os

from analysis_village.numucc_1p0pi.makedf.makedf import *

_sel = os.environ.get("GENIE_KNOB_GROUP", "").strip()
if _sel.lower() == "slim":
    DFS, ARGS, NAMES = build_genie_slim_config_sel_all()
else:
    DFS, ARGS, NAMES = build_genie_knobgroup_config_sel_all(group_filter=_sel if _sel else None)
