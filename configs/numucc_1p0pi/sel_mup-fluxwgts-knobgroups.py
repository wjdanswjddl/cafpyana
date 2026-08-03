# BNB flux multisim: all knobs in one pass (makedf.bnbsyst.regen_systematics).
# HDF keys evt, hdr — same layout as sel_mup-g4wgts.py.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS, ARGS, NAMES = build_flux_knobgroup_config()
