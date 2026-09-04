# BNB flux multisim at loose ``sel_all``: all knobs in one pass (makedf.bnbsyst.regen_systematics).
# HDF keys evt, trk, hdr — same layout as sel_all-g4wgts.py.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS, ARGS, NAMES = build_flux_knobgroup_config_sel_all()
