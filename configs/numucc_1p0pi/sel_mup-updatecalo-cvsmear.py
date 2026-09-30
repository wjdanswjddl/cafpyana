# Product B (sel_mup) dE/dx smear unisim: one evt table holding the CV-redo
# selection (dedx_var=0) and the smear13 / smear26 selections (dedx_var=13/26),
# all from the same trkdf in the same job.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_dedxsmear, make_hdrdf, make_mcnudf]
ARGS = [{}, {}, {}]
NAMES = ["evt", "hdr", "mcnu"]
