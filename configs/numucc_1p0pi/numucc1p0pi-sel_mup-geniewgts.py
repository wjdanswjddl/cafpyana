# slimmed GENIE wgt for selected events and all MC nu
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_wgts_genie, make_mcnudf_wgts_genie, make_hdrdf]
NAMES = ["evt", "mcnu", "hdr"]
