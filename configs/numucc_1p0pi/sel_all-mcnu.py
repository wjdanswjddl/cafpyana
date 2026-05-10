# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all, make_mcnudf, make_hdrdf]
NAMES = ["evt", "mcnu", "hdr"]
