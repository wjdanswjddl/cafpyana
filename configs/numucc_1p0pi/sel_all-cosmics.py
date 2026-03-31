# selected event rates for MC files
# use for event selection comparison
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_cosmics, make_trkdf, make_mcnudf, make_hdrdf]
NAMES = ["evt", "trk", "mcnu", "hdr"]
