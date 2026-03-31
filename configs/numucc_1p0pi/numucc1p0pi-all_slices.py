# for unfolding
# save reco slices after event selection
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf, make_mcnudf, make_hdrdf]
NAMES = ["evt", "mcnu", "hdr"]
