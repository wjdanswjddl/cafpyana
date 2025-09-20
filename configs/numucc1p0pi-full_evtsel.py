# for unfolding
# save reco slices after event selection
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf_final, make_mcnudf, make_hdrdf]
NAMES = ["evt", "mcnu", "hdr"]
