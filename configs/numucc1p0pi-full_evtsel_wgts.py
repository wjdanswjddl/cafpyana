# for unfolding
# save reco slices after event selection
# save all MC nu & syst wgts
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf_final_wgts, make_hdrdf]
NAMES = ["evt", "hdr"]
