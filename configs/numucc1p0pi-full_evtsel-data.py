# for unfolding
# save reco slices after event selection
# save all MC nu & syst wgts

from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf_final, make_trkdf, make_hdrdf, make_potdf_bnb]
NAMES = ["evt", "trk", "hdr", "bnbpot"]
