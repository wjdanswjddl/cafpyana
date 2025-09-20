# for unfolding for off-beam data
# save all slices, tracks (off-beam has no genie MC info)
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf, make_trkdf, make_hdrdf]
NAMES = ["evt", "trk", "hdr"]
