# for full event selection studies for data
# save all slices, tracks
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf, make_trkdf, make_hdrdf, make_potdf_bnb]
NAMES = ["evt", "trk", "hdr", "bnbpot"]
