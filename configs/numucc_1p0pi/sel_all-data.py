# all events for data files
# use for event selection comparison
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all, make_trkdf, make_hdrdf, make_potdf_bnb]
NAMES = ["evt", "trk", "hdr", "bnbpot"]
