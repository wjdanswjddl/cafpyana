# for full event selection studies for off-beam data
# save all slices, tracks (off-beam has no genie MC info)
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf, make_trkdf, make_hdrdf]
NAMES = ["evt", "trk", "hdr"]
