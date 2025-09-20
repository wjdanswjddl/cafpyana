# for full event selection studies
# save all slices, tracks, and MC nu (no syst wgts)
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf, make_mcnudf, make_hdrdf, make_trkdf]
NAMES = ["evt", "mcnu", "hdr", "trk"]
