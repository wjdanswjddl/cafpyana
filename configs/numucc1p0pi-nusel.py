# for full event selection studies
# save all slices, tracks, and MC nu (no syst wgts)
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf_nu, make_hdrdf, make_trkdf]
NAMES = ["evt", "hdr", "trk"]
