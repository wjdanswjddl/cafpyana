# for full event selection studies
# save all slices, tracks, and MC nu (no syst wgts)
# save just evts, for investigating early evt sel stages with full stats
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf, make_hdrdf]
NAMES = ["evt", "hdr"]
