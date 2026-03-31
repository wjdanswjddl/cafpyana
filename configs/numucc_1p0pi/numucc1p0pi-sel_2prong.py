# 2-prong selected events
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_2prong, make_hdrdf, make_metadf]
ARGS = [{}, {}, {}]
NAMES = ["evt", "hdr", "meta"]
