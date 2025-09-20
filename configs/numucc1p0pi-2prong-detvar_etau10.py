# save events after 2-prong slice selection (before quality cuts on the prongs)
# the two tracks will be saved under "t1" and "t2"
from analysis_village.numucc1p0pi.makedf.makedf import *

DFS = [make_numucc1p0pi_evtdf_2prong_etau10, make_hdrdf]
NAMES = ["evt", "hdr"]
