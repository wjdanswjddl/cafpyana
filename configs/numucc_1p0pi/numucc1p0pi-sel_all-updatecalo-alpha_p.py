# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all_updatecalo_alpha_p, make_trkdf_updatecalo_alpha_p, make_hdrdf]
NAMES = ["evt", "trk", "hdr"]
