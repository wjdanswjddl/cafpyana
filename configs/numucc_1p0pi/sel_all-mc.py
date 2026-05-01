# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all, make_trkdf, make_trkhitdf_plane0, make_trkhitdf_plane1, make_trkhitdf_plane2, make_hdrdf]
ARGS = [{}, {}, {}, {}, {}, {}, {}, {}]
NAMES = ["evt", "trk", "hit0", "hit1", "hit2", "hdr"]
