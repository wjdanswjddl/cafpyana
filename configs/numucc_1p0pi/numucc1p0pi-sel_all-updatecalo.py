# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all_updatecalo, make_trkdf_updatecalo, make_metadf, make_hdrdf, make_mcnudf]
NAMES = ["evt", "trk", "meta", "hdr", "mcnu"]
