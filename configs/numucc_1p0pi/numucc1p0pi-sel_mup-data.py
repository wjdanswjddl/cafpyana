# selected event rates for data files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup, make_hdrdf, make_potdf_bnb, make_triggerdf]
ARGS = [{}, {}, {}, {}]
NAMES = ["evt", "hdr", "bnbpot", "trigger"]
