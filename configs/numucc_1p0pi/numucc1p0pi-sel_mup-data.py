# selected event rates for data files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup, make_hdrdf, make_triggerdf, make_potdf_bnb]
ARGS = [{}, {}, {}, {}]
NAMES = ["evt", "hdr", "trigger", "bnbpot"]

