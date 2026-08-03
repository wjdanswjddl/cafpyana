# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_2prong, make_hdrdf, make_potdf_bnb]
ARGS = [{}, {}, {}]
NAMES = ["evt", "hdr", "bnbpot"]

