# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *
from makedf.geniesyst import gen1_systematics

DFS = [make_pandora_evtdf_mup_genieslimwgts, make_hdrdf, make_metadf, make_mcnudf_genieslimwgts]
ARGS = [{}, {}, {}, {"genie_systematics": gen1_systematics}]
NAMES = ["evt", "hdr", "meta", "mcnu"]
