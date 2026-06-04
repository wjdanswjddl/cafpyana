# selected event rates for MC files
# use for detector variation samples
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_metadf, make_pandora_evtdf_mup_updateefield]
ARGS = [{}, {"updateefield": True}]
NAMES = ["meta", "evt"]
