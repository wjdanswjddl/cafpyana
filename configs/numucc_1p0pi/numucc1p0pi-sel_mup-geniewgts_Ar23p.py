# all nu & selected evts w/ GENIE weights, for getting GENIE cov matrices
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_wgts_ar23p, make_mcnudf_ar23p, make_hdrdf]
ARGS = [{}, {}, {}]
NAMES = ["evt", "mcnu", "hdr"]
