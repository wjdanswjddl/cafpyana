# all nu & selected evts w/ GENIE weights, for getting GENIE cov matrices
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_wgts_MEC, make_mcnudf_MEC, make_hdrdf]
NAMES = ["evt", "mcnu", "hdr"]
