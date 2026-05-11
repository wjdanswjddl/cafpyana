# Geant4 multisim at mup selection (HDF keys evt, hdr).
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_mc_multisim, make_hdrdf]
ARGS = [
    dict(wgt_types=["g4"], slim=False, multisim_nuniv=1000),
    {},
]
NAMES = ["evt", "hdr"]
