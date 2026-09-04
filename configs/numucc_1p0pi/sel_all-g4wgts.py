# Geant4 multisim at loose ``sel_all`` (HDF keys evt, trk, hdr).
# Layout matches ``sel_all-mc`` + weights for ``syst_multisim_chunk.py --input-stage sel_all``.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_all_mc_multisim, make_trkdf, make_hdrdf]
ARGS = [
    dict(wgt_types=["g4"], slim=False, multisim_nuniv=1000),
    {},
    {},
]
NAMES = ["evt", "trk", "hdr"]
