# MC statistical uncertainty: Poisson universes per neutrino row (computed in-memory;
# no CAF globalTree weights). HDF keys evt, hdr — same layout as sel_mup-g4wgts.py.
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_mc_multisim, make_hdrdf]
ARGS = [
    dict(wgt_types=["mcstat"], slim=False, multisim_nuniv=100),
    {},
]
NAMES = ["evt", "hdr"]
