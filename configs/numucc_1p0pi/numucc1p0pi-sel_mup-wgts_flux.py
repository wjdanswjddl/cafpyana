# selected event rates & slimmed wgts (GENIE, Flux, G4) for uncertainties on event rates for MC files
from analysis_village.numucc_1p0pi.makedf.makedf import *

DFS = [make_pandora_evtdf_mup_wgts_flux, make_hdrdf]
NAMES = ["evt", "hdr"]
