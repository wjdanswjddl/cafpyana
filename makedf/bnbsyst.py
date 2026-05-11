from . import getsyst
import pandas as pd

# Regen systematic variations
regen_systematics = [
    'expskin_Flux',
    'horncurrent_Flux',
    'kminus_Flux',
    'kplus_Flux',
    'kzero_Flux',
    'nucleoninexsec_Flux',
    'nucleonqexsec_Flux',
    'nucleontotxsec_Flux',
    'piminus_Flux',
    'pioninexsec_Flux',
    'pionqexsec_Flux',
    'piontotxsec_Flux',
    'piplus_Flux'
]

bnb_systematics_beam = [
    'expskin_Flux',
    'horncurrent_Flux',
]

bnb_systematics_hadron = [
    'kminus_Flux',
    'kplus_Flux',
    'kzero_Flux',
    'piminus_Flux',
    'piplus_Flux'
]

bnb_systematics_xsec = [
    'pioninexsec_Flux',
    'pionqexsec_Flux',
    'piontotxsec_Flux',
    'nucleoninexsec_Flux',
    'nucleonqexsec_Flux',
    'nucleontotxsec_Flux',
]

# Flux multisim bundles for df configs (see analysis_village ... build_flux_knobgroup_config).
# Values are (systematics list, multisim_nuniv) matching sel_mup-wgts_flux_* defaults.
BNB_FLUX_GROUPS = {
    "beam": (bnb_systematics_beam, 200),
    "hadron": (bnb_systematics_hadron, 200),
    "xsec": (bnb_systematics_xsec, 1000),
}


def bnbsyst(f, nuind, multisim_nuniv=1000, slim=False, systematics=None):
    if systematics is None:
        systematics = regen_systematics

    bnbwgtdf = getsyst.getsyst(f, systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim, slimname="Flux")

    if slim:  # keep only the multiplied "Flux.univ_" columns
        flux_cols = [c for c in bnbwgtdf.columns if c[0] == "Flux"]
        bnbwgtdf = bnbwgtdf[flux_cols]
        
    return bnbwgtdf

