from . import getsyst

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

def bnbsyst(f, nuind):
    return getsyst.getsyst(f, regen_systematics, nuind)

