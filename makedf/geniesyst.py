from . import getsyst
import pandas as pd

# Regen systematic variations
regen_systematics = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    "GENIEReWeight_SBN_v1_multisigma_RPA_CCQE",
    "GENIEReWeight_SBN_v1_multisigma_CoulombCCQE",
    # "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",

    # MEC
    "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",

    # RES
    # "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1eta",
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",

    # Non-Res
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi",
    "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi",

    # DIS
    # "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    # "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",

    # NCEL
    # "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
]

regen_systematics_multisims = [
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
]

def geniesyst(f, nuind, genie_multisim_nuniv=100):
    geniewgtdf = getsyst.getsyst(f, regen_systematics, nuind)
    geniemswgtdf = getsyst.getsyst(f, regen_systematics_multisims, nuind)

    genie_cols = pd.MultiIndex.from_product(
        [["GENIE"], [f"univ_{i}" for i in range(genie_multisim_nuniv)]],
    )
    genie_wgt = pd.DataFrame(
        1.0,
        index=geniewgtdf.index,
        columns=genie_cols,
    )
    for syst in regen_systematics_multisims:
        genie_wgt *= geniemswgtdf[syst].to_numpy()
    geniewgtdf = pd.concat([geniewgtdf, genie_wgt], axis=1)

    return geniewgtdf
