from . import getsyst
import pandas as pd

regen_systematics = [
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    'GENIEReWeight_SBN_v1_multisim_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisim_CoulombCCQE',
    'GENIEReWeight_SBN_v1_multisim_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisim_NormNCMEC',
    'GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1gamma',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1eta',
    'GENIEReWeight_SBN_v1_multisim_COHVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse',

    # CCQE
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
    "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',
    'GENIEReWeight_SBN_v1_multisigma_MFP_N',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_N',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_N',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_N',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_N',

    # Zexp
    'ZExpPCAWeighter_SBNnusyst_b1',
    'ZExpPCAWeighter_SBNnusyst_b2', 
    'ZExpPCAWeighter_SBNnusyst_b3',
    'ZExpPCAWeighter_SBNnusyst_b4'

   # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", 

    # MEC
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", 

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", 
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", 
]

regen_systematics_sbnd = [
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    'GENIEReWeight_SBN_v1_multisim_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisim_CoulombCCQE',
    'GENIEReWeight_SBN_v1_multisim_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisim_NormNCMEC',
    'GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1gamma',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1eta',
    'GENIEReWeight_SBN_v1_multisim_COHVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse',

    # CCQE
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
    "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',
    'GENIEReWeight_SBN_v1_multisigma_MFP_N',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_N',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_N',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_N',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_N',

    # Zexp
    'ZExpPCAWeighter_SBNnusyst_b1',
    'ZExpPCAWeighter_SBNnusyst_b2', 
    'ZExpPCAWeighter_SBNnusyst_b3',
    'ZExpPCAWeighter_SBNnusyst_b4'

   # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", 

    # MEC
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", 

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", 
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", 
]

regen_systematics_sbnd_multisims = [
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    'GENIEReWeight_SBN_v1_multisim_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisim_CoulombCCQE',
    'GENIEReWeight_SBN_v1_multisim_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisim_NormNCMEC',
    'GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1gamma',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1eta',
    'GENIEReWeight_SBN_v1_multisim_COHVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse',
    'GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse',
]

regen_systematics_sbnd_multisigma = [
    # CCQE
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
    "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',
    'GENIEReWeight_SBN_v1_multisigma_MFP_N',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_N',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_N',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_N',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_N',

    # Zexp
    'ZExpPCAWeighter_SBNnusyst_b1',
    'ZExpPCAWeighter_SBNnusyst_b2', 
    'ZExpPCAWeighter_SBNnusyst_b3',
    'ZExpPCAWeighter_SBNnusyst_b4'
]

regen_systematics_sbnd_morph = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", 

    # MEC
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", 

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", 
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", 
]


def geniesyst(f, nuind):
    return getsyst.getsyst(f, regen_systematics, nuind)

def geniesyst_icarus(f, nuind):
    return getsyst.getsyst(f, regen_systematics, nuind)

def geniesyst_sbnd(f, nuind, genie_multisim_nuniv=100):
    geniewgtdf = getsyst.getsyst(f, regen_systematics_sbnd, nuind)
    geniemswgtdf = getsyst.getsyst(f, regen_systematics_sbnd_multisims, nuind)

    genie_cols = pd.MultiIndex.from_product(
        [["GENIE"], [f"univ_{i}" for i in range(genie_multisim_nuniv)]],
    )
    genie_wgt = pd.DataFrame(
        1.0,
        index=geniewgtdf.index,
        columns=genie_cols,
    )
    for syst in regen_systematics_sbnd_multisims:
        genie_wgt *= geniemswgtdf[syst].to_numpy()
    geniewgtdf = pd.concat([geniewgtdf, genie_wgt], axis=1)

    return geniewgtdf
