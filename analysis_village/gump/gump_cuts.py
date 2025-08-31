# Standard library imports
import os
import sys

# Third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

# Local imports
import analysis_village.gump.kinematics
from makedf.util import *


# Fiducial volume cuts for SBND and ICARUS
SBNDFVCuts = {
    "x": {"min": -199.15 + 10, "max": 199.15 - 10},
    "y": {"min": -200. + 10, "max": 200. - 10},
    "z": {"min": 0.0 + 10, "max": 500. - 50}
}

ICARUSFVCuts = {
    "C0": {
        "x": {"min": -358.49 + 10, "max": -61.94 - 10},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    },
    "C1": {
        "x": {"min": 358.49, "max": 61.94},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    }
}

def fv_cut(df, det):
    if det == "ICARUS":
        return (((df.x < ICARUSFVCuts['C0']['x']['max']) & (df.x > ICARUSFVCuts['C0']['x']['min'])) |\
                ((df.x < ICARUSFVCuts['C1']['x']['max']) & (df.x > ICARUSFVCuts['C1']['x']['min']))) &\
                 (df.y < ICARUSFVCuts['C0']['y']['max']) & (df.y > ICARUSFVCuts['C0']['y']['min']) &\
                 (df.z < ICARUSFVCuts['C0']['z']['max']) & (df.z > ICARUSFVCuts['C0']['z']['min'])

    elif det == "SBND":
        return (df.x < SBNDFVCuts['x']['max']) & (df.x > SBNDFVCuts['x']['min']) &\
               (df.y < SBNDFVCuts['y']['max']) & (df.y > SBNDFVCuts['y']['min']) &\
               (df.z < SBNDFVCuts['z']['max']) & (df.z > SBNDFVCuts['z']['min'])

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def cosmic_cut(df):
    return (df.nu_score > 0.5)

def twoprong_cut(df):
    return (np.isnan(df.other_shw_length) & np.isnan(df.other_trk_length))

def pid_cut(mu_chi2_mu_cand, mu_chi2_prot_cand, prot_chi2_mu_cand,
            prot_chi2_prot_cand, mu_len):

    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 25, 100, 50
    mu_cut = (mu_chi2_mu_cand < MUSEL_MUSCORE_TH) & \
             (prot_chi2_mu_cand > MUSEL_PSCORE_TH) & \
             (mu_len > MUSEL_LEN_TH)

    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut = (mu_chi2_prot_cand > PSEL_MUSCORE_TH) & \
            (prot_chi2_prot_cand < PSEL_PSCORE_TH)

    return mu_cut & p_cut

def stub_cut(df):
    cut = (df.has_stub == 0)
    return cut

mode_list = [0, 10, 1, 2, 3]
mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

def breakdown_mode(var, df):
    """Break down variable by interaction mode."""
    ret = [var[df.genie_mode == i] for i in mode_list]
    ret.append(var[sum([df.genie_mode == i for i in mode_list]) == 0])
    return ret

def tmatchdf(df, mcdf):
    # filter the columns from the mcdf we want
    # map old name to new name
    tosave = {
      "pdg": "pdg",
      "is_sig": "is_sig",
      "genie_mode": "genie_mode"
    }

    def savecol(c):
        return c[0] in list(tosave.keys()) 

    mcdf_cols = [c for c in mcdf.columns if savecol(c)]
    mcdf = mcdf[mcdf_cols]

    # Get the column depth matching on both sides
    mcdf.columns = pd.MultiIndex.from_tuples([(tosave[c[0]], c[1]) if c[0] in tosave else (c[0], c[1]) for c in mcdf.columns])
    df.columns = pd.MultiIndex.from_tuples([(c, "") for c in df.columns])

    tmatch_df = pd.merge(df, mcdf, how="left", left_on=["__ntuple", "entry", "tmatch_idx"], right_index=True) 

    # Fill nan's for not-matching columns
    for c in tosave.values():
        tmatch_df[(c, "")] = tmatch_df[(c, "")].fillna(0.)

    # Systematic weights are all 1 for no truth match
    tmatch_df.fillna(1, inplace=True)

    return tmatch_df
