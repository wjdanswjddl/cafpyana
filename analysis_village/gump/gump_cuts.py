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

def pid_cut(df):
    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 25, 100, 50
    mu_cut = (df.mu.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & \
             (df.mu.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH) & \
             (df.mu.trk.len > MUSEL_LEN_TH)

    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut = (df.p.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & \
            (df.p.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH)

    return mu_cut & p_cut

def stub_cut(df):
    cut = (df.pass_proton_stub == 0)
    return cut

def is_cosmic(df):
    """Return mask for cosmic events."""
    return (df.slc.truth.pdg == -1)

def is_FV(df, det):
    """Return mask for events in fiducial volume."""
    return fv_cut(df.slc.vertex, det)

def is_numu(df):
    """Return mask for numu events."""
    return (np.abs(df.slc.truth.pdg) == 14)

def is_CC(df):
    """Return mask for CC events."""
    return (df.slc.truth.iscc == 1)

def is_NC(df):
    """Return mask for NC events."""
    return (df.slc.truth.iscc == 0)

def is_1p0pi(df):
    """Return mask for 1mu, 1p, 0pi events."""
    return (df.slc.truth.nmu_27MeV == 1) & (df.slc.truth.np_50MeV == 1) & (df.slc.truth.npi_30MeV == 0) & (df.slc.truth.npi0 == 0)

def is_signal(df, det):
    """Return mask for signal events."""
    return is_numu(df) & is_CC(df) & is_1p0pi(df) & is_FV(df, det)

def is_outFV(df, det):
    """Return mask for signal events outside FV."""
    return is_numu(df) & is_CC(df) & is_1p0pi(df) & np.invert(is_FV(df, det))

def is_othernumuCC(df, det):
    """Return mask for other numu CC events in FV."""
    return is_numu(df) & is_CC(df) & np.invert(is_1p0pi(df)) & is_FV(df, det)

mode_list = [0, 10, 1, 2, 3]
mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

def breakdown_mode(var, df):
    """Break down variable by interaction mode."""
    ret = [var[df.genie_mode == i] for i in mode_list]
    ret.append(var[sum([df.genie_mode == i for i in mode_list]) == 0])
    return ret

def breakdown_top(var, df, det):
    """Break down variable by topological category."""
    ret = [var[is_signal(df, det)],
            var[is_othernumuCC(df, det)],
            var[is_NC(df)],
            var[is_outFV(df, det)],
            var[is_cosmic(df)],
            var[np.invert(is_signal(df, det) | is_othernumuCC(df, det) | is_NC(df) | is_outFV(df, det) | is_cosmic(df))]
            ]
    return ret
