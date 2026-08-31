import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
from makedf.util import *


# def InFV(data): # cm
#     xmin = -190.
#     ymin = -190.
#     zmin = 10.
#     xmax = 190.
#     ymax =  190.
#     zmax =  450.
#     return (np.abs(data.x) > 10) & (np.abs(data.x) < 190) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)


# # --- FV recommendation from TPC studies by Sungbin ---
# def InFV_nohiyz(data):
#     xmin = 10.
#     xmax = 190.
#     zmin = 10.
#     zmax = 450.
#     ymax_highz = 100.
#     pass_xz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
#     pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
#     return pass_xz & pass_y

# def InFV_nohiyz_trk(data):
#     xmax = 190.
#     zmin = 10.
#     zmax = 450.
#     ymax_highz = 100.
#     pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
#     pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
#     return pass_xz & pass_y


# === event breakdown ===
def IsNu(df):
    return (np.abs(df.pdg) == 14) | (np.abs(df.pdg) == 12)

def IsCosmic(df):
    return ~IsNu(df)

def IsNuOutFV(df):
    return IsNu(df) & ~InFV(df.position, det="SBND")

def IsNuInFV(df):
    return IsNu(df) & InFV(df.position, det="SBND")

def IsNuInFV_NuOther(df):
    return IsNuInFV(df) & (df.pdg != 14)

def IsNuInFV_NumuNC(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 0)

# ---- numu CC in FV, breakdown in topology
# TODO: upper limit on muon energy!
def Is_1p0pi(df):
    return (df.nmu_220MeVc == 1) & (df.np_300MeVc == 1) & (df.npi_70MeVc == 0) & (df.npi0 == 0) &\
        (np.sqrt(df.mu.genp.x**2 + df.mu.genp.y**2 + df.mu.genp.z**2) < 1) &\
        (np.sqrt(df.p.genp.x**2 + df.p.genp.y**2 + df.p.genp.z**2) < 1)

def Is_Np0pi(df):
    return (df.nmu_220MeVc == 1) & (df.np_300MeVc > 1) & (df.npi_70MeVc == 0) & (df.npi0 == 0) 

def IsNuInFV_NumuCC_Other(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              ~Is_1p0pi(df) & ~Is_Np0pi(df)

def IsNuInFV_NumuCC_Np0pi(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              Is_Np0pi(df)

def IsNuInFV_NumuCC_1p0pi(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              Is_1p0pi(df)

# --- numu CC in FV, breakdown in genie mode
def IsNuInFV_NumuCC_QE(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              (df.genie_mode == 0)

def IsNuInFV_NumuCC_MEC(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              (df.genie_mode == 10)

def IsNuInFV_NumuCC_RES(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              (df.genie_mode == 1)

def IsNuInFV_NumuCC_DIS(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              (df.genie_mode == 2)

# def IsNuInFV_NumuCC_COH(df):
#     return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
#               (df.genie_mode == 3)

def IsNuInFV_NumuCC_OtherMode(df):
    return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
              ~(df.genie_mode == 0) & ~(df.genie_mode == 1) & ~(df.genie_mode == 2) & ~(df.genie_mode == 10)


# --- for event topoloy breakdown ---
def IsNu(df):
    return ~df.pdg.isna()


def IsSignal(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position, det="SBND")
    is_1mu1p0pi = (df.nmu_220MeVc == 1) & (df.npi_70MeVc == 0) & (df.np_300MeVc == 1) & (df.npi0 == 0) & (df.mu.totp < 1) & (df.p.totp < 1) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNp0pi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position, det="SBND")
    is_1mu1p0pi = (df.nmu_220MeVc == 1) & (df.npi_70MeVc == 0) & (df.np_300MeVc > 1) & (df.npi0 == 0) & (df.mu.totp < 1) & (df.p.totp < 1) #& (df.mu.genE > 0.25) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNcpi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position, det="SBND")
    is_1mu1p0pi = (df.nmu_220MeVc == 1) & (df.npi_70MeVc > 0) & (df.npi0 == 0) #& (df.mu.genE > 0.25) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def get_int_category(df):
    # cut_notnu = ~IsNu(df)
    # cut_nu_outfv = IsNu(df) & ~InFV(df.position)
    # cut_signal = IsSignal(df)
    # cut_1muNp0pi = Is1muNp0pi(df)
    # cut_1muNcpi = Is1muNcpi(df)

    cut_cosmic = IsCosmic(df)
    cut_nu_outfv = IsNuOutFV(df)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(df)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(df)
    cut_nu_infv_numu_cc_other = IsNuInFV_NumuCC_Other(df)
    cut_nu_infv_numu_cc_np0pi = IsNuInFV_NumuCC_Np0pi(df)
    cut_nu_infv_numu_cc_1p0pi = IsNuInFV_NumuCC_1p0pi(df)

    # assert there's no overlap between the categories, AND that all categories are covered just in case i messed something up...
    assert (cut_cosmic & cut_nu_outfv & cut_nu_infv_nu_other & cut_nu_infv_numu_nc & cut_nu_infv_numu_cc_other & cut_nu_infv_numu_cc_np0pi & cut_nu_infv_numu_cc_1p0pi).sum() == 0
    assert (cut_cosmic | cut_nu_outfv | cut_nu_infv_nu_other | cut_nu_infv_numu_nc | cut_nu_infv_numu_cc_other | cut_nu_infv_numu_cc_np0pi | cut_nu_infv_numu_cc_1p0pi).sum() == len(df)

    # category 1 NEEDS TO BE THE SIGNAL MODE
    nuint_categ = pd.Series(10, index=df.index)
    nuint_categ[cut_cosmic] = -1  # not nu
    nuint_categ[cut_nu_outfv] = 0  # nu out of FV
    nuint_categ[cut_nu_infv_numu_cc_1p0pi] = 1    # nu in FV, signal
    nuint_categ[cut_nu_infv_numu_cc_np0pi] = 2  # 1mu, 0cpi, Np, 0pi0
    nuint_categ[cut_nu_infv_numu_cc_other] = 3  # 1mu, Ncpi, 0pi0
    nuint_categ[cut_nu_infv_numu_nc] = 4  # nu in FV, numu NC
    nuint_categ[cut_nu_infv_nu_other] = 5  # nu in FV, other

    return nuint_categ


def get_genie_category(df):
    cut_cosmic = IsCosmic(df)
    cut_nu_outfv = IsNuOutFV(df)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(df)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(df)
    # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
    cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode(df)
    cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS(df)
    cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES(df)
    cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC(df)
    cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE(df)

    assert (cut_cosmic & cut_nu_outfv & cut_nu_infv_nu_other & cut_nu_infv_numu_nc & cut_nu_infv_numu_othermode & cut_nu_infv_numu_cc_dis & cut_nu_infv_numu_cc_res & cut_nu_infv_numu_cc_me & cut_nu_infv_numu_cc_qe).sum() == 0
    assert (cut_cosmic | cut_nu_outfv | cut_nu_infv_nu_other | cut_nu_infv_numu_nc | cut_nu_infv_numu_othermode | cut_nu_infv_numu_cc_dis | cut_nu_infv_numu_cc_res | cut_nu_infv_numu_cc_me | cut_nu_infv_numu_cc_qe).sum() == len(df)

    genie_categ = pd.Series(10, index=df.index)
    genie_categ[cut_cosmic] = -1  # not nu
    genie_categ[cut_nu_outfv] = 0  # nu out of FV
    genie_categ[cut_nu_infv_numu_cc_qe] = 1  # nu in FV, QE
    genie_categ[cut_nu_infv_numu_cc_me] = 2  # nu in FV, MEC
    genie_categ[cut_nu_infv_numu_cc_res] = 3  # nu in FV, RES
    genie_categ[cut_nu_infv_numu_cc_dis] = 4  # nu in FV, DIS
    genie_categ[cut_nu_infv_numu_othermode] = 5  # nu in FV, other mode
    genie_categ[cut_nu_infv_numu_nc] = 6  # nu in FV, numu NC
    genie_categ[cut_nu_infv_nu_other] = 7  # nu in FV, other

    return genie_categ



# --- for plotting ---
# nu / cosmic breakdown
nu_cosmics_labels = ["Cosmic", r"Out-FV $\nu$", r"FV $\nu$"]
nu_cosmics_colors = ["gray", "C0", "C1"]

# signal / backgroundtopology breakdown
# the signal mode code MUST be the first item in the list for all the code below to work
topology_list = [1, 2, 3, 4, 5, 0, -1]
# mode_labels = [ r"FV other $\nu$", r"FV $\nu_{\mu}$ NC",
#                 r"FV $\nu_{\mu}$ CC Other", r"FV $\nu_{\mu}$ CC Np0$\pi$", r"FV $\nu_{\mu}$ CC 1p0$\pi$",
#                 r"Out-FV $\nu$", "Cosmic"]
# mode_colors = ["crimson", "darkgreen", 
#                 "coral", "darkslateblue", "mediumslateblue", "sienna", "gray"]
topology_labels = [r"FV $\nu_{\mu}$ CC 1p0$\pi$", r"FV $\nu_{\mu}$ CC Np0$\pi$", r"FV $\nu_{\mu}$ CC Other", r"FV $\nu$ NC", r"FV other $\nu$", r"Out-FV $\nu$", "Cosmic"]
topology_colors = ["mediumslateblue", "darkslateblue", "coral", "darkgreen", "crimson", "sienna", "gray"] 


# --- GENIE interaction mode breakdown ---
# genie_mode_list = [0, 10, 1, 2, 3]
genie_mode_list = [1, 2, 3, 4, 5, 6, 7, 0, -1]
genie_mode_labels = [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC MEC', r'$\nu_{\mu}$ CC RES', r'$\nu_{\mu}$ CC SIS/DIS', r'$\nu_{\mu}$ CC COH', 
                     r"$\nu$ NC", r"FV other $\nu$", r"Out-FV $\nu$", "Cosmic"]
genie_mode_colors = ["#9b5580", "#390C1E", "#2c7c94", "#D88A3B", "#BFB17C", 
                     "darkgreen", "crimson", "sienna","gray"] 