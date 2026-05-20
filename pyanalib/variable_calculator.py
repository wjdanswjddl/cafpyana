import numpy as np
import pandas as pd
from makedf.util import *
from makedf.constants import *

from pyanalib.pandas_helpers import multicol_add, pad_column_name

_CC1PI0PI_TKI_VAR_NAMES = (
    "del_alpha",
    "del_phi",
    "del_Tp",
    "del_p",
    "del_Tp_x",
    "del_Tp_y",
)


def get_cc1p0pi_tki(mudf, pdf, P_mu_col, P_p_col):
    """
    Calculate TKI variables for numu CC 1p0pi selected events

    Inputs:
    - mudf : pandas df with muon information. Must contain P_mu_col and dir
    - pdf : pandas df with leading-proton information. Must contain P_p_col and dir
    - P_mu_col : tuple of str -- column name in mudf that holds the absolute muon momentum.
    - P_p_col : tuple of str -- column name in pdf that holds the absolute proton momentum.

    Returns: 
    - A dictionary with one entry per TKI observable:
       - del_alpha: angle between transverse momentum of muon and transverse momentum imbalance
       - del_phi:   angle between transverse momentum of muon and transverse momentum of proton
       - del_Tp:    magnitude of the transverse momentum imbalance
       - del_Tp_x:  (p̂_ν × p̂_T^μ) · δp⃗_T  with p̂_ν along +z
       - del_Tp_y:  -p̂_T^μ · δp⃗_T
       - del_p:     magnitude of the 3D imbalance

    Notes
    -----
    - The masses and nuclear constants (MUON_MASS, PROTON_MASS, MASS_A, MASS_Ap) are imported from makedf.constants

    """

    mu_p = mudf[P_mu_col]
    mu_p_x = mu_p * mudf["dir"]["x"]
    mu_p_y = mu_p * mudf["dir"]["y"]
    mu_p_z = mu_p * mudf["dir"]["z"]
    mu_phi_x = mu_p_x/mag2d(mu_p_x, mu_p_y)
    mu_phi_y = mu_p_y/mag2d(mu_p_x, mu_p_y)

    p_p = pdf[P_p_col]
    p_p_x = p_p * pdf["dir"]["x"]
    p_p_y = p_p * pdf["dir"]["y"]
    p_p_z = p_p * pdf["dir"]["z"]
    p_phi_x = p_p_x/mag2d(p_p_x, p_p_y)
    p_phi_y = p_p_y/mag2d(p_p_x, p_p_y)

    mu_Tp_x = mudf["dir"]["x"] * mu_p
    mu_Tp_y = mudf["dir"]["y"] * mu_p
    mu_Tp_z = mudf["dir"]["z"] * mu_p
    mu_Tp = mag2d(mu_Tp_x, mu_Tp_y)

    p_Tp_x = pdf["dir"]["x"] * p_p
    p_Tp_y = pdf["dir"]["y"] * p_p
    p_Tp_z = pdf["dir"]["z"] * p_p
    p_Tp = mag2d(p_Tp_x, p_Tp_y)

    _del_Tp_x = mu_Tp_x + p_Tp_x
    _del_Tp_y = mu_Tp_y + p_Tp_y
    del_Tp = mag2d(_del_Tp_x, _del_Tp_y)

    del_alpha = np.arccos(-(mu_Tp_x*_del_Tp_x + mu_Tp_y*_del_Tp_y)/(mu_Tp*del_Tp))
    del_phi = np.arccos(-(mu_Tp_x*p_Tp_x + mu_Tp_y*p_Tp_y)/(mu_Tp*p_Tp))

    mu_E = mag2d(mu_p, MUON_MASS)
    p_E = mag2d(p_p, PROTON_MASS)

    # R = MASS_A + mu_p_z + p_p_z - mu_E - p_E
    # del_Lp = 0.5*R - mag2d(MASS_Ap, del_Tp)**2/(2*R)
    e_cal = mu_E - MUON_MASS + p_E - PROTON_MASS + 0.0309 # https://link.springer.com/article/10.1140/epjc/s10052-019-6750-3
    del_Lp = mu_p_z + p_p_z - e_cal
    del_p = mag2d(del_Tp, del_Lp)

    # δp_{T,x} = (p̂_ν × p̂_T^μ) · δp⃗_T,  δp_{T,y} = -p̂_T^μ · δp⃗_T  (p̂_ν = +z)
    del_Tp_x = -mu_phi_y * _del_Tp_x + mu_phi_x * _del_Tp_y
    del_Tp_y = -(mu_phi_x * _del_Tp_x + mu_phi_y * _del_Tp_y)

    return {
        "del_alpha": del_alpha * 180/np.pi,
        "del_phi": del_phi * 180/np.pi,
        "del_Tp": del_Tp,
        "del_Tp_x": del_Tp_x,
        "del_Tp_y": del_Tp_y,
        "del_p": del_p,
    }


def add_reco_cc1p0pi_tki_evtdf(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Attach reco CC1pi TKI columns ``del_*`` onto ``evtdf`` via :func:`get_cc1p0pi_tki`.

    Uses ``mu.pfp.trk`` / ``p.pfp.trk`` range momentum ``P`` / ``p_muon`` and ``p_proton``.
    Skips if ``del_Tp`` already appears at column level 0.
    """
    if evtdf is None or len(evtdf) == 0:
        return evtdf
    try:
        if "del_Tp" in evtdf.columns.get_level_values(0):
            return evtdf
    except Exception:
        pass
    slc_mudf = evtdf.mu.pfp.trk
    slc_pdf = evtdf.p.pfp.trk
    slc_P_mu_col = pad_column_name(("P", "p_muon"), slc_mudf)
    slc_P_p_col = pad_column_name(("P", "p_proton"), slc_pdf)
    tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
    for var_name in _CC1PI0PI_TKI_VAR_NAMES:
        evtdf = multicol_add(evtdf, tki_reco[var_name].rename(var_name))
    return evtdf


def add_mc_cc1p0pi_tki_mcnu(mc_nu_df: pd.DataFrame) -> pd.DataFrame:
    """Attach MC-truth CC1pi TKI ``del_*`` onto ``mcnu`` via :func:`get_cc1p0pi_tki`.

    Resolves muon / proton blocks as:

    * ``mc.mu`` / ``mc.p`` when truth columns live under a leading ``mc`` level (e.g. GENIE
      chunk map after ``mcnu`` column prefixing).
    * Top-level ``mu`` / ``p`` when those exist **and** the ``mc`` slice has no ``mu``
      (mixed layout: GENIE weights under ``mc``, truth lepton / proton still at top level).
    """
    if mc_nu_df is None or len(mc_nu_df) == 0:
        return mc_nu_df
    try:
        if "del_Tp" in mc_nu_df.columns.get_level_values(0):
            return mc_nu_df
    except Exception:
        pass
    if not isinstance(mc_nu_df.columns, pd.MultiIndex):
        return mc_nu_df

    levels0 = set(mc_nu_df.columns.get_level_values(0))
    mc_mudf = None
    mc_pdf = None
    # Prefer truth nested under mc when that subtree actually contains mu / p.
    if "mc" in levels0:
        mc_blk = mc_nu_df["mc"]
        sub0 = set(mc_blk.columns.get_level_values(0))
        if "mu" in sub0 and "p" in sub0:
            mc_mudf = mc_blk["mu"]
            mc_pdf = mc_blk["p"]
    # Mixed mcnu: only GENIE (or other) blocks use ``mc``; lepton / proton stay at top level.
    if mc_mudf is None and "mu" in levels0 and "p" in levels0:
        mc_mudf = mc_nu_df["mu"]
        mc_pdf = mc_nu_df["p"]
    if mc_mudf is None or mc_pdf is None:
        return mc_nu_df

    try:
        mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
        mc_P_p_col = pad_column_name(("totp",), mc_pdf)
        tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return mc_nu_df
    for var_name in _CC1PI0PI_TKI_VAR_NAMES:
        mc_nu_df = multicol_add(mc_nu_df, tki_mc[var_name].rename("{}".format(var_name)))
    return mc_nu_df


def add_truth_cc1p0pi_tki_evtdf(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Truth-level TKI on ``evtdf`` under ``('mc', '<var>', '', ...)`` via :func:`get_cc1p0pi_tki`."""
    if evtdf is None or len(evtdf) == 0:
        return evtdf
    nlevel = evtdf.columns.nlevels
    key_mc_delTp = ("mc", "del_Tp") + ("",) * (nlevel - 2)
    try:
        if key_mc_delTp in evtdf.columns:
            return evtdf
    except Exception:
        pass
    try:
        slc_mudf = evtdf.mu.pfp.trk.truth.p
        slc_pdf = evtdf.p.pfp.trk.truth.p
    except (KeyError, AttributeError, TypeError):
        return evtdf
    slc_P_mu_col = pad_column_name(("totp",), slc_mudf)
    slc_P_p_col = pad_column_name(("totp",), slc_pdf)
    tki_truth = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
    for var_name in _CC1PI0PI_TKI_VAR_NAMES:
        tname = ("mc", var_name) + ("",) * (nlevel - 2)
        evtdf = multicol_add(evtdf, tki_truth[var_name].rename(tname))
    return evtdf