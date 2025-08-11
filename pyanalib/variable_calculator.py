import numpy as np
import pandas as pd
from makedf.util import *
from makedf.constants import *

def get_tki(mudf, pdf, P_mu_col, P_p_col):

    # Caculate transverse kinematics
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

    mu_Tp_x = mu_phi_y*mu_p_x - mu_phi_x*mu_p_y
    mu_Tp_y = mu_phi_x*mu_p_x - mu_phi_y*mu_p_y
    mu_Tp = mag2d(mu_Tp_x, mu_Tp_y)

    p_Tp_x = mu_phi_y*p_p_x - mu_phi_x*p_p_y
    p_Tp_y = mu_phi_x*p_p_x - mu_phi_y*p_p_y
    p_Tp = mag2d(p_Tp_x, p_Tp_y)

    del_Tp_x = mu_Tp_x + p_Tp_x
    del_Tp_y = mu_Tp_y + p_Tp_y
    del_Tp = mag2d(del_Tp_x, del_Tp_y)

    del_alpha = np.arccos(-(mu_Tp_x*del_Tp_x + mu_Tp_y*del_Tp_y)/(mu_Tp*del_Tp))
    del_phi = np.arccos(-(mu_Tp_x*p_Tp_x + mu_Tp_y*p_Tp_y)/(mu_Tp*p_Tp))

    mu_E = mag2d(mu_p, MUON_MASS)
    p_E = mag2d(p_p, PROTON_MASS)

    R = MASS_A + mu_p_z + p_p_z - mu_E - p_E
    del_Lp = 0.5*R - mag2d(MASS_Ap, del_Tp)**2/(2*R)
    del_p = mag2d(del_Tp, del_Lp)

    return {
        "del_alpha": del_alpha,
        "del_phi": del_phi,
        "del_Tp": del_Tp,
        "del_p": del_p,
    }