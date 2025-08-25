from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from analysis_village.gump.kinematics import *

def pass_slc_with_n_pfps(df, n = 2):
    group_levels = ['entry', 'rec.slc..index']

    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with at least 2 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def tki(dirs, range_P_muon, range_P_proton):
    if(range_P_muon.iloc[0] < range_P_muon.iloc[1]):
    mass_0 = PDG["muon"][2]
    mass_1 = PDG["proton"][2]
    p_0 = range_P_muon.iloc[0]
    p_1 = range_P_proton.iloc[1]
    # swap ordering if chi2 ratio bad


def measure_tki(group):
    dirs = group[('pfp','trk','dir')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_proton = group[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]
    return tki(dirs, range_P_muon, range_P_proton)

def InFV_nohiyz(data):
    xmin = 10.
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def InFV_nohiyz_trk(data):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def add_contained_col(df):
    contained = InFV_nohiyz_trk(df.pfp.trk.start) & InFV_nohiyz_trk(df.pfp.trk.end)
    df['contained'] = contained

def Signal(df): # definition
    is_fv = InFV_nohiyz(df.position)
    return is_fv

def make_pandora_no_cuts_df(f):
    pandora_df = make_pandora_df(f)

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    pandora_df = pass_slc_with_n_pfps(pandora_df)

    # ----- loose PID for candidates ----
    pandora_df[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = pandora_df.pfp.trk.chi2pid.I2.chi2_muon/pandora_df.pfp.trk.chi2pid.I2.chi2_proton
    mu_df = pandora_df[(pandora_df.pfp.trackScore > 0.5)].sort_values(by=("pfp", "trk", "chi2pid", "I2", "mu_over_p", ""), ascending=False).groupby(level=[0,1]).head(1)
    idx_mu = mu_df.index
   
    idx_pfps = pandora_df.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = pandora_df.loc[idx_not_mu]

    p_df = notmudf[(notmudf.pfp.trackScore > 0.5)].sort_values(by=("pfp", "trk", "chi2pid", "I2", "mu_over_p", ""), ascending=False).groupby(level=[0,1]).tail(1)
    idx_p = p_df.index

    # mu_df.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mu_df.columns])
    # p_df.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in p_df.columns])
    # comb = multicol_merge(mu_df.droplevel(-1), p_df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

    # mu_df = comb.mu
    # p_df = comb.p

    # note if there are any other track/s
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = pandora_df.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.pfp.trackScore < 0.5]
    other_shw_length = othershwdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    # longest other track
    othertrkdf = otherdf[otherdf.pfp.trackScore > 0.5]
    other_trk_length = othertrkdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")

    range_p_mu_series = mu_df.pfp.trk.rangeP.p_muon
    range_p_mu_series = range_p_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True) 

    if pandora_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry','rec.slc..index']
        )
        del_p_series = pd.Series(dtype='float', name='del_p', index=empty_index)
    else:
        del_p_series = pandora_df.groupby(['entry','rec.slc..index']).apply(measure_tki)

    #del_p, del_Tp, del_phi, del_alpha, mu_E, p_E = transverse_kinematics(mu_df.pfp.trk.rangeP.p_muon, mu_df.pfp.trk.dir, p_df.pfp.trk.rangeP.p_proton, p_df.pfp.trk.dir)
    #nu_E = neutrino_energy(mu_df.pfp.trk.rangeP.p_muon, mu_df.pfp.trk.dir, p_df.pfp.trk.rangeP.p_proton, p_df.pfp.trk.dir)

    #del_p = del_p.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #del_Tp = del_Tp.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #del_phi = del_phi.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #del_alpha = del_alpha.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #mu_E = mu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #p_E = p_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    #nu_E = nu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    ######## (9) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)


    ## (10) creat a slice-based reco df
    slcdf = pd.DataFrame({
        # 'nu_E': nu_E,
        # 'mu_E': mu_E,
        # 'p_E': p_E,
        # 'del_p': del_p,
        # 'del_Tp': del_Tp,
        # 'del_phi': del_phi,
        'del_p': del_p_series,
        'range_p_mu_series': range_p_mu_series,
        'tmatch_idx': tmatch_idx_series
    })

    return slcdf

def make_gump_nudf(f):

    nudf = make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)
    wgtdf = pd.concat([bnbsyst.bnbsyst(f, nudf.ind), geniesyst.geniesyst_sbnd(f, nudf.ind)], axis=1)

    is_fv = InFV_nohiyz_trk(nudf.position)
    is_signal = Signal(nudf)
    is_cc = nudf.iscc
    genie_mode = nudf.genie_mode
    w = nudf.w

    try :
        nuint_categ = pd.Series(8, index=nudf.index)
    except Exception as e:
        print(f"Error init nuint_categ")
        return

    nuint_categ[~is_fv] = -1  # Out of FV
    nuint_categ[is_fv & is_signal] = 1 # Signal
    nuint_categ[is_fv & ~is_cc & ~is_signal] = 0  # NC
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 3)] = 2  # Non-signal CCCOH
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 0)] = 3  # CCQE
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 10)] = 4  # 2p2h
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w < 1.4) | (genie_mode == 1))] = 5  # RES
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w > 2.0) | (genie_mode == 2))] = 6  # DIS

    nudf['nuint_categ'] = nuint_categ

    mu_p_series = magdf(nudf.mu.genp)
    #cpi_p_series = magdf(nudf.cpi.genp)

    mu_cos_theta_series = nudf.mu.genp.z / mu_p_series
    #cpi_cos_theta_series = nudf.cpi.genp.z / cpi_p_series

    this_nudf = pd.DataFrame({
        'true_p_mu': mu_p_series,
        #'true_p_pi': cpi_p_series,
        'true_cos_theta_mu': mu_cos_theta_series,
        #'true_cos_theta_pi': cpi_cos_theta_series,
        #'true_t': true_t_series,
        'nuint_categ': nuint_categ
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
