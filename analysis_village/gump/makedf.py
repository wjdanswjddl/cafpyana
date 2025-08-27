from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from analysis_village.gump.kinematics import *
from analysis_village.gump.gump_cuts import *

def pass_slc_with_n_pfps(df, n = 2):
    group_levels = ['entry', 'rec.slc..index']

    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with at least 2 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def measure_nu_E(group):
    dirs = group[('pfp','trk','dir')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_proton = group[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]
    mu_over_p = group[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]
    
    ret = transverse_kinematics(range_P_muon.iloc[0], dirs.iloc[0], range_P_proton.iloc[1], dirs.iloc[1])
    if(mu_over_p.iloc[0] < mu_over_p.iloc[1]):
        return neutrino_energy(range_P_muon.iloc[0], dirs.iloc[0], range_P_proton.iloc[1], dirs.iloc[1])
    else:
        return neutrino_energy(range_P_muon.iloc[1], dirs.iloc[1], range_P_proton.iloc[0], dirs.iloc[0])

def measure_tki(group):
    dirs = group[('pfp','trk','dir')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    #dir_muon = 
    range_P_proton = group[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]
    #dir_proton = 
    mu_over_p = group[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]
    
    if(mu_over_p.iloc[0] < mu_over_p.iloc[1]):
        return transverse_kinematics(range_P_muon.iloc[0], dirs.iloc[0], range_P_proton.iloc[1], dirs.iloc[1])
    else:
        return transverse_kinematics(range_P_muon.iloc[1], dirs.iloc[1], range_P_proton.iloc[0], dirs.iloc[0])


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

def make_pandora_gump_df(f):

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    pandora_df = make_pandora_df(f)
    stub_df = make_stubs(f)

    pandora_df = multicol_merge(pandora_df, stub_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    pandora_df[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = pandora_df.pfp.trk.chi2pid.I2.chi2_muon/pandora_df.pfp.trk.chi2pid.I2.chi2_proton
    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = pandora_df.pfp[(pandora_df.pfp.trackScore> 0.5)].sort_values([("trk", "chi2pid","I2","mu_over_p", "")]).groupby(level=[0, 1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    pandora_df = multicol_merge(pandora_df, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index
        
    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = pandora_df.pfp.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = pandora_df.pfp.loc[idx_not_mu]
    pdf = notmudf[(notmudf.trackScore > 0.5)].sort_values([("trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    pandora_df = multicol_merge(pandora_df, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_p = pdf.index
    
    # note if there are any other track/showers
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = pandora_df.pfp.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.trackScore < 0.5]
    other_shw_length = othershwdf.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    pandora_df = multicol_add(pandora_df, other_shw_length)
    # longest other track
    othertrkdf = otherdf[otherdf.trackScore > 0.5]
    other_trk_length = othertrkdf.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")
    pandora_df = multicol_add(pandora_df, other_trk_length)

    # FV cut
    pandora_df = pandora_df[fv_cut(pandora_df.slc.vertex, DETECTOR)]
    # Cosmic cut
    pandora_df.slc = pandora_df.slc[cosmic_cut(pandora_df.slc)]
    # Two prong cut
    pandora_df = pandora_df[twoprong_cut(pandora_df)] 
    # PID cut
    pandora_df = pandora_df[pid_cut(pandora_df)]
    # Stub cut
    pandora_df = pandora_df[stub_cut(pandora_df)]


    if pandora_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry','rec.slc..index']
        )
        del_p = pd.Series(dtype='float', name='del_p', index=empty_index)
        del_Tp = pd.Series(dtype='float', name='del_Tp', index=empty_index)
        del_phi = pd.Series(dtype='float', name='del_phi', index=empty_index)
        del_alpha = pd.Series(dtype='float', name='del_alpha', index=empty_index)
        mu_E = pd.Series(dtype='float', name='mu_E', index=empty_index)
        p_E = pd.Series(dtype='float', name='p_E', index=empty_index)
        nu_E = pd.Series(dtype='float', name='nu_E', index=empty_index)
    else:
        tki = transverse_kinematics(pandora_df.mu.trk.rangeP.p_muon, pandora_df.mu.trk.dir, pandora_df.p.trk.rangeP.p_proton, pandora_df.p.trk.dir)
        nu_E = neutrino_energy(pandora_df.mu.trk.rangeP.p_muon, pandora_df.mu.trk.dir, pandora_df.p.trk.rangeP.p_proton, pandora_df.p.trk.dir)
        del_p = tki['del_p']
        del_Tp = tki['del_Tp']
        del_phi = tki['del_phi']
        del_alpha = tki['del_alpha']
        mu_E = tki['mu_E']
        p_E = tki['p_E']

        nu_E = nu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        nu_E = nu_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        mu_E = mu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        mu_E = mu_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        p_E = p_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        p_E = p_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_p = del_p.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_p = del_p.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_Tp = del_Tp.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_Tp = del_Tp.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_phi = del_phi.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_phi = del_phi.reset_index(level='rec.slc.reco.stub..index', drop=True)


    ######## (9) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx

    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.stub..index', drop=True)


    ## (10) create a slice-based reco df
    slcdf = pd.DataFrame({
        'nu_E': nu_E,
        'mu_E': mu_E,
        'p_E': p_E,
        'del_p': del_p,
        'del_Tp': del_Tp,
        'del_phi': del_phi,
        'tmatch_idx': tmatch_idx_series
    })

    return slcdf

def make_pandora_no_cuts_df(f):

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    pandora_df = make_pandora_df(f)
    stub_df = make_stubs(f)

    pandora_df = multicol_merge(pandora_df, stub_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    pandora_df[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = pandora_df.pfp.trk.chi2pid.I2.chi2_muon/pandora_df.pfp.trk.chi2pid.I2.chi2_proton
    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = pandora_df.pfp[(pandora_df.pfp.trackScore> 0.5)].sort_values([("trk", "chi2pid","I2","mu_over_p", "")]).groupby(level=[0, 1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    pandora_df = multicol_merge(pandora_df, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index
        
    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = pandora_df.pfp.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = pandora_df.pfp.loc[idx_not_mu]
    pdf = notmudf[(notmudf.trackScore > 0.5)].sort_values([("trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    pandora_df = multicol_merge(pandora_df, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_p = pdf.index
    
    # note if there are any other track/showers
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = pandora_df.pfp.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.trackScore < 0.5]
    other_shw_length = othershwdf.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    pandora_df = multicol_add(pandora_df, other_shw_length)
    # longest other track
    othertrkdf = otherdf[otherdf.trackScore > 0.5]
    other_trk_length = othertrkdf.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")
    pandora_df = multicol_add(pandora_df, other_trk_length)

    if pandora_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry','rec.slc..index']
        )
        del_p = pd.Series(dtype='float', name='del_p', index=empty_index)
        del_Tp = pd.Series(dtype='float', name='del_Tp', index=empty_index)
        del_phi = pd.Series(dtype='float', name='del_phi', index=empty_index)
        del_alpha = pd.Series(dtype='float', name='del_alpha', index=empty_index)
        mu_E = pd.Series(dtype='float', name='mu_E', index=empty_index)
        p_E = pd.Series(dtype='float', name='p_E', index=empty_index)
        nu_E = pd.Series(dtype='float', name='nu_E', index=empty_index)
    else:
        tki = transverse_kinematics(pandora_df.mu.trk.rangeP.p_muon, pandora_df.mu.trk.dir, pandora_df.p.trk.rangeP.p_proton, pandora_df.p.trk.dir)
        nu_E = neutrino_energy(pandora_df.mu.trk.rangeP.p_muon, pandora_df.mu.trk.dir, pandora_df.p.trk.rangeP.p_proton, pandora_df.p.trk.dir)
        del_p = tki['del_p']
        del_Tp = tki['del_Tp']
        del_phi = tki['del_phi']
        del_alpha = tki['del_alpha']
        mu_E = tki['mu_E']
        p_E = tki['p_E']

        nu_E = nu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        nu_E = nu_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        mu_E = mu_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        mu_E = mu_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        p_E = p_E.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        p_E = p_E.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_p = del_p.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_p = del_p.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_Tp = del_Tp.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_Tp = del_Tp.reset_index(level='rec.slc.reco.stub..index', drop=True)

        del_phi = del_phi.reset_index(level='rec.slc.reco.pfp..index', drop=True)
        del_phi = del_phi.reset_index(level='rec.slc.reco.stub..index', drop=True)


    ######## (9) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx

    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.stub..index', drop=True)


    ## (10) create a slice-based reco df
    slcdf = pd.DataFrame({
        'nu_E': nu_E,
        'mu_E': mu_E,
        'p_E': p_E,
        'del_p': del_p,
        'del_Tp': del_Tp,
        'del_phi': del_phi,
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

    muon_p_series = magdf(nudf.mu.genp)
    proton_p_series = magdf(nudf.p.genp)
    
    true_tki = transverse_kinematics(muon_p_series, nudf.mu.genp, proton_p_series, nudf.mu.genp)
    true_del_p = true_tki['del_p']
    true_del_Tp = true_tki['del_Tp']
    true_del_phi = true_tki['del_phi']
    true_del_alpha = true_tki['del_alpha']
    true_mu_E = true_tki['mu_E']
    true_p_E = true_tki['p_E']

    true_nu_E = neutrino_energy(muon_p_series, nudf.mu.genp, proton_p_series, nudf.mu.genp)

    this_nudf = pd.DataFrame({
        'true_nu_E': true_nu_E,
        'true_mu_E': true_mu_E,
        'true_p_E': true_p_E,
        'true_del_p': true_del_p,
        'true_del_Tp': true_del_Tp,
        'true_del_phi': true_del_phi,
        'true_del_p': true_del_p,
        'nuint_categ': nuint_categ
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
