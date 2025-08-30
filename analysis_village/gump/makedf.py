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

def make_pandora_gump_df(f):

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    pandora_df = make_pandora_df(f)
    stub_df = make_stubs(f, det=DETECTOR)

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
    pandora_df = pandora_df[pid_cut(pandora_df.mu.trk.chi2pid.I2.chi2_muon, 
                                    pandora_df.mu.trk.chi2pid.I2.chi2_proton, 
                                    pandora_df.p.trk.chi2pid.I2.chi2_muon, 
                                    pandora_df.p.trk.chi2pid.I2.chi2_proton,
                                    pandora_df.mu.trk.len)]
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

    slcdf = make_slcdf(f)
    trkdf = make_trkdf(f, False)
    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    trkdf = trkdf[trkdf.pfp.dist_to_vertex < 10]
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.pfp.trk.chi2pid.I2.chi2_muon/trkdf.pfp.trk.chi2pid.I2.chi2_proton
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    # track containment
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = (InFV(trkdf.pfp.trk.start, 0, det=DETECTOR)) & (InFV(trkdf.pfp.trk.end, 0, det=DETECTOR))

    # reco momentum -- range for contained, MCS for exiting
    trkdf[("pfp", "trk", "P", "p_muon", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_muon", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_muon", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_muon","", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_muon", "", "")]

    trkdf[("pfp", "trk", "P", "p_pion", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_pion", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_pion", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_pion", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_pion", "", "")]

    trkdf[("pfp", "trk", "P", "p_proton", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_proton", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_proton", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_proton", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_proton", "", "")]

    # opening angles
    trkdf[("pfp", "trk", "cos", "x", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "x", "", "")] = (trkdf.pfp.trk.end.x-trkdf.pfp.trk.start.x)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "cos", "y", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "y", "", "")] = (trkdf.pfp.trk.end.y-trkdf.pfp.trk.start.y)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "cos", "z", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "z", "", "")] = (trkdf.pfp.trk.end.z-trkdf.pfp.trk.start.z)/trkdf.pfp.trk.len
    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = trkdf[(trkdf.pfp.trackScore> 0.5)].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid","I2","mu_over_p", "")]).groupby(level=[0, 1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    pd.set_option('display.max_rows', None)
    idx_mu = mudf.index

    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = trkdf.pfp.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = trkdf.loc[idx_not_mu]
    pdf = notmudf[(notmudf.pfp.trackScore > 0.5)].sort_values(notmudf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    slcdf = multicol_merge(slcdf, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_p = pdf.index

    # note if there are any other track/showers
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = trkdf.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.pfp.trackScore < 0.5]
    other_shw_length = othershwdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    slcdf = multicol_add(slcdf, other_shw_length)
    # longest other track
    othertrkdf = otherdf[otherdf.pfp.trackScore > 0.5]
    other_trk_length = othertrkdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")
    slcdf = multicol_add(slcdf, other_trk_length)

    if slcdf.empty:
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
        has_stub = pd.Series(dtype='float', name='has_stub', index=empty_index)
    else:
        pd.set_option('display.max_rows', None)
        tki = transverse_kinematics(slcdf.mu.pfp.trk.P.p_muon, slcdf.mu.pfp.trk.cos, slcdf.p.pfp.trk.P.p_proton, slcdf.p.pfp.trk.cos)
        nu_E = neutrino_energy(slcdf.mu.pfp.trk.P.p_muon, slcdf.mu.pfp.trk.cos, slcdf.p.pfp.trk.P.p_proton, slcdf.p.pfp.trk.cos)
        del_p = tki['del_p']
        del_Tp = tki['del_Tp']
        del_phi = tki['del_phi']
        del_alpha = tki['del_alpha']
        mu_E = tki['mu_E']
        p_E = tki['p_E']

    ######## (9) - c: slc.tmatch.idx for truth matching
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan
    tmatch_idx_series = slcdf.slc.tmatch.idx
    slc_vtx = slcdf.slc.vertex
    nu_score = slcdf.slc.nu_score
    other_shw_length = slcdf.other_shw_length
    other_trk_length = slcdf.other_trk_length
    mu_chi2_of_mu_cand = slcdf.mu.pfp.trk.chi2pid.I2.chi2_muon
    prot_chi2_of_mu_cand = slcdf.mu.pfp.trk.chi2pid.I2.chi2_proton
    mu_chi2_of_prot_cand = slcdf.p.pfp.trk.chi2pid.I2.chi2_muon
    prot_chi2_of_prot_cand = slcdf.p.pfp.trk.chi2pid.I2.chi2_proton
    mu_len = slcdf.mu.pfp.trk.len

    stubdf = make_stubs(f)
    has_any_stub_series = stubdf.groupby('rec.slc..index')['pass_proton_stub'].transform('any')
    slc_has_stub_series = pd.Series(index=mu_len.index, dtype=bool)

    for k in slc_has_stub_series.keys():
        try:
            threek = k + (0,)
            slc_has_stub_series[k] = has_any_stub_series[threek]
        except KeyError:
            slc_has_stub_series[k] = False

    ## (10) create a slice-based reco df
    slcdf = pd.DataFrame({
        'other_shw_length': other_shw_length,
        'other_trk_length': other_trk_length,
        'slc_vtx_x': slc_vtx.x,
        'slc_vtx_y': slc_vtx.y,
        'slc_vtx_z': slc_vtx.z,
        'nu_score': nu_score,
        'mu_chi2_of_mu_cand': mu_chi2_of_mu_cand,
        'mu_chi2_of_prot_cand': mu_chi2_of_prot_cand,
        'prot_chi2_of_mu_cand': prot_chi2_of_mu_cand,
        'prot_chi2_of_prot_cand': prot_chi2_of_prot_cand,
        'mu_len': mu_len,
        'nu_E': nu_E,
        'mu_E': mu_E,
        'p_E': p_E,
        'del_p': del_p,
        'del_Tp': del_Tp,
        'del_phi': del_phi,
        'tmatch_idx': tmatch_idx_series,
        'has_stub': slc_has_stub_series
    })

    return slcdf

def make_gump_nudf(f):
    nudf = make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)
    wgtdf = pd.concat([bnbsyst.bnbsyst(f, nudf.ind), geniesyst.geniesyst_sbnd(f, nudf.ind)], axis=1)
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    is_fv = fv_cut(nudf.position, DETECTOR)
    is_cc = nudf.iscc
    is_nc = (nudf.iscc == 0)
    is_cosmic = (nudf.pdg == -1)
    genie_mode = nudf.genie_mode
    pdg = nudf.pdg
    is_1p0pi = (nudf.nmu_27MeV == 1) & (nudf.np_50MeV == 1) & (nudf.npi_30MeV == 0) & (nudf.npi0 == 0)
    is_numu = (nudf.pdg == 14)
    is_other_numucc = (is_numu & is_cc & (is_1p0pi == 0) & is_fv)
    is_sig = is_fv & is_1p0pi & is_numu & is_cc
    w = nudf.w

    nudf['nuint_categ'] = genie_mode 

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
        'genie_mode': genie_mode, 
        'is_sig': is_sig, 
        'pdg': pdg 
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
