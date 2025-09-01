from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from analysis_village.gump.kinematics import *
from analysis_village.gump.gump_cuts import *

# to do: make_pandora_with_cuts using correct formatting 
# and can be turned into ttree for PROfit

def make_pandora_no_cuts_df(f):

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    slcdf = make_slcdf(f)
    StartingRows = len(slcdf)

    trkdf = make_trkdf(f, False)
    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    trkdf = trkdf[trkdf.pfp.dist_to_vertex < 10]
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.pfp.trk.chi2pid.I2.chi2_muon/trkdf.pfp.trk.chi2pid.I2.chi2_proton

    # track containment
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = fv_cut(trkdf.pfp.trk.start, DETECTOR) & fv_cut(trkdf.pfp.trk.end, DETECTOR)

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
        print("found empty slice!")
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
        is_contained = pd.Series(dtype='float', name='is_contained', index=empty_index)
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
        is_contained = slcdf.p.pfp.trk.is_contained & slcdf.mu.pfp.trk.is_contained

    ######## (9) - c: slc.tmatch.idx for truth matching
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan
    tmatch_idx_series = slcdf.slc.tmatch.idx

    slc_vtx = slcdf.slc.vertex
    nu_score = slcdf.slc.nu_score
    is_clear_cosmic = slcdf.slc.is_clear_cosmic
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
        'is_clear_cosmic': is_clear_cosmic,
        'nu_score': nu_score,
        'is_contained': is_contained,
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

    EndingRows = len(slcdf)
    if(StartingRows != EndingRows):
        print("You lost a row somewhere! This could impact efficiency calculations!")
        print("Starting Rows:", StartingRows)
        print("Ending Rows:", EndingRows)
        sys.exit()
    return slcdf

def make_gump_nudf(f):
    nudf = make_mcdf(f)
    StartingLength = len(nudf)
    nudf["ind"] = nudf.index.get_level_values(1)
    # wgtdf = pd.concat([bnbsyst.bnbsyst(f, nudf.ind), geniesyst.geniesyst_sbnd(f, nudf.ind)], axis=1)
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
    nmu = nudf.nmu_27MeV
    np = nudf.np_50MeV
    npi = nudf.npi_30MeV
    npi0 = nudf.npi0
    nn = nudf.nn_0MeV
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
        'is_fv': is_fv, 
        'is_1p0pi': is_1p0pi, 
        'is_numu': is_numu, 
        'is_cc': is_cc, 
        'is_sig': is_sig, 
        'pdg': pdg,
        'nmu': nmu,
        'nn': nn,
        'np': np,
        'npi': npi,
        'npi0': npi0 
    })

    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    EndingLength = len(this_nudf)
    if StartingLength != EndingLength:
        print("StartingLength: ",StartingLength)
        print("EndingLength: ",EndingLength)
        sys.exit()
    return this_nudf
