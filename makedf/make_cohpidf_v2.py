from .makedf import *
from pyanalib.pandas_helpers import *
from .util import *
## -- truth level flags
def Signal(df): # definition
    is_fv = InFV(df.position, inzback = 0, det = "SBND")
    is_numu = (df.pdg == 14) | (df.pdg == -14)
    is_cc = (df.iscc == 1)
    is_coh = (df.genie_mode == 3)
    is_1pi0p = (df.nmu_40MeV == 1) & (df.npi_30MeV == 1) & (df.np_50MeV == 0) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_1pi0p & is_coh

def CCCOH(df):
    is_cc = df.iscc
    genie_mode = df.genie_mode
    return is_cc & (genie_mode == 3)

## -- reco level flags
def pass_slc_with_n_pfps(df, n = 2):
    group_levels = ['entry', 'rec.slc..index']
    
    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with at least 2 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def Avg(df, pid, drop_0=True):  # average score of 3 planes, exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    #average = df[[("chi2pid", "I0", "chi2_"+pid), ("chi2pid", "I1", "chi2_"+pid), ("chi2pid", "I2", "chi2_"+pid)]].mean(skipna=drop_0, axis=1)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def add_contained_col(df):
    containd = InFV(df.pfp.trk.start) & InFV(df.pfp.trk.end)
    df[('pfp', 'containd', '', '', '', '')] = containd

def reco_t(dir_x, dir_y, dir_z, range_P_muon, range_P_pion):
    # -- assume first particle is muon and the other is pion
    mass_0 = PDG["muon"][2]
    mass_1 = PDG["pipm"][2]
    p_0 = range_P_muon.iloc[0]
    p_1 = range_P_pion.iloc[1]
    # -- if second track is longer, swap the mass assumption
    if(range_P_muon.iloc[0] > range_P_muon.iloc[1]):
        mass_0 = PDG["pipm"][2]
        mass_1 = PDG["muon"][2]
        p_0 = range_P_pion.iloc[0]
        p_1 = range_P_muon.iloc[1]
    E_0 = np.sqrt(mass_0**2 + p_0**2)
    E_1 = np.sqrt(mass_1**2 + p_1**2)

    # -- each term
    px_sq = np.power(p_0 * dir_x.iloc[0] + p_1 * dir_x.iloc[1], 2.)
    py_sq = np.power(p_0 * dir_y.iloc[0] + p_1 * dir_y.iloc[1], 2.)
    pz_sq = np.power(E_0 + E_1 - p_0 * dir_z.iloc[0] - p_1 * dir_z.iloc[1], 2.)
    abs_t = px_sq + py_sq + pz_sq
    
    #print(abs_t)
    return abs_t

def measure_reco_t(group):
    dir_x = group[('pfp', 'trk', 'dir', 'x', '', '')]
    dir_y = group[('pfp', 'trk', 'dir', 'y', '', '')]
    dir_z = group[('pfp', 'trk', 'dir', 'z', '', '')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_pion = group[('pfp', 'trk', 'rangeP', 'p_pion', '', '')]

    # Call reco_t function
    return reco_t(dir_x, dir_y, dir_z, range_P_muon, range_P_pion)

def opening_angle(dir_x, dir_y, dir_z):
    this_cos_theta = dir_x.iloc[0] * dir_x.iloc[1] + dir_y.iloc[0] * dir_y.iloc[1] + dir_z.iloc[0] * dir_z.iloc[1]
    return this_cos_theta

def measure_opening_angle(group):
    dir_x = group[('pfp', 'trk', 'dir', 'x', '', '')]
    dir_y = group[('pfp', 'trk', 'dir', 'y', '', '')]
    dir_z = group[('pfp', 'trk', 'dir', 'z', '', '')]

    # Call reco_t function
    return opening_angle(dir_x, dir_y, dir_z)

def beam_totp_angle(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass):
    if n_trk_mupid != 2:
        return -999.  
    dir_x = dir_x[mu_pid_pass]
    dir_y = dir_y[mu_pid_pass]
    dir_z = dir_z[mu_pid_pass]
    range_P_muon = range_P_muon[mu_pid_pass]
    range_P_pion = range_P_pion[mu_pid_pass]
    if(range_P_muon.size != 2):
        print("error, dir_x.len != 2")
        return -888.
    
    # -- assume first particle is muon and the other is pion
    p_0 = range_P_muon.iloc[0]
    p_1 = range_P_pion.iloc[1]
    # -- if second track is longer, swap the mass assumption
    if(range_P_muon.iloc[0] > range_P_muon.iloc[1]):
        p_0 = range_P_pion.iloc[0]
        p_1 = range_P_muon.iloc[1]

    totpx = p_0 * dir_x.iloc[0] + p_1 * dir_x.iloc[1]
    totpy = p_0 * dir_y.iloc[0] + p_1 * dir_y.iloc[1]
    totpz = p_0 * dir_z.iloc[0] + p_1 * dir_z.iloc[1]

    totp_cos = totpz / np.power(np.power(totpx, 2.) + np.power(totpy, 2.) + np.power(totpz, 2.) , 0.5)
    return totp_cos
    
def measure_beam_totp_angle(group):
    n_trk_mupid = group[('n_trk_mupid', '', '')].iloc[0]
    dir_x = group[('trk', 'dir', 'x')]
    dir_y = group[('trk', 'dir', 'y')]
    dir_z = group[('trk', 'dir', 'z')]
    range_P_muon = group[('trk', 'rangeP', 'p_muon')]
    range_P_pion = group[('trk', 'rangeP', 'p_pion')]
    mu_pid_pass = group[('trk', 'mu_pid_pass', '')]

    # Call reco_t function
    return beam_totp_angle(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass)


def true_t(this_nuint_categ, E_nu, E_mu, px_mu, py_mu, pz_mu, E_pi, px_pi, py_pi, pz_pi):
    #print(this_nuint_categ)
    if this_nuint_categ.iloc[0] != 1:
        return -999.
    else:
        t = np.power(E_nu.iloc[0] - E_mu.iloc[0] - E_pi.iloc[0], 2.) - np.power(E_nu.iloc[0] - pz_mu.iloc[0] - pz_pi.iloc[0], 2.)
        - np.power(py_mu.iloc[0] - py_pi.iloc[0], 2.) - np.power(px_mu.iloc[0] - px_pi.iloc[0], 2.)
        t = np.fabs(t)
        return t

def get_true_t(group):
    this_nuint_categ = group[('nuint_categ', '', '')]
    E_nu = group[('E', '', '')]

    E_mu = group[('mu', 'genE', '')]
    px_mu = group[('mu', 'genp', 'x')]
    py_mu = group[('mu', 'genp', 'y')]
    pz_mu = group[('mu', 'genp', 'z')]

    E_pi = group[('cpi', 'genE', '')]
    px_pi = group[('cpi', 'genp', 'x')]
    py_pi = group[('cpi', 'genp', 'y')]
    pz_pi = group[('cpi', 'genp', 'z')]    

    return true_t(this_nuint_categ, E_nu, E_mu, px_mu, py_mu, pz_mu, E_pi, px_pi, py_pi, pz_pi)

## -- data fram maker for cohpi analysis
def make_cohpidf_v2(f):
    
    pandora_df = make_pandora_df(f)
    
    #### (1) FV cut
    pandora_df = pandora_df[InFV(df = pandora_df.slc.vertex, inzback = 0, det = "SBND")]

    #### (2) Not clear cosmic cut
    pandora_df = pandora_df[pandora_df.slc.is_clear_cosmic == 0]
    #pandora_df = pandora_df[pandora_df[('slc', 'is_clear_cosmic', '', '', '', '')] == 0]
    
    #### (3) Track multiplicity cut
    ######## (3) - a: keep only pfp objects with length > 4 cm and dist_to_vertex < 6 cm
    pandora_df = pandora_df[(pandora_df.pfp.trk.len > 4.) & (pandora_df.pfp.dist_to_vertex < 6.)]
    ######## (3) - b: keep onlhy slices with exactly two pfp object passing the requirement
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (4) The two trakcs should have track score > 0.5
    pandora_df = pandora_df[pandora_df.pfp.trackScore > 0.5]
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (5) Tho two tracks should satisfy chi2 pid cut 
    pandora_df = pandora_df[(pandora_df.pfp.trk.chi2pid.I2.chi2_muon < 25.) & (pandora_df.pfp.trk.chi2pid.I2.chi2_proton > 100.)]
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (6) Nuscore > 0.65
    pandora_df = pandora_df[pandora_df.slc.nu_score > 0.65]

    #### (7) Cosine of the opening angle between the two track should be greater than 0.2
    opening_angle_series = pandora_df.groupby(['entry', 'rec.slc..index']).apply(measure_opening_angle)
    if pandora_df.empty:
        pandora_df[('slc', 'opening_angle', '', '', '', '')] = pd.Series(dtype='float')
    else:
        pandora_df[('slc', 'opening_angle', '', '', '', '')] = opening_angle_series
    pandora_df = pandora_df[pandora_df.slc.opening_angle > 0.2]

    #### (8) collect only slc variables of interest
    ######## (8) - a: reco momentum/angle of the two tracks
    long_trk_df = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    short_trk_df = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)

    range_p_mu_series = long_trk_df.pfp.trk.rangeP.p_muon
    range_p_mu_series = range_p_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    range_p_pi_series = short_trk_df.pfp.trk.rangeP.p_pion
    range_p_pi_series = range_p_pi_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    cos_theta_mu_series = long_trk_df.pfp.trk.dir.z
    cos_theta_mu_series = cos_theta_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    cos_theta_pi_series = short_trk_df.pfp.trk.dir.z
    cos_theta_pi_series = cos_theta_pi_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    ######## (8) - b: reco |t|
    if pandora_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry', 'rec.slc..index']
        )
        reco_t_series = pd.Series(dtype='float', name='reco_t', index=empty_index)
        #reco_t_series = pd.Series(dtype='float', name='reco_t')
    else:
        reco_t_series = pandora_df.groupby(['entry', 'rec.slc..index']).apply(measure_reco_t)
    #reco_t_series = reco_t_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    ######## (8) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    
    ## (9) creat a slice-based reco df
    slcdf = pd.DataFrame({
        'range_p_mu': range_p_mu_series,
        'range_p_pi': range_p_pi_series,
        'cos_theta_mu': cos_theta_mu_series,
        'cos_theta_pi': cos_theta_pi_series,
        'reco_t': reco_t_series,
        'tmatch_idx': tmatch_idx_series
    })

    #print(slcdf.columns)
    return slcdf
    
def make_cohpi_nudf(f):
    #nudf = make_mcdf(f)
    nudf = make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)
    wgtdf = geniesyst.geniesyst_sbnd(f, nudf.ind)
    
    is_fv = InFV(df = nudf.position, inzback = 0, det = "SBND")
    is_signal = Signal(nudf)
    is_cc = nudf.iscc
    genie_mode = nudf.genie_mode
    w = nudf.w

    try :
        nuint_categ = pd.Series(8, index=nudf.index)
        #print(f"done init nuint_categ")
    except Exception as e:
        print(f"Error init nuint_categ")
        return

    nuint_categ[~is_fv] = -1  # Out of FV
    nuint_categ[is_fv & ~is_cc] = 0  # NC
    nuint_categ[is_fv & is_cc & is_signal] = 1  # Signal
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 3)] = 2  # Non-signal CCCOH
    nuint_categ[is_fv & is_cc & (genie_mode == 0)] = 3  # CCQE
    nuint_categ[is_fv & is_cc & (genie_mode == 10)] = 4  # 2p2h
    nuint_categ[is_fv & is_cc & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w < 1.4) | (genie_mode == 1))] = 5  # RES
    nuint_categ[is_fv & is_cc & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w > 2.0) | (genie_mode == 2))] = 6  # DIS
    nuint_categ[is_fv & is_cc & ((1.4 < w) & (w < 2.0) & (genie_mode != 1) & (genie_mode != 2) & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10))] = 7  # INEL

    #nudf = pd.DataFrame(index=nudf.index)
    #print(nuint_categ)
    nudf['nuint_categ'] = nuint_categ
    
    true_t_series = nudf.groupby(['entry', 'rec.mc.nu..index']).apply(get_true_t)

    mu_p_series = magdf(nudf.mu.genp)
    cpi_p_series = magdf(nudf.cpi.genp)

    mu_cos_theta_series = nudf.mu.genp.z / mu_p_series
    cpi_cos_theta_series = nudf.cpi.genp.z / cpi_p_series

    this_nudf = pd.DataFrame({
        'p_mu': mu_p_series,
        'p_pi': cpi_p_series,
        'cos_theta_mu': mu_cos_theta_series,
        'cos_theta_pi': cpi_cos_theta_series,
        'true_t': true_t_series,
        'nuint_categ': nuint_categ
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])
    
    print(this_nudf.columns)
    print(wgtdf.columns)
    
    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
