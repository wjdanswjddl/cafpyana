from .makedf import *

## -- private branch definitions
pfp_trk_branches = [
    "rec.slc.reco.pfp.trk.start.x", "rec.slc.reco.pfp.trk.start.y", "rec.slc.reco.pfp.trk.start.z",
    "rec.slc.reco.pfp.trk.end.x", "rec.slc.reco.pfp.trk.end.y", "rec.slc.reco.pfp.trk.end.z",
    "rec.slc.reco.pfp.trk.dir.x", "rec.slc.reco.pfp.trk.dir.y", "rec.slc.reco.pfp.trk.dir.z",
    "rec.slc.reco.pfp.trk.phi", "rec.slc.reco.pfp.trk.costh",
    "rec.slc.reco.pfp.trk.len",
    "rec.slc.reco.pfp.trk.rangeP.p_muon",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_muon",
    "rec.slc.reco.pfp.trk.mcsP.bwdP_muon",
    "rec.slc.reco.pfp.trk.mcsP.is_bwd_muon",
    "rec.slc.reco.pfp.trk.rangeP.p_pion",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_pion",
    "rec.slc.reco.pfp.trk.rangeP.p_proton",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_proton",
    "rec.slc.reco.pfp.trk.bestplane",
]

pfp_trk_mc_branches_names = [
    "interaction_id",
    "parent",
    "pdg",
    "G4ID",
    "end_process",
    "start_process",
    "startE",
    "start.x", "start.y", "start.z",
    "startp.x", "startp.y", "startp.z",
    "end.x", "end.y", "end.z",
    "endp.x", "endp.y", "endp.z",
    "genp.x", "genp.y", "genp.z",
    "genE",
    "length",
    "cont_tpc",
]
pfp_trk_mc_branches = ["rec.slc.reco.pfp.trk.truth.p." + n for n in pfp_trk_mc_branches_names]

pfp_trk_chi2_branches = [
    "rec.slc.reco.pfp.trk.chi2pid.2.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_proton",
    "rec.slc.reco.pfp.trk.chi2pid.1.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_proton",
    "rec.slc.reco.pfp.trk.chi2pid.0.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_proton",
]

pandora_branches = ["rec.slc.reco.pfp.trackScore"]
cnn_branches = [ "rec.slc.reco.pfp.cnnscore.michel", "rec.slc.reco.pfp.cnnscore.endmichel", "rec.slc.reco.pfp.cnnscore.nclusters",
                "rec.slc.reco.pfp.cnnscore.noise", "rec.slc.reco.pfp.cnnscore.shower", "rec.slc.reco.pfp.cnnscore.track"]


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
def dist_pfptrk_vertex(df):
    this_vertex_x = df[('slc', 'vertex', 'x')]
    this_vertex_y = df[('slc', 'vertex', 'y')]
    this_vertex_z = df[('slc', 'vertex', 'z')]

    this_pfp_start_x = df[('trk', 'start', 'x')]
    this_pfp_start_y = df[('trk', 'start', 'y')]
    this_pfp_start_z = df[('trk', 'start', 'z')]

    this_dist = np.sqrt(
        (this_vertex_x - this_pfp_start_x) ** 2 +
        (this_vertex_y - this_pfp_start_y) ** 2 +
        (this_vertex_z - this_pfp_start_z) ** 2
    )

    return this_dist

def Avg(df, pid, drop_0=True):  # average score of 3 planes, exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    #average = df[[("chi2pid", "I0", "chi2_"+pid), ("chi2pid", "I1", "chi2_"+pid), ("chi2pid", "I2", "chi2_"+pid)]].mean(skipna=drop_0, axis=1)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def reco_t(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass):
    #print("reco_t")
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
    n_trk_mupid = group[('n_trk_mupid', '', '')].iloc[0]
    dir_x = group[('trk', 'dir', 'x')]
    dir_y = group[('trk', 'dir', 'y')]
    dir_z = group[('trk', 'dir', 'z')]
    range_P_muon = group[('trk', 'rangeP', 'p_muon')]
    range_P_pion = group[('trk', 'rangeP', 'p_pion')]
    mu_pid_pass = group[('trk', 'mu_pid_pass', '')]

    # Call reco_t function
    return reco_t(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass)

def opening_angle(n_trk_mupid, dir_x, dir_y, dir_z, mu_pid_pass):
    #print("opening_angle")
    if n_trk_mupid != 2:
        return -999.
    dir_x = dir_x[mu_pid_pass]
    dir_y = dir_y[mu_pid_pass]
    dir_z = dir_z[mu_pid_pass]
    if(dir_x.size != 2):
        print("error, dir_x.len != 2")
        return -888.
    
    this_cos_theta = dir_x.iloc[0] * dir_x.iloc[1] + dir_y.iloc[0] * dir_y.iloc[1] + dir_z.iloc[0] * dir_z.iloc[1]
    return this_cos_theta

def measure_opening_angle(group):
    n_trk_mupid = group[('n_trk_mupid', '', '')].iloc[0]
    dir_x = group[('trk', 'dir', 'x')]
    dir_y = group[('trk', 'dir', 'y')]
    dir_z = group[('trk', 'dir', 'z')]
    mu_pid_pass = group[('trk', 'mu_pid_pass', '')]

    # Call reco_t function
    return opening_angle(n_trk_mupid, dir_x, dir_y, dir_z, mu_pid_pass)

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

## -- data fram maker for cohpi analysis
def make_cohpidf(f):
    
    ## 1) Truth df
    nudf = make_mcdf(f)
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

    nudf = pd.DataFrame(index=nudf.index)

    nudf['nuint_categ'] = nuint_categ

    ## 2) slc df
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf.loc[np.invert(slcdf[("rec","slc","tmatch","eff")] > 0.5) & (slcdf[("rec","slc","tmatch","idx")] >= 0), ("rec","slc","tmatch","idx")] = np.nan
    
    matchdf = multicol_merge(slcdf.reset_index(), nudf.reset_index(),
                            left_on=[("entry", "",""), ("rec", "slc","tmatch", "idx")],
                            right_on=[("entry", "",""), ("rec.mc.nu..index", "","")], 
                            how="left") ## -- save all sllices
    
    matchdf = matchdf.set_index(["entry", "rec.slc..index"], verify_integrity=True)
    
    #### 2 - 1) add pfptrack-related columns
    pfptrkdf = loadbranches(f["recTree"], pfp_trk_branches)
    pfptrkdf = pfptrkdf.rec.slc.reco.pfp
    pfptrkchi2df = loadbranches(f["recTree"], pfp_trk_chi2_branches)
    pfptrkchi2df = pfptrkchi2df.rec.slc.reco.pfp.trk
    pfptrkdf = pfptrkdf.join(pfptrkchi2df)

    pfptruthdf = loadbranches(f["recTree"], pfp_trk_mc_branches)
    pfptruthdf = pfptruthdf.rec.slc.reco.pfp.trk.truth
    pfpdf = pd.merge(pfptrkdf, pfptruthdf, left_index=True, right_index=True, how="inner")

    pandoradf = loadbranches(f["recTree"], pandora_branches)
    pandoradf = pandoradf.rec.slc
    cnniddf = loadbranches(f["recTree"], cnn_branches)
    cnniddf = cnniddf.rec.slc.reco
    scoresdf = pd.merge(pandoradf, cnniddf, left_index=True, right_index=True, how="inner")
    pfpdf = pd.merge(pfpdf, scoresdf, left_index=True, right_index=True, how="inner")

    #### 2 - 2) define reco-level event selection variables
    ###### -- FV
    is_reco_fv = InFV(matchdf.rec.slc.vertex, inzback = 0, det = "SBND")
    matchdf[('rec', 'is_reco_fv', '', '')] = is_reco_fv

    ###### -- multiplicity of pfp tracks with length < 4 cm cut
    cut_trk_len = pfpdf.trk.len > 4.
    pfpdf[('trk', 'len', 'pass')] = cut_trk_len
    n_trk_df = cut_trk_len.reset_index(name='len')
    all_combinations = (
        n_trk_df[['entry', 'rec.slc..index']].drop_duplicates().set_index(['entry', 'rec.slc..index'])
    )
    n_trk_df = (
        n_trk_df[n_trk_df['len'] == True]
        .groupby(['entry', 'rec.slc..index'])
        .size()
        .reindex(all_combinations.index, fill_value=0)
    )
    matchdf[('rec', 'n_trk_4cm', '', '')] = n_trk_df

    ###### -- multiplicity of pfp track distance (vertex, trk starting point) cut
    masterdf = pd.merge(matchdf.rec, pfpdf, left_index=True, right_index=True, how="inner")
    this_df_series = dist_pfptrk_vertex(masterdf)
    masterdf['dist_pfptrk_vertex'] = this_df_series
    cut_vtx_dist = masterdf.dist_pfptrk_vertex < 6.
    cut_vtx_dist = cut_trk_len & cut_vtx_dist
    pfpdf[('trk', 'vtxdist', 'pass')] = cut_vtx_dist
    n_pass_vtxdist = cut_vtx_dist.reset_index(name='vtxdist')
    n_pass_vtxdist = (
        n_pass_vtxdist[n_pass_vtxdist['vtxdist'] == True]
	.groupby(['entry', 'rec.slc..index'])
        .size()
        .reindex(all_combinations.index, fill_value=0)
    )
    matchdf[('rec', 'n_trk_vtxdist', '', '')] = n_pass_vtxdist

    ###### -- multiplicity of pfp track chi2 pid cut
    cut_pidscore = (Avg(pfpdf, "muon", drop_0=True) < 25) & (Avg(pfpdf, "proton", drop_0=True) > 100)
    cut_pidscore = cut_pidscore & cut_trk_len & cut_vtx_dist
    pfpdf[('trk', 'mu_pid_pass', '')] = cut_pidscore
    n_pass_mupid_df = cut_pidscore.reset_index(name='pidscore')
    n_pass_mupid_df = (
        n_pass_mupid_df[n_pass_mupid_df['pidscore'] == True]
        .groupby(['entry', 'rec.slc..index'])
        .size()
        .reindex(all_combinations.index, fill_value=0)
    )
    matchdf[('rec', 'n_trk_mupid', '', '')] = n_pass_mupid_df

    ###### -- reco t
    masterdf = pd.merge(matchdf.rec, pfpdf, left_index=True, right_index=True, how="inner")
    reco_t_series = masterdf.groupby(['entry', 'rec.slc..index']).apply(measure_reco_t)
    reco_t_df = reco_t_series.to_frame(name='reco_t_value')
    reco_t_df.index.set_names(['entry', 'rec.slc..index'], inplace=True)
    matchdf[('rec', 'reco_t', '', '')] = reco_t_df

    ###### -- two track opening angle
    opening_angle_series = masterdf.groupby(['entry', 'rec.slc..index']).apply(measure_opening_angle)
    opening_angle_df = opening_angle_series.to_frame(name='reco_opening_angle')
    opening_angle_df.index.set_names(['entry', 'rec.slc..index'], inplace=True)
    matchdf[('rec', 'opening_angle', '', '')] = opening_angle_df

    ###### -- two track momentum sum angle
    beam_totp_angle_series = masterdf.groupby(['entry', 'rec.slc..index']).apply(measure_beam_totp_angle)
    beam_totp_angle_df = beam_totp_angle_series.to_frame(name='reco_beam_totp_angle')
    beam_totp_angle_df.index.set_names(['entry', 'rec.slc..index'], inplace=True)
    matchdf[('rec', 'beam_totp_angle', '', '')] = beam_totp_angle_df

    ###### -- proton stub multiplicity - use only collection plane (plane == 2)
    stubdf = make_stubs(f, det="SBND")
    stubdf = stubdf[stubdf.plane == 2]
    cut_stub_proton = stubdf.pass_proton_stub
    n_stub_proton_df = cut_stub_proton.reset_index(name='pass_proton_stub')
    all_combinations = (
        n_stub_proton_df[['entry', 'rec.slc..index']].drop_duplicates().set_index(['entry', 'rec.slc..index'])
    )
    n_stub_proton_df = (
        n_stub_proton_df[n_stub_proton_df['pass_proton_stub'] == True]
        .groupby(['entry', 'rec.slc..index'])
        .size()
        .reindex(all_combinations.index, fill_value=0)
    )
    n_stub_series = pd.Series(-1, index=matchdf.index)
    n_stub_series.update(n_stub_proton_df)
    matchdf[('rec', 'n_stub_proton', '', '')] = n_stub_series

    return matchdf

