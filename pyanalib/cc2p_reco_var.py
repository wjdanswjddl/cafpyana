import numpy as np

PDG = {
    "muon": [13, "muon", 0.105,],
    "proton": [2212, "proton", 0.938272,],
    "neutron": [2112, "neutron", 0.9395654,],
    "pizero": [111, "pizero", 0.1349768],
    "pipm": [211, "piplus", 0.13957039],
    "argon": [1000180400, "argon", (18*0.938272 + 22*0.9395654)],
    "gamma": [22, "gamma", 0 ],
    "lambda": [3122, "lambda", 1.115683],
    "kaon_p": [321, "kaon_p",  0.493677],
    "sigma_p": [3222, "sigma_p", 1.18936],
    "kaon_0": [311, "kaon_0", 0.497648],
    "sigma_0": [3212, "sigma_0", 1.19246],
    "lambda_p_c": [4122, "lambda_p_c", 2.28646],
    "sigma_pp_c": [4222, "sigma_pp_c", 2.45397],
    "electron": [11, "electron", 0.510998950],
    "sigma_p_c": [4212, "sigma_p_c", 2.4529],
}

def InFV_trk(data): # cm
    xmin = -190.
    ymin = -190.
    zmin = 10.
    xmax = 190.
    ymax =  190.
    zmax =  490.
    return (data.x > xmin) & (data.x < xmax) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)


def InFV(data): # cm
    xmin = -190.
    ymin = -190.
    zmin = 10.
    xmax = 190.
    ymax =  190.
    zmax =  450.
    return (np.abs(data.x) > 10) & (np.abs(data.x) < 190) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)

## -- truth level flags
def Signal(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_2p0pi = (df.nmu_40MeV == 1) & (df.npi_30MeV == 0) & (df.np_50MeV == 2) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_2p0pi

## -- reco level flags
def pass_slc_with_n_pfps(df, n = 3):
    group_levels = ['entry', 'rec.slc..index']
    
    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with exactly 3 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]
    
    filtered_df = get_n_recopid_per_slc(filtered_df)
    
    filtered_df = filtered_df[ (filtered_df["muon_counter"] == 1) & (filtered_df["proton_counter"] == 2) & (filtered_df["pion_counter"] == 0)]    

    return filtered_df

def Avg(df, pid, drop_0=True):  # exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def get_pid_result(row):
    
    chi2_muon = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon', '')]
    chi2_proton = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton', '')]
    len = row[('pfp', 'trk', 'len', '', '', '')]

    if chi2_muon < 25. and chi2_proton > 100.:
        return 13  # muon
    elif chi2_proton < 90.:
        return 2212  # proton
    else:
        return 211  # charged pion   
        
def add_n_slice_col(reco_df):
    df_reset = reco_df.reset_index()
    slc_counts = (
        df_reset[['entry', 'rec.slc..index']]
        .drop_duplicates()
        .groupby(['entry'])
        .size()
        .reset_index(name='n_slc_per_entry')
    )

    slc_counts.columns = pd.MultiIndex.from_tuples([
        ('entry', '', '', '', '', ''),
        ('slc', 'n_slc_per_entry', '', '', '', '')
    ])
    df_reset = df_reset.merge(slc_counts, on=[('entry', '', '', '', '', '')])
    df_reset = df_reset.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"], verify_integrity=True)
    return df_reset     

def get_n_recopid_per_slc(df):
    
    pid_series = df.pfp.trk.reco_pid
    this_df = pid_series.reset_index()

    muons = this_df[this_df["reco_pid"] == 13]
    protons = this_df[this_df["reco_pid"] == 2212]
    pions = this_df[this_df["reco_pid"] == 211]    

    muon_counts = muons.groupby(["entry", "rec.slc..index"]).size().rename("n_mu")
    proton_counts = protons.groupby(["entry", "rec.slc..index"]).size().rename("n_proton")
    pion_counts = pions.groupby(["entry", "rec.slc..index"]).size().rename("n_pion")    

    this_df = this_df.merge(muon_counts, on=["entry", "rec.slc..index"], how="left")
    this_df = this_df.merge(proton_counts, on=["entry", "rec.slc..index"], how="left")
    this_df = this_df.merge(pion_counts, on=["entry", "rec.slc..index"], how="left")    

    this_df["n_mu"] = this_df["n_mu"].fillna(0).astype(int)
    this_df["n_proton"] = this_df["n_proton"].fillna(0).astype(int)
    this_df["n_pion"] = this_df["n_pion"].fillna(0).astype(int)    

    this_df.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"], inplace=True)

    df[('muon_counter', '', '', '', '', '')] = this_df.n_mu
    df[('proton_counter', '', '', '', '', '')] = this_df.n_proton
    df[('pion_counter', '', '', '', '', '')] = this_df.n_pion    

    return df

def add_contained_col(df):
    contained = InFV(df.pfp.trk.start) & InFV(df.pfp.trk.end)
    df[('pfp', 'contained', '', '', '', '')] = contained
    
def reco_deltapt(muon_dir_x, muon_dir_y, muon_dir_z, range_P_muon, 
                 proton_dir_x, proton_dir_y, proton_dir_z, range_P_proton):

    # # -- each term
    px_sq = np.power(muon_dir_x * range_P_muon + proton_dir_x * range_P_proton, 2.)
    py_sq = np.power(muon_dir_y * range_P_muon + proton_dir_y * range_P_proton, 2.)
    deltapt = np.sqrt(px_sq + py_sq)
    
    #print(deltapt)
    return deltapt

def measure_reco_deltapt(group):
    
    muons = group[group.pfp.trk.reco_pid == 13]
    protons = group[group.pfp.trk.reco_pid == 2212]    
    
    muon_dir_x = muons[('pfp', 'trk', 'dir', 'x', '', '')]
    muon_dir_y = muons[('pfp', 'trk', 'dir', 'y', '', '')]
    muon_dir_z = muons[('pfp', 'trk', 'dir', 'z', '', '')]
    muon_range = muons[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    
    proton_dir_x = protons[('pfp', 'trk', 'dir', 'x', '', '')]
    proton_dir_y = protons[('pfp', 'trk', 'dir', 'y', '', '')]
    proton_dir_z = protons[('pfp', 'trk', 'dir', 'z', '', '')]
    proton_range = protons[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]    
    
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_proton = group[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]

    # Call reco_deltapt function
    return reco_deltapt(muon_dir_x, muon_dir_y, muon_dir_z, muon_range, proton_dir_x, proton_dir_y, proton_dir_z, proton_range)  