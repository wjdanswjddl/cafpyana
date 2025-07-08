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