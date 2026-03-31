import pandas as pd
import numpy as np
from os import path, makedirs
from datetime import datetime
from tqdm import tqdm

# local imports
import sys
sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/cafpyana_2026Jan17/cafpyana') # absolute path for running on EAF
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.constants import *
from analysis_village.numucc_1p0pi.files_config import *
from pyanalib.covariance import *
from makedf.mcstat import get_MCstat_unc

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


file_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09"

n_max_concat = 1000
# -- Offbeam Data
offbeam_keys2load = ['hdr', 'evt', 'trk']

tot_gates = 0
n_data_dict = {}
dirs = ["x", "y", "z"] 
for stages in ["all", "clear_cosmic", "vertex_in_fv"]:
    n_data_dict[stages] = {}
    for this_dir in dirs:
        n_data_dict[stages][this_dir] = []

for tag in tqdm(generate_tags("ad")):
    offbeam_dfs = load_and_concat_mc_dfs(
        file_dir=file_dir,
        chunk_tags=[tag],
        df_tag="_all",
        keys2load=offbeam_keys2load,
        n_max_concat=n_max_concat,
        sub_dir="data",
        sample_dir="OffBeam"
    )
    offbeam_hdr_df = offbeam_dfs['hdr']
    offbeam_evt_df = offbeam_dfs['evt']
    offbeam_trk_df = offbeam_dfs['trk']

    offbeam_gates = offbeam_hdr_df.noffbeambnb.sum()
    offbeam_gates = offbeam_hdr_df[offbeam_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
    print("offbeam cosmics data gates: {:.2e}".format(offbeam_gates))
    tot_gates += offbeam_gates

    offbeam_trk_df_longest = offbeam_trk_df.sort_values(by=[('pfp', 'trk', 'len')], ascending=False)
    offbeam_trk_df_longest = offbeam_trk_df_longest.groupby(level=[0,1,2]).head(1)

    bins = np.linspace(-1,1,51)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    for which_dir in ["x", "y", "z"]:
        var_data = offbeam_trk_df_longest.pfp.trk.dir[which_dir]
        n_data, _ = np.histogram(var_data, bins=bins)
        n_data_dict["all"][which_dir].append(n_data)

    offbeam_trk_df_longest = offbeam_trk_df_longest[offbeam_evt_df.slc.is_clear_cosmic == 0]
    for which_dir in ["x", "y", "z"]:
        var_data = offbeam_trk_df_longest.pfp.trk.dir[which_dir]
        n_data, _ = np.histogram(var_data, bins=bins)
        n_data_dict["clear_cosmic"][which_dir].append(n_data)

    offbeam_trk_df_longest = offbeam_trk_df_longest[InFV(offbeam_trk_df_longest.pfp.trk.start, det="SBND_nohighyz")]
    for which_dir in ["x", "y", "z"]:
        var_data = offbeam_trk_df_longest.pfp.trk.dir[which_dir]
        n_data, _ = np.histogram(var_data, bins=bins)
        n_data_dict["vertex_in_fv"][which_dir].append(n_data)

#save as pickle
n_data_dict["tot_gates"] = tot_gates
with open(file_dir + "/offbeam_data_distribution_dict.pkl", "wb") as f:
    pickle.dump(n_data_dict, f)

## -- Intime MC
intime_keys2load = ['hdr', 'evt', 'trk']

tot_gates = 0
n_mc_dict = {}
dirs = ["x", "y", "z"] 
for stages in ["all", "clear_cosmic", "vertex_in_fv"]:
    n_mc_dict[stages] = {}
    for this_dir in dirs:
        n_mc_dict[stages][this_dir] = []

for tag in tqdm(generate_tags("au")):
    intime_dfs = load_and_concat_mc_dfs(
        file_dir=file_dir,
        chunk_tags=[tag],
        df_tag="_all",
        keys2load=intime_keys2load,
        n_max_concat=n_max_concat,
        sub_dir="MC",
        sample_dir="intime"
    )
    intime_hdr_df = intime_dfs['hdr']
    intime_evt_df = intime_dfs['evt']
    intime_trk_df = intime_dfs['trk']

    intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['ngenevt'].sum()
    print("intime cosmics data gates: {:.2e}".format(intime_gates))
    tot_gates += intime_gates

    intime_trk_df_longest = intime_trk_df.sort_values(by=[('pfp', 'trk', 'len')], ascending=False)
    intime_trk_df_longest = intime_trk_df_longest.groupby(level=[0,1,2]).head(1)
    intime_trk_df_longest = intime_trk_df_longest[intime_evt_df.slc.is_clear_cosmic == 0]
    intime_trk_df_longest = intime_trk_df_longest[InFV(intime_trk_df_longest.pfp.trk.start, det="SBND_nohighyz")]

    bins = np.linspace(-1,1,51)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    for which_dir in ["x", "y", "z"]:
        var_mc = intime_trk_df_longest.pfp.trk.dir[which_dir]
        n_mc, _ = np.histogram(var_mc, bins=bins)
        n_mc_dict["all"][which_dir].append(n_mc)

    intime_trk_df_longest = intime_trk_df_longest[intime_evt_df.slc.is_clear_cosmic == 0]
    for which_dir in ["x", "y", "z"]:
        var_mc = intime_trk_df_longest.pfp.trk.dir[which_dir]
        n_mc, _ = np.histogram(var_mc, bins=bins)
        n_mc_dict["clear_cosmic"][which_dir].append(n_mc)

    intime_trk_df_longest = intime_trk_df_longest[InFV(intime_trk_df_longest.pfp.trk.start, det="SBND_nohighyz")]
    for which_dir in ["x", "y", "z"]:
        var_mc = intime_trk_df_longest.pfp.trk.dir[which_dir]
        n_mc, _ = np.histogram(var_mc, bins=bins)
        n_mc_dict["vertex_in_fv"][which_dir].append(n_mc)

n_mc_dict["tot_gates"] = tot_gates
with open(file_dir + "/intime_mc_distribution_dict.pkl", "wb") as f:
    pickle.dump(n_mc_dict, f)
