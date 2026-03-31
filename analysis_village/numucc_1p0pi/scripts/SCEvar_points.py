import pandas as pd
import numpy as np
import sys
from os import path, makedirs
from datetime import datetime

# local imports
# sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/cafpyana_2026Jan17/cafpyana') # absolute path for running on EAF
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.categories import *
from makedf.geniesyst import *
from analysis_village.numucc_1p0pi.utils import *

import matplotlib.pyplot as plt 
from matplotlib.patches import Patch

plt.style.use("presentation.mplstyle")

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
import multiprocessing as mp
import pickle

# load dfs
# CV
ret_dfs = load_and_concat_mc_dfs(
    file_dir="/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10",
    sub_dir="",
    sample_dir="",
    df_tag="",
    chunk_tags=["SystVar_CV_sel_2prong"],
    keys2load=['meta', 'evt'],
    n_max_concat=999
)

cv_mc_df = ret_dfs["evt"]
cv_meta_df = ret_dfs["meta"]

cv_mc_df.loc[:,'topo_categ'] = get_topo_category(cv_mc_df)
cv_mc_df.loc[:,'genie_categ'] = get_genie_category(cv_mc_df)

# SCE vars
ret_0xSCE_dfs = load_and_concat_mc_dfs(
    file_dir="/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10",
    sub_dir="",
    sample_dir="",
    df_tag="",
    chunk_tags=["SystVar_0xSCE_sel_2prong"],
    keys2load=['meta', 'evt'],
    n_max_concat=999
)

_0xSCE_mc_df = ret_0xSCE_dfs["evt"]
_0xSCE_meta_df = ret_0xSCE_dfs["meta"]

_0xSCE_mc_df.loc[:,'topo_categ'] = get_topo_category(_0xSCE_mc_df)
_0xSCE_mc_df.loc[:,'genie_categ'] = get_genie_category(_0xSCE_mc_df)

ret_2xSCE_dfs = load_and_concat_mc_dfs(
    file_dir="/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10",
    sub_dir="",
    sample_dir="",
    df_tag="",
    chunk_tags=["SystVar_2xSCE_sel_2prong"],
    keys2load=['meta', 'evt'],
    n_max_concat=999
)

_2xSCE_mc_df = ret_2xSCE_dfs["evt"]
_2xSCE_meta_df = ret_2xSCE_dfs["meta"]

_2xSCE_mc_df.loc[:,'topo_categ'] = get_topo_category(_2xSCE_mc_df)
_2xSCE_mc_df.loc[:,'genie_categ'] = get_genie_category(_2xSCE_mc_df)


# select events that exist in all detvar samples
cv_idx_df = cv_meta_df.reset_index().set_index(["run", "subrun", "evt", "E"])
meta_0xSCE_idx_df = ret_0xSCE_dfs["meta"].reset_index().set_index(["run", "subrun", "evt", "E"])
meta_2xSCE_idx_df = ret_2xSCE_dfs["meta"].reset_index().set_index(["run", "subrun", "evt", "E"])

cv_idx_list = cv_idx_df.index
common_idx_list = [idx for idx in cv_idx_list if idx in meta_0xSCE_idx_df.index and idx in meta_2xSCE_idx_df.index]
# common_idx_list = [idx for idx in meta_0xSCE_idx_df.index if idx in meta_2xSCE_idx_df.index]

meta_0xSCE_idx_df = meta_0xSCE_idx_df.loc[common_idx_list]
ret_0xSCE_dfs["meta"] = meta_0xSCE_idx_df
meta_2xSCE_idx_df = meta_2xSCE_idx_df.loc[common_idx_list]
ret_2xSCE_dfs["meta"] = meta_2xSCE_idx_df

print(len(meta_0xSCE_idx_df.index))
print(len(meta_2xSCE_idx_df.index))
print(len(common_idx_list))

new_dfs = []
for detvar_dfs in [ret_dfs, ret_0xSCE_dfs, ret_2xSCE_dfs]:
    for k in tqdm(detvar_dfs.keys()):
        if "meta" in k:
            continue

        detvar_dfs[k]["E"] = detvar_dfs[k].mc.E.copy()
        this_meta_df = detvar_dfs["meta"].copy()
        this_meta_df = this_meta_df.reset_index().set_index(["run", "subrun", "evt", "E"])
        this_meta_df = this_meta_df.loc[common_idx_list]
        this_common_idx = this_meta_df.reset_index().set_index(["__ntuple", "entry", "E"]).index
        this_sel_df = detvar_dfs[k].reset_index().set_index(["__ntuple", "entry", "E"]) # .loc[this_common_idx]

        this_sel_idx = this_sel_df.index
        this_sel_common_idx = [idx for idx in this_sel_idx if idx in this_common_idx]
        this_sel_df = this_sel_df.loc[this_sel_common_idx]
        this_sel_df = this_sel_df.reset_index().set_index(["__ntuple", "entry", "rec.slc..index"])

        # mc_tot_pot = detvar_dfs["meta"]['pot'].sum()
        # print("mc_tot_pot: %.3e" %(mc_tot_pot))
        # data_tot_pot = 4.57e18
        # mc_pot_scale = data_tot_pot / mc_tot_pot
        # this_sel_df["pot_weight"] = mc_pot_scale * np.ones(len(this_sel_df))

        detvar_dfs[k] = this_sel_df.groupby(level=[0,1]).head(1)

        # plt.hist(this_sel_df.mc.E, bins=bins, histtype="step", label=k)

    new_dfs.append(detvar_dfs)


# get start positions
nevts = len(ret_dfs["meta"])
print("nevts: ", nevts)

start_x_list = []
start_y_list = []
start_z_list = []

# Worker function for multiprocess; needs to be at top level for Pool
def process_event(eidx):
    eidx = eidx + 10000
    try:
        idx0 = ret_dfs["meta"].reset_index().iloc[eidx]["__ntuple"].values[0]
        idx1 = ret_dfs["meta"].reset_index().iloc[eidx]["entry"].values[0]
        Eval_cv = ret_dfs["evt"].loc[idx0, idx1].E
        start_x = ret_dfs["evt"].loc[idx0, idx1].trk1.pfp.trk.start.x
        start_y = ret_dfs["evt"].loc[idx0, idx1].trk1.pfp.trk.start.y
        start_z = ret_dfs["evt"].loc[idx0, idx1].trk1.pfp.trk.start.z

        idx0_2x = ret_2xSCE_dfs["meta"].iloc[eidx]["__ntuple"].values[0]
        idx1_2x = ret_2xSCE_dfs["meta"].iloc[eidx]["entry"].values[0]
        Evar_2x = ret_2xSCE_dfs["evt"].loc[idx0_2x, idx1_2x].E
        start_x_2x = ret_2xSCE_dfs["evt"].loc[idx0_2x, idx1_2x].trk1.pfp.trk.start.x
        start_y_2x = ret_2xSCE_dfs["evt"].loc[idx0_2x, idx1_2x].trk1.pfp.trk.start.y
        start_z_2x = ret_2xSCE_dfs["evt"].loc[idx0_2x, idx1_2x].trk1.pfp.trk.start.z

        idx0_0x = ret_0xSCE_dfs["meta"].iloc[eidx]["__ntuple"].values[0]
        idx1_0x = ret_0xSCE_dfs["meta"].iloc[eidx]["entry"].values[0]
        Evar_0x = ret_0xSCE_dfs["evt"].loc[idx0_0x, idx1_0x].E
        start_x_0x = ret_0xSCE_dfs["evt"].loc[idx0_0x, idx1_0x].trk1.pfp.trk.start.x
        start_y_0x = ret_0xSCE_dfs["evt"].loc[idx0_0x, idx1_0x].trk1.pfp.trk.start.y
        start_z_0x = ret_0xSCE_dfs["evt"].loc[idx0_0x, idx1_0x].trk1.pfp.trk.start.z

        assert Eval_cv.values[0] == Evar_2x.values[0] == Evar_0x.values[0]
        return (
            [start_x.values[0], start_x_2x.values[0], start_x_0x.values[0]],
            [start_y.values[0], start_y_2x.values[0], start_y_0x.values[0]],
            [start_z.values[0], start_z_2x.values[0], start_z_0x.values[0]],
            start_x.values[0]  # For skipping duplicates in post-process
        )
    except Exception:
        return None

from tqdm import tqdm

# Use multiprocessing Pool to process events in parallel
results = []
with mp.Pool(processes=20) as pool:
    for res in tqdm(pool.imap(process_event, range(30000)), total=30000):
        if res is not None:
            results.append(res)

# Post-process: filter unique start_x values in order of processed events.
seen = set()
for sx, sy, sz, val_x in results:
    if val_x in seen:
        continue
    seen.add(val_x)
    start_x_list.append(sx)
    start_y_list.append(sy)
    start_z_list.append(sz)


save_dict = {
    "x": start_x_list,
    "y": start_y_list,
    "z": start_z_list,
}


with open("start_pos_list-SCEvar-2.pkl", "wb") as f:
    pickle.dump(save_dict, f)
