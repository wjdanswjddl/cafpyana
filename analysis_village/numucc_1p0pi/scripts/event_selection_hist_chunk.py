# Event Selection

import pandas as pd
import numpy as np
import sys
from os import path, makedirs
from datetime import datetime
import pickle

# local imports
# sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana') # absolute path for running on EAF
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.categories import *
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.files_config import *
from analysis_village.numucc_1p0pi.makedf.selections import *
from pyanalib.split_df_helpers import *
from pyanalib.pandas_helpers import *
from pyanalib.covariance import *

import matplotlib.pyplot as plt 
from matplotlib.patches import Patch
from matplotlib.colors import LogNorm

plt.style.use("presentation.mplstyle")

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
# turn off SettingWithCopyWarning:o
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--df_dir", type=str)
parser.add_argument("--filename", type=str, required=True)
parser.add_argument("--file_type", type=str, default="mc", choices=["mc", "data", "offbeam", "intime", "dirt"])
parser.add_argument("--save_content", action="store_true")
parser.add_argument("--save_name", type=str)
args = parser.parse_args()

file_type = args.file_type


def get_n_breakdown(df, var_config, plot_type="topology", file_type="mc"):

    vardf, _ = get_clipped_evts(df, var_config.var_evt_reco_col, var_config.bins)

    if file_type == "mc":
        if plot_type == "topology":
            cuts = get_topo_category(df, ret_cuts=True)
        elif plot_type == "genie":
            cuts = get_genie_category(df, ret_cuts=True)
        elif plot_type == "pdg":
            cuts = get_pdg_category(df, ret_cuts=True)
        else:
            raise ValueError("Invalid plot_type: %s, please choose between [topology, genie, or pdg]" % plot_type)

        var_categ = [vardf[i] for i in cuts]

        n_breakdown = []
        for vc in var_categ:
            n_vc, _ = np.histogram(vc, bins=var_config.bins)
            n_breakdown.append(n_vc)
        n_breakdown = np.array(n_breakdown)

    else:
        n_v, _ = np.histogram(vardf, bins=var_config.bins)
        n_breakdown = np.array([n_v])

    return n_breakdown



# plot breakdown type
plot_type = "topology"

today_str = datetime.now().strftime("%Y%m%d")

# save contents for plots
save_content = args.save_content
# save_content_base_dir = "/exp/sbnd/data/users/munjung/xsec/RESULTS/numuCC_1p0pi"
# save_content_dir = path.join(save_content_base_dir, f"event_selection-{today_str}")
# save_content_dir = args.save_dir
# if save_content:
#     if not path.exists(save_content_dir):
#         makedirs(save_content_dir)
#     print("saving nevts in ", save_content_dir)


# ===== load dfs =====
# df_dir_mc = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_01_102644__sel_all-mc-BNB_cosmics"
# filename_str_mc = "sel_all-mc-BNB_cosmics_100"
# keys2load_mc = ['evt', 'trk', 'hit0', 'hit1', 'hit2', 'hdr']
keys2load = ['evt', 'trk', 'hdr']
# loaded_dfs = dfs_from_dir(search_dir=args.df_dir, filename_str=args.filename, keys2load=keys2load, n_max_concat=999)
filename = path.join(args.df_dir, args.filename)
loaded_dfs = load_dfs(filename, keys2load, n_max_concat=999)
hdr_df, evt_df, trk_df = loaded_dfs["hdr"], loaded_dfs["evt"], loaded_dfs["trk"]

# TODO: move to dfmaker
# ===== calculate avg chi2 for pid =====
chimu_avg = avg_chi2(trk_df, "chi2_muon")
trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu_avg
chip_avg = avg_chi2(trk_df, "chi2_proton")
trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip_avg
# ===== get mcs range diff =====
def get_mcs_range_diff(trks):
    mcs_range_diff = (trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon
    trks[("pfp", "trk", "mcs_range_diff", "", "", "")] = mcs_range_diff
    return trks

trk_df = get_mcs_range_diff(trk_df)

# ===== exposure calculation =====
print(f"file_type: {file_type}")
if file_type == "mc":
    tot_pot = hdr_df['pot'].sum()
    tot_gates = hdr_df.nbnbinfo.sum()

elif file_type == "data":
    tot_pot = hdr_df['pot'].sum()
    tot_gates = np.nan

elif file_type == "offbeam":
    tot_pot = np.nan
    tot_gates = hdr_df[hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()

elif file_type == "intime":
    tot_pot = np.nan
    tot_gates = hdr_df[hdr_df['first_in_subrun'] == 1]['ngenevt'].sum()

elif file_type == "dirt":
    tot_pot = hdr_df['pot'].sum()
    tot_gates = np.nan

else:
    raise ValueError(f"Invalid file_type: {file_type}, please choose between [mc, data, offbeam, intime, dirt]")
print(f"tot_pot: {tot_pot:.3e}", f"tot_gates: {tot_gates:.3e}")

# ===== save breakdowns for summary plots =====
breakdown_dict = {"topology": {}, "genie": {}}

# save dfs per stage, df_dict used for efficiency plot
stage_names = []
df_dict = {} 

# save histogram contents
df_hist_contents = {"tot_pot": tot_pot, "tot_gates": tot_gates}

# ===== all reconstructed slices =====
stage_key = "allreco"
stage_names.append(stage_key)
df_dict[stage_key] = evt_df
# =====================================


# ===== not clear cosmic =====
stage_key = "is_clear_cosmic"
stage_names.append(stage_key)
evt_df = cut_clear_cosmic(evt_df)
df_dict[stage_key] = evt_df
# =====================================


# ===== vertex in fiducial volume =====
stage_key = "vertex_in_fv"
stage_names.append(stage_key)
evt_df = cut_vertex_in_fv(evt_df, det=DETECTOR)
df_dict[stage_key] = evt_df
# =====================================

# ****** nu-score plot ******
var_config = VariableConfig.nu_score()
n_mc_breakdown = get_n_breakdown(evt_df, var_config, plot_type="topology", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ----- get valid tracks
trk_df = get_valid_trks(trk_df)
trk_df = match_trkdf_to_slcdf(trk_df, evt_df)
evt_df = get_trk_info(evt_df, trk_df, SAVE_NTRKS)

# ****** n tracks plot ******
var_config = VariableConfig.n_trks()
n_mc_breakdown = get_n_breakdown(evt_df, var_config, plot_type="topology", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ===== 2prong cut =====
stage_key = "2prong"
stage_names.append(stage_key)
evt_df = cut_2prong(evt_df)
df_dict[stage_key] = evt_df
# =====================================


# ===== 2prong contained =====
stage_key = "2prong-contained"
stage_names.append(stage_key)
evt_df = cut_2prong_contained(evt_df, det=DETECTOR)
df_dict[stage_key] = evt_df
# =====================================

# ----- take both tracks for plotting
def get_both_trks(df):
    return pd.concat([df.trk1, df.trk2])

# ****** track score plot ******
var_config = VariableConfig.track_score()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ===== 2prong track score cut =====
stage_key = "2prong-trackscore"
stage_names.append(stage_key)
evt_df = cut_2prong_trackscore(evt_df, TRACKSCORE_TH)
df_dict[stage_key] = evt_df
# =====================================

# ****** vtx dist plot ******
var_config = VariableConfig.vtx_dist()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ===== 2prong vtx dist cut =====
stage_key = "2prong-vtxdist"
stage_names.append(stage_key)
evt_df = cut_2prong_vtxdist(evt_df, VTXDIST_TH)
df_dict[stage_key] = evt_df
# =====================================

# ****** trk len plot ******
var_config = VariableConfig.trk_len()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ****** mcs range diff plot ******
var_config = VariableConfig.mcs_range_diff()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************


# ****** chi2_mu plot ******
var_config = VariableConfig.chi2_mu()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ****** chi2_proton plot ******
var_config = VariableConfig.chi2_proton()
n_mc_breakdown = get_n_breakdown(get_both_trks(evt_df), var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ****** chi2_mu vs. chi2_proton plot ******
var_config_1 = VariableConfig.chi2_mu()
var_config_2 = VariableConfig.chi2_proton()

n, bins1, bins2 = np.histogram2d(get_both_trks(evt_df)[var_config_1.var_evt_reco_col], 
           get_both_trks(evt_df)[var_config_2.var_evt_reco_col], 
           bins=[var_config_1.bins, var_config_2.bins])

df_hist_contents["_".join([stage_key, "chi2","2D"])] = {"n": n, "bins1": bins1, "bins2": bins2}
# *********************************

#  ----- PID selection definitions -----
def len_cut(trks):
    trks = trks[trks.pfp.trk.len > 50]
    return trks

def mcs_range_diff_cut(trks):
    trks = trks[np.abs(trks.pfp.trk.mcs_range_diff) < QUAL_TH]
    return trks

def is_mu_candidate(trks):
    trks = len_cut(trks)
    trks = mcs_range_diff_cut(trks)
    chimu_avg = trks.pfp.trk.chi2pid.avg.chi2_muon
    chip_avg = trks.pfp.trk.chi2pid.avg.chi2_proton
    chi2_cut = (chimu_avg > 0) & (chimu_avg < MU_CHI2MU_TH) & (chip_avg > MU_CHI2P_TH) 
    trks = trks[chi2_cut]
    return trks

def is_not_mu_candidate(trks):
    nlevels = len(evt_df.index.names)
    mcs_range_diff = np.abs((trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon)
    chimu_avg = trks.pfp.trk.chi2pid.avg.chi2_muon
    chip_avg = trks.pfp.trk.chi2pid.avg.chi2_proton
    mu_cut = (chimu_avg > 0) & (chimu_avg < MU_CHI2MU_TH) & \
            (chip_avg > MU_CHI2P_TH) & \
            (trks.pfp.trk.len > MU_LEN_TH) & \
            (mcs_range_diff < QUAL_TH)
    not_mu_candidate = pd.concat([trks[~mu_cut], trks[mu_cut].groupby(level=list(range(nlevels))).nth(1)])
    return not_mu_candidate


mus = is_mu_candidate(trk_df)

notmus = is_not_mu_candidate(trk_df)

# ****** chi2_mu plot ******
var_config = VariableConfig.chi2_mu()
n_mc_breakdown = get_n_breakdown(notmus, var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join(["not_mu", var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ****** chi2_proton plot ******
var_config = VariableConfig.chi2_proton()
n_mc_breakdown = get_n_breakdown(notmus, var_config, plot_type="pdg", file_type=file_type)
df_hist_contents["_".join(["not_mu", var_config.var_save_name])] = {"n": n_mc_breakdown, "bins": var_config.bins}
# *********************************

# ****** chi2_mu vs. chi2_proton plot ******
var_config_1 = VariableConfig.chi2_mu()
var_config_2 = VariableConfig.chi2_proton()

n, bins1, bins2 = np.histogram2d(notmus[var_config_1.var_evt_reco_col], 
           notmus[var_config_2.var_evt_reco_col], 
           bins=[var_config_1.bins, var_config_2.bins])

df_hist_contents["_".join(["not_mu", "chi2","2D"])] = {"n": n, "bins1": bins1, "bins2": bins2}
# *********************************

# ----- get mu & p candidates
evt_df = get_mu_p_candidate(evt_df, 
                       mu_chi2mu_th=MU_CHI2MU_TH, mu_chi2p_th=MU_CHI2P_TH, mu_len_th=MU_LEN_TH, qual_th=QUAL_TH,
                       p_chi2mu_th=-1, p_chi2p_th=P_CHI2P_TH, p_len_th=P_LEN_TH)


# ===== 2prong muX cut =====
stage_key = "2prong-muX"
stage_names.append(stage_key)

evt_df = cut_has_mu(evt_df)
evt_df = cut_mu_kinematics(evt_df, mu_Plo_th=MU_PLO_TH, mu_Phi_th=MU_PHI_TH)

df_dict[stage_key] = evt_df
# =====================================

# ===== 2prong mup cut =====
stage_key = "2prong-mup"
stage_names.append(stage_key)

evt_df = cut_has_p(evt_df)
evt_df = cut_p_kinematics(evt_df, p_Plo_th=P_PLO_TH, p_Phi_th=P_PHI_TH)

df_dict[stage_key] = evt_df
# =====================================

# save content to file
with open(args.save_name+".pkl", "wb") as f:
    pickle.dump(df_hist_contents, f)
print("DONE")