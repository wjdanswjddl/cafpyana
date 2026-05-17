# Event Selection

import pandas as pd
import numpy as np
import sys
from os import path, makedirs
from datetime import datetime
import pickle

# local imports
# sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana')  # repo root when cwd is not cafpyana
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

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--save_fig", action="store_true")
parser.add_argument("--save_content", action="store_true")
args = parser.parse_args()


# ===== run / save configs =====
show_plot = False

# plot breakdown type
plot_type = "topology"

# save plots
save_fig = args.save_fig
save_fig_base_dir = "/exp/sbnd/data/users/munjung/xsec/PLOTS/numuCC_1p0pi"
today_str = datetime.now().strftime("%Y%m%d")
save_fig_dir = path.join(save_fig_base_dir, f"event_selection-{today_str}")
if save_fig:
    if not path.exists(save_fig_dir):
        makedirs(save_fig_dir)
    print("saving plots in ", save_fig_dir)


# save contents for plots
save_content = args.save_content
save_content_base_dir = "/exp/sbnd/data/users/munjung/xsec/RESULTS/numuCC_1p0pi"
save_content_dir = path.join(save_content_base_dir, f"event_selection-{today_str}")
if save_content:
    if not path.exists(save_content_dir):
        makedirs(save_content_dir)
    print("saving nevts in ", save_content_dir)


# ===== load dfs =====
df_dir_mc = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_01_102644__sel_all-mc-BNB_cosmics"
filename_str_mc = "sel_all-mc-BNB_cosmics_100"
keys2load_mc = ['evt', 'trk', 'hit0', 'hit1', 'hit2', 'hdr']
df_mc = dfs_from_dir(search_dir=df_dir_mc, filename_str=filename_str_mc, keys2load=keys2load_mc, n_max_concat=999)

df_dir_data = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_01_102054__sel_all-data-BNB_cosmics"
filename_str_data = "sel_all-data-BNB_cosmics_10"
keys2load_data = ['evt', 'trk', 'hit0', 'hit1', 'hit2', 'hdr', 'bnbpot', 'trigger']
df_data = dfs_from_dir(search_dir=df_dir_data, filename_str=filename_str_data, keys2load=keys2load_data, n_max_concat=999)

df_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/lowE"
keys2load_lowE = ['hdr', 'evt', 'trk']
df_lowE = dfs_from_dir(search_dir=df_dir, filename_str="aa_all", keys2load=keys2load_lowE, n_max_concat=999)

df_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/OffBeam"
keys2load_offbeam = ['hdr', 'evt', 'trk']
df_offbeam = dfs_from_dir(search_dir=df_dir, filename_str="aa_all", keys2load=keys2load_offbeam, n_max_concat=999)

df_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/MC/intime"
keys2load_intime = ['hdr', 'evt', 'trk']
df_intime = dfs_from_dir(search_dir=df_dir, filename_str="aa_all", keys2load=keys2load_intime, n_max_concat=999)

mc_hdr_df, mc_evt_df, mc_trk_df = df_mc["hdr"], df_mc["evt"], df_mc["trk"]
data_hdr_df, data_evt_df, data_trk_df = df_data["hdr"], df_data["evt"], df_data["trk"]
dirt_hdr_df, dirt_evt_df, dirt_trk_df = df_lowE["hdr"], df_lowE["evt"], df_lowE["trk"]
offbeam_hdr_df, offbeam_evt_df, offbeam_trk_df = df_offbeam["hdr"], df_offbeam["evt"], df_offbeam["trk"]
intime_hdr_df, intime_evt_df, intime_trk_df = df_intime["hdr"], df_intime["evt"], df_intime["trk"]

# TODO: move to dfmaker
# ===== calculate avg chi2 for pid =====
chimu_avg = avg_chi2(mc_trk_df, "chi2_muon")
mc_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu_avg
chip_avg = avg_chi2(mc_trk_df, "chi2_proton")
mc_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip_avg

chimu_avg = avg_chi2(data_trk_df, "chi2_muon")
data_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu_avg
chip_avg = avg_chi2(data_trk_df, "chi2_proton")
data_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip_avg

chimu_avg = avg_chi2(intime_trk_df, "chi2_muon")
intime_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu_avg
chip_avg = avg_chi2(intime_trk_df, "chi2_proton")
intime_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip_avg

chimu_avg = avg_chi2(dirt_trk_df, "chi2_muon")
dirt_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_muon", "")] = chimu_avg
chip_avg = avg_chi2(dirt_trk_df, "chi2_proton")
dirt_trk_df[("pfp", "trk", "chi2pid", "avg", "chi2_proton", "")] = chip_avg

# ===== get mcs range diff =====
def get_mcs_range_diff(trks):
    mcs_range_diff = (trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon
    trks[("pfp", "trk", "mcs_range_diff", "", "", "")] = mcs_range_diff
    return trks

mc_trk_df, data_trk_df, intime_trk_df, offbeam_trk_df, dirt_trk_df = (
    get_mcs_range_diff(df)
    for df in (mc_trk_df, data_trk_df, intime_trk_df, offbeam_trk_df, dirt_trk_df)
)

# ===== exposure calculation =====
data_tot_pot = data_hdr_df['pot'].sum()
pot_str = get_pot_str(data_tot_pot)
data_evt_df["pot_weight"] = 1.0  # np.ones(len(data_evt_df))
data_trk_df["pot_weight"] = 1.0  # np.ones(len(data_trk_df))
data_gates = data_hdr_df.nbnbinfo.sum()
print(f"data_tot_pot: {data_tot_pot:.3e}", f"data tot gates : {data_gates:.3e}")

mc_tot_pot = mc_hdr_df['pot'].sum()
mc_pot_scale = data_tot_pot / mc_tot_pot if mc_tot_pot > 0 else 0
mc_evt_df["pot_weight"] = mc_pot_scale
mc_trk_df["pot_weight"] = mc_pot_scale
print(f"mc_tot_pot: {mc_tot_pot:.3e}", f"mc_pot_scale: {mc_pot_scale:.3e}")

dirt_tot_pot = dirt_hdr_df['pot'].sum()
dirt_pot_scale = data_tot_pot / dirt_tot_pot if dirt_tot_pot > 0 else 0
dirt_evt_df["pot_weight"] = dirt_pot_scale
dirt_trk_df["pot_weight"] = dirt_pot_scale
print(f"dirt_tot_pot: {dirt_tot_pot:.3e}", f"dirt_pot_scale: {dirt_pot_scale:.3e}")

# -- cosmics
f = 0.08
offbeam_gates = offbeam_hdr_df.loc[offbeam_hdr_df['first_in_subrun'] == 1, 'noffbeambnb'].sum()
scale_offbeam_to_lightdata = (1-f) * data_gates / offbeam_gates if offbeam_gates > 0 else 0
offbeam_evt_df["gates_weight"] = scale_offbeam_to_lightdata
offbeam_evt_df["pot_weight"] = scale_offbeam_to_lightdata
offbeam_trk_df["pot_weight"] = scale_offbeam_to_lightdata
print(f"offbeam cosmics data gates: {offbeam_gates:.2e}", f"goal scale: {scale_offbeam_to_lightdata:.2f}")

intime_gates = intime_hdr_df.loc[intime_hdr_df['first_in_subrun'] == 1, 'ngenevt'].sum()
scale_intime_to_lightdata = (1-f) * data_gates / intime_gates if intime_gates > 0 else 0
intime_evt_df["gates_weight"] = scale_intime_to_lightdata
intime_evt_df["pot_weight"] = scale_intime_to_lightdata
intime_trk_df["pot_weight"] = scale_intime_to_lightdata
print(f"intime cosmics data gates: {intime_gates:.2e}", f"goal scale: {scale_intime_to_lightdata:.2f}")

pot_str = get_pot_str(data_tot_pot)
plot_labels_bar = ["Events (POT={})".format(pot_str), "", ""]
plot_labels_hist = ["", "Events (POT={})".format(pot_str), ""]

mc_df     = mc_evt_df
data_df   = data_evt_df
intime_df = intime_evt_df
offbeam_df = offbeam_evt_df
dirt_df   = dirt_evt_df

# ===== save breakdowns for summary plots =====
breakdown_dict = {"topology": {}, "genie": {}}

# save dfs per stage, df_dict used for efficiency plot
stage_names = []
df_dict = {} 
df_dict_data = {}
df_dict_intime = {} 
df_dict_offbeam = {}
df_dict_dirt = {}

# save histogram contents
df_hist_contents = {}


# ===== all reconstructed slices =====
stage_key = "allreco"
stage_names.append(stage_key)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ===== not clear cosmic =====
stage_key = "is_clear_cosmic"
stage_names.append(stage_key)

mc_df, data_df, intime_df, dirt_df = (
    cut_clear_cosmic(df)
    for df in (mc_df, data_df, intime_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ===== vertex in fiducial volume =====
stage_key = "vertex_in_fv"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_vertex_in_fv(df, det=DETECTOR)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ****** nu-score plot ******
var_config = VariableConfig.nu_score()
plot_labels = ["Neutrino Score", "Events (POT={})".format(pot_str), ""]

save_name = save_fig_dir + "/selected-{}_{}.png".format(var_config.var_save_name, plot_type)
ret_hist = overlay_hists(plot_type,
                         mc_df=mc_df, 
                         data_df=data_df, 
                         intime_df=offbeam_df, 
                         dirt_df=dirt_df,
                         ratio=True,
                         ax_ylim_ratio=1.8,
                         vline=[[0.45, 1]],
                         var_config=var_config,
                         plot_labels=plot_labels,
                         save_fig=save_fig, 
                         save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ===== nu-score cut =====
stage_key = "nu_score"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_nu_score(df, NU_SCORE_TH)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ----- get valid tracks
mc_trk_df, data_trk_df, intime_trk_df, offbeam_trk_df, dirt_trk_df = (
    get_valid_trks(df)
    for df in (mc_trk_df, data_trk_df, intime_trk_df, offbeam_trk_df, dirt_trk_df)
)

mc_trk_df, data_trk_df, intime_trk_df, offbeam_trk_df, dirt_trk_df = (
    match_trkdf_to_slcdf(trk_df, df)
    for trk_df, df in ((mc_trk_df, mc_df), (data_trk_df, data_df), (intime_trk_df, intime_df), (offbeam_trk_df, offbeam_df), (dirt_trk_df, dirt_df))
)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    get_trk_info(df, trk_df, SAVE_NTRKS)
    for df, trk_df in ((mc_df, mc_trk_df), (data_df, data_trk_df), (intime_df, intime_trk_df), (offbeam_df, offbeam_trk_df), (dirt_df, dirt_trk_df))
)


# ****** n tracks plot ******
var_config = VariableConfig.n_trks()
plot_labels = ["Number of tracks", "Events (POT={})".format(pot_str), ""]

save_name = save_fig_dir + "/selected-{}_{}.png".format(var_config.var_save_name, plot_type)
ret_hist = overlay_hists(plot_type,
                                mc_df=mc_df, 
                                data_df=data_df, 
                                intime_df=offbeam_df, 
                                dirt_df=dirt_df,
                                var_config=var_config,
                                plot_labels=plot_labels,
                                save_fig=save_fig, 
                                save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ===== 2prong cut =====
stage_key = "2prong"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_2prong(df)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ===== 2prong contained =====
stage_key = "2prong-contained"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_2prong_contained(df, det=DETECTOR)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================


# ----- take both tracks for plotting
def get_both_trks(df):
    return pd.concat([df.trk1, df.trk2])

# ****** track score plot ******
var_config = VariableConfig.track_score()
plot_labels = [var_config.var_labels[0], "Tracks / Bin (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                         mc_df=get_both_trks(mc_df), 
                         data_df=get_both_trks(data_df), 
                         intime_df=get_both_trks(offbeam_df), 
                         dirt_df=get_both_trks(dirt_df),
                         ratio=True,
                         ax_ylim_ratio=1.8,
                         vline=[[0.5, 1]],
                         var_config=var_config,
                         plot_labels=plot_labels,
                         save_fig=save_fig, 
                         save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ===== 2prong track score cut =====
stage_key = "2prong-trackscore"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_2prong_trackscore(df, TRACKSCORE_TH)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================

# ****** vtx dist plot ******
var_config = VariableConfig.vtx_dist()

plot_labels = [var_config.var_labels[0], "Tracks / Bin (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                         mc_df=get_both_trks(mc_df), 
                         data_df=get_both_trks(data_df), 
                         intime_df=get_both_trks(offbeam_df), 
                         dirt_df=get_both_trks(dirt_df),
                         ratio=True,
                         vline=[[VTXDIST_TH, 0]],
                         var_config=var_config,
                         plot_labels=plot_labels,
                         save_fig=save_fig, 
                         save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ===== 2prong vtx dist cut =====
stage_key = "2prong-vtxdist"
stage_names.append(stage_key)

mc_df, data_df, intime_df, offbeam_df, dirt_df = (
    cut_2prong_vtxdist(df, VTXDIST_TH)
    for df in (mc_df, data_df, intime_df, offbeam_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, offbeam_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_offbeam, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================

# ****** trk len plot ******
var_config = VariableConfig.trk_len()
plot_labels = [var_config.var_labels[0], "Tracks / Bin (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                            mc_df=get_both_trks(mc_df), 
                            data_df=get_both_trks(data_df), 
                            intime_df=get_both_trks(offbeam_df), 
                            dirt_df=get_both_trks(dirt_df),
                            ratio=True,
                            vline=[[50, 1]],
                            var_config=var_config,
                            plot_labels=plot_labels,
                            save_fig=save_fig, 
                            save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ****** mcs range diff plot ******
var_config = VariableConfig.mcs_range_diff()
plot_labels = [var_config.var_labels[0], "Tracks / Bin (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                                mc_df=get_both_trks(mc_df), 
                                data_df=get_both_trks(data_df), 
                                intime_df=get_both_trks(offbeam_df), 
                                dirt_df=get_both_trks(dirt_df),
                                ratio=True,
                                vline=[[-QUAL_TH,0], [QUAL_TH,1]],
                                var_config=var_config,
                                plot_labels=plot_labels,
                                save_fig=save_fig, 
                                save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************


# ****** chi2_mu plot ******
var_config = VariableConfig.chi2_mu()
plot_labels = [var_config.var_labels[0], "Tracks/ Bin  (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                                mc_df=get_both_trks(mc_df), 
                                data_df=get_both_trks(data_df), 
                                intime_df=get_both_trks(offbeam_df), 
                                dirt_df=get_both_trks(dirt_df),
                                ratio=True,
                                vline=[[MU_CHI2MU_TH, 0]],
                                var_config=var_config,
                                plot_labels=plot_labels,
                                save_fig=save_fig, 
                                save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ****** chi2_proton plot ******
var_config = VariableConfig.chi2_proton()
plot_labels = [var_config.var_labels[0], "Tracks / Bin (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                                mc_df=get_both_trks(mc_df), 
                                data_df=get_both_trks(data_df), 
                                intime_df=get_both_trks(offbeam_df), 
                                dirt_df=get_both_trks(dirt_df),
                                ratio=True,
                                vline=[[MU_CHI2P_TH, 1]],
                                ax_ylim_ratio=1.8,
                                var_config=var_config,
                                plot_labels=plot_labels,
                                save_fig=save_fig, 
                                save_name=save_name)
df_hist_contents["_".join([stage_key, var_config.var_save_name])] = ret_hist
# *********************************

# ****** chi2_mu vs. chi2_proton plot ******

# MC
var_config_1 = VariableConfig.chi2_mu()
var_config_2 = VariableConfig.chi2_proton()

n, bins1, bins2, _ =plt.hist2d(mc_trk_df[var_config_1.var_evt_reco_col], 
           mc_trk_df[var_config_2.var_evt_reco_col], 
           bins=[var_config_1.bins, var_config_2.bins],
           cmap="GnBu", norm=LogNorm())
plt.colorbar(label="Tracks")
plt.xlabel("Muon $\\chi^2$")
plt.ylabel("Proton $\\chi^2$")
plt.text(0.025, 1.07, r"$\mathbf{SBND}$ Preliminary    $\mathbf{SBND}$ Simulation", transform=plt.gca().transAxes, fontsize=14, color='gray', ha='left', va='top')

df_hist_contents["_".join([stage_key, "chi2","2D","MC"])] = {"n": n, "bins1": bins1, "bins2": bins2}

if save_fig:
    save_name = save_fig_dir + "/chi2muon_vs_chi2proton-MC.pdf"
    plt.savefig(save_name, bbox_inches="tight")
plt.show()

# Data
n, bins, _, _ = plt.hist2d(data_trk_df[var_config_1.var_evt_reco_col], 
           data_trk_df[var_config_2.var_evt_reco_col], 
           bins=[var_config_1.bins, var_config_2.bins],
           cmap="GnBu", norm=LogNorm())
plt.colorbar(label="Tracks")
plt.xlabel("Muon $\\chi^2$")
plt.ylabel("Proton $\\chi^2$")
plt.text(0.025, 1.07, r"$\mathbf{SBND}$ Preliminary    $\mathbf{SBND}$ Data", transform=plt.gca().transAxes, fontsize=14, color='gray', ha='left', va='top')

df_hist_contents["_".join([stage_key, "chi2","2D","Data"])] = {"n": n, "bins1": bins1, "bins2": bins2}

if save_fig:
    save_name = save_fig_dir + "/chi2muon_vs_chi2proton-data.pdf"
    plt.savefig(save_name, bbox_inches="tight")
plt.show()

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
    nlevels = len(mc_df.index.names)
    mcs_range_diff = np.abs((trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon)
    chimu_avg = trks.pfp.trk.chi2pid.avg.chi2_muon
    chip_avg = trks.pfp.trk.chi2pid.avg.chi2_proton
    mu_cut = (chimu_avg > 0) & (chimu_avg < MU_CHI2MU_TH) & \
            (chip_avg > MU_CHI2P_TH) & \
            (trks.pfp.trk.len > MU_LEN_TH) & \
            (mcs_range_diff < QUAL_TH)
    not_mu_candidate = pd.concat([trks[~mu_cut], trks[mu_cut].groupby(level=list(range(nlevels))).nth(1)])
    return not_mu_candidate


mc_mus = is_mu_candidate(mc_trk_df)
intime_mus = is_mu_candidate(intime_trk_df)
dirt_mus = is_mu_candidate(dirt_trk_df)
data_mus = is_mu_candidate(data_trk_df)

mc_notmus = is_not_mu_candidate(mc_trk_df)
data_notmus = is_not_mu_candidate(data_trk_df)
intime_notmus = is_not_mu_candidate(intime_trk_df)
offbeam_notmus = is_not_mu_candidate(offbeam_trk_df)
dirt_notmus = is_not_mu_candidate(dirt_trk_df)

# ****** chi2_mu plot ******
var_config = VariableConfig.chi2_mu()
plot_labels = [var_config.var_labels[0], "Events (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                            mc_df=mc_notmus, 
                            data_df=data_notmus, 
                            intime_df=offbeam_notmus, 
                            dirt_df=dirt_notmus,
                            ratio=True,
                            # vline=[MU_CHI2MU_TH],
                            var_config=var_config,
                            plot_labels=plot_labels,
                            save_fig=save_fig, 
                            save_name=save_name)

df_hist_contents["_".join(["not_mu", var_config.var_save_name])] = ret_hist
# *********************************

# ****** chi2_proton plot ******
var_config = VariableConfig.chi2_proton()
plot_labels = [var_config.var_labels[0], "Events (POT={})".format(pot_str), ""]
save_name = save_fig_dir + "/2prong-{}_{}.png".format(var_config.var_save_name, "pdg")
ret_hist = overlay_hists("pdg",
                          mc_df=mc_notmus, 
                          data_df=data_notmus, 
                          intime_df=offbeam_notmus, 
                          dirt_df=dirt_notmus,
                          ratio=True,
                          vline=[[MU_CHI2P_TH, 1]],
                          var_config=var_config,
                          plot_labels=plot_labels,
                          save_fig=save_fig, 
                          save_name=save_name)

df_hist_contents["_".join(["not_mu", var_config.var_save_name])] = ret_hist
# *********************************

# ----- get mu & p candidates
mc_df, data_df, intime_df, dirt_df = (
    get_mu_p_candidate(df, 
                       mu_chi2mu_th=MU_CHI2MU_TH, mu_chi2p_th=MU_CHI2P_TH, mu_len_th=MU_LEN_TH, qual_th=QUAL_TH,
                       p_chi2mu_th=-1, p_chi2p_th=P_CHI2P_TH, p_len_th=P_LEN_TH)
    for df in (mc_df, data_df, intime_df, dirt_df)
)

# ===== 2prong muX cut =====
stage_key = "2prong-muX"
stage_names.append(stage_key)

mc_df, data_df, intime_df, dirt_df = (
    cut_has_mu(df)
    for df in (mc_df, data_df, intime_df, dirt_df)
)
mc_df, data_df, intime_df, dirt_df = (
    cut_mu_kinematics(df, mu_Plo_th=MU_PLO_TH, mu_Phi_th=MU_PHI_TH)
    for df in (mc_df, data_df, intime_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================

# ===== 2prong mup cut =====
stage_key = "2prong-mup"
stage_names.append(stage_key)

mc_df, data_df, intime_df, dirt_df = (
    cut_has_p(df)
    for df in (mc_df, data_df, intime_df, dirt_df)
)
mc_df, data_df, intime_df, dirt_df = (
    cut_p_kinematics(df, p_Plo_th=P_PLO_TH, p_Phi_th=P_PHI_TH)
    for df in (mc_df, data_df, intime_df, dirt_df)
)

for d, d_dict in zip(
    (mc_df, data_df, intime_df, dirt_df),
    (df_dict, df_dict_data, df_dict_intime, df_dict_dirt)
    ):
    d_dict[stage_key] = d
# =====================================

print("DONE")