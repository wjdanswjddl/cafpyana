
from os import path, makedirs
from datetime import datetime
from functools import partial
import pickle
import argparse

import numpy as np
import pandas as pd

import sys
# sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana') # absolute path for running on EAF
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.files_config import *
plt.style.use("presentation.mplstyle")

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

parser = argparse.ArgumentParser(description='Event selection settings')
parser.add_argument('--chunk_idx', type=int, default=0, help='Chunk index for data')
parser.add_argument('--n_time_splits', type=int, default=15, help='Number of time splits for data')
args = parser.parse_args()

print("Processing n_time_splits: ", args.n_time_splits)
print("Processing chunk_idx: ", args.chunk_idx)


def get_syst_unc(var_config):
    date_str = "20260220"
    mcstat_syst = np.load(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-{date_str}/mcstat_syst_dict.npz", allow_pickle=True)
    mcstat_syst =dict(mcstat_syst)[var_config.var_save_name].item()['MCstat']['cov_frac']

    g4_syst = np.load(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-{date_str}/g4_syst_dict.npz", allow_pickle=True)
    g4_syst =dict(g4_syst)[var_config.var_save_name].item()['G4']['cov_frac']

    flux_syst = np.load(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-{date_str}/flux_syst_dict.npz", allow_pickle=True)
    flux_syst =dict(flux_syst)[var_config.var_save_name].item()['flux']['cov_frac']

    date_str = "20260222"
    cosmics_syst = np.load(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-{date_str}/cosmics_syst_dict.npz", allow_pickle=True)
    cosmics_syst =dict(cosmics_syst)[var_config.var_save_name].item()['Cosmics']['cov_frac']

    date_str = "20260219"
    genie_syst = pickle.load(open(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/cov_mat_dict-{date_str}.pkl", "rb"))
    genie_syst = genie_syst[var_config.var_save_name]['genie'] 

    # detvar_syst = pickle.load(open("/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10/nevts/det_unc_dict-20260216.pkl", "rb"))
    # detvar_syst = detvar_syst[var_config.var_save_name]['detvar']

    # detvar_syst = pickle.load(open("/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10/nevts/det_unc_dict-20260216.pkl", "rb"))
    # detvar_syst = np.sqrt(detvar_syst[var_config.var_save_name]['ccal']**2 \
    #     + detvar_syst[var_config.var_save_name]['alpha']**2 \
    #     + detvar_syst[var_config.var_save_name]['beta']**2 \
    #     + detvar_syst[var_config.var_save_name]['R']**2) / 2.


    # flat uncertainties
    pot_frac_unc = 0.02
    ntargets_frac_unc = 0.01

    # flat uncertainties
    frac_uncert_total = np.zeros(len(var_config.bin_centers))
    cov_total = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))
    systs      = [mcstat_syst, genie_syst, flux_syst, g4_syst, cosmics_syst]
    for syst in systs:
        cov_total += syst
        syst_uncert = np.sqrt(np.diag(syst))
        frac_uncert_total += syst_uncert ** 2

    flat_systs = [pot_frac_unc, ntargets_frac_unc]
    for syst in flat_systs:
        cov_total += np.diag(syst * np.ones(len(var_config.bin_centers)) ** 2)
        syst_uncert = syst * np.ones(len(var_config.bin_centers))
        frac_uncert_total += syst_uncert ** 2

    # frac_uncert_total += detvar_syst ** 2

    frac_uncert_total = np.sqrt(frac_uncert_total)
    # syst = frac_uncert_total

    return cov_total, frac_uncert_total


save_fig = True

today_str = datetime.now().strftime("%Y%m%d")
save_fig_dir = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20/chunk{args.chunk_idx}")
save_fig_dir_perTPC = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-perTPC/chunk{args.chunk_idx}")
save_fig_dir_TPC1 = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-TPC1/chunk{args.chunk_idx}")
save_fig_dir_TPC2 = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-TPC2/chunk{args.chunk_idx}")
save_fig_dir_fwd = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-fwd_muons/chunk{args.chunk_idx}")
save_fig_dir_bwd = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-bwd_muons/chunk{args.chunk_idx}")
save_fig_dir_crosser_fwd = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-crosser_muons_fwd/chunk{args.chunk_idx}")
save_fig_dir_crosser_bwd = path.join(save_fig_base_dir, f"selected_events-data-2prong-1e20-crosser_muons_bwd/chunk{args.chunk_idx}")

if save_fig:
    if not path.exists(save_fig_dir):
        makedirs(save_fig_dir)
    print("saving plots in ", save_fig_dir)

    if not path.exists(save_fig_dir_fwd):
        makedirs(save_fig_dir_fwd)
    print("saving plots in ", save_fig_dir_fwd)

    if not path.exists(save_fig_dir_bwd):
        makedirs(save_fig_dir_bwd)
    print("saving plots in ", save_fig_dir_bwd)

    if not path.exists(save_fig_dir_crosser_fwd):
        makedirs(save_fig_dir_crosser_fwd)
    print("saving plots in ", save_fig_dir_crosser_fwd)

    if not path.exists(save_fig_dir_crosser_bwd):
        makedirs(save_fig_dir_crosser_bwd)
    print("saving plots in ", save_fig_dir_crosser_bwd)

    if not path.exists(save_fig_dir_perTPC):
        makedirs(save_fig_dir_perTPC)
    print("saving plots in ", save_fig_dir_perTPC)

    if not path.exists(save_fig_dir_TPC1):
        makedirs(save_fig_dir_TPC1)
    print("saving plots in ", save_fig_dir_TPC1)

    if not path.exists(save_fig_dir_TPC2):
        makedirs(save_fig_dir_TPC2)
    print("saving plots in ", save_fig_dir_TPC2)

# dfs = get_ana_dfs(option="selected_events")
# mc_evt_df     = dfs["mc"]
# mc_hdr_df     = dfs["mc_hdr"]
# data_evt_df_   = dfs["data"]
# data_hdr_df_   = dfs["data_hdr"]
# intime_evt_df = dfs["intime"]
# intime_hdr_df = dfs["intime_hdr"]
# dirt_evt_df   = dfs["dirt"]
# dirt_hdr_df   = dfs["dirt_hdr"]
# pot_label     = dfs["pot_label"]


from pyanalib.split_df_helpers_new import *

df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_05_164044__sel_2prong-data-BNB_cosmics"
keys2load_data = ['hdr', 'evt']
df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_2prong-data-BNB_cosmics", keys2load=keys2load_data, n_max_concat=999)
data_evt_df_ = df_data['evt']
data_hdr_df_ = df_data['hdr']

df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_05_233424__sel_2prong-data-OffBeamLight"
keys2load_data = ['hdr', 'evt']
df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_2prong-data-OffBeamLight", keys2load=keys2load_data, n_max_concat=999)
intime_evt_df = df_data['evt']
intime_hdr_df = df_data['hdr']

df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_05_231714__sel_2prong-mc-BNB_cosmics"
keys2load_data = ['hdr', 'evt']
df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_2prong-mc-BNB_cosmics", keys2load=keys2load_data, n_max_concat=999)
mc_evt_df = df_data['evt']
mc_hdr_df = df_data['hdr']


mc_evt_df[('trk1', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['trk1', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['trk1', 'pfp', 'trk', 'dir', 'y', '', '']))
mc_evt_df[('trk2', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['trk2', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['trk2', 'pfp', 'trk', 'dir', 'y', '', '']))

data_evt_df_['trk1', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['trk1', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['trk1', 'pfp', 'trk', 'dir', 'y', '', '']))
data_evt_df_['trk2', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['trk2', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['trk2', 'pfp', 'trk', 'dir', 'y', '', '']))

intime_evt_df[('trk1', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['trk1', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['trk1', 'pfp', 'trk', 'dir', 'y', '', '']))
intime_evt_df[('trk2', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['trk2', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['trk2', 'pfp', 'trk', 'dir', 'y', '', '']))

# dirt_evt_df[('trk1', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(dirt_evt_df['trk1', 'pfp', 'trk', 'dir', 'x', '', ''], dirt_evt_df['trk1', 'pfp', 'trk', 'dir', 'y', '', '']))
# dirt_evt_df[('p', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(dirt_evt_df['p', 'pfp', 'trk', 'dir', 'x', '', ''], dirt_evt_df['p', 'pfp', 'trk', 'dir', 'y', '', '']))

# TODO: for the integrated plot to work
mc_evt_df.loc[mc_evt_df.mc.iscc.isna(), ("mc","iscc")] = 999
data_evt_df_["mc", "iscc"] = 999


# Split data_hdr into 20 time-ordered chunks (sorted by run, then evt within run).
# np.array_split sizes differ by at most 1 when len % 20 != 0.
_n_time_splits = args.n_time_splits
_sorted = data_hdr_df_.sort_values(["run", "evt"], kind="mergesort")
data_hdr_df_splits = [
    _sorted.iloc[idx]
    for idx in np.array_split(np.arange(len(_sorted)), _n_time_splits)
]
for i in range(_n_time_splits):
    data_hdr_df = data_hdr_df_splits[i]


chunk_idx = args.chunk_idx
data_hdr_df = data_hdr_df_splits[chunk_idx]
evt_idxs = data_evt_df_.index
common_idxs = data_evt_df_.reset_index(level=[2]).index.intersection(data_hdr_df.index)
data_evt_df = data_evt_df_.reset_index(level=[2]).loc[common_idxs].reset_index().set_index(["__ntuple", "entry", "rec.slc..index"])

data_tot_pot = data_hdr_df['pot'].sum()
data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
print("data_tot_pot: %.3e" %(data_tot_pot))
pot_str = get_pot_str(data_tot_pot)
pot_label = f"Events / Bin (POT={pot_str})"
data_gates = data_hdr_df.nbnbinfo.sum()
print("data tot gates : %.3e" %(data_gates))

mc_tot_pot = mc_hdr_df['pot'].sum()
mc_pot_scale = data_tot_pot / mc_tot_pot
print("mc_pot_scale: %.3e" %(mc_pot_scale))
mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))

intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
f = 0.0753
print("DATA GATES: ", data_gates)
print("INTIME GATES: ", intime_gates)
scale_intime_to_lightdata = (1-f)*data_gates/intime_gates
print("intime data scale: {:.2f}".format(scale_intime_to_lightdata))
intime_evt_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))


# ==== cuts ====
def perTPC_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=0) & InFV(df["trk1"].pfp.trk.end, det="SBND_TPC1", incathode=0) & InFV(df["trk2"].pfp.trk.end, det="SBND_TPC1", incathode=0)
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=0) & InFV(df["trk1"].pfp.trk.end, det="SBND_TPC2", incathode=0) & InFV(df["trk2"].pfp.trk.end, det="SBND_TPC2", incathode=0)
    perTPC_cut = in_TPC1_cut | in_TPC2_cut
    return perTPC_cut

mc_evt_df_perTPC = mc_evt_df[perTPC_cut(mc_evt_df)]
data_evt_df_perTPC = data_evt_df[perTPC_cut(data_evt_df)]
intime_evt_df_perTPC = intime_evt_df[perTPC_cut(intime_evt_df)]

def inTPC1_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=0) # & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=2.5) & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=2.5)
    return in_TPC1_cut

mc_evt_df_TPC1 = mc_evt_df[inTPC1_cut(mc_evt_df)]
data_evt_df_TPC1 = data_evt_df[inTPC1_cut(data_evt_df)]
intime_evt_df_TPC1 = intime_evt_df[inTPC1_cut(intime_evt_df)]


def inTPC2_cut(df):
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=0) #& InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=2.5) & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=2.5)
    return in_TPC2_cut

mc_evt_df_TPC2 = mc_evt_df[inTPC2_cut(mc_evt_df)]
data_evt_df_TPC2 = data_evt_df[inTPC2_cut(data_evt_df)]
intime_evt_df_TPC2 = intime_evt_df[inTPC2_cut(intime_evt_df)]


def fwd_muons_cut(df):
    cut = df["trk1"].pfp.trk.dir.z > 0
    return cut

mc_evt_df_fwd_muons = mc_evt_df[fwd_muons_cut(mc_evt_df)]
data_evt_df_fwd_muons = data_evt_df[fwd_muons_cut(data_evt_df)]
intime_evt_df_fwd_muons = intime_evt_df[fwd_muons_cut(intime_evt_df)]

def bwd_muons_cut(df):
    cut = df["trk1"].pfp.trk.dir.z < 0
    return cut

mc_evt_df_bwd_muons = mc_evt_df[bwd_muons_cut(mc_evt_df)]
data_evt_df_bwd_muons = data_evt_df[bwd_muons_cut(data_evt_df)]
intime_evt_df_bwd_muons = intime_evt_df[bwd_muons_cut(intime_evt_df)]

def crosser_muons_cut(df):
    cut = df["trk1"].pfp.trk.end.x * df["trk1"].pfp.trk.start.x < 0
    return cut

mc_evt_df_crosser_muons_fwd = mc_evt_df[crosser_muons_cut(mc_evt_df) & fwd_muons_cut(mc_evt_df)]
data_evt_df_crosser_muons_fwd = data_evt_df[crosser_muons_cut(data_evt_df) & fwd_muons_cut(data_evt_df)]
intime_evt_df_crosser_muons_fwd = intime_evt_df[crosser_muons_cut(intime_evt_df) & fwd_muons_cut(intime_evt_df)]

mc_evt_df_crosser_muons_bwd = mc_evt_df[crosser_muons_cut(mc_evt_df) & bwd_muons_cut(mc_evt_df)]
data_evt_df_crosser_muons_bwd = data_evt_df[crosser_muons_cut(data_evt_df) & bwd_muons_cut(data_evt_df)]
intime_evt_df_crosser_muons_bwd = intime_evt_df[crosser_muons_cut(intime_evt_df) & bwd_muons_cut(intime_evt_df)]


# ===- plotter per cut ====

eps = 1e-8
ratio = True
approval = "internal"
textloc = [0.03, 0.55]
ax_ylim_ratio = 1.9

data_vs_mc_plotter = partial(
    overlay_hists,
    mc_df=mc_evt_df,
    data_df=data_evt_df,
    intime_df=intime_evt_df,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_perTPC = partial(
    overlay_hists,
    mc_df=mc_evt_df_perTPC,
    data_df=data_evt_df_perTPC,
    intime_df=intime_evt_df_perTPC,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_TPC1 = partial(
    overlay_hists,
    mc_df=mc_evt_df_TPC1,
    data_df=data_evt_df_TPC1,
    intime_df=intime_evt_df_TPC1,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_TPC2 = partial(
    overlay_hists,
    mc_df=mc_evt_df_TPC2,
    data_df=data_evt_df_TPC2,
    intime_df=intime_evt_df_TPC2,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_fwd_muons = partial(
    overlay_hists,
    mc_df=mc_evt_df_fwd_muons,
    data_df=data_evt_df_fwd_muons,
    intime_df=intime_evt_df_fwd_muons,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_bwd_muons = partial(
    overlay_hists,
    mc_df=mc_evt_df_bwd_muons,
    data_df=data_evt_df_bwd_muons,
    intime_df=intime_evt_df_bwd_muons,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)


data_vs_mc_plotter_crosser_muons_fwd = partial(
    overlay_hists,
    mc_df=mc_evt_df_crosser_muons_fwd,
    data_df=data_evt_df_crosser_muons_fwd,
    intime_df=intime_evt_df_crosser_muons_fwd,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    plot=False,
    save_fig=save_fig, 
)

data_vs_mc_plotter_crosser_muons_bwd = partial(
    overlay_hists,
    mc_df=mc_evt_df_crosser_muons_bwd,
    data_df=data_evt_df_crosser_muons_bwd,
    intime_df=intime_evt_df_crosser_muons_bwd,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    plot=False,
    save_fig=save_fig, 
)


# ==== plots ====

# approved vars
var_configs = [
    VariableConfig.all_events(),
    # VariableConfig.muon_momentum(),
    # VariableConfig.muon_direction(),
    # VariableConfig.proton_momentum(),
    # VariableConfig.proton_direction(),
    #VariableConfig.tki_del_Tp(),
    # VariableConfig.tki_del_Tp_x(),
    # VariableConfig.tki_del_Tp_y(),
    # VariableConfig.tki_del_p(),
    # VariableConfig.tki_del_alpha(),
    # VariableConfig.tki_del_phi()
    ]

for var_config in var_configs:
#
   cov, frac_uncert = get_syst_unc(var_config)
#    plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
#
   for breakdown_type in ["topology"]:
#        ret = data_vs_mc_plotter(breakdown_type=breakdown_type,
#                                var_config=var_config,
#                                plot_labels=plot_labels_hist,
#                                syst=cov,
#                                textchi2=True,
#                                save_name=path.join(save_fig_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)))
#
#        ret = data_vs_mc_plotter_perTPC(breakdown_type=breakdown_type,
#                                var_config=var_config,
#                                plot_labels=plot_labels_hist,
#                                syst=cov,
#                                textchi2=True,
#                                save_name=path.join(save_fig_dir_perTPC, "{}_{}".format(var_config.var_save_name, breakdown_type)))
#
#        ret = data_vs_mc_plotter_TPC1(breakdown_type=breakdown_type,
#                                var_config=var_config,
#                                plot_labels=plot_labels_hist,
#                                syst=cov,
#                                textchi2=True,
#                                save_name=path.join(save_fig_dir_TPC1, "{}_{}".format(var_config.var_save_name, breakdown_type)))
#
#        ret = data_vs_mc_plotter_TPC2(breakdown_type=breakdown_type,
#                                var_config=var_config,
#                                plot_labels=plot_labels_hist,
#                                syst=cov,
#                                textchi2=True,
#                                save_name=path.join(save_fig_dir_TPC2, "{}_{}".format(var_config.var_save_name, breakdown_type)))
#
#        # ret = data_vs_mc_plotter(breakdown_type=breakdown_type,
#        #                         var_config=var_config,
#        #                         plot_labels=plot_labels_hist,
#        #                         syst=cov,
#        #                         syst_decomp=True,
#        #                         textchi2=True,
#        #                         save_name=path.join(save_fig_dir, "{}_{}-syst-decomp".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, "Forward Muons"]
        # ret = data_vs_mc_plotter_fwd_muons(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_fwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, "Backward Muons"]
        # ret = data_vs_mc_plotter_bwd_muons(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_bwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))


        plot_labels_hist = [var_config.var_labels[1], pot_label, "Crosser Muons (Forward)"]
        ret = data_vs_mc_plotter_crosser_muons_fwd(breakdown_type=breakdown_type,
                                var_config=var_config,
                                plot_labels=plot_labels_hist,
                                syst=cov,
                                textchi2=True,
                                save_name=path.join(save_fig_dir_crosser_fwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        plot_labels_hist = [var_config.var_labels[1], pot_label, "Crosser Muons (Backward)"]
        ret = data_vs_mc_plotter_crosser_muons_bwd(breakdown_type=breakdown_type,
                                var_config=var_config,
                                plot_labels=plot_labels_hist,
                                syst=cov,
                                textchi2=True,
                                save_name=path.join(save_fig_dir_crosser_bwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))


# more vars

var_configs = [
    # VariableConfig.muon_direction_phi(),
    # VariableConfig.proton_direction_phi(),
    VariableConfig.trk1_direction_phi(),
    VariableConfig.trk2_direction_phi(),
    # VariableConfig.vertex_x(),
    # VariableConfig.muon_end_x(),
    # VariableConfig.muon_direction_x(),
    # VariableConfig.muon_direction_y(),
    # VariableConfig.proton_direction_x(),
    # VariableConfig.proton_direction_y(),
    # VariableConfig.vertex_y(),
    # VariableConfig.vertex_z(),
    # VariableConfig.muon_end_y(),
    # VariableConfig.muon_end_z()
    ]

for var_config in var_configs:


    for breakdown_type in ["topology"]:
    # for breakdown_type in ["topology", "genie", "genie_sb"]:
        frac_unc, cov = get_frac_unc(mc_evt_df, intime_evt_df, intime_evt_df, var_config)

        # plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
        # ret = data_vs_mc_plotter(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
        # frac_unc, cov = get_frac_unc(mc_evt_df_perTPC, intime_evt_df_perTPC, intime_evt_df_perTPC, var_config)
        # ret = data_vs_mc_plotter_perTPC(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_perTPC, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
        # frac_unc, cov = get_frac_unc(mc_evt_df_TPC1, intime_evt_df_TPC1, intime_evt_df_TPC1, var_config)
        # ret = data_vs_mc_plotter_TPC1(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_TPC1, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
        # frac_unc, cov = get_frac_unc(mc_evt_df_TPC2, intime_evt_df_TPC2, intime_evt_df_TPC2, var_config)
        # ret = data_vs_mc_plotter_TPC2(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_TPC2, "{}_{}".format(var_config.var_save_name, breakdown_type)))


        # plot_labels_hist = [var_config.var_labels[1], pot_label, "Forward Muons"]
        # ret = data_vs_mc_plotter_fwd_muons(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_fwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        # plot_labels_hist = [var_config.var_labels[1], pot_label, "Backward Muons"]
        # ret = data_vs_mc_plotter_bwd_muons(breakdown_type=breakdown_type,
        #                         var_config=var_config,
        #                         plot_labels=plot_labels_hist,
        #                         syst=cov,
        #                         textchi2=True,
        #                         save_name=path.join(save_fig_dir_bwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        plot_labels_hist = [var_config.var_labels[1], pot_label, "Crosser Muons (Forward)"]
        ret = data_vs_mc_plotter_crosser_muons_fwd(breakdown_type=breakdown_type,
                                var_config=var_config,
                                plot_labels=plot_labels_hist,
                                syst=cov,
                                textchi2=True,
                                save_name=path.join(save_fig_dir_crosser_fwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))

        plot_labels_hist = [var_config.var_labels[1], pot_label, "Crosser Muons (Backward)"]
        ret = data_vs_mc_plotter_crosser_muons_bwd(breakdown_type=breakdown_type,
                                var_config=var_config,
                                plot_labels=plot_labels_hist,
                                syst=cov,
                                textchi2=True,
                                save_name=path.join(save_fig_dir_crosser_bwd, "{}_{}".format(var_config.var_save_name, breakdown_type)))
