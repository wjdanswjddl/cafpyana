
from os import path, makedirs

_SYST_RESULTS_DIR = None
_GENIE_COV_PKL = None
from datetime import datetime
from functools import partial
import pickle
import argparse
import json

import numpy as np
import pandas as pd

import sys
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana')  # repo root when cwd is not cafpyana
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    FINAL_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.files_config import *
plt.style.use("presentation.mplstyle")

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser(description='Event selection settings')
parser.add_argument(
    '--chunk_idx',
    type=int,
    default=0,
    help='Legacy name: exposure-batch index for data (time-ordered slice; see exposure_access).',
)
parser.add_argument(
    '--exposure-batch-index',
    type=int,
    dest='chunk_idx',
    default=argparse.SUPPRESS,
    help='Preferred alias for --chunk_idx (one exposure batch).',
)
parser.add_argument(
    '--n_time_splits',
    type=int,
    default=15,
    help='Number of time-ordered exposure batches to split data into',
)
parser.add_argument(
    '--do_octant_plots',
    action='store_true',
    help='If set, make SBND octant diagnostic + per-octant overlay plots (default: off).',
)
parser.add_argument(
    '--do_quadrant_plots',
    action='store_true',
    help='If set, make SBND quadrant diagnostic + per-quadrant overlay plots (E/W and Top/Bottom; N/S combined). Default: off.',
)
parser.add_argument(
    '--syst-results-dir',
    type=str,
    default=None,
    help='Folder with mcstat_syst_dict.npz, g4_syst_dict.npz, flux_syst_dict.npz, cosmics_syst_dict.npz '
    '(e.g. syst_multisim_aggregate.py / get_systematics_multisim.py / get_systematics_mcstat_flux_g4.py output). '
    'Default: hardcoded dated paths.',
)
parser.add_argument(
    '--genie-cov-pkl',
    type=str,
    default=None,
    help='Override path to GENIE covariance pickle (cov_mat_dict-*.pkl).',
)
args = parser.parse_args()
_SYST_RESULTS_DIR = args.syst_results_dir
_GENIE_COV_PKL = args.genie_cov_pkl

print("Processing n_time_splits: ", args.n_time_splits)
print("Processing chunk_idx: ", args.chunk_idx)
print("Processing do_octant_plots: ", args.do_octant_plots)
print("Processing do_quadrant_plots: ", args.do_quadrant_plots)


def _zero_cov(var_config):
    n = len(var_config.bin_centers)
    return np.zeros((n, n))


def get_syst_unc(var_config):
    plots_base = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi"

    def _cov_frac_from_npz(npz_obj, var_sn, inner_key):
        z = dict(npz_obj)
        if var_sn not in z:
            return _zero_cov(var_config)
        ret = z[var_sn].item()[inner_key]
        return ret["cov_frac"]

    if _SYST_RESULTS_DIR is not None:
        bundle = _SYST_RESULTS_DIR
        mcstat_npz = np.load(path.join(bundle, "mcstat_syst_dict.npz"), allow_pickle=True)
        g4_npz = np.load(path.join(bundle, "g4_syst_dict.npz"), allow_pickle=True)
        flux_npz = np.load(path.join(bundle, "flux_syst_dict.npz"), allow_pickle=True)
        cosmics_npz = np.load(path.join(bundle, "cosmics_syst_dict.npz"), allow_pickle=True)
    else:
        date_str = "20260220"
        mcstat_npz = np.load(path.join(plots_base, f"systematics-{date_str}", "mcstat_syst_dict.npz"), allow_pickle=True)
        g4_npz = np.load(path.join(plots_base, f"systematics-{date_str}", "g4_syst_dict.npz"), allow_pickle=True)
        flux_npz = np.load(path.join(plots_base, f"systematics-{date_str}", "flux_syst_dict.npz"), allow_pickle=True)
        date_str = "20260222"
        cosmics_npz = np.load(path.join(plots_base, f"systematics-{date_str}", "cosmics_syst_dict.npz"), allow_pickle=True)

    vn = var_config.var_save_name
    mcstat_syst = _cov_frac_from_npz(mcstat_npz, vn, "MCstat")
    g4_syst = _cov_frac_from_npz(g4_npz, vn, "G4")
    flux_syst = _cov_frac_from_npz(flux_npz, vn, "flux")
    cosmics_syst = _cov_frac_from_npz(cosmics_npz, vn, "Cosmics")

    if _GENIE_COV_PKL is not None:
        genie_path = _GENIE_COV_PKL
    else:
        genie_path = path.join(plots_base, "cov_mat_dict-20260219.pkl")
    genie_blob = pickle.load(open(genie_path, "rb"))
    try:
        genie_syst = genie_blob[var_config.var_save_name]["genie"]
    except KeyError:
        genie_syst = _zero_cov(var_config) 

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

print("[save_fig] is set to ", save_fig)
today_str = datetime.now().strftime("%Y%m%d")
save_fig_dir = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}/chunk{args.chunk_idx}")
save_fig_dir_perTPC = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-perTPC/chunk{args.chunk_idx}")
save_fig_dir_inTPC1 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-contTPC1/chunk{args.chunk_idx}")
save_fig_dir_inTPC2 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-contTPC2/chunk{args.chunk_idx}")
# save_fig_dir_TPC1 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20{today_str}-TPC1/chunk{args.chunk_idx}")
# save_fig_dir_TPC2 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-TPC2/chunk{args.chunk_idx}")
save_fig_dir_fwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-fwd_muons/chunk{args.chunk_idx}")
save_fig_dir_bwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-bwd_muons/chunk{args.chunk_idx}")
save_fig_dir_crosser_fwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-crosser_muons_fwd/chunk{args.chunk_idx}")
save_fig_dir_crosser_bwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-{today_str}-crosser_muons_bwd/chunk{args.chunk_idx}")

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

    if not path.exists(save_fig_dir_inTPC1):
        makedirs(save_fig_dir_inTPC1)
    print("saving plots in ", save_fig_dir_inTPC1)

    if not path.exists(save_fig_dir_inTPC2):
        makedirs(save_fig_dir_inTPC2)
    print("saving plots in ", save_fig_dir_inTPC2)

dfs = get_ana_dfs(option="selected_events")
mc_evt_df     = dfs["mc"]
mc_hdr_df     = dfs["mc_hdr"]
data_evt_df_   = dfs["data"]
data_hdr_df_   = dfs["data_hdr"]
intime_evt_df = dfs["intime"]
intime_hdr_df = dfs["intime_hdr"]
dirt_evt_df   = dfs["dirt"]
dirt_hdr_df   = dfs["dirt_hdr"]
pot_label     = dfs["pot_label"]


# from pyanalib.split_df_helpers_new import *

# df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_05_164044__sel_mup-data-BNB_cosmics"
# keys2load_data = ['hdr', 'evt']
# df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_mup-data-BNB_cosmics", keys2load=keys2load_data, n_max_concat=999)
# data_evt_df_ = df_data['evt']
# data_hdr_df_ = df_data['hdr']

# df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_05_233424__sel_mup-data-OffBeamLight"
# keys2load_data = ['hdr', 'evt']
# df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_mup-data-OffBeamLight", keys2load=keys2load_data, n_max_concat=999)
# intime_evt_df = df_data['evt']
# intime_hdr_df = df_data['hdr']

# df_dir = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_06_042101__sel_mup-mc-BNB_cosmics"
# keys2load_data = ['hdr', 'evt']
# df_data = dfs_from_dir(search_dir=df_dir, filename_str="sel_mup-mc-BNB_cosmics", keys2load=keys2load_data, n_max_concat=999)
# mc_evt_df = df_data['evt']
# mc_hdr_df = df_data['hdr']


mc_evt_df[('mu', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['mu', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
mc_evt_df[('p', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['p', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['p', 'pfp', 'trk', 'dir', 'y', '', '']))

data_evt_df_['mu', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['mu', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
data_evt_df_['p', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['p', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['p', 'pfp', 'trk', 'dir', 'y', '', '']))

intime_evt_df[('mu', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['mu', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
intime_evt_df[('p', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['p', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['p', 'pfp', 'trk', 'dir', 'y', '', '']))

# dirt_evt_df[('mu', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(dirt_evt_df['mu', 'pfp', 'trk', 'dir', 'x', '', ''], dirt_evt_df['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
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
perTPC_inset = 5
def perTPC_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df["mu"].pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df["p"].pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df["mu"].pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df["p"].pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)
    perTPC_cut = in_TPC1_cut | in_TPC2_cut
    return perTPC_cut

mc_evt_df_perTPC = mc_evt_df[perTPC_cut(mc_evt_df)]
data_evt_df_perTPC = data_evt_df[perTPC_cut(data_evt_df)]
intime_evt_df_perTPC = intime_evt_df[perTPC_cut(intime_evt_df)]

def inTPC1_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)
    return in_TPC1_cut

mc_evt_df_inTPC1 = mc_evt_df[inTPC1_cut(mc_evt_df)]
data_evt_df_inTPC1 = data_evt_df[inTPC1_cut(data_evt_df)]
intime_evt_df_inTPC1 = intime_evt_df[inTPC1_cut(intime_evt_df)]


def inTPC2_cut(df):
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)
    return in_TPC2_cut

mc_evt_df_inTPC2 = mc_evt_df[inTPC2_cut(mc_evt_df)]
data_evt_df_inTPC2 = data_evt_df[inTPC2_cut(data_evt_df)]
intime_evt_df_inTPC2 = intime_evt_df[inTPC2_cut(intime_evt_df)]


# ============================================================
# Inspect Distributions in Detector Quadrants and Octants
#   - SBND volume octants: TPC E/W (x), N/S (z), top/bottom (y)
#     E/W: negative x → East (x < x0); positive x → West (x >= x0)
#     N/S: lower z → North (z < z0); higher z → South (z >= z0)
# ============================================================

def _weighted_counts(series, weights):
    # deterministic label ordering if series is categorical
    if isinstance(series.dtype, pd.CategoricalDtype):
        cats = series.dtype.categories
        out = pd.Series(0.0, index=cats, dtype=float)
        grp = pd.DataFrame({"k": series, "w": weights}).groupby("k")["w"].sum()
        out.loc[grp.index] = grp.values
        return out
    return pd.DataFrame({"k": series, "w": weights}).groupby("k")["w"].sum().sort_index()


def _add_sbnd_octant_labels(df, x0, y0, z0):
    x = df.slc.vertex.x.astype(float)
    y = df.slc.vertex.y.astype(float)
    z = df.slc.vertex.z.astype(float)

    ew = np.where(x < x0, "E", "W")  # SBND: negative x is East
    ns = np.where(z < z0, "N", "S")  # SBND: lower z is North
    tb = np.where(y >= y0, "Top", "Bottom")
    octant = np.char.add(np.char.add(np.char.add(ew, "-"), np.char.add(ns, "-")), tb)

    ew = pd.Categorical(ew, categories=["W", "E"], ordered=True)
    ns = pd.Categorical(ns, categories=["S", "N"], ordered=True)
    tb = pd.Categorical(tb, categories=["Bottom", "Top"], ordered=True)
    octant = pd.Categorical(
        octant,
        categories=[
            "W-S-Bottom", "W-S-Top", "W-N-Bottom", "W-N-Top",
            "E-S-Bottom", "E-S-Top", "E-N-Bottom", "E-N-Top",
        ],
        ordered=True,
    )

    out = df.copy()
    out[("sbnd", "octant", "ew", "", "", "", "")] = ew
    out[("sbnd", "octant", "ns", "", "", "", "")] = ns
    out[("sbnd", "octant", "tb", "", "", "", "")] = tb
    out[("sbnd", "octant", "octant", "", "", "", "")] = octant
    return out


def _add_sbnd_quadrant_labels(df, x0, y0):
    x = df.slc.vertex.x.astype(float)
    y = df.slc.vertex.y.astype(float)

    # SBND: negative x is East; lower y is Bottom
    ew = np.where(x < x0, "E", "W")
    tb = np.where(y >= y0, "Top", "Bottom")
    quad = np.char.add(np.char.add(ew, "-"), tb)

    ew = pd.Categorical(ew, categories=["E", "W"], ordered=True)
    tb = pd.Categorical(tb, categories=["Bottom", "Top"], ordered=True)
    quad = pd.Categorical(
        quad,
        categories=["E-Bottom", "E-Top", "W-Bottom", "W-Top"],
        ordered=True,
    )

    out = df.copy()
    out[("sbnd", "quad", "ew", "", "", "", "")] = ew
    out[("sbnd", "quad", "tb", "", "", "", "")] = tb
    out[("sbnd", "quad", "quad", "", "", "", "")] = quad
    return out


def _plot_octant_bars(df, title, save_base):
    w = df["pot_weight"].astype(float)
    ew = df["sbnd"].octant.ew
    ns = df["sbnd"].octant.ns
    tb = df["sbnd"].octant.tb
    oc = df["sbnd"].octant.octant

    counts_ew = _weighted_counts(ew, w)
    counts_ns = _weighted_counts(ns, w)
    counts_tb = _weighted_counts(tb, w)
    counts_oc = _weighted_counts(oc, w)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()

    for ax, counts, xlabel in [
        (axs[0], counts_ew, "TPC E/W (x split)"),
        (axs[1], counts_ns, "N/S (z split)"),
        (axs[2], counts_tb, "Top/Bottom (y split)"),
        (axs[3], counts_oc, "Octant (E/W, N/S, Top/Bottom)"),
    ]:
        xs = np.arange(len(counts.index))
        ax.bar(xs, counts.values, color="C0", alpha=0.85)
        ax.set_xticks(xs)
        ax.set_xticklabels([str(x) for x in counts.index], rotation=30, ha="right")
        ax.set_ylabel("Weighted Events")
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="y", alpha=0.25)

    fig.suptitle(title, fontsize=18)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    if save_fig:
        fig.savefig(save_base + fig_ext, bbox_inches="tight", dpi=dpi)
    plt.close(fig)

    return {"ew": counts_ew, "ns": counts_ns, "tb": counts_tb, "octant": counts_oc}


def _plot_quadrant_bars(df, title, save_base):
    w = df["pot_weight"].astype(float)
    ew = df["sbnd"].quad.ew
    tb = df["sbnd"].quad.tb
    qd = df["sbnd"].quad.quad

    counts_ew = _weighted_counts(ew, w)
    counts_tb = _weighted_counts(tb, w)
    counts_qd = _weighted_counts(qd, w)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, counts, xlabel in [
        (axs[0], counts_ew, "TPC E/W (x split)"),
        (axs[1], counts_tb, "Top/Bottom (y split)"),
        (axs[2], counts_qd, "Quadrant (E/W, Top/Bottom)"),
    ]:
        xs = np.arange(len(counts.index))
        ax.bar(xs, counts.values, color="C0", alpha=0.85)
        ax.set_xticks(xs)
        ax.set_xticklabels([str(x) for x in counts.index], rotation=30, ha="right")
        ax.set_ylabel("Weighted Events")
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="y", alpha=0.25)

    fig.suptitle(title, fontsize=16)
    fig.tight_layout(rect=[0, 0.02, 1, 0.92])
    if save_fig:
        fig.savefig(save_base + fig_ext, bbox_inches="tight", dpi=dpi)
    plt.close(fig)

    return {"ew": counts_ew, "tb": counts_tb, "quad": counts_qd}


# Use a consistent split point (x0,y0,z0) from the selected DATA sample.
# This avoids per-sample midpoints and makes E/W, N/S, Top/Bottom comparable.
_vx = data_evt_df.slc.vertex.x.astype(float).to_numpy()
_vy = data_evt_df.slc.vertex.y.astype(float).to_numpy()
_vz = data_evt_df.slc.vertex.z.astype(float).to_numpy()
_x0 = 0 #0.5 * (np.nanmin(_vx) + np.nanmax(_vx))
_y0 = 0 #0.5 * (np.nanmin(_vy) + np.nanmax(_vy))
_z0 = 250 #0.5 * (np.nanmin(_vz) + np.nanmax(_vz))

if args.do_octant_plots or args.do_quadrant_plots:
    print_sbnd_octant_vertex_ranges(_x0, _y0, _z0)

    for _name, _df in [
        ("data", data_evt_df),
        ("mc", mc_evt_df),
        ("intime", intime_evt_df),
        ("data_perTPC", data_evt_df_perTPC),
        ("mc_perTPC", mc_evt_df_perTPC),
        ("intime_perTPC", intime_evt_df_perTPC),
        ("data_inTPC1", data_evt_df_inTPC1),
        ("mc_inTPC1", mc_evt_df_inTPC1),
        ("intime_inTPC1", intime_evt_df_inTPC1),
        ("data_inTPC2", data_evt_df_inTPC2),
        ("mc_inTPC2", mc_evt_df_inTPC2),
        ("intime_inTPC2", intime_evt_df_inTPC2),
    ]:
        _labeled = _add_sbnd_quadrant_labels(_df, _x0, _y0) if args.do_quadrant_plots else _df
        _labeled = _add_sbnd_octant_labels(_labeled, _x0, _y0, _z0) if args.do_octant_plots else _labeled
        locals()[_name] = _labeled

    # overwrite the primary dfs used later (keeps downstream code untouched)
    data_evt_df = locals()["data"]
    mc_evt_df = locals()["mc"]
    intime_evt_df = locals()["intime"]
    data_evt_df_perTPC = locals()["data_perTPC"]
    mc_evt_df_perTPC = locals()["mc_perTPC"]
    intime_evt_df_perTPC = locals()["intime_perTPC"]
    data_evt_df_inTPC1 = locals()["data_inTPC1"]
    mc_evt_df_inTPC1 = locals()["mc_inTPC1"]
    intime_evt_df_inTPC1 = locals()["intime_inTPC1"]
    data_evt_df_inTPC2 = locals()["data_inTPC2"]
    mc_evt_df_inTPC2 = locals()["mc_inTPC2"]
    intime_evt_df_inTPC2 = locals()["intime_inTPC2"]

    if args.do_octant_plots:
        # quick diagnostic plots (saved next to the nominal outputs)
        _plot_octant_bars(
            data_evt_df,
            title=f"SBND octants (data, selected events) — split@({(_x0):.1f},{(_y0):.1f},{(_z0):.1f})",
            save_base=path.join(save_fig_dir, "sbnd_octants_data"),
        )
        _plot_octant_bars(
            mc_evt_df,
            title="SBND octants (MC, selected events)",
            save_base=path.join(save_fig_dir, "sbnd_octants_mc"),
        )
        _plot_octant_bars(
            intime_evt_df,
            title="SBND octants (off-beam/intime, selected events)",
            save_base=path.join(save_fig_dir, "sbnd_octants_intime"),
        )

    if args.do_quadrant_plots:
        _plot_quadrant_bars(
            data_evt_df,
            title=f"SBND quadrants (data, selected events) — split@({(_x0):.1f},{(_y0):.1f})",
            save_base=path.join(save_fig_dir, "sbnd_quadrants_data"),
        )
        _plot_quadrant_bars(
            mc_evt_df,
            title="SBND quadrants (MC, selected events)",
            save_base=path.join(save_fig_dir, "sbnd_quadrants_mc"),
        )
        _plot_quadrant_bars(
            intime_evt_df,
            title="SBND quadrants (off-beam/intime, selected events)",
            save_base=path.join(save_fig_dir, "sbnd_quadrants_intime"),
        )


def fwd_muons_cut(df):
    cut = df["mu"].pfp.trk.dir.z > 0
    return cut

mc_evt_df_fwd_muons = mc_evt_df[fwd_muons_cut(mc_evt_df)]
data_evt_df_fwd_muons = data_evt_df[fwd_muons_cut(data_evt_df)]
intime_evt_df_fwd_muons = intime_evt_df[fwd_muons_cut(intime_evt_df)]

def bwd_muons_cut(df):
    cut = df["mu"].pfp.trk.dir.z < 0
    return cut

mc_evt_df_bwd_muons = mc_evt_df[bwd_muons_cut(mc_evt_df)]
data_evt_df_bwd_muons = data_evt_df[bwd_muons_cut(data_evt_df)]
intime_evt_df_bwd_muons = intime_evt_df[bwd_muons_cut(intime_evt_df)]

def crosser_muons_cut(df):
    cut = df["mu"].pfp.trk.end.x * df["mu"].pfp.trk.start.x < 0
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

data_vs_mc_plotter_inTPC1= partial(
    overlay_hists,
    mc_df=mc_evt_df_inTPC1,
    data_df=data_evt_df_inTPC1,
    intime_df=intime_evt_df_inTPC1,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

data_vs_mc_plotter_inTPC2 = partial(
    overlay_hists,
    mc_df=mc_evt_df_inTPC2,
    data_df=data_evt_df_inTPC2,
    intime_df=intime_evt_df_inTPC2,
    dirt_df = None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig, 
    plot=False
)

# data_vs_mc_plotter_TPC1 = partial(
#     overlay_hists,
#     mc_df=mc_evt_df_TPC1,
#     data_df=data_evt_df_TPC1,
#     intime_df=intime_evt_df_TPC1,
#     dirt_df = None,
#     ax_ylim_ratio=ax_ylim_ratio,
#     ratio=ratio,
#     textloc=textloc,
#     approval=approval,
#     save_fig=save_fig, 
#     plot=False
# )

# data_vs_mc_plotter_TPC2 = partial(
#     overlay_hists,
#     mc_df=mc_evt_df_TPC2,
#     data_df=data_evt_df_TPC2,
#     intime_df=intime_evt_df_TPC2,
#     dirt_df = None,
#     ax_ylim_ratio=ax_ylim_ratio,
#     ratio=ratio,
#     textloc=textloc,
#     approval=approval,
#     save_fig=save_fig, 
#     plot=False
# )

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

chi2_records = []

def _record_chi2(ret, var_name, breakdown_type, cut_label):
    if ret.get("chi2_val") is not None:
        chi2_records.append({
            "var_name": var_name,
            "breakdown_type": breakdown_type,
            "cut_label": cut_label,
            "chi2": float(ret["chi2_val"]),
            "p_val": float(ret["p_val"]),
            "ndof": int(ret["ndof"]),
        })

def _plot_pull(ret, var_config, breakdown_type, save_dir):
    if ret.get("chi2_pull") is None:
        return
    plot_chi2_pull(
        ret["chi2_pull"],
        var_config,
        chi2_val=ret["chi2_val"],
        p_val=ret["p_val"],
        ndof=ret["ndof"],
        plot=False,
        save_fig=save_fig,
        save_name=path.join(save_dir, "{}_{}_chi2pull".format(var_config.var_save_name, breakdown_type)),
    )

# approved vars (shared tuple with cumulative / syst helpers)
var_configs = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)

_cut_variants_main = [
    ("nominal",  data_vs_mc_plotter,       save_fig_dir),
    ("perTPC",   data_vs_mc_plotter_perTPC, save_fig_dir_perTPC),
    ("inTPC1",   data_vs_mc_plotter_inTPC1, save_fig_dir_inTPC1),
    ("inTPC2",   data_vs_mc_plotter_inTPC2, save_fig_dir_inTPC2),
]

for var_config in var_configs:
   cov, frac_uncert = get_syst_unc(var_config)

   for breakdown_type in ["topology"]:
        plot_labels_hist = [var_config.var_labels[1], pot_label, ""]

        for cut_label, plotter, this_save_dir in _cut_variants_main:
            ret = plotter(breakdown_type=breakdown_type,
                          var_config=var_config,
                          plot_labels=plot_labels_hist,
                          syst=cov,
                          textchi2=True,
                          save_name=path.join(this_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)))
            _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label)
            _plot_pull(ret, var_config, breakdown_type, this_save_dir)


if args.do_octant_plots:
    # ==== Octant Plots ====
    var_configs_topology = list(var_configs)

    _SBND_OCTANT_LABELS = (
        "W-S-Bottom", "W-S-Top", "W-N-Bottom", "W-N-Top",
        "E-S-Bottom", "E-S-Top", "E-N-Bottom", "E-N-Top",
    )

    def _octant_rowmask(df, oct_label):
        return df["sbnd"].octant.octant.astype(str) == oct_label

    _octant_volume_variants = [
        ("nominal", save_fig_dir, data_evt_df, mc_evt_df, intime_evt_df),
        ("perTPC", save_fig_dir_perTPC, data_evt_df_perTPC, mc_evt_df_perTPC, intime_evt_df_perTPC),
        ("inTPC1", save_fig_dir_inTPC1, data_evt_df_inTPC1, mc_evt_df_inTPC1, intime_evt_df_inTPC1),
        ("inTPC2", save_fig_dir_inTPC2, data_evt_df_inTPC2, mc_evt_df_inTPC2, intime_evt_df_inTPC2),
    ]

    for vol_tag, save_root, d_full, m_full, i_full in _octant_volume_variants:
        oct_root = path.join(save_root, "by_octant")

        for oct_label in _SBND_OCTANT_LABELS:
            oct_slug = oct_label.replace("-", "_")
            oct_save_dir = path.join(oct_root, oct_slug)
            if save_fig and not path.exists(oct_save_dir):
                makedirs(oct_save_dir)

            d_o = d_full.loc[_octant_rowmask(d_full, oct_label)]
            m_o = m_full.loc[_octant_rowmask(m_full, oct_label)]
            i_o = i_full.loc[_octant_rowmask(i_full, oct_label)]
            if len(d_o) < 1 and len(m_o) < 1:
                continue

            plotter_oct = partial(
                overlay_hists,
                mc_df=m_o,
                data_df=d_o,
                intime_df=i_o,
                dirt_df=None,
                ax_ylim_ratio=ax_ylim_ratio,
                ratio=ratio,
                textloc=textloc,
                approval=approval,
                save_fig=save_fig,
                plot=False,
            )
            cut_label_oct = "octant_{}_{}".format(vol_tag, oct_slug)

            for var_config in var_configs_topology:
                cov, frac_uncert = get_syst_unc(var_config)
                for breakdown_type in ["topology"]:
                    plot_labels_hist = [
                        var_config.var_labels[1],
                        pot_label,
                        "{} · Octant {}".format(vol_tag, oct_label),
                    ]
                    ret = plotter_oct(
                        breakdown_type=breakdown_type,
                        var_config=var_config,
                        plot_labels=plot_labels_hist,
                        syst=cov,
                        textchi2=True,
                        save_name=path.join(oct_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)),
                    )
                    _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label_oct)
                    _plot_pull(ret, var_config, breakdown_type, oct_save_dir)

if args.do_quadrant_plots:
    # ==== Quadrant Plots (N/S combined) ====
    var_configs_topology_quad = list(var_configs_topology)

    _SBND_QUADRANT_LABELS = (
        "E-Bottom", "E-Top",
        "W-Bottom", "W-Top",
    )

    def _quadrant_rowmask(df, quad_label):
        return df["sbnd"].quad.quad.astype(str) == quad_label

    _quad_volume_variants = [
        ("nominal", save_fig_dir, data_evt_df, mc_evt_df, intime_evt_df),
        ("perTPC", save_fig_dir_perTPC, data_evt_df_perTPC, mc_evt_df_perTPC, intime_evt_df_perTPC),
        ("inTPC1", save_fig_dir_inTPC1, data_evt_df_inTPC1, mc_evt_df_inTPC1, intime_evt_df_inTPC1),
        ("inTPC2", save_fig_dir_inTPC2, data_evt_df_inTPC2, mc_evt_df_inTPC2, intime_evt_df_inTPC2),
    ]

    for vol_tag, save_root, d_full, m_full, i_full in _quad_volume_variants:
        quad_root = path.join(save_root, "by_quadrant")

        for quad_label in _SBND_QUADRANT_LABELS:
            quad_slug = quad_label.replace("-", "_")
            quad_save_dir = path.join(quad_root, quad_slug)
            if save_fig and not path.exists(quad_save_dir):
                makedirs(quad_save_dir)

            d_q = d_full.loc[_quadrant_rowmask(d_full, quad_label)]
            m_q = m_full.loc[_quadrant_rowmask(m_full, quad_label)]
            i_q = i_full.loc[_quadrant_rowmask(i_full, quad_label)]
            if len(d_q) < 1 and len(m_q) < 1:
                continue

            plotter_quad = partial(
                overlay_hists,
                mc_df=m_q,
                data_df=d_q,
                intime_df=i_q,
                dirt_df=None,
                ax_ylim_ratio=ax_ylim_ratio,
                ratio=ratio,
                textloc=textloc,
                approval=approval,
                save_fig=save_fig,
                plot=False,
            )
            cut_label_quad = "quadrant_{}_{}".format(vol_tag, quad_slug)

            for var_config in var_configs_topology_quad:
                cov, frac_uncert = get_syst_unc(var_config)
                for breakdown_type in ["topology"]:
                    plot_labels_hist = [
                        var_config.var_labels[1],
                        pot_label,
                        "{} · Quadrant {}".format(vol_tag, quad_label),
                    ]
                    ret = plotter_quad(
                        breakdown_type=breakdown_type,
                        var_config=var_config,
                        plot_labels=plot_labels_hist,
                        syst=cov,
                        textchi2=True,
                        save_name=path.join(quad_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)),
                    )
                    _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label_quad)
                    _plot_pull(ret, var_config, breakdown_type, quad_save_dir)


# More final-sample evt vars (same tuple as ``final_selected_evt_vars``; appended to
# PER_EVT_PLOTS in scripts/syst_detvar_chunk.py for detvar unisim histograms).
_FINAL_EVT_VARS_EXTRA = FINAL_SELECTED_EVT_VARIABLE_CONFIGS

var_configs = []
_seen_final_evt = set()
for _vc in _FINAL_EVT_VARS_EXTRA:
    if _vc.var_save_name not in _seen_final_evt:
        var_configs.append(_vc)
        _seen_final_evt.add(_vc.var_save_name)

for var_config in var_configs:

    for breakdown_type in ["topology"]:
    # for breakdown_type in ["topology", "genie", "genie_sb"]:
        plot_labels_hist = [var_config.var_labels[1], pot_label, ""]

        _cut_variants_phi = [
            ("nominal",      data_vs_mc_plotter,              save_fig_dir,          mc_evt_df,          intime_evt_df),
            ("perTPC",       data_vs_mc_plotter_perTPC,        save_fig_dir_perTPC,   mc_evt_df_perTPC,   intime_evt_df_perTPC),
            ("inTPC1",       data_vs_mc_plotter_inTPC1,        save_fig_dir_inTPC1,   mc_evt_df_inTPC1,   intime_evt_df_inTPC1),
            ("inTPC2",       data_vs_mc_plotter_inTPC2,        save_fig_dir_inTPC2,   mc_evt_df_inTPC2,   intime_evt_df_inTPC2),
            # ("crosser_fwd",  data_vs_mc_plotter_crosser_muons_fwd, save_fig_dir_crosser_fwd, mc_evt_df_crosser_muons_fwd, intime_evt_df_crosser_muons_fwd),
            # ("crosser_bwd",  data_vs_mc_plotter_crosser_muons_bwd, save_fig_dir_crosser_bwd, mc_evt_df_crosser_muons_bwd, intime_evt_df_crosser_muons_bwd),
        ]

        for cut_label, plotter, this_save_dir, this_mc_df, this_intime_df in _cut_variants_phi:
            frac_unc, cov = get_frac_unc(this_mc_df, this_intime_df, this_intime_df, var_config)
            _plot_labels = plot_labels_hist
            if "crosser" in cut_label:
                direction = "Forward" if "fwd" in cut_label else "Backward"
                _plot_labels = [var_config.var_labels[1], pot_label, f"Crosser Muons ({direction})"]
            ret = plotter(breakdown_type=breakdown_type,
                          var_config=var_config,
                          plot_labels=_plot_labels,
                          syst=cov,
                          textchi2=True,
                          save_name=path.join(this_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)))
            _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label)
            _plot_pull(ret, var_config, breakdown_type, this_save_dir)

        print("Saved plots for", var_config.var_save_name)


if args.do_octant_plots:
    # ==== Octant Plots ====
    var_configs_phi = list(var_configs)

    for vol_tag, save_root, d_full, m_full, i_full in _octant_volume_variants:
        oct_root = path.join(save_root, "by_octant")

        for oct_label in _SBND_OCTANT_LABELS:
            oct_slug = oct_label.replace("-", "_")
            oct_save_dir = path.join(oct_root, oct_slug)
            if save_fig and not path.exists(oct_save_dir):
                makedirs(oct_save_dir)

            d_o = d_full.loc[_octant_rowmask(d_full, oct_label)]
            m_o = m_full.loc[_octant_rowmask(m_full, oct_label)]
            i_o = i_full.loc[_octant_rowmask(i_full, oct_label)]
            if len(d_o) < 1 and len(m_o) < 1:
                continue

            plotter_oct = partial(
                overlay_hists,
                mc_df=m_o,
                data_df=d_o,
                intime_df=i_o,
                dirt_df=None,
                ax_ylim_ratio=ax_ylim_ratio,
                ratio=ratio,
                textloc=textloc,
                approval=approval,
                save_fig=save_fig,
                plot=False,
            )
            cut_label_oct_phi = "octant_phi_{}_{}".format(vol_tag, oct_slug)

            for var_config in var_configs_phi:
                for breakdown_type in ["topology"]:
                    frac_unc, cov = get_frac_unc(m_o, i_o, i_o, var_config)
                    plot_labels_hist = [
                        var_config.var_labels[1],
                        pot_label,
                        "{} · Octant {}".format(vol_tag, oct_label),
                    ]
                    ret = plotter_oct(
                        breakdown_type=breakdown_type,
                        var_config=var_config,
                        plot_labels=plot_labels_hist,
                        syst=cov,
                        textchi2=True,
                        save_name=path.join(oct_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)),
                    )
                    _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label_oct_phi)
                    _plot_pull(ret, var_config, breakdown_type, oct_save_dir)

if args.do_quadrant_plots:
    # ==== Quadrant Plots (N/S combined) ====
    var_configs_phi_quad = list(var_configs_phi)

    _SBND_QUADRANT_LABELS = (
        "E-Bottom", "E-Top",
        "W-Bottom", "W-Top",
    )

    def _quadrant_rowmask(df, quad_label):
        return df["sbnd"].quad.quad.astype(str) == quad_label

    _quad_volume_variants = [
        ("nominal", save_fig_dir, data_evt_df, mc_evt_df, intime_evt_df),
        ("perTPC", save_fig_dir_perTPC, data_evt_df_perTPC, mc_evt_df_perTPC, intime_evt_df_perTPC),
        ("inTPC1", save_fig_dir_inTPC1, data_evt_df_inTPC1, mc_evt_df_inTPC1, intime_evt_df_inTPC1),
        ("inTPC2", save_fig_dir_inTPC2, data_evt_df_inTPC2, mc_evt_df_inTPC2, intime_evt_df_inTPC2),
    ]

    for vol_tag, save_root, d_full, m_full, i_full in _quad_volume_variants:
        quad_root = path.join(save_root, "by_quadrant")

        for quad_label in _SBND_QUADRANT_LABELS:
            quad_slug = quad_label.replace("-", "_")
            quad_save_dir = path.join(quad_root, quad_slug)
            if save_fig and not path.exists(quad_save_dir):
                makedirs(quad_save_dir)

            d_q = d_full.loc[_quadrant_rowmask(d_full, quad_label)]
            m_q = m_full.loc[_quadrant_rowmask(m_full, quad_label)]
            i_q = i_full.loc[_quadrant_rowmask(i_full, quad_label)]
            if len(d_q) < 1 and len(m_q) < 1:
                continue

            plotter_quad = partial(
                overlay_hists,
                mc_df=m_q,
                data_df=d_q,
                intime_df=i_q,
                dirt_df=None,
                ax_ylim_ratio=ax_ylim_ratio,
                ratio=ratio,
                textloc=textloc,
                approval=approval,
                save_fig=save_fig,
                plot=False,
            )
            cut_label_quad_phi = "quadrant_phi_{}_{}".format(vol_tag, quad_slug)

            for var_config in var_configs_phi_quad:
                for breakdown_type in ["topology"]:
                    frac_unc, cov = get_frac_unc(m_q, i_q, i_q, var_config)
                    plot_labels_hist = [
                        var_config.var_labels[1],
                        pot_label,
                        "{} · Quadrant {}".format(vol_tag, quad_label),
                    ]
                    ret = plotter_quad(
                        breakdown_type=breakdown_type,
                        var_config=var_config,
                        plot_labels=plot_labels_hist,
                        syst=cov,
                        textchi2=True,
                        save_name=path.join(quad_save_dir, "{}_{}".format(var_config.var_save_name, breakdown_type)),
                    )
                    _record_chi2(ret, var_config.var_save_name, breakdown_type, cut_label_quad_phi)
                    _plot_pull(ret, var_config, breakdown_type, quad_save_dir)


# ==== save chi2 records ====
chi2_save_path = path.join(save_fig_dir, "chi2_records.json")
with open(chi2_save_path, "w") as f:
    json.dump(chi2_records, f, indent=2)
print(f"chi2 records saved to {chi2_save_path} ({len(chi2_records)} entries)")
