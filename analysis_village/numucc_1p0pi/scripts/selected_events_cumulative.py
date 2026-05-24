
# selected_events_cumulative.py
#
# Cumulative variant of selected_events.py (Stage 3 data access: Gen 1 cumulative batches).
# For each exposure-batch index K, this script uses ALL time-ordered batches 0..K (inclusive)
# instead of only batch K. Iterating K = 0, 1, ..., n_time_splits-1 produces monotonically
# increasing integrated data POT. Terminology: see analysis_village.numucc_1p0pi.exposure_access.
#
# This script processes ALL requested batches in a single Python invocation:
#   - get_ana_dfs(...) is called ONCE (the dominant I/O cost).
#   - get_syst_unc(var_config) is memoized per var_config.
#   - The per-chunk work (cumulative POT, weights, cuts, plots) is the only
#     thing that runs inside the loop.

import os
from os import path, makedirs
from datetime import datetime
from functools import partial
import pickle
import argparse
import json
import time

import numpy as np
import pandas as pd

import sys
# sys.path.append('../../../')
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana')  # repo root when cwd is not cafpyana
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.dataset_locations import SELECTED_EVENTS_GLOBS, sorted_glob
plt.style.use("presentation.mplstyle")

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser(description='Event selection settings (cumulative, single-process loop)')
parser.add_argument('--n_time_splits', type=int, default=15,
                    help='Number of time-ordered exposure batches for data')
parser.add_argument('--chunk_idxs', type=int, nargs='+', default=None,
                    help='Legacy: cumulative batch indices K (each run uses batches 0..K). '
                         'Default: all 0..n_time_splits-1.')
parser.add_argument(
    '--exposure-batch-indices',
    type=int,
    nargs='+',
    default=None,
    dest='chunk_idxs',
    help='Preferred alias for --chunk_idxs.',
)
parser.add_argument('--chunk_idx', type=int, default=None,
                    help='Legacy single index; equivalent to --chunk_idxs <chunk_idx>.')
parser.add_argument(
    '--exposure-batch-index',
    type=int,
    default=None,
    dest='chunk_idx',
    help='Preferred alias for --chunk_idx.',
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
args = parser.parse_args()

# resolve chunk index list
if args.chunk_idxs is not None:
    chunk_idxs_to_run = list(args.chunk_idxs)
elif args.chunk_idx is not None:
    chunk_idxs_to_run = [args.chunk_idx]
else:
    chunk_idxs_to_run = list(range(args.n_time_splits))

for _idx in chunk_idxs_to_run:
    if _idx < 0 or _idx >= args.n_time_splits:
        raise ValueError(
            f"chunk_idx={_idx} is out of range for n_time_splits={args.n_time_splits}"
        )

print("Processing n_time_splits: ", args.n_time_splits)
print("Processing cumulative chunk_idxs: ", chunk_idxs_to_run)
print("Processing do_octant_plots: ", args.do_octant_plots)
print("Processing do_quadrant_plots: ", args.do_quadrant_plots)


# ============================================================
# Systematic uncertainty loader (memoized per var_config name)
# ============================================================

_syst_cache = {}

def _load_syst_unc_uncached(var_config):
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

    pot_frac_unc = 0.02
    ntargets_frac_unc = 0.01

    frac_uncert_total = np.zeros(len(var_config.bin_centers))
    cov_total = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))
    systs = [mcstat_syst, genie_syst, flux_syst, g4_syst, cosmics_syst]
    for syst in systs:
        cov_total += syst
        syst_uncert = np.sqrt(np.diag(syst))
        frac_uncert_total += syst_uncert ** 2

    flat_systs = [pot_frac_unc, ntargets_frac_unc]
    for syst in flat_systs:
        cov_total += np.diag(syst * np.ones(len(var_config.bin_centers)) ** 2)
        syst_uncert = syst * np.ones(len(var_config.bin_centers))
        frac_uncert_total += syst_uncert ** 2

    frac_uncert_total = np.sqrt(frac_uncert_total)
    return cov_total, frac_uncert_total


def get_syst_unc(var_config):
    key = var_config.var_save_name
    if key not in _syst_cache:
        _syst_cache[key] = _load_syst_unc_uncached(var_config)
    return _syst_cache[key]


save_fig = True
today_str = datetime.now().strftime("%Y%m%d")

# Zero-pad chunk index so that lexicographic sorting on directory names
# (cum00, cum01, ..., cum14) matches numeric order. This makes GIF assembly
# from `cum*/<var>_<bd>.<ext>` work with a plain glob+sort.
_cum_pad_width = max(2, len(str(max(0, args.n_time_splits - 1))))
def _cum_tag(chunk_idx):
    return f"cum{chunk_idx:0{_cum_pad_width}d}"


# ============================================================
# Octant / Quadrant utility functions (defined once)
#   - SBND volume octants: TPC E/W (x), N/S (z), top/bottom (y)
#     E/W: negative x → East (x < x0); positive x → West (x >= x0)
#     N/S: lower z → South (z < z0); higher z → North (z >= z0)
# ============================================================

def _weighted_counts(series, weights):
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

    ew = np.where(x < x0, "E", "W")
    ns = np.where(z < z0, "S", "N")
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


_SBND_OCTANT_LABELS = (
    "W-S-Bottom", "W-S-Top", "W-N-Bottom", "W-N-Top",
    "E-S-Bottom", "E-S-Top", "E-N-Bottom", "E-N-Top",
)
_SBND_QUADRANT_LABELS = (
    "E-Bottom", "E-Top",
    "W-Bottom", "W-Top",
)

def _octant_rowmask(df, oct_label):
    return df["sbnd"].octant.octant.astype(str) == oct_label

def _quadrant_rowmask(df, quad_label):
    return df["sbnd"].quad.quad.astype(str) == quad_label


# ============================================================
# Cuts (functions over dfs; pure, no I/O)
# ============================================================

perTPC_inset = 5

def perTPC_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df["mu"].pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df["p"].pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df["mu"].pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df["p"].pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)
    return in_TPC1_cut | in_TPC2_cut

def inTPC1_cut(df):
    return InFV(df.slc.vertex, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset) & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)

def inTPC2_cut(df):
    return InFV(df.slc.vertex, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset) & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)

def fwd_muons_cut(df):
    return df["mu"].pfp.trk.dir.z > 0

def bwd_muons_cut(df):
    return df["mu"].pfp.trk.dir.z < 0

def crosser_muons_cut(df):
    return df["mu"].pfp.trk.end.x * df["mu"].pfp.trk.start.x < 0


# ============================================================
# Variable configurations (built once)
# ============================================================

# Primary kinematics (same as ``selected_events.py`` / ``CORE_SELECTED_EVT_VARIABLE_CONFIGS``)
var_configs_main = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)

# Φ / vertex / endpoint extras require derived columns below; includes core + merged finals
var_configs_phi = with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))


# ============================================================
# Heavy I/O: load all dfs ONCE
# ============================================================

t_load_start = time.time()

def _concat_hdf_splits_from_files(files, keys2load, n_max_splits=999):
    """Load + concat split-key HDF shards while keeping __ntuple globally unique.

    This mirrors the non-legacy directory-based workflow used by notebooks/scripts
    via pyanalib.split_df_helpers_new.dfs_from_dir, but takes an explicit file list.
    """
    if not files:
        return {k: pd.DataFrame() for k in keys2load}

    df_lists = {k: [] for k in keys2load}
    ntuple_offset = np.int64(0)

    for f in files:
        try:
            this = load_dfs(f, keys2load, n_max_concat=n_max_splits)
        except Exception as e:
            print(f"[I/O] WARNING: failed to load {f}: {e}")
            continue

        # Dense remapping prevents ntuple_offset from exploding with raw index magnitude.
        ref_df = this[keys2load[0]]
        if isinstance(ref_df.index, pd.MultiIndex):
            raw_ntuple_vals = ref_df.index.get_level_values(0)
        else:
            raw_ntuple_vals = ref_df.index
        unique_ntuples = np.array(sorted(raw_ntuple_vals.unique()))
        ntuple_remap = {old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)}
        n_unique = np.int64(len(unique_ntuples))

        for k in keys2load:
            df = this[k]
            if isinstance(df.index, pd.MultiIndex):
                names = df.index.names
                idx_loc = names.index("__ntuple") if "__ntuple" in names else 0
                new_tuples = []
                for tup in df.index:
                    tup = list(tup)
                    tup[idx_loc] = ntuple_remap[tup[idx_loc]]
                    new_tuples.append(tuple(tup))
                df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
            else:
                if df.index.name == "__ntuple":
                    df.index = df.index.map(ntuple_remap)
            df_lists[k].append(df)

        ntuple_offset += n_unique

    return {k: pd.concat(df_lists[k], axis=0, sort=False) if df_lists[k] else pd.DataFrame()
            for k in keys2load}


def _load_selected_events_sample(sample, keys2load=("hdr", "evt")):
    if sample not in SELECTED_EVENTS_GLOBS:
        raise KeyError(f"unknown sample {sample!r}; expected one of {tuple(SELECTED_EVENTS_GLOBS)}")
    pattern = SELECTED_EVENTS_GLOBS[sample]
    files = sorted_glob(pattern)
    if not files:
        print(f"[I/O] WARNING: sample={sample} matched 0 files for glob: {pattern}")
    else:
        print(f"[I/O] sample={sample} matched {len(files)} files")
    return _concat_hdf_splits_from_files(files, list(keys2load), n_max_splits=999)


mc_dfs = _load_selected_events_sample("mc", keys2load=("hdr", "evt"))
data_dfs = _load_selected_events_sample("data", keys2load=("hdr", "evt"))
intime_dfs = _load_selected_events_sample("intime", keys2load=("hdr", "evt"))
dirt_dfs = _load_selected_events_sample("dirt", keys2load=("hdr", "evt"))

mc_hdr_df, mc_evt_df = mc_dfs["hdr"], mc_dfs["evt"]
data_hdr_df_, data_evt_df_ = data_dfs["hdr"], data_dfs["evt"]
intime_hdr_df, intime_evt_df = intime_dfs["hdr"], intime_dfs["evt"]
dirt_hdr_df, dirt_evt_df = dirt_dfs["hdr"], dirt_dfs["evt"]

print(f"[I/O] selected-events dfs loaded in {time.time() - t_load_start:.1f}s")

# Phi columns (added once on full dfs)
mc_evt_df[('mu', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['mu', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
mc_evt_df[('p', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(mc_evt_df['p', 'pfp', 'trk', 'dir', 'x', '', ''], mc_evt_df['p', 'pfp', 'trk', 'dir', 'y', '', '']))

data_evt_df_['mu', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['mu', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
data_evt_df_['p', 'pfp', 'trk', 'phi', '', '', ''] = np.degrees(np.arctan2(data_evt_df_['p', 'pfp', 'trk', 'dir', 'x', '', ''], data_evt_df_['p', 'pfp', 'trk', 'dir', 'y', '', '']))

intime_evt_df[('mu', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['mu', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['mu', 'pfp', 'trk', 'dir', 'y', '', '']))
intime_evt_df[('p', 'pfp', 'trk', 'phi', '', '', '')] = np.degrees(np.arctan2(intime_evt_df['p', 'pfp', 'trk', 'dir', 'x', '', ''], intime_evt_df['p', 'pfp', 'trk', 'dir', 'y', '', '']))

# TODO: for the integrated plot to work
mc_evt_df.loc[mc_evt_df.mc.iscc.isna(), ("mc","iscc")] = 999
data_evt_df_["mc", "iscc"] = 999

# Time-ordered chunk splits (computed once)
_n_time_splits = args.n_time_splits
_sorted = data_hdr_df_.sort_values(["run", "evt"], kind="mergesort")
data_hdr_df_splits = [
    _sorted.iloc[idx]
    for idx in np.array_split(np.arange(len(_sorted)), _n_time_splits)
]


# ============================================================
# Per-chunk processing (everything below depends on chunk_idx)
# ============================================================

def process_chunk(chunk_idx):
    t0 = time.time()
    print(f"\n=========== cumulative chunk_idx={chunk_idx} (chunks 0..{chunk_idx}) ===========")
    # Belt-and-suspenders: drop any matplotlib figures lingering from a previous
    # chunk so a long sweep can't run out of memory or silently fail to save.
    plt.close('all')
    cum_tag = _cum_tag(chunk_idx)

    # output dirs
    save_fig_dir = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}/{cum_tag}")
    save_fig_dir_perTPC = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-perTPC/{cum_tag}")
    save_fig_dir_inTPC1 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-contTPC1/{cum_tag}")
    save_fig_dir_inTPC2 = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-contTPC2/{cum_tag}")
    save_fig_dir_fwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-fwd_muons/{cum_tag}")
    save_fig_dir_bwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-bwd_muons/{cum_tag}")
    save_fig_dir_crosser_fwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-crosser_muons_fwd/{cum_tag}")
    save_fig_dir_crosser_bwd = path.join(save_fig_base_dir, f"selected_events-data-mup-1e20-cumulative-{today_str}-crosser_muons_bwd/{cum_tag}")

    if save_fig:
        for d in [save_fig_dir, save_fig_dir_fwd, save_fig_dir_bwd,
                  save_fig_dir_crosser_fwd, save_fig_dir_crosser_bwd,
                  save_fig_dir_perTPC, save_fig_dir_inTPC1, save_fig_dir_inTPC2]:
            if not path.exists(d):
                makedirs(d)
        print("saving plots under: ", save_fig_dir)

    # cumulative data selection: chunks 0..chunk_idx
    data_hdr_df = pd.concat(data_hdr_df_splits[: chunk_idx + 1])
    print(f"  using {chunk_idx + 1}/{_n_time_splits} chunks "
          f"({len(data_hdr_df)} hdr rows of {len(_sorted)} total)")

    common_idxs = data_evt_df_.reset_index(level=[2]).index.intersection(data_hdr_df.index)
    data_evt_df = (data_evt_df_.reset_index(level=[2])
                                .loc[common_idxs]
                                .reset_index()
                                .set_index(["__ntuple", "entry", "rec.slc..index"]))

    # POT / scaling
    data_tot_pot = data_hdr_df['pot'].sum()
    data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
    pot_str = get_pot_str(data_tot_pot)
    pot_label = f"Events / Bin (POT={pot_str})"
    data_gates = data_hdr_df.nbnbinfo.sum()
    print("  data_tot_pot: %.3e" % data_tot_pot)
    print("  data tot gates : %.3e" % data_gates)

    mc_tot_pot = mc_hdr_df['pot'].sum()
    mc_pot_scale = data_tot_pot / mc_tot_pot
    print("  mc_pot_scale: %.3e" % mc_pot_scale)
    mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))

    intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
    f = 0.0753
    scale_intime_to_lightdata = (1 - f) * data_gates / intime_gates
    print("  intime data scale: {:.2f}".format(scale_intime_to_lightdata))
    intime_evt_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
    intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))

    # cuts (re-derive subsets per chunk so they pick up the latest pot_weight)
    mc_evt_df_perTPC = mc_evt_df[perTPC_cut(mc_evt_df)]
    data_evt_df_perTPC = data_evt_df[perTPC_cut(data_evt_df)]
    intime_evt_df_perTPC = intime_evt_df[perTPC_cut(intime_evt_df)]

    mc_evt_df_inTPC1 = mc_evt_df[inTPC1_cut(mc_evt_df)]
    data_evt_df_inTPC1 = data_evt_df[inTPC1_cut(data_evt_df)]
    intime_evt_df_inTPC1 = intime_evt_df[inTPC1_cut(intime_evt_df)]

    mc_evt_df_inTPC2 = mc_evt_df[inTPC2_cut(mc_evt_df)]
    data_evt_df_inTPC2 = data_evt_df[inTPC2_cut(data_evt_df)]
    intime_evt_df_inTPC2 = intime_evt_df[inTPC2_cut(intime_evt_df)]

    # Optional octant/quadrant decoration. Use a dict to track (re-)labeled frames.
    if args.do_octant_plots or args.do_quadrant_plots:
        _x0, _y0, _z0 = 0, 0, 250
        print_sbnd_octant_vertex_ranges(_x0, _y0, _z0)

        chunk_dfs = {
            "data": data_evt_df, "mc": mc_evt_df, "intime": intime_evt_df,
            "data_perTPC": data_evt_df_perTPC, "mc_perTPC": mc_evt_df_perTPC, "intime_perTPC": intime_evt_df_perTPC,
            "data_inTPC1": data_evt_df_inTPC1, "mc_inTPC1": mc_evt_df_inTPC1, "intime_inTPC1": intime_evt_df_inTPC1,
            "data_inTPC2": data_evt_df_inTPC2, "mc_inTPC2": mc_evt_df_inTPC2, "intime_inTPC2": intime_evt_df_inTPC2,
        }
        for _name, _df in list(chunk_dfs.items()):
            _labeled = _add_sbnd_quadrant_labels(_df, _x0, _y0) if args.do_quadrant_plots else _df
            _labeled = _add_sbnd_octant_labels(_labeled, _x0, _y0, _z0) if args.do_octant_plots else _labeled
            chunk_dfs[_name] = _labeled

        data_evt_df = chunk_dfs["data"]
        mc_evt_df_loc = chunk_dfs["mc"]
        intime_evt_df_loc = chunk_dfs["intime"]
        data_evt_df_perTPC = chunk_dfs["data_perTPC"]
        mc_evt_df_perTPC = chunk_dfs["mc_perTPC"]
        intime_evt_df_perTPC = chunk_dfs["intime_perTPC"]
        data_evt_df_inTPC1 = chunk_dfs["data_inTPC1"]
        mc_evt_df_inTPC1 = chunk_dfs["mc_inTPC1"]
        intime_evt_df_inTPC1 = chunk_dfs["intime_inTPC1"]
        data_evt_df_inTPC2 = chunk_dfs["data_inTPC2"]
        mc_evt_df_inTPC2 = chunk_dfs["mc_inTPC2"]
        intime_evt_df_inTPC2 = chunk_dfs["intime_inTPC2"]

        if args.do_octant_plots:
            _plot_octant_bars(
                data_evt_df,
                title=f"SBND octants (data, cum0..{chunk_idx}) — split@({_x0:.1f},{_y0:.1f},{_z0:.1f})",
                save_base=path.join(save_fig_dir, "sbnd_octants_data"),
            )
            _plot_octant_bars(
                mc_evt_df_loc,
                title="SBND octants (MC, selected events)",
                save_base=path.join(save_fig_dir, "sbnd_octants_mc"),
            )
            _plot_octant_bars(
                intime_evt_df_loc,
                title="SBND octants (off-beam/intime, selected events)",
                save_base=path.join(save_fig_dir, "sbnd_octants_intime"),
            )

        if args.do_quadrant_plots:
            _plot_quadrant_bars(
                data_evt_df,
                title=f"SBND quadrants (data, cum0..{chunk_idx}) — split@({_x0:.1f},{_y0:.1f})",
                save_base=path.join(save_fig_dir, "sbnd_quadrants_data"),
            )
            _plot_quadrant_bars(
                mc_evt_df_loc,
                title="SBND quadrants (MC, selected events)",
                save_base=path.join(save_fig_dir, "sbnd_quadrants_mc"),
            )
            _plot_quadrant_bars(
                intime_evt_df_loc,
                title="SBND quadrants (off-beam/intime, selected events)",
                save_base=path.join(save_fig_dir, "sbnd_quadrants_intime"),
            )
    else:
        mc_evt_df_loc = mc_evt_df
        intime_evt_df_loc = intime_evt_df

    # forward / backward / crosser muons (built fresh for this chunk's cumulative data)
    mc_evt_df_fwd_muons = mc_evt_df_loc[fwd_muons_cut(mc_evt_df_loc)]
    data_evt_df_fwd_muons = data_evt_df[fwd_muons_cut(data_evt_df)]
    intime_evt_df_fwd_muons = intime_evt_df_loc[fwd_muons_cut(intime_evt_df_loc)]

    mc_evt_df_bwd_muons = mc_evt_df_loc[bwd_muons_cut(mc_evt_df_loc)]
    data_evt_df_bwd_muons = data_evt_df[bwd_muons_cut(data_evt_df)]
    intime_evt_df_bwd_muons = intime_evt_df_loc[bwd_muons_cut(intime_evt_df_loc)]

    mc_evt_df_crosser_muons_fwd = mc_evt_df_loc[crosser_muons_cut(mc_evt_df_loc) & fwd_muons_cut(mc_evt_df_loc)]
    data_evt_df_crosser_muons_fwd = data_evt_df[crosser_muons_cut(data_evt_df) & fwd_muons_cut(data_evt_df)]
    intime_evt_df_crosser_muons_fwd = intime_evt_df_loc[crosser_muons_cut(intime_evt_df_loc) & fwd_muons_cut(intime_evt_df_loc)]

    mc_evt_df_crosser_muons_bwd = mc_evt_df_loc[crosser_muons_cut(mc_evt_df_loc) & bwd_muons_cut(mc_evt_df_loc)]
    data_evt_df_crosser_muons_bwd = data_evt_df[crosser_muons_cut(data_evt_df) & bwd_muons_cut(data_evt_df)]
    intime_evt_df_crosser_muons_bwd = intime_evt_df_loc[crosser_muons_cut(intime_evt_df_loc) & bwd_muons_cut(intime_evt_df_loc)]

    # ===- plotter per cut ====
    eps = 1e-8
    ratio = True
    approval = "internal"
    textloc = [0.03, 0.55]
    ax_ylim_ratio = 1.9

    common_plotter_kwargs = dict(
        dirt_df=None,
        ax_ylim_ratio=ax_ylim_ratio,
        ratio=ratio,
        textloc=textloc,
        approval=approval,
        save_fig=save_fig,
        plot=False,
    )

    data_vs_mc_plotter = partial(
        overlay_hists, mc_df=mc_evt_df_loc, data_df=data_evt_df, intime_df=intime_evt_df_loc, **common_plotter_kwargs)
    data_vs_mc_plotter_perTPC = partial(
        overlay_hists, mc_df=mc_evt_df_perTPC, data_df=data_evt_df_perTPC, intime_df=intime_evt_df_perTPC, **common_plotter_kwargs)
    data_vs_mc_plotter_inTPC1 = partial(
        overlay_hists, mc_df=mc_evt_df_inTPC1, data_df=data_evt_df_inTPC1, intime_df=intime_evt_df_inTPC1, **common_plotter_kwargs)
    data_vs_mc_plotter_inTPC2 = partial(
        overlay_hists, mc_df=mc_evt_df_inTPC2, data_df=data_evt_df_inTPC2, intime_df=intime_evt_df_inTPC2, **common_plotter_kwargs)
    data_vs_mc_plotter_fwd_muons = partial(
        overlay_hists, mc_df=mc_evt_df_fwd_muons, data_df=data_evt_df_fwd_muons, intime_df=intime_evt_df_fwd_muons, **common_plotter_kwargs)
    data_vs_mc_plotter_bwd_muons = partial(
        overlay_hists, mc_df=mc_evt_df_bwd_muons, data_df=data_evt_df_bwd_muons, intime_df=intime_evt_df_bwd_muons, **common_plotter_kwargs)
    data_vs_mc_plotter_crosser_muons_fwd = partial(
        overlay_hists, mc_df=mc_evt_df_crosser_muons_fwd, data_df=data_evt_df_crosser_muons_fwd, intime_df=intime_evt_df_crosser_muons_fwd, **common_plotter_kwargs)
    data_vs_mc_plotter_crosser_muons_bwd = partial(
        overlay_hists, mc_df=mc_evt_df_crosser_muons_bwd, data_df=data_evt_df_crosser_muons_bwd, intime_df=intime_evt_df_crosser_muons_bwd, **common_plotter_kwargs)

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
                "cum_chunk_idx": chunk_idx,
                "n_chunks_included": chunk_idx + 1,
                "n_time_splits": _n_time_splits,
                "data_tot_pot": float(data_tot_pot),
                "data_gates": float(data_gates),
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

    _cut_variants_main = [
        ("nominal",  data_vs_mc_plotter,        save_fig_dir),
        ("perTPC",   data_vs_mc_plotter_perTPC, save_fig_dir_perTPC),
        ("inTPC1",   data_vs_mc_plotter_inTPC1, save_fig_dir_inTPC1),
        ("inTPC2",   data_vs_mc_plotter_inTPC2, save_fig_dir_inTPC2),
    ]

    for var_config in var_configs_main:
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

    # ==== Octant Plots (main vars) ====
    if args.do_octant_plots:
        _octant_volume_variants = [
            ("nominal", save_fig_dir, data_evt_df, mc_evt_df_loc, intime_evt_df_loc),
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
                    overlay_hists, mc_df=m_o, data_df=d_o, intime_df=i_o, **common_plotter_kwargs)
                cut_label_oct = "octant_{}_{}".format(vol_tag, oct_slug)
                for var_config in var_configs_main:
                    cov, frac_uncert = get_syst_unc(var_config)
                    for breakdown_type in ["topology"]:
                        plot_labels_hist = [
                            var_config.var_labels[1],
                            pot_label,
                            "{} · Octant {} · cum0..{}".format(vol_tag, oct_label, chunk_idx),
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

    # ==== Quadrant Plots (main vars; N/S combined) ====
    if args.do_quadrant_plots:
        _quad_volume_variants = [
            ("nominal", save_fig_dir, data_evt_df, mc_evt_df_loc, intime_evt_df_loc),
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
                    overlay_hists, mc_df=m_q, data_df=d_q, intime_df=i_q, **common_plotter_kwargs)
                cut_label_quad = "quadrant_{}_{}".format(vol_tag, quad_slug)
                for var_config in var_configs_main:
                    cov, frac_uncert = get_syst_unc(var_config)
                    for breakdown_type in ["topology"]:
                        plot_labels_hist = [
                            var_config.var_labels[1],
                            pot_label,
                            "{} · Quadrant {} · cum0..{}".format(vol_tag, quad_label, chunk_idx),
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

    # ==== more vars (vertex_x/y/z): use get_frac_unc per cut variant ====

    _cut_variants_phi = [
        ("nominal",      data_vs_mc_plotter,         save_fig_dir,        mc_evt_df_loc,      intime_evt_df_loc),
        ("perTPC",       data_vs_mc_plotter_perTPC,  save_fig_dir_perTPC, mc_evt_df_perTPC,   intime_evt_df_perTPC),
        ("inTPC1",       data_vs_mc_plotter_inTPC1,  save_fig_dir_inTPC1, mc_evt_df_inTPC1,   intime_evt_df_inTPC1),
        ("inTPC2",       data_vs_mc_plotter_inTPC2,  save_fig_dir_inTPC2, mc_evt_df_inTPC2,   intime_evt_df_inTPC2),
        # ("crosser_fwd",  data_vs_mc_plotter_crosser_muons_fwd, save_fig_dir_crosser_fwd, mc_evt_df_crosser_muons_fwd, intime_evt_df_crosser_muons_fwd),
        # ("crosser_bwd",  data_vs_mc_plotter_crosser_muons_bwd, save_fig_dir_crosser_bwd, mc_evt_df_crosser_muons_bwd, intime_evt_df_crosser_muons_bwd),
    ]

    for var_config in var_configs_phi:
        for breakdown_type in ["topology"]:
            plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
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
            print("  saved plots for", var_config.var_save_name)

    # ==== Octant plots (phi vars) ====
    if args.do_octant_plots:
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
                    overlay_hists, mc_df=m_o, data_df=d_o, intime_df=i_o, **common_plotter_kwargs)
                cut_label_oct_phi = "octant_phi_{}_{}".format(vol_tag, oct_slug)
                for var_config in var_configs_phi:
                    for breakdown_type in ["topology"]:
                        frac_unc, cov = get_frac_unc(m_o, i_o, i_o, var_config)
                        plot_labels_hist = [
                            var_config.var_labels[1],
                            pot_label,
                            "{} · Octant {} · cum0..{}".format(vol_tag, oct_label, chunk_idx),
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

    # ==== Quadrant plots (phi vars; N/S combined) ====
    if args.do_quadrant_plots:
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
                    overlay_hists, mc_df=m_q, data_df=d_q, intime_df=i_q, **common_plotter_kwargs)
                cut_label_quad_phi = "quadrant_phi_{}_{}".format(vol_tag, quad_slug)
                for var_config in var_configs_phi:
                    for breakdown_type in ["topology"]:
                        frac_unc, cov = get_frac_unc(m_q, i_q, i_q, var_config)
                        plot_labels_hist = [
                            var_config.var_labels[1],
                            pot_label,
                            "{} · Quadrant {} · cum0..{}".format(vol_tag, quad_label, chunk_idx),
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

    # ==== save chi2 records (per cumulative chunk) ====
    chi2_save_path = path.join(save_fig_dir, "chi2_records.json")
    with open(chi2_save_path, "w") as f:
        json.dump(chi2_records, f, indent=2)
    print(f"  chi2 records saved to {chi2_save_path} ({len(chi2_records)} entries)")

    # Make sure no figures are still open between chunks.
    plt.close('all')

    # ==== verify on-disk frames per cut variant (helps catch silent failures) ====
    saved_frames = {}
    for tag, d in [
        ("nominal", save_fig_dir),
        ("perTPC", save_fig_dir_perTPC),
        ("inTPC1", save_fig_dir_inTPC1),
        ("inTPC2", save_fig_dir_inTPC2),
        ("fwd", save_fig_dir_fwd),
        ("bwd", save_fig_dir_bwd),
        ("crosser_fwd", save_fig_dir_crosser_fwd),
        ("crosser_bwd", save_fig_dir_crosser_bwd),
    ]:
        if path.isdir(d):
            try:
                files = [
                    f for f in os.listdir(d)
                    if f.lower().endswith((".png", ".pdf", ".jpg", ".jpeg", ".svg"))
                ]
            except OSError:
                files = []
            saved_frames[tag] = len(files)
        else:
            saved_frames[tag] = 0
    summary = ", ".join(f"{k}={v}" for k, v in saved_frames.items())
    print(f"  saved frames per cut [chunk={chunk_idx}]: {summary}")

    elapsed = time.time() - t0
    print(f"  cum chunk {chunk_idx} done in {elapsed:.1f}s")
    return chi2_records, {"chunk_idx": chunk_idx, "cum_tag": cum_tag,
                          "save_fig_dir": save_fig_dir, "saved_frames": saved_frames,
                          "data_tot_pot": float(data_tot_pot),
                          "data_gates": float(data_gates)}


# ============================================================
# Drive the loop
# ============================================================

all_chi2_records = []
all_manifest_entries = []
t_loop_start = time.time()
for _chunk_idx in chunk_idxs_to_run:
    _chi2_recs, _manifest = process_chunk(_chunk_idx)
    all_chi2_records.extend(_chi2_recs)
    all_manifest_entries.append(_manifest)
print(f"\n[loop] processed {len(chunk_idxs_to_run)} cumulative chunks in {time.time() - t_loop_start:.1f}s")

# Combined outputs across all processed chunks
if chunk_idxs_to_run:
    combined_dir = path.join(
        save_fig_base_dir,
        f"selected_events-data-mup-1e20-cumulative-{today_str}",
    )
    if save_fig and not path.exists(combined_dir):
        makedirs(combined_dir)

    # Combined chi2 records
    combined_path = path.join(combined_dir, "chi2_records_all.json")
    with open(combined_path, "w") as f:
        json.dump(all_chi2_records, f, indent=2)
    print(f"combined chi2 records ({len(all_chi2_records)} entries) saved to {combined_path}")

    # Frame manifest: for each cum chunk, where its plots live and how many were
    # actually written. This is exactly what `make_gifs.ipynb` (or any other
    # GIF-builder) needs to assemble per-variable animations across cum0..N.
    manifest_path = path.join(combined_dir, "frame_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump({
            "today": today_str,
            "n_time_splits": args.n_time_splits,
            "chunk_idxs_run": chunk_idxs_to_run,
            "cum_pad_width": _cum_pad_width,
            "chunks": all_manifest_entries,
        }, f, indent=2)
    print(f"frame manifest saved to {manifest_path}")

    # Print a one-line "GIF hint" for each cut/variable pair under the nominal
    # output: this is the glob you'd hand to ImageMagick / ffmpeg.
    if all_manifest_entries:
        print("\nGIF assembly hints (frames listed in cum order):")
        for var_save_name in sorted({r["var_name"] for r in all_chi2_records}):
            for cut_label in sorted({r["cut_label"] for r in all_chi2_records
                                     if r["var_name"] == var_save_name and "_" not in r["cut_label"][:7]}):
                # only show top-level cut variants in the hint
                hint = path.join(
                    combined_dir + ("" if cut_label == "nominal" else f"-{cut_label}"),
                    f"{_cum_tag(0)[:3]}*",
                    f"{var_save_name}_topology{fig_ext}",
                )
                print(f"  - {var_save_name} [{cut_label}]: {hint}")
