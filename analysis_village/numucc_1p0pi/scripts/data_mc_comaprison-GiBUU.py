# # Data vs GiBUU MC — final selection, per-TPC, beam-quality data
# 
# Overlaid histograms for **fully selected, per-TPC contained** events.
# 
# - **Data**: `evt_good` from `beam_data_1e20_qualitycut.df` (`beam_quality.ipynb`)
# - **MC**: GiBUU `dfs_from_dir` on `sel_mup-mc-GiBUU/merged_perTPC`
# - **MC weights**: `pot_weight = (data_tot_pot / mc_tot_pot) × mc.genweight`
# - **Intime cosmics**: `dfs_from_dir` on `sel_mup-mc-Intime/merged_perTPC` (also loads OffBeamLight for reference)
# - **Plotting**: `overlay_hists` (topology + GENIE sideband breakdowns), same as `data_mc_comparison.ipynb`
# - **Systematics bands**: summed category fractional covariance from `systematics-summary.ipynb` export (`CategorySummary/category_syst_summary.npz`; **`total_rate`** = GENIE rate)
# 


# %%
from os import path, makedirs
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd

import sys
sys.path.append('/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana')
from pyanalib.split_df_helpers import load_dfs
from pyanalib.split_df_helpers_new import dfs_from_dir
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.files_config import *
plt.style.use("presentation.mplstyle")

# %%
# --- paths ---
DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"
QUALITY_DF_DIR = path.join(DFS_ROOT, "2026_05_16_230705__sel_mup-data-1e20/merged_perTPC")
QUALITY_DF_PATH = path.join(QUALITY_DF_DIR, "beam_data_1e20_qualitycut.df")

DIR_MC = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_05_18_104112__sel_mup-mc-GiBUU/merged_perTPC"
DIR_INTIME = path.join(DFS_ROOT, "2026_05_11_083624__sel_mup-mc-Intime/merged_perTPC")
DIR_OFFBEAM = path.join(DFS_ROOT, "2026_05_11_040015__sel_mup-data-OffBeamLight/merged_perTPC")

KEYS2LOAD = ["hdr", "evt"]
N_MAX_CONCAT = 999

FOM_CUT = 0.98
MIN_RUN_DURATION_MIN = 20.0

save_fig = True
today_str = datetime.now().strftime("%Y%m%d")
save_fig_dir = path.join(
    save_fig_base_dir,
    f"data_mc_comparison_gibuu_{today_str}",
)
if save_fig and not path.exists(save_fig_dir):
    makedirs(save_fig_dir)
print("quality-cut df:", QUALITY_DF_PATH)
print("saving plots in", save_fig_dir)

# %% [markdown]
# ## Load MC, intime MC, and OffBeamLight (`dfs_from_dir`)

# %%
from pyanalib.split_df_helpers_new import dfs_from_dir

mc_dfs = dfs_from_dir(
    DIR_MC,
    filename_str="sel_mup",
    keys2load=KEYS2LOAD,
    n_max_concat=N_MAX_CONCAT,
)
intime_dfs = dfs_from_dir(
    DIR_INTIME,
    filename_str="sel_mup-mc-Intime",
    keys2load=KEYS2LOAD,
    n_max_concat=N_MAX_CONCAT,
)
offbeam_dfs = dfs_from_dir(
    DIR_OFFBEAM,
    filename_str="sel_mup-data-OffBeamLight",
    keys2load=KEYS2LOAD,
    n_max_concat=N_MAX_CONCAT,
)

mc_evt_df = mc_dfs["evt"]
mc_hdr_df = mc_dfs["hdr"]
intime_evt_df = intime_dfs["evt"]
intime_hdr_df = intime_dfs["hdr"]
offbeam_evt_df = offbeam_dfs["evt"]
offbeam_hdr_df = offbeam_dfs["hdr"]

print(f"mc:      evt={len(mc_evt_df):,}  hdr={len(mc_hdr_df):,}")
print(f"intime:  evt={len(intime_evt_df):,}  hdr={len(intime_hdr_df):,}")
print(f"offbeam: evt={len(offbeam_evt_df):,}  hdr={len(offbeam_hdr_df):,}")

mc_evt_df.loc[mc_evt_df.mc.iscc.isna(), ("mc", "iscc")] = 999


# %% [markdown]
# ## Load beam-quality-filtered data

# %%
quality_dfs = load_dfs(
    QUALITY_DF_PATH,
    keys2load=["hdr", "trigger", "evt_good"],
    n_max_concat=1,
)
data_evt_df = quality_dfs["evt_good"]
data_hdr_df = quality_dfs["hdr"].join(quality_dfs["trigger"])

data_evt_df[("mc", "iscc")] = 999
print(f"evt_good rows: {len(data_evt_df):,}")

# %%
data_tot_pot = data_hdr_df["pot"].sum() * 0.98
data_gates = data_hdr_df.nbnbinfo.sum()

# %%
pot_str = get_pot_str(data_tot_pot)
pot_label = f"Events / Bin (POT={pot_str})"

data_evt_df["pot_weight"] = np.ones(len(data_evt_df))

mc_tot_pot = mc_hdr_df["pot"].sum()
mc_pot_scale = data_tot_pot / mc_tot_pot
print(f"mc_pot_scale: {mc_pot_scale:.3e}")
gw = np.asarray(mc_evt_df["mc"]["genweight"], dtype=float)
mc_evt_df["pot_weight"] = mc_pot_scale * gw

intime_gates = offbeam_hdr_df[offbeam_hdr_df["first_in_subrun"] == 1]["noffbeambnb"].sum()
f = 0.0753
scale_intime_to_lightdata = (1 - f) * data_gates / intime_gates
print(f"intime data scale: {scale_intime_to_lightdata:.2f}")
intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))


# %% [markdown]
# ## Per-TPC contained cut (MC / intime)
# 
# Data `evt_good` is already from `merged_perTPC` (full selection + per-TPC). Apply the same
# containment cut to MC and intime MC (`dfs_from_dir` prints a reminder to recheck FV / TKI).

# %%
perTPC_inset = 10

def perTPC_cut(df):
    in_TPC1 = (
        InFV(df.slc.vertex, det="SBND_TPC1", incathode=perTPC_inset)
        & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)
        & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=perTPC_inset)
    )
    in_TPC2 = (
        InFV(df.slc.vertex, det="SBND_TPC2", incathode=perTPC_inset)
        & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)
        & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=perTPC_inset)
    )
    return in_TPC1 | in_TPC2

mc_evt_df = mc_evt_df.loc[perTPC_cut(mc_evt_df)]
intime_evt_df = intime_evt_df.loc[perTPC_cut(intime_evt_df)]

print(f"perTPC: data={len(data_evt_df):,} mc={len(mc_evt_df):,} intime={len(intime_evt_df):,}")

# %%
for df in (mc_evt_df, data_evt_df, intime_evt_df):
    df[("mu", "pfp", "trk", "phi", "", "", "")] = np.degrees(
        np.arctan2(
            df[("mu", "pfp", "trk", "dir", "x", "", "")],
            df[("mu", "pfp", "trk", "dir", "y", "", "")],
        )
    )
    df[("p", "pfp", "trk", "phi", "", "", "")] = np.degrees(
        np.arctan2(
            df[("p", "pfp", "trk", "dir", "x", "", "")],
            df[("p", "pfp", "trk", "dir", "y", "", "")],
        )
    )

# %%
from pyanalib.variable_calculator import get_cc1p0pi_tki
from pyanalib.pandas_helpers import pad_column_name

def evt_df_fixed(df):
    slc_mudf = df.mu.pfp.trk
    slc_pdf = df.p.pfp.trk
    tki_reco = get_cc1p0pi_tki(
        slc_mudf,
        slc_pdf,
        pad_column_name(('P', 'p_muon'), slc_mudf),
        pad_column_name(('P', 'p_proton'), slc_pdf),
    )
    df['del_Tp_x'] = tki_reco['del_Tp_x']
    df['del_Tp_y'] = tki_reco['del_Tp_y']

    mc_mudf = df.mu.pfp.trk.truth.p
    mc_pdf = df.p.pfp.trk.truth.p
    tki_mc = get_cc1p0pi_tki(
        mc_mudf,
        mc_pdf,
        pad_column_name(('totp',), mc_mudf),
        pad_column_name(('totp',), mc_pdf),
    )
    df['mc_del_Tp_x'] = tki_mc['del_Tp_x']
    df['mc_del_Tp_y'] = tki_mc['del_Tp_y']
    df[('mc', 'del_Tp_x')] = tki_mc['del_Tp_x']
    df[('mc', 'del_Tp_y')] = tki_mc['del_Tp_y']

    n_before = len(df)
    df = df[np.abs(df.slc.vertex.x) > 10]
    if 'topo_categ' not in df.columns:
        df = df.copy()
        df.loc[:, 'topo_categ'] = get_topo_category(df)
    return df, n_before


mc_evt_df, n_evt_before_fv = evt_df_fixed(mc_evt_df)
data_evt_df, n_evt_before_fv = evt_df_fixed(data_evt_df)
intime_evt_df, n_evt_before_fv = evt_df_fixed(intime_evt_df)
offbeam_evt_df, n_evt_before_fv = evt_df_fixed(offbeam_evt_df)

# %%
from pathlib import Path

from analysis_village.numucc_1p0pi.dataset_locations import PLOTS_BASE
from analysis_village.numucc_1p0pi.syst_disk_layout import category_summary_npz_path

# Same tree as systematics-summary.ipynb (export cell -> CategorySummary/category_syst_summary.npz).
SYST_DISK_ROOT = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final"
CATEGORY_SUMMARY_OUT = (
    Path(PLOTS_BASE) / "syst_uncertainty_breakdown" / "final_selected" / "category_syst_summary.npz"
)
CATEGORY_SUMMARY_NPZ = category_summary_npz_path(SYST_DISK_ROOT)
if not path.isfile(CATEGORY_SUMMARY_NPZ) and CATEGORY_SUMMARY_OUT.is_file():
    CATEGORY_SUMMARY_NPZ = str(CATEGORY_SUMMARY_OUT)

import os

os.environ["NUMUCC_SYST_DISK_ROOT"] = SYST_DISK_ROOT
OVERLAY_SYST_KIND = "rate"  # total_rate in category summary (GENIE rate, not xsec)

print("SYST_DISK_ROOT =", SYST_DISK_ROOT)
print("category summary:", CATEGORY_SUMMARY_NPZ, "exists =", path.isfile(CATEGORY_SUMMARY_NPZ))
if not path.isfile(CATEGORY_SUMMARY_NPZ):
    print(
        "WARNING: run the export cell in systematics-summary.ipynb "
        "to write category_syst_summary.npz before plotting syst bands."
    )

# %%
ratio = True
approval = "internal"
textloc = [0.03, 0.55]
ax_ylim_ratio = 1.9

data_vs_mc_plotter = partial(
    overlay_hists,
    mc_df=mc_evt_df,
    data_df=data_evt_df,
    intime_df=intime_evt_df,
    dirt_df=None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig,
    syst=None,
    syst_kind=OVERLAY_SYST_KIND,
    syst_disk_root=SYST_DISK_ROOT,
    category_syst_summary_path=CATEGORY_SUMMARY_NPZ,
    load_syst_from_summary=True,
)

# %% [markdown]
# ## Overlay histograms

# %%
var_configs = [
    VariableConfig.all_events(),
    VariableConfig.muon_momentum(),
    VariableConfig.muon_direction(),
    VariableConfig.proton_momentum(),
    VariableConfig.proton_direction(),
    VariableConfig.tki_del_Tp(),
    VariableConfig.tki_del_Tp_x(),
    VariableConfig.tki_del_Tp_y(),
    VariableConfig.tki_del_p(),
    VariableConfig.tki_del_alpha(),
    VariableConfig.tki_del_phi(),
]

for var_config in var_configs:
    for breakdown_type in ["topology", "genie_sb"]:
        plot_labels_hist = [var_config.var_labels[1], "Events / Bin", ""]
        ret = data_vs_mc_plotter(
            breakdown_type=breakdown_type,
            var_config=var_config,
            plot_labels=plot_labels_hist,
            textchi2=True,
            save_name=path.join(
                save_fig_dir,
                f"{var_config.var_save_name}_{breakdown_type}",
            ),
        )

# %% [markdown]
# ## Overlay with background GENIE uncertainty
# 
# Same overlays as above, with an additional hatched band for **GENIE background-rate** uncertainty from `systematics-genie-SB.ipynb` (`genie_bkgd_rate` in `cov_mat_dict.pkl`). The band width is \(\sqrt{\mathrm{diag}(C_\mathrm{frac})}\times N_\mathrm{bkgd}\) per bin (non-signal topologies), drawn at the total MC prediction in **dark orange** (`+++` hatch) over the gray total-systematics band.

# %%
from analysis_village.numucc_1p0pi.syst_disk_layout import FILE_GENIE, SUB_GENIE, category_out_dir
from analysis_village.numucc_1p0pi.utils import resolve_genie_sb_cov_mat_pkl

GENIE_SB_SYST_DISK_ROOT = path.join(save_fig_base_dir, "systematics-notebook-genie-SB-integrated")
GENIE_SB_COV_MAT_PKL = resolve_genie_sb_cov_mat_pkl(
    path.join(category_out_dir(GENIE_SB_SYST_DISK_ROOT, SUB_GENIE), FILE_GENIE)
)
print("GENIE SB cov_mat_dict:", GENIE_SB_COV_MAT_PKL, "exists =", path.isfile(GENIE_SB_COV_MAT_PKL))

data_vs_mc_plotter_bkgd = partial(
    overlay_hists,
    mc_df=mc_evt_df,
    data_df=data_evt_df,
    intime_df=intime_evt_df,
    dirt_df=None,
    ax_ylim_ratio=ax_ylim_ratio,
    ratio=ratio,
    textloc=textloc,
    approval=approval,
    save_fig=save_fig,
    syst=None,
    syst_kind=OVERLAY_SYST_KIND,
    syst_disk_root=SYST_DISK_ROOT,
    category_syst_summary_path=CATEGORY_SUMMARY_NPZ,
    load_syst_from_summary=True,
    show_bkgd_syst_band=True,
    genie_sb_cov_mat_pkl=GENIE_SB_COV_MAT_PKL,
    bkgd_syst_band_color="darkorange",
)

# %%
for var_config in var_configs:
    for breakdown_type in ["topology", "genie_sb"]:
        plot_labels_hist = [var_config.var_labels[1], "Events / Bin", ""]
        ret = data_vs_mc_plotter_bkgd(
            breakdown_type=breakdown_type,
            var_config=var_config,
            plot_labels=plot_labels_hist,
            textchi2=True,
            save_name=path.join(
                save_fig_dir,
                f"{var_config.var_save_name}_{breakdown_type}_bkgd_syst",
            ),
        )

# %%



