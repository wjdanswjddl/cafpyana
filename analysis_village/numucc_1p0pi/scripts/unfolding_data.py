#!/usr/bin/env python3
"""Wiener-SVD unfolding on beam-quality data vs GENIE CCQE MC (script version of unfolding-data.ipynb).

Loads ``evt_good`` from ``beam_data_1e20_qualitycut.df``, GENIE CCQE MC via ``dfs_from_dir``,
applies per-TPC containment, converts event rates to cross section, unfolds with category
``total_xsec`` systematics, and optionally overlays GiBUU truth folded with GENIE ``AddSmear``.

Writes ``unfolding_results.npz`` under the output directory when ``--save-result`` is set.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from datetime import datetime
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from analysis_village.numucc_1p0pi.constants import M_AR, N_A, RHO  # noqa: E402
from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_GLOBS, PLOTS_BASE  # noqa: E402
from analysis_village.numucc_1p0pi.files_config import save_fig_base_dir  # noqa: E402
from analysis_village.numucc_1p0pi.syst_disk_layout import category_summary_npz_path  # noqa: E402
from analysis_village.unfolding.wienersvd import WienerSVD  # noqa: E402
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    InFV,
    cov_from_fraccov,
    dpi,
    fig_ext,
    get_category_summary_syst_unc,
    get_chi2,
    get_integrated_flux,
    get_pot_str,
    get_response_matrix,
    get_topo_category,
    overlay_hists,
    plot_heatmap,
    plot_unfolded_result,
    signal_cut,
    signal_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
from analysis_village.flux.raytrace_volume_defs import (  # noqa: E402
    FV_SPLIT_TRUNCY_BOXES,
    RAYTRACE_VOLUME_LABEL,
)
from pyanalib.pandas_helpers import pad_column_name  # noqa: E402
from pyanalib.split_df_helpers import load_dfs  # noqa: E402
from pyanalib.split_df_helpers_new import dfs_from_dir  # noqa: E402
from pyanalib.variable_calculator import get_cc1p0pi_tki  # noqa: E402

warnings.filterwarnings("ignore", category=PerformanceWarning)

DEFAULT_DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"
DEFAULT_QUALITY_DF = (
    f"{DEFAULT_DFS_ROOT}/2026_05_16_230705__sel_mup-data-1e20/merged_perTPC/"
    "beam_data_1e20_qualitycut.df"
)
DEFAULT_GIBUU_DF_DIR = (
    f"{DEFAULT_DFS_ROOT}/2026_05_18_104112__sel_mup-mc-GiBUU/merged_perTPC"
)
DEFAULT_FLUX_FILE = "/exp/sbnd/data/users/munjung/flux/SBND_gsimple_raytrace/Gen1.root"
DEFAULT_SYST_DISK_ROOT = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final"

VARIABLE_BUILDERS = {
    "muon_momentum": VariableConfig.muon_momentum,
    "muon_direction": VariableConfig.muon_direction,
    "proton_momentum": VariableConfig.proton_momentum,
    "proton_direction": VariableConfig.proton_direction,
    "opening_angle": VariableConfig.opening_angle,
    "tki_del_Tp": VariableConfig.tki_del_Tp,
    "tki_del_Tp_x": VariableConfig.tki_del_Tp_x,
    "tki_del_Tp_y": VariableConfig.tki_del_Tp_y,
    "tki_del_p": VariableConfig.tki_del_p,
    "tki_del_alpha": VariableConfig.tki_del_alpha,
    "tki_del_phi": VariableConfig.tki_del_phi,
}

DEFAULT_VARIABLES = [
    "muon_momentum",
    "muon_direction",
    "proton_momentum",
    "proton_direction",
    "tki_del_Tp",
    "tki_del_p",
    "tki_del_alpha",
    "tki_del_phi",
]

_CAT_SYST_LABELS = {
    "flux": "Flux",
    "g4": "G4",
    "mcstat": "MC stat.",
    "detector": "Detector",
    "cosmics": "Cosmics",
    "genie_rate": "GENIE (rate)",
    "genie_xsec": "GENIE (xsec)",
    "pot": "Exposure",
    "ntargets": "Targets",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--quality-df-path", default=DEFAULT_QUALITY_DF)
    p.add_argument(
        "--genie-dir",
        default=None,
        help="GENIE CCQE merged_perTPC dir (default: dirname of GENIE_GROUP_GLOBS['CCQE'])",
    )
    p.add_argument("--dfs-root", default=DEFAULT_DFS_ROOT)
    p.add_argument("--flux-file", default=DEFAULT_FLUX_FILE)
    p.add_argument("--syst-disk-root", default=DEFAULT_SYST_DISK_ROOT)
    p.add_argument("--category-summary-npz", default=None)
    p.add_argument("--syst-kind", default="xsec", choices=("xsec", "rate"))
    p.add_argument("--output-dir", default=None, help="Plot/NPZ output dir (default: dated under save_fig_base_dir)")
    p.add_argument("--save-result", action="store_true", default=True)
    p.add_argument("--no-save-result", action="store_false", dest="save_result")
    p.add_argument("--save-fig", action="store_true", default=None,
                   help="Save plots (default: same as --save-result)")
    p.add_argument("--no-plots", action="store_true", help="Skip all matplotlib figures")
    p.add_argument("--plot-flux", action="store_true", help="Show/save integrated-flux diagnostic plot")
    p.add_argument("--variables", default=",".join(DEFAULT_VARIABLES),
                   help="Comma-separated var_save_name list")
    p.add_argument("--load-gibuu", action="store_true", help="Load GiBUU mcnu and add comparison plots/NPZ fields")
    p.add_argument("--gibuu-df-dir", default=DEFAULT_GIBUU_DF_DIR)
    p.add_argument("--data-pot-scale", type=float, default=0.9822,
                   help="Multiply data hdr POT sum by this factor")
    p.add_argument("--c-type", type=int, default=2)
    p.add_argument("--norm-type", type=float, default=0.0)
    p.add_argument("--covariance-scale", type=float, default=1.1,
                   help="Multiplier on syst covariance passed to WienerSVD")
    p.add_argument("--per-tpc-inset", type=float, default=10.0)
    return p.parse_args(argv)


def resolve_category_summary_npz(args: argparse.Namespace) -> str:
    if args.category_summary_npz:
        return args.category_summary_npz
    npz = category_summary_npz_path(args.syst_disk_root)
    fallback = (
        Path(PLOTS_BASE) / "syst_uncertainty_breakdown" / "final_selected" / "category_syst_summary.npz"
    )
    if not os.path.isfile(npz) and fallback.is_file():
        return str(fallback)
    return npz


def per_tpc_cut(df: pd.DataFrame, inset: float) -> pd.Series:
    in_tpc1 = (
        InFV(df.slc.vertex, det="SBND_TPC1", incathode=inset)
        & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=inset)
        & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=inset)
    )
    in_tpc2 = (
        InFV(df.slc.vertex, det="SBND_TPC2", incathode=inset)
        & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=inset)
        & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=inset)
    )
    return in_tpc1 | in_tpc2


def evt_df_fixed(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    slc_mudf = df.mu.pfp.trk
    slc_pdf = df.p.pfp.trk
    tki_reco = get_cc1p0pi_tki(
        slc_mudf,
        slc_pdf,
        pad_column_name(("P", "p_muon"), slc_mudf),
        pad_column_name(("P", "p_proton"), slc_pdf),
    )
    df["del_Tp_x"] = tki_reco["del_Tp_x"]
    df["del_Tp_y"] = tki_reco["del_Tp_y"]

    mc_mudf = df.mu.pfp.trk.truth.p
    mc_pdf = df.p.pfp.trk.truth.p
    tki_mc = get_cc1p0pi_tki(
        mc_mudf,
        mc_pdf,
        pad_column_name(("totp",), mc_mudf),
        pad_column_name(("totp",), mc_pdf),
    )
    df["mc_del_Tp_x"] = tki_mc["del_Tp_x"]
    df["mc_del_Tp_y"] = tki_mc["del_Tp_y"]
    df[("mc", "del_Tp_x")] = tki_mc["del_Tp_x"]
    df[("mc", "del_Tp_y")] = tki_mc["del_Tp_y"]

    n_before = len(df)
    df = df[np.abs(df.slc.vertex.x) > 10]
    if "topo_categ" not in df.columns:
        df = df.copy()
        df.loc[:, "topo_categ"] = get_topo_category(df)
    return df, n_before


def compute_xsec_unit(data_tot_pot: float, flux_file: str, plot_flux: bool) -> float:
    print(f"Fiducial volume (FV_split_truncY): {RAYTRACE_VOLUME_LABEL['FV_split_truncY']}")
    v_sbnd = 0.0
    for i, box in enumerate(FV_SPLIT_TRUNCY_BOXES):
        dx = box["x_range"][1] - box["x_range"][0]
        dy = box["y_range"][1] - box["y_range"][0]
        dz = box["z_range"][1] - box["z_range"][0]
        v_box = dx * dy * dz
        v_sbnd += v_box
        print(
            f"  slab {i + 1}: x={box['x_range']} y={box['y_range']} z={box['z_range']} "
            f"-> {v_box:.4e} cm³"
        )
    print(f"V_SBND = {v_sbnd:.6e} cm³")

    integrated_flux_per_pot = get_integrated_flux(flux_file, plot=plot_flux)
    integrated_flux = integrated_flux_per_pot * data_tot_pot
    print(f"Integrated flux × POT = {integrated_flux:.6e} ν/cm²")
    n_targets = (RHO * v_sbnd / M_AR) * N_A
    print(f"# argon target nuclei: {n_targets:.4e}")
    xsec_unit = 1.0 / (integrated_flux * n_targets)
    print(f"xsec_unit = {xsec_unit:.6e} cm²/nucleon")
    return xsec_unit


def pack_unfold_results(unfold: dict, var_config: VariableConfig) -> dict:
    bins = var_config.bins
    bin_widths = np.diff(bins)
    if len(bins) == 2:
        bin_widths = np.array([1.0])
    unfolded = np.asarray(unfold["unfold"], dtype=float)
    return {
        "bins": np.asarray(bins, dtype=float),
        "bin_centers": np.asarray(var_config.bin_centers, dtype=float),
        "bin_widths": bin_widths,
        "unfold": unfolded,
        "unfold_per_bin_width": unfolded / bin_widths,
        "stat_err": np.sqrt(np.maximum(np.diag(unfold["StatUnfoldCov"]), 0.0)),
        "syst_err": np.sqrt(np.maximum(np.diag(unfold["SystUnfoldCov"]), 0.0)),
        "total_err": np.sqrt(np.maximum(np.diag(unfold["UnfoldCov"]), 0.0)),
        "stat_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["StatUnfoldCov"]), 0.0)) / bin_widths,
        "syst_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["SystUnfoldCov"]), 0.0)) / bin_widths,
        "total_err_per_bin_width": np.sqrt(np.maximum(np.diag(unfold["UnfoldCov"]), 0.0)) / bin_widths,
        "AddSmear": np.asarray(unfold["AddSmear"], dtype=float),
        "UnfoldCov": np.asarray(unfold["UnfoldCov"], dtype=float),
        "StatUnfoldCov": np.asarray(unfold["StatUnfoldCov"], dtype=float),
        "SystUnfoldCov": np.asarray(unfold["SystUnfoldCov"], dtype=float),
    }


def mcnu_truth_xsec_model(mcnu_df: pd.DataFrame, var_config: VariableConfig, xsec_unit: float) -> np.ndarray:
    from analysis_village.numucc_1p0pi.utils import get_clipped_evts

    nudf_signal = signal_cut(mcnu_df[mcnu_df.topo_categ == 1])
    var_allmc, wgt_allmc = get_clipped_evts(
        nudf_signal, var_config.var_nu_col, var_config.bins, var_save_name=var_config.var_save_name
    )
    nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=var_config.bins)
    return np.asarray(nevts_allmc, dtype=float) * xsec_unit


def get_syst_unc(
    var_config: VariableConfig,
    data_evt_df: pd.DataFrame,
    *,
    syst_kind: str,
    syst_disk_root: str,
    category_summary_npz: str,
    scale: float = 1.0,
    plot: bool = False,
    save_fig: bool = False,
    save_name: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    from analysis_village.numucc_1p0pi.syst_category_summary import (
        TOTAL_RATE,
        TOTAL_XSEC,
        load_category_syst_summary,
    )

    frac_uncert_total, frac_cov_matrix_total = get_category_summary_syst_unc(
        var_config,
        syst_kind=syst_kind,
        syst_disk_root=syst_disk_root,
        category_syst_summary_path=category_summary_npz,
    )
    total_key = TOTAL_XSEC if syst_kind == "xsec" else TOTAL_RATE

    if plot:
        summary = load_category_syst_summary(category_summary_npz)
        pack = summary["by_var"][var_config.var_save_name]
        for cat_key, cat_blk in pack["categories"].items():
            w = np.asarray(cat_blk["frac_unc_pct"], dtype=float) / 100.0
            label = _CAT_SYST_LABELS.get(cat_key, cat_key)
            plt.hist(
                var_config.bin_centers,
                bins=var_config.bins,
                weights=w,
                histtype="step",
                linewidth=2,
                label=label,
            )
        tot_w = np.asarray(pack[total_key]["frac_unc_pct"], dtype=float) / 100.0
        plt.hist(
            var_config.bin_centers,
            bins=var_config.bins,
            weights=tot_w,
            histtype="step",
            linewidth=2,
            color="k",
            label="Total Syst.",
        )
        n_data, _ = np.histogram(data_evt_df[var_config.var_evt_reco_col], bins=var_config.bins)
        n_data_err = np.sqrt(n_data * scale)
        data_stat_uncert = n_data_err / (n_data * scale)
        plt.hist(
            var_config.bin_centers,
            bins=var_config.bins,
            weights=data_stat_uncert,
            histtype="step",
            linewidth=2,
            label="Data stat.",
            linestyle="--",
            color="k",
        )
        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.ylim(0, max(frac_uncert_total) * 1.4)
        plt.xlabel(var_config.var_labels[1])
        plt.ylabel("Fractional uncertainty")
        plt.legend(fontsize=11, ncol=3, loc="upper center")
        plt.grid(which="major", linestyle="-", linewidth=0.7, alpha=0.7)
        plt.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
        plt.minorticks_on()
        if save_fig and save_name:
            plt.savefig(save_name + fig_ext, bbox_inches="tight", dpi=dpi)
        plt.show()

    return frac_uncert_total, frac_cov_matrix_total


def load_genie_mc(genie_dir: str, n_max_concat: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genie_dfs = dfs_from_dir(
        genie_dir,
        filename_str="sel_mup-wgts_genie_CCQE",
        keys2load=["hdr", "evt", "mcnu"],
        n_max_concat=n_max_concat,
    )
    mc_evt_df = genie_dfs["evt"]
    mc_hdr_df = genie_dfs["hdr"]
    mc_nu_df = genie_dfs["mcnu"]
    print(
        f"GENIE CCQE: evt={len(mc_evt_df):,}  hdr={len(mc_hdr_df):,}  mcnu={len(mc_nu_df):,}"
    )
    return mc_evt_df, mc_hdr_df, mc_nu_df


def load_beam_data(quality_df_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    quality_dfs = load_dfs(
        quality_df_path,
        keys2load=["hdr", "trigger", "evt_good"],
        n_max_concat=1,
    )
    data_evt_df = quality_dfs["evt_good"]
    data_hdr_df = quality_dfs["hdr"].join(quality_dfs["trigger"])
    data_evt_df[("mc", "iscc")] = 999
    print(f"evt_good rows: {len(data_evt_df):,}")
    return data_evt_df, data_hdr_df


def apply_pot_weights(
    data_evt_df: pd.DataFrame,
    data_hdr_df: pd.DataFrame,
    mc_evt_df: pd.DataFrame,
    mc_hdr_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    data_pot_scale: float,
) -> tuple[float, float, float, str]:
    data_tot_pot = data_hdr_df["pot"].sum() * data_pot_scale
    data_gates = data_hdr_df.nbnbinfo.sum()
    pot_str = get_pot_str(data_tot_pot)
    pot_label = f"Events / Bin (POT={pot_str})"
    data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
    print(f"data_tot_pot: {data_tot_pot:.3e}")
    print(f"data tot gates: {data_gates:.3e}")

    mc_tot_pot = mc_hdr_df["pot"].sum()
    mc_pot_scale = data_tot_pot / mc_tot_pot
    print(f"mc_tot_pot: {mc_tot_pot:.3e}")
    print(f"mc_pot_scale: {mc_pot_scale:.3e}")
    mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
    mc_nu_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_nu_df))

    if "topo_categ" not in mc_evt_df.columns:
        mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
    if "topo_categ" not in mc_nu_df.columns:
        mc_nu_df.loc[:, "topo_categ"] = get_topo_category(mc_nu_df)

    return data_tot_pot, mc_tot_pot, mc_pot_scale, pot_label


def load_gibuu_mcnu(
    gibuu_df_dir: str,
    data_tot_pot: float,
    n_max_concat: int,
) -> tuple[pd.DataFrame | None, float | None]:
    try:
        gibuu_dfs = dfs_from_dir(
            gibuu_df_dir,
            filename_str="sel_mup",
            keys2load=["hdr", "mcnu"],
            n_max_concat=n_max_concat,
        )
        gibuu_hdr_df = gibuu_dfs["hdr"]
        gibuu_mcnu_df = gibuu_dfs["mcnu"]
        gibuu_tot_pot = gibuu_hdr_df["pot"].sum()
        gibuu_pot_scale = data_tot_pot / gibuu_tot_pot
        gw_nu = np.asarray(gibuu_mcnu_df["mc"]["genweight"], dtype=float)
        gibuu_mcnu_df["pot_weight"] = gibuu_pot_scale * gw_nu
        gibuu_mcnu_df.loc[gibuu_mcnu_df["mc"]["iscc"].isna(), ("mc", "iscc")] = 999
        gibuu_mcnu_df.loc[:, "topo_categ"] = get_topo_category(gibuu_mcnu_df)
        print(
            "GiBUU: loaded %d mcnu rows; gibuu_tot_pot=%.3e  gibuu_pot_scale=%.3e"
            % (len(gibuu_mcnu_df), gibuu_tot_pot, gibuu_pot_scale)
        )
        return gibuu_mcnu_df, gibuu_pot_scale
    except Exception as exc:
        print("GiBUU load failed; comparison disabled:", exc)
        return None, None


def save_unfolding_npz(
    path_out: str,
    unfold_results_by_var: dict,
    *,
    manifest: dict,
    xsec_unit: float,
    data_tot_pot: float,
) -> None:
    npz_payload = {
        "manifest": np.array([manifest], dtype=object),
        "xsec_unit": np.float64(xsec_unit),
        "data_tot_pot": np.float64(data_tot_pot),
    }
    for vn, pack in unfold_results_by_var.items():
        for key, val in pack.items():
            npz_payload[f"{vn}::{key}"] = val
    np.savez_compressed(path_out, **npz_payload)
    print("saved unfolding results:", path_out)


def run_unfolding(
    var_configs: list[VariableConfig],
    *,
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    data_evt_df: pd.DataFrame,
    pot_label: str,
    xsec_unit: float,
    unfolding_plotter,
    breakdown_type: str,
    save_fig_dir: str,
    save_fig: bool,
    save_plots: bool,
    syst_kind: str,
    syst_disk_root: str,
    category_summary_npz: str,
    c_type: int,
    norm_type: float,
    cov_scale: float,
) -> tuple[dict, dict, dict]:
    unfold_cache: dict = {}
    unfold_results_by_var: dict = {}

    for var_config in var_configs:
        print("\n===", var_config.var_save_name, "===")
        _, syst_cov_matrix = get_syst_unc(
            var_config,
            data_evt_df,
            syst_kind=syst_kind,
            syst_disk_root=syst_disk_root,
            category_summary_npz=category_summary_npz,
            plot=False,
        )

        plot_labels_hist = [var_config.var_labels[1], pot_label, ""]
        ret = unfolding_plotter(
            var_config=var_config,
            plot_labels=plot_labels_hist,
            syst=syst_cov_matrix,
            save_name=os.path.join(save_fig_dir, f"{var_config.var_save_name}_{breakdown_type}"),
        )

        ret_signal_hists = signal_hists(
            mc_evt_df, mc_nu_df, var_config, mode="unfold", return_data=True, plot=False
        )

        if len(var_config.bins) == 2:
            reco_vs_true = np.array([[1.0]])
        else:
            reco_vs_true, _, _ = np.histogram2d(
                ret_signal_hists["var_sel_truth"],
                ret_signal_hists["var_sel_reco"],
                weights=ret_signal_hists["wgt_sel_truth"],
                bins=var_config.bins,
            )

        if save_plots:
            save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-reco_vs_true"
            plot_heatmap(
                reco_vs_true,
                var_config.bins,
                plot_labels=[var_config.var_labels[2], var_config.var_labels[1], "Smearing"],
                verbose=True,
                save_fig=save_fig,
                save_name=save_fig_name,
            )

        eff = ret_signal_hists["nevts_sel_truth"] / ret_signal_hists["nevts_allmc"]
        print("efficiency (sel_truth/allmc):", eff)

        response = get_response_matrix(reco_vs_true, eff)

        if save_plots:
            save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-response_matrix"
            plot_heatmap(
                response,
                var_config.bins,
                plot_labels=[var_config.var_labels[2], var_config.var_labels[1], "Response"],
                save_fig=save_fig,
                verbose=True,
                save_name=save_fig_name,
            )

        nevts_sel_data = ret["total_data"] - ret["total_mc_bkgd"]
        measured = nevts_sel_data * xsec_unit
        model = ret_signal_hists["nevts_allmc"] * xsec_unit
        covariance = (
            cov_from_fraccov(syst_cov_matrix, ret_signal_hists["nevts_sel_reco"]) * xsec_unit**2
        )
        unfold = WienerSVD(
            response,
            model,
            measured,
            covariance * cov_scale,
            c_type,
            norm_type,
            stat_scaling=xsec_unit,
        )

        if save_plots:
            save_name = f"{save_fig_dir}/{var_config.var_save_name}-unfolded_event_rates-DATA"
            plot_unfolded_result(
                unfold,
                measured,
                {"GENIE": model},
                var_config,
                xsec_unit=xsec_unit,
                save_fig=save_fig,
                save_name=save_name,
                data=True,
                closure_test=False,
            )
            save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-data-add_smear"
            plot_heatmap(
                unfold["AddSmear"],
                var_config.bins,
                plot_labels=[var_config.var_labels[2], var_config.var_labels[1], "$A_c$"],
                save_fig=save_fig,
                save_name=save_fig_name,
            )

        unfold_cache[var_config.var_save_name] = {
            "unfold": unfold,
            "model": model,
            "measured": measured,
            "syst_cov_matrix": syst_cov_matrix,
        }
        unfold_results_by_var[var_config.var_save_name] = {
            **pack_unfold_results(unfold, var_config),
            "measured": np.asarray(measured, dtype=float),
            "model_genie": np.asarray(model, dtype=float),
            "xsec_unit": np.float64(xsec_unit),
        }

    return unfold_cache, unfold_results_by_var


def run_gibuu_comparison(
    var_configs: list[VariableConfig],
    *,
    gibuu_mcnu_df: pd.DataFrame,
    unfold_cache: dict,
    unfold_results_by_var: dict,
    xsec_unit: float,
    save_fig_dir: str,
    save_fig: bool,
    save_plots: bool,
    update_npz_fields: bool,
) -> None:
    for var_config in var_configs:
        cache = unfold_cache[var_config.var_save_name]
        unfold = cache["unfold"]
        model = cache["model"]
        measured = cache["measured"]

        model_gibuu = mcnu_truth_xsec_model(gibuu_mcnu_df, var_config, xsec_unit)
        if update_npz_fields:
            unfold_results_by_var[var_config.var_save_name]["model_gibuu"] = np.asarray(
                model_gibuu, dtype=float
            )
            unfold_results_by_var[var_config.var_save_name]["model_gibuu_folded"] = (
                unfold["AddSmear"] @ model_gibuu
            )

        if not save_plots:
            continue

        plot_unfolded_result(
            unfold,
            measured,
            {"GiBUU": model_gibuu},
            var_config,
            xsec_unit=xsec_unit,
            save_fig=save_fig,
            save_name=f"{save_fig_dir}/{var_config.var_save_name}-unfolded_event_rates-DATA-GiBUU",
            data=True,
            closure_test=False,
        )
        plot_unfolded_result(
            unfold,
            measured,
            {"GENIE": model, "GiBUU": model_gibuu},
            var_config,
            xsec_unit=xsec_unit,
            save_fig=save_fig,
            save_name=f"{save_fig_dir}/{var_config.var_save_name}-unfolded_event_rates-DATA-GENIE-GiBUU",
            data=True,
            closure_test=False,
        )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.no_plots:
        plt.ioff()
        import matplotlib
        matplotlib.use("Agg")

    save_fig = args.save_result if args.save_fig is None else args.save_fig
    save_plots = save_fig and not args.no_plots
    today_str = datetime.now().strftime("%Y%m%d")
    save_fig_dir = args.output_dir or os.path.join(
        save_fig_base_dir, f"unfolding-data-{today_str}"
    )
    if args.save_result and not os.path.isdir(save_fig_dir):
        os.makedirs(save_fig_dir, exist_ok=True)
    unfold_results_npz = os.path.join(save_fig_dir, "unfolding_results.npz")

    genie_dir = args.genie_dir or os.path.dirname(GENIE_GROUP_GLOBS["CCQE"])
    category_summary_npz = resolve_category_summary_npz(args)
    os.environ["NUMUCC_SYST_DISK_ROOT"] = args.syst_disk_root

    print("output dir:", save_fig_dir)
    print("unfolding results:", unfold_results_npz)
    print("quality-cut df:", args.quality_df_path)
    print("GENIE CCQE dir:", genie_dir)
    print("category summary:", category_summary_npz, "exists =", os.path.isfile(category_summary_npz))

    if os.path.isfile(category_summary_npz):
        from analysis_village.numucc_1p0pi.syst_category_summary import load_category_syst_summary

        _syst_summary = load_category_syst_summary(category_summary_npz)
        _manifest = _syst_summary.get("manifest") or {}
        print("category summary schema:", _manifest.get("schema", "?"))
        print("flat terms [%]:", _manifest.get("flat_frac_unc_pct", {}))
        print("variables in NPZ:", len(_syst_summary["by_var"]))
    else:
        print("WARNING: category summary missing — re-run export in systematics-summary.ipynb")

    mc_evt_df, mc_hdr_df, mc_nu_df = load_genie_mc(genie_dir, n_max_concat=999)
    data_evt_df, data_hdr_df = load_beam_data(args.quality_df_path)
    data_tot_pot, mc_tot_pot, _, pot_label = apply_pot_weights(
        data_evt_df,
        data_hdr_df,
        mc_evt_df,
        mc_hdr_df,
        mc_nu_df,
        args.data_pot_scale,
    )

    inset = args.per_tpc_inset
    data_evt_df = data_evt_df.loc[per_tpc_cut(data_evt_df, inset)]
    mc_evt_df = mc_evt_df.loc[per_tpc_cut(mc_evt_df, inset)]
    print(f"perTPC: data={len(data_evt_df):,} mc={len(mc_evt_df):,}")

    mc_evt_df, _ = evt_df_fixed(mc_evt_df)
    data_evt_df, _ = evt_df_fixed(data_evt_df)

    xsec_unit = compute_xsec_unit(data_tot_pot, args.flux_file, args.plot_flux)

    var_names = [v.strip() for v in args.variables.split(",") if v.strip()]
    unknown = [v for v in var_names if v not in VARIABLE_BUILDERS]
    if unknown:
        raise ValueError(f"Unknown variables: {unknown}; valid: {sorted(VARIABLE_BUILDERS)}")
    var_configs = [VARIABLE_BUILDERS[v]() for v in var_names]

    mc_evt_df.loc[mc_evt_df.mc.iscc.isna(), ("mc", "iscc")] = 999
    data_evt_df[("mc", "iscc")] = 999

    breakdown_type = "topology"
    unfolding_plotter = partial(
        overlay_hists,
        breakdown_type=breakdown_type,
        mc_df=mc_evt_df,
        data_df=data_evt_df,
        intime_df=None,
        ax_ylim_ratio=1.6,
        ratio=True,
        textloc=[0.05, 0.55],
        approval="internal",
        save_fig=save_fig,
        syst_kind=args.syst_kind,
        syst_disk_root=args.syst_disk_root,
        category_syst_summary_path=category_summary_npz,
        load_syst_from_summary=True,
    )

    unfold_cache, unfold_results_by_var = run_unfolding(
        var_configs,
        mc_evt_df=mc_evt_df,
        mc_nu_df=mc_nu_df,
        data_evt_df=data_evt_df,
        pot_label=pot_label,
        xsec_unit=xsec_unit,
        unfolding_plotter=unfolding_plotter,
        breakdown_type=breakdown_type,
        save_fig_dir=save_fig_dir,
        save_fig=save_fig,
        save_plots=save_plots,
        syst_kind=args.syst_kind,
        syst_disk_root=args.syst_disk_root,
        category_summary_npz=category_summary_npz,
        c_type=args.c_type,
        norm_type=args.norm_type,
        cov_scale=args.covariance_scale,
    )

    gibuu_mcnu_df = None
    gibuu_pot_scale = None
    if args.load_gibuu:
        gibuu_mcnu_df, gibuu_pot_scale = load_gibuu_mcnu(
            args.gibuu_df_dir, data_tot_pot, n_max_concat=999
        )

    if gibuu_mcnu_df is not None:
        run_gibuu_comparison(
            var_configs,
            gibuu_mcnu_df=gibuu_mcnu_df,
            unfold_cache=unfold_cache,
            unfold_results_by_var=unfold_results_by_var,
            xsec_unit=xsec_unit,
            save_fig_dir=save_fig_dir,
            save_fig=save_fig,
            save_plots=save_plots,
            update_npz_fields=args.save_result,
        )

    if args.save_result and unfold_results_by_var:
        manifest = {
            "schema": "numucc1p0pi_unfolding_results_v1",
            "date": today_str,
            "quality_df": args.quality_df_path,
            "genie_dir": genie_dir,
            "gibuu_df_dir": args.gibuu_df_dir if args.load_gibuu else None,
            "data_tot_pot": float(data_tot_pot),
            "mc_tot_pot": float(mc_tot_pot),
            "gibuu_pot_scale": float(gibuu_pot_scale) if gibuu_pot_scale is not None else None,
            "xsec_unit": float(xsec_unit),
            "flux_file": args.flux_file,
            "flux_fv": "FV_split_truncY",
            "variables": list(unfold_results_by_var.keys()),
            "C_type": args.c_type,
            "Norm_type": args.norm_type,
        }
        save_unfolding_npz(
            unfold_results_npz,
            unfold_results_by_var,
            manifest=manifest,
            xsec_unit=xsec_unit,
            data_tot_pot=data_tot_pot,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
