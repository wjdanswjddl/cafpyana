#!/usr/bin/env python
"""Legacy orchestrator for MCstat / bundled Flux / bundled G4 (+ cosmics).

**Prefer for new work**

- Monolithic **MCstat + Flux + G4** with two metrics (total rate + background-subtracted):
  ``get_systematics_multisim.py``.
- Cosmics unisim alone: ``get_systematics_cosmics.py``.
- **Chunked** neutrino multisim (HDF splits / grid jobs): ``syst_multisim_chunk.py``
  then ``syst_multisim_aggregate.py``, or ``run_syst_multisim_chunked.sh``.

This script remains useful when you need **exactly** the older single-metric
(``bkgd_subtract=True`` only) covariance workflow or the convenience wrapper
``--chunks-dir`` → ``syst_multisim_aggregate.py``.
"""

import argparse
import subprocess
import sys
import pandas as pd
import numpy as np
from os import path, makedirs
from datetime import datetime
from tqdm import tqdm

# local imports
sys.path.append('../../../')
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    build_var_configs,
    drop_bad_g4_weights,
    save_legacy_category_npzs,
)
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.constants import *
from analysis_village.numucc_1p0pi.files_config import *
from pyanalib.covariance import *

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

save_result = True
save_fig = save_result

# Set in ``main()`` before processing (plot helpers read this global).
save_fig_dir = path.join(save_fig_base_dir, "systematics-placeholder")


def _syst_plot_tag(syst_name):
    if isinstance(syst_name, tuple):
        return str(syst_name[1])
    return str(syst_name)


# ===== functions to process systematics =====
def process_systematics(mc_evt_df, var_config, syst_name, syst_dict):
    univ_events, cv_events = get_univ_rates(evtdf=mc_evt_df,
                                            var_config=var_config,
                                            n_univ=100,
                                            bkgd_subtract=True,
                                            syst_name=syst_name)
                                        
    ret = get_covariance_matrix(univ_events, cv_events)

    tag = _syst_plot_tag(syst_name)
    save_fig_name = "{}/{}-{}-universes".format(save_fig_dir, var_config.var_save_name, tag)
    plot_univ_hists(univ_events, cv_events, syst_name, var_config, save_fig=save_fig, save_name=save_fig_name)

    # frac_unc = np.sqrt(np.diag(ret["cov_frac"]))
    # plot_frac_unc([frac_unc], var_config)

    for matrix_type in ["cov", "cov_frac", "corr"]:
        print(f"plotting {matrix_type} matrix for {syst_name}")
        save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-{tag}-{matrix_type}.pdf"
        plot_heatmap(ret[matrix_type], 
                    var_config.bins, 
                    plot_labels=[var_config.var_labels[1], var_config.var_labels[1], f"{matrix_type.capitalize()}"],
                    save_fig=save_fig, save_name=save_fig_name)

    if isinstance(syst_name, str):
        syst_dict[syst_name][var_config.var_save_name] = ret
    else:
        syst_dict[syst_name[1]][var_config.var_save_name] = ret
    return syst_dict

def process_systematics_cosmics(nu_mc_df, offbeam_data_df, intime_mc_df, var_config, syst_dict):
    univ_events, cv_events = get_univ_rates(evtdf=nu_mc_df,
                                            var_config=var_config,
                                            n_univ=1,
                                            bkgd_subtract=True,
                                            syst_name=("mc", "Flux"))

    data_events, _ = np.histogram(offbeam_data_df[var_config.var_evt_reco_col], bins=var_config.bins)
    mc_events, _   = np.histogram(intime_mc_df[var_config.var_evt_reco_col], weights=intime_mc_df["pot_scale"], bins=var_config.bins)

    if var_config.var_save_name == "integrated":
        data_events = np.array([len(offbeam_data_df)])
        mc_events = np.array([len(intime_mc_df)])
        print(data_events, mc_events)

    # TODO: fill 0 with 1
    mc_events[mc_events == 0] = 1
    data_events[data_events == 0] = 1

    fig, ax = plt.subplots()
    plt.hist(var_config.bin_centers, weights=data_events, bins=var_config.bins, histtype="step", color="red", label="Data")
    plt.hist(var_config.bin_centers, weights=mc_events, bins=var_config.bins, histtype="step", color="black", label="MC")
    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel("Events / Bin")
    plt.legend()
    save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-cosmics-universes.pdf"
    plt.savefig(save_fig_name, bbox_inches="tight", dpi=dpi)
    plt.close()

    # treat data as a unisim systematic universe
    # the MC-data difference is the variation in event rate
    syst_name = "cosmics"
    univ_events =np.array([cv_events + (data_events - mc_events)])
    ret = get_covariance_matrix(univ_events, cv_events)

    # plot_univ_hists(univ_events, cv_events, syst_name, var_config)
    # frac_unc = np.sqrt(np.diag(ret["cov_frac"]))
    # plot_frac_unc([frac_unc], var_config)

    matrix_type = "cov"
    save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name}-{matrix_type}.pdf"
    plot_heatmap(ret[matrix_type], 
                var_config.bins, 
                plot_labels=[var_config.var_labels[1], var_config.var_labels[1], "Covariance"],
                save_fig=save_fig, save_name=save_fig_name)

    syst_dict[syst_name][var_config.var_save_name] = ret
    return syst_dict

def save_syst_dict(syst_dict, save_filename):
    print("saving dict with keys: ", syst_dict.keys())
    print("for systs: ", syst_dict[list(syst_dict.keys())[0]].keys())
    print("saving syst_dict as npz in %s" % (save_filename))
    np.savez(save_filename, **syst_dict)


def parse_args():
    p = argparse.ArgumentParser(
        description="Flux / G4 / MCstat (+ cosmics) covariance matrices for plotting pipelines.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--mc-df-stage",
        choices=("final", "sel_all"),
        default="final",
        help='"final" uses tight-selection MC CAFs (default chunk tags); '
        '"sel_all" uses loose dfs + wgts (same convention as event_selection MC).',
    )
    p.add_argument(
        "--var-set",
        choices=("final", "intermediate", "both"),
        default="final",
        help="Which VariableConfig set to propagate through uncertainties.",
    )
    p.add_argument(
        "--out-tag",
        default=None,
        help="Output subdirectory name under plots base (default: today + stage tags).",
    )
    p.add_argument("--no-plots", action="store_true", help="Skip figures (covariance only).")
    p.add_argument("--no-legacy-npz", action="store_true", help="Skip per-category NPZs used by selected_events.py.")
    p.add_argument(
        "--chunks-dir",
        default=None,
        help="If set, skip monolithic MC load and run syst_multisim_aggregate.py on nu__*.pkl chunks only.",
    )
    p.add_argument(
        "--skip-cosmics",
        action="store_true",
        help="With --chunks-dir: skip cosmics block in aggregate (nu uncertainties only).",
    )
    return p.parse_args()


def main():
    global save_fig_dir, save_fig

    args = parse_args()
    today_str = datetime.now().strftime("%Y%m%d")
    if args.out_tag:
        sub = args.out_tag
    else:
        mc_part = "mc_final" if args.mc_df_stage == "final" else "mc_selall"
        sub = "{}_{}_{}".format(today_str, mc_part, args.var_set)
    save_fig_dir = path.join(save_fig_base_dir, "systematics-{}".format(sub))

    if save_fig and not args.no_plots:
        if not path.exists(save_fig_dir):
            makedirs(save_fig_dir)
        print("saving plots in ", save_fig_dir)
    elif args.no_plots:
        save_fig = False
        if not path.exists(save_fig_dir):
            makedirs(save_fig_dir)

    if args.chunks_dir:
        agg_py = path.join(path.dirname(path.abspath(__file__)), "syst_multisim_aggregate.py")
        cmd = [
            sys.executable,
            agg_py,
            "--chunks_dir",
            args.chunks_dir,
            "--syst-disk-root",
            save_fig_dir,
            "--mc-df-stage",
            args.mc_df_stage,
            "--var-set",
            args.var_set,
        ]
        if args.skip_cosmics:
            cmd.append("--skip-cosmics")
        if args.no_plots:
            cmd.append("--no-plots")
        if args.no_legacy_npz:
            cmd.append("--no-legacy-npz")
        print("[get_systematics] chunked mode:", " ".join(cmd))
        subprocess.check_call(cmd)
        return

    save_filename = path.join(save_fig_dir, "syst_dict.npz")

    var_configs = build_var_configs(args.var_set)

    if args.mc_df_stage == "final":
        ret = get_ana_dfs(option="systs", systs_mc_df_tag="", systs_chunk_tags=None)
    else:
        ret = get_ana_dfs(
            option="systs",
            systs_mc_df_tag="-sel_all-wgts",
            systs_chunk_tags=generate_tags("ah")[1:],
        )

    mc_evt_df = ret["evt"]
    mc_evt_df = drop_bad_g4_weights(mc_evt_df)

    syst_names = ["MCstat", "Flux", "G4"]
    syst_dict = {name: {} for name in syst_names}

    print("Processing neutrino systematics...", "n_evts=", len(mc_evt_df))
    for syst_name_outer in tqdm(syst_names):
        syst_key = ("mc", syst_name_outer) if syst_name_outer in ("Flux", "G4") else syst_name_outer
        for var_config in tqdm(var_configs):
            syst_dict = process_systematics(mc_evt_df, var_config, syst_key, syst_dict)
        save_syst_dict(syst_dict, save_filename)

    print("Processing cosmics systematics...")
    dfs = get_ana_dfs(option="cosmics_systs")
    offbeam_data_df = dfs["data"]
    intime_mc_df = dfs["mc"]

    syst_dict["cosmics"] = {}
    for var_config in tqdm(var_configs):
        syst_dict = process_systematics_cosmics(
            mc_evt_df, offbeam_data_df, intime_mc_df, var_config, syst_dict
        )
        save_syst_dict(syst_dict, save_filename)

    if not args.no_legacy_npz:
        save_legacy_category_npzs(syst_dict, save_fig_dir)
        print("Wrote syst_disk_layout NPZs under", save_fig_dir)


if __name__ == "__main__":
    main()
