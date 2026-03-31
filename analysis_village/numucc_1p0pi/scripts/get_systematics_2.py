import pandas as pd
import numpy as np
from os import path, makedirs
from datetime import datetime

# local imports
import sys
sys.path.append('../../../')
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from pyanalib.covariance import *
from makedf.mcstat import get_MCstat_unc
from makedf.geniesyst import *

# turn off PerformanceWarning 
# triggered by mismatched column levels
import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


def get_systematics(mc_evt_df, mc_nu_df, var_config, syst_name, syst_type="genie", plot=False):
    # ===== get universes =====
    # Note: these are background subtracted
    matrices = {}
    n_univ = len([c for c in mc_evt_df[syst_name].columns if "univ" in c[0]])
    for cov_type in ["xsec", "rate"]:
        univ_events, cv_events = get_univ_rates(cov_type, mc_evt_df, mc_nu_df, var_config, syst_name, n_univ)
        ret = get_covariance_matrix(univ_events, cv_events)
        save_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name[1]}_{cov_type}-univ_hists"
        plot_univ_hists(univ_events, cv_events, syst_name, var_config, plot=plot, save_fig=save_fig, save_name=save_name)
        matrices[cov_type] = ret

        # ===== plot covariance matrices =====
        matrix_type = "cov"
        plot_labels = [var_config.var_labels[1], var_config.var_labels[1], "Covariance"]
        save_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name[1]}_{cov_type}-{matrix_type}"
        plot_heatmap(ret[matrix_type], 
                    var_config.bins, 
                    plot_labels=plot_labels,
                    plot=plot,
                    save_fig=save_fig, save_name=save_name)

        matrix_type = "cov_frac"
        plot_labels = [var_config.var_labels[1], var_config.var_labels[1], "Fractional Covariance"]
        save_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name[1]}_{cov_type}-{matrix_type}"
        plot_heatmap(ret[matrix_type], 
                    var_config.bins, 
                    plot_labels=plot_labels,
                    plot=plot,
                    save_fig=save_fig, save_name=save_name)

        matrix_type = "corr"
        plot_labels = [var_config.var_labels[1], var_config.var_labels[1], "Correlation"]
        save_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name[1]}_{cov_type}-{matrix_type}"
        plot_heatmap(ret[matrix_type], 
                    var_config.bins, 
                    plot_labels=plot_labels,
                    plot=plot,
                    save_fig=save_fig, save_name=save_name)

    return matrices


if __name__ == "__main__":
    # ===== save configs =====
    today_str = datetime.now().strftime("%Y%m%d")
    save_result = True
    save_fig = True

    # syst_type = "genie-CCQE"
    # syst_list = qe_genie_systematics
    # df_tag="_geniewgts_CCQE"
    # subdir = "genie_wgts-CCQE"

    # syst_type = "genie-MEC"
    # syst_list = mec_genie_systematics
    # df_tag = ""
    # subdir = "genie_wgts-MEC"

    syst_lists = [ar23p_genie_systematics]
    for gidx, genie_tag in enumerate(["Ar23p"]):

        syst_type = f"genie-{genie_tag}"
        syst_list = syst_lists[gidx]
        df_tag = ""
        subdir = f"genie_wgts-{genie_tag}"

        save_fig_base_dir = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi"
        save_fig_dir = path.join(save_fig_base_dir, "systematics-{}-{}".format(syst_type, today_str))
        if save_fig:
            if not path.exists(save_fig_dir):
                makedirs(save_fig_dir)
            print("saving plots in ", save_fig_dir)

        # ===== files to process =====
        file_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09"
        n_max_concat = 3
        mc_keys2load = ['hdr', 'mcnu', 'evt'] 

        concat_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=["batch1_"+t for t in generate_tags("ah")],
            df_tag=df_tag,
            keys2load=mc_keys2load,
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="BNB_cosmics/"+subdir
        )

        mc_hdr_df = concat_dfs['hdr']
        mc_nu_df  = concat_dfs['mcnu']
        mc_evt_df = concat_dfs['evt']

        # ===== total pot =====
        mc_tot_pot = mc_hdr_df['pot'].sum()
        print("mc_tot_pot: %.3e" %(mc_tot_pot))
        mc_pot_scale = 1.0
        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
        mc_nu_df["pot_weight"]  = mc_pot_scale * np.ones(len(mc_nu_df))

        # ===== variables to process =====
        var_configs = [VariableConfig.all_events(),
                    VariableConfig.muon_momentum(),
                    VariableConfig.muon_direction(),
                    VariableConfig.proton_momentum(),
                    VariableConfig.proton_direction(),
                    #    VariableConfig.opening_angle(),
                    VariableConfig.tki_del_alpha(),
                    VariableConfig.tki_del_phi(),
                    VariableConfig.tki_del_Tp(),
                    VariableConfig.tki_del_p(),
                    VariableConfig.tki_del_Tp_x(),
                    VariableConfig.tki_del_Tp_y()]

        # ===== systs to process =====
        syst_names = [("mc", syst_list[sidx]) for sidx in range(len(syst_list))]

        # ===== process systematics =====
        syst_dict = {} 

        for var_config in var_configs:
            print(f"Processing {var_config.var_save_name}...")
            syst_dict[var_config.var_save_name] = {}
            for syst_name in syst_names:
                print(f"Processing {syst_name[1]}...")
                matrices = get_systematics(mc_evt_df, mc_nu_df, var_config, syst_name)

                if save_result:
                    syst_dict[var_config.var_save_name][syst_name[1]] = matrices

        # ===== save results =====
        # print what we're saving
        print("saving dict with keys: ", syst_dict.keys())
        print("for systs: ", syst_dict[list(syst_dict.keys())[0]].keys())
        save_filename = f"{save_fig_dir}/{syst_type}_syst_dict.npz"
        print("saving syst_dict as npz in %s" % (save_filename))
        np.savez(save_filename, **syst_dict)
