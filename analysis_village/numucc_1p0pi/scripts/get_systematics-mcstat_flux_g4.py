# Get the total systematics for Flux, G4, MCStat

import pandas as pd
import numpy as np
from os import path, makedirs
from datetime import datetime
from tqdm import tqdm

# local imports
import sys
sys.path.append('../../../')
from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from analysis_village.numucc_1p0pi.constants import *
from analysis_village.numucc_1p0pi.files_config import *
from pyanalib.covariance import *
from makedf.mcstat import get_MCstat_unc

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


# ===== save configs =====
save_result = True
save_fig = save_result

today_str = datetime.now().strftime("%Y%m%d")
save_fig_dir = path.join(save_fig_base_dir, "systematics-{}".format(today_str))

if save_fig:
    if not path.exists(save_fig_dir):
        makedirs(save_fig_dir)
    print("saving plots in ", save_fig_dir)

# ===== functions to process systematics =====
def process_systematics(mc_evt_df, var_config, syst_name, syst_dict):
    univ_events, cv_events = get_univ_rates(evtdf=mc_evt_df,
                                            var_config=var_config,
                                            n_univ=100,
                                            bkgd_subtract=True,
                                            syst_name=syst_name)
                                        
    ret = get_covariance_matrix(univ_events, cv_events)

    save_fig_name = "{}/{}-{}-universes".format(save_fig_dir, var_config.var_save_name, syst_name)
    plot_univ_hists(univ_events, cv_events, syst_name, var_config, save_fig=save_fig, save_name=save_fig_name)

    # frac_unc = np.sqrt(np.diag(ret["cov_frac"]))
    # plot_frac_unc([frac_unc], var_config)

    for matrix_type in ["cov", "cov_frac", "corr"]:
        print(f"plotting {matrix_type} matrix for {syst_name}")
        save_fig_name = f"{save_fig_dir}/{var_config.var_save_name}-{syst_name}-{matrix_type}.pdf"
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


if __name__ == "__main__":

    save_filename = f"{save_fig_dir}/syst_dict.npz"
    
    # ===== variables to process =====
    var_configs = [VariableConfig.all_events(),
                    VariableConfig.vertex_x(),
                    VariableConfig.vertex_y(),
                    VariableConfig.vertex_z(),
                    VariableConfig.muon_momentum(),
                    VariableConfig.muon_direction(),
                    VariableConfig.muon_direction_x(),
                    VariableConfig.muon_direction_y(),
                    VariableConfig.proton_momentum(),
                    VariableConfig.proton_direction(),
                    VariableConfig.proton_direction_x(),
                    VariableConfig.proton_direction_y(),
                    VariableConfig.opening_angle(),
                    VariableConfig.tki_del_alpha(),
                    VariableConfig.tki_del_phi(),
                    VariableConfig.tki_del_Tp(),
                    VariableConfig.tki_del_p(),
                    VariableConfig.tki_del_Tp_x(),
                    VariableConfig.tki_del_Tp_y()]

    # ===== NEUTRINO DATA =====
    print("Processing neutrino systematics...")
    # ===== load samples =====
    ret = get_ana_dfs(option="systs")
    mc_hdr_df = ret['hdr']
    mc_evt_df = ret['evt']

    # TODO
    bad_idx = []
    for i in range(100):
        var = mc_evt_df.mc.G4["univ_{}".format(i)]
        if len(var[var > 1e3]) > 0:
            idx = mc_evt_df[var > 1e3].index
            bad_idx.append(idx)
            print(f"dropping {len(bad_idx)} events with G4 weight > 1e3")
            mc_evt_df = mc_evt_df.drop(idx[0])


    #  ===== process systematics =====
    syst_names = ["MCstat", "Flux", "G4"]
    syst_dict = {syst_name: {} for syst_name in syst_names}

    for syst_name in tqdm(syst_names):
        for var_config in tqdm(var_configs):
            if syst_name == "Flux" or syst_name == "G4":
                syst_name = ("mc", syst_name)
            syst_dict = process_systematics(mc_evt_df, var_config, syst_name, syst_dict)
            save_syst_dict(syst_dict, save_filename)



    # ===== COSMIC DATA =====
    print("Processing cosmics systematics...")
    # ===== load samples =====
    dfs = get_ana_dfs(option="cosmics_systs")
    offbeam_data_df = dfs["data"]
    intime_mc_df = dfs["mc"]

    # ===== process systematics =====
    syst_dict["cosmics"] = {}
    for var_config in tqdm(var_configs):
        syst_dict = process_systematics_cosmics(mc_evt_df, offbeam_data_df, intime_mc_df, var_config, syst_dict)
        save_syst_dict(syst_dict, save_filename)
