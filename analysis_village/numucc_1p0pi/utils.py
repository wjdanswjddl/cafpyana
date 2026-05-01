import numpy as np
import pandas as pd
from tqdm import tqdm
import string
import statsmodels.api as sm
import pickle

import sys
sys.path.append('../../')
from pyanalib.split_df_helpers import *
from pyanalib.stat_helpers import *
from pyanalib.covariance import *

from makedf.constants import *
from analysis_village.unfolding.wienersvd import *
from analysis_village.numucc_1p0pi.categories import *
from analysis_village.numucc_1p0pi.constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
from matplotlib.legend import Legend
plt.style.use("presentation.mplstyle")
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

pdg_labels = [r"$\mu^{\pm}$", r"$p$", r"$\pi^{\pm}$", r"Other"]
pdg_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]

dpi = 300
fig_ext = ".png"

# ======= util to load systematic uncertainties ======
def get_syst_unc(var_config, plot=False, save_fig=False, save_name=None):
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

    detector_syst = np.load(f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics_studies_detvar-20260309/detector_syst_dict.npz", allow_pickle=True)
    detector_syst =dict(detector_syst)['detector'].item()[var_config.var_save_name]['cov_frac']

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
    frac_cov_matrix_total = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))
    systs      = [mcstat_syst, genie_syst, flux_syst, g4_syst, cosmics_syst, detector_syst]
    syst_names = ["MC stat.", "GENIE", "Flux", "G4", "Cosmics" , "Detector"]

    for syst_name, syst in zip(syst_names, systs):
        syst_uncert = np.sqrt(np.diag(syst))

        if syst_name == "Cosmics":
            syst_uncert = np.max(syst_uncert) * np.ones(len(var_config.bin_centers))

        frac_uncert_total += syst_uncert ** 2
        frac_cov_matrix_total += syst
        if plot:
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=syst_uncert,   histtype="step", linewidth=2, label=syst_name)

    flat_systs = [pot_frac_unc, ntargets_frac_unc]
    flat_syst_names = ["POT", "Ntargets"]
    for syst_name, syst in zip(flat_syst_names, flat_systs):
        syst_uncert = syst * np.ones(len(var_config.bin_centers))
        frac_uncert_total += syst_uncert ** 2
        frac_cov_matrix_total += np.diag(syst_uncert ** 2)
        if plot:
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=syst_uncert,   histtype="step", linewidth=2, label=syst_name)

    # frac_uncert_total += detvar_syst ** 2

    frac_uncert_total = np.sqrt(frac_uncert_total)
    syst = frac_uncert_total

    if plot:
        plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_total,    histtype="step", linewidth=2, color="k",  label="Total")

        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.ylim(0, max(frac_uncert_total) * 1.4)

        plt.xlabel(var_config.var_labels[1])
        plt.ylabel("Uncertainty [%]")
        plt.legend(fontsize=11, ncol=3, loc="upper center")

        plt.grid(which='major', linestyle='-', linewidth=0.7, alpha=0.7)
        plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
        plt.minorticks_on()

        if save_fig:
            plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

        if not plot:
            plt.close()
        else:
            plt.show();

    return syst, frac_cov_matrix_total


def get_frac_unc(mc_evt_df=None, intime_evt_df=None, offbeam_evt_df=None, var_config=None):

    # nu mc unc
    var = mc_evt_df[var_config.var_evt_reco_col]
    var, weights = get_clipped_evts(mc_evt_df, var_config.var_evt_reco_col, var_config.bins)

    # frac_unc = np.zeros(len(var_config.bin_centers))
    n_cv, _ = np.histogram(var, bins=var_config.bins)
    n_univs = []
    tot_cov_mat = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))
    for syst in ["G4", "GENIE", "Flux"]:
    # for syst in ["GENIE", "Flux"]:
        for i in range(100):
            weights = mc_evt_df.mc[syst]["univ_{}".format(i)] * 2
            # set nan to 1
            weights = weights.fillna(1)
            # clip at 5
            weights = np.clip(weights, 0, 5)
            n_univ, _ = np.histogram(var, bins=var_config.bins, weights=weights)
            n_univs.append(n_univ)

        cov_mat = np.cov(np.array(n_univs).T)
        tot_cov_mat += cov_mat
        nu_frac_unc = np.sqrt(np.diag(cov_mat / n_cv**2))
        nu_frac_unc = np.where(n_cv == 0, 0, nu_frac_unc)

    tot_cov_mat = fraccov_from_cov(tot_cov_mat, n_cv)

    # cosmic mc unc is the difference between the intime mc and offbeam data
    intime_var, intime_weights = get_clipped_evts(intime_evt_df, var_config.var_evt_reco_col, var_config.bins)
    offbeam_var, offbeam_weights = get_clipped_evts(offbeam_evt_df, var_config.var_evt_reco_col, var_config.bins)
    n_intime, _ = np.histogram(intime_var, bins=var_config.bins)
    n_offbeam, _ = np.histogram(offbeam_var, bins=var_config.bins)
    cosmic_mc_unc = np.sqrt(np.diag(np.cov(np.array([n_intime, n_offbeam]).T) / n_intime**2))
    cosmic_mc_unc = np.where(n_intime == 0, 0, cosmic_mc_unc)
    # also multiply to the cosmic component of nu mc
    var_numc_cosmics = var[IsCosmic(mc_evt_df)]
    n_numc_cosmics, _ = np.histogram(var_numc_cosmics, bins=var_config.bins)

    cosmic_cov = np.diag((cosmic_mc_unc * (n_offbeam+n_numc_cosmics) / n_cv) ** 2)
    tot_cov_mat += cosmic_cov


    tot_frac_unc = (nu_frac_unc * n_cv + cosmic_mc_unc * n_intime + cosmic_mc_unc * n_numc_cosmics) / (n_cv + n_intime)
    # tot_frac_unc = (nu_frac_unc * n_cv + cosmic_mc_unc * n_intime) / (n_cv + n_intime)
    return tot_frac_unc, tot_cov_mat


#  ====== calculation functions ======
def get_pot_str(tot_pot):
    pot_str = "{:.2e}".format(tot_pot).replace("e+0", "e").replace("e+", "e")\
        .replace("e", " $\\times 10^{") + "}$"
    pot_str = pot_str.replace(".00", "")
    pot_str = pot_str.replace("$\\times 10^{", "$\\times 10^{")  
    pot_str = pot_str.replace("}$", "}") 
    pot_str = pot_str if pot_str.endswith("$") else pot_str + "$"
    pot_str = pot_str.replace(" $", "$")
    # print(pot_str)
    return pot_str


def generate_tags(end_tag=""):
    tags = []
    for first in string.ascii_lowercase:
        for second in string.ascii_lowercase:
            tag = first + second
            if tag == end_tag:
                break
            tags.append(tag)
        if tag == end_tag:
            break
    return tags


def get_clipped_evts(df, var_col, bins, verbose=False):
    var = df[var_col]
    var = np.clip(var, bins[0], bins[-1] - EPSILON)

    if 'pot_weight' in df.columns:
        weights = df.loc[:, 'pot_weight']
    else:
        if verbose:
            print("No pot_weight column found, return 1 as pot scale (expected for data)")
        weights = np.ones_like(var)
    return var, weights

def get_eff_err(success,total):  # success/total
    err = [[],[]]
    eff = success/total
    for i in range(len(success)):
        this_success = success[i]
        this_tot = total[i]
        interval = sm.stats.proportion_confint(this_success,this_tot,method='wilson')
        err[0].append(abs(eff[i]-interval[0]))
        err[1].append(abs(eff[i]-interval[1]))
    return err


def get_univ_rates(cov_type="rate", 
                    syst_type="GENIE",
                    evtdf=None, 
                    nudf=None, 
                    var_config=None, 
                    syst_name="", 
                    n_univ=100, 
                    bkgd_subtract=True,
                    xsec_unit=0,
                    plot=False):
    """
    for the GENIE uncertainty on the xsec measurement
    """

    if cov_type == "xsec":
        print("getting {} universes for {} uncertainty on the xsec".format(n_univ, syst_name))
        scale_factor = xsec_unit
        if xsec_unit == 0:
            print("pass xsec_unit as an argument to get_univ_rates")

    elif cov_type == "rate":
        print("getting {} universes for {} uncertainty on the event rate".format(n_univ, syst_name))
        scale_factor = 1.0

    else:
        raise ValueError("Invalid covariance type: {}, choose in [xsec, rate]".format(cov_type))

    bins = var_config.bins

    evtdf_signal = evtdf[evtdf.topo_categ == 1]
    # reco variable histogram, topology breakdown
    evtdf_div_topo = [evtdf[evtdf.topo_categ == mode]for mode in topology_list]

    if nudf is not None:
        nudf_signal = nudf[nudf.topo_categ == 1]

    ret = signal_hists(evtdf, nudf, var_config, return_data=True, plot=plot)

    univ_events = []
    univ_effs   = []
    univ_smears = []

    for uidx in tqdm(range(n_univ), desc="Getting universes"):
        univ_col = f"univ_{uidx}"

        # ---- uncertainty on the signal rate ----
        # only consider effect on the response matrix for the signal channel
        if cov_type == "xsec" and syst_type == "GENIE":
            # smearing matrix
            # handle case where there's a single bin, in which case there's no smearing
            if len(bins) == 2:
                reco_vs_true = np.array([[1.0]])
            else:
                reco_vs_true, _, _ = np.histogram2d(ret["var_sel_truth"], 
                                                    ret["var_sel_reco"], 
                                                    weights=ret["wgt_sel_truth"]*evtdf_signal[syst_name][univ_col],
                                                    bins=bins)
            univ_smears.append(reco_vs_true)

            # efficiency
            signal_allmc_univ, _ = np.histogram(ret["var_allmc"],
                                               weights=ret["wgt_allmc"]*nudf_signal[syst_name][univ_col],
                                               bins=bins)
            signal_sel_univ, _ = np.histogram(ret["var_sel_truth"],
                                               weights=ret["wgt_sel_truth"]*evtdf_signal[syst_name][univ_col],
                                               bins=bins)
            eff = signal_sel_univ / signal_allmc_univ
            univ_effs.append(eff)

            response_univ = get_response_matrix(reco_vs_true, eff)
            signal_univ = response_univ @ ret["nevts_allmc"] # note that we multiply the CV signal rate!
            # signal_univ = signal_cv

        elif cov_type == "rate":
            signal_univ, _ = np.histogram(ret["var_sel_reco"], 
                                          weights=ret["wgt_sel_reco"]*evtdf_signal[syst_name][univ_col],
                                          bins=bins)

        else:
            raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

        # TODO: this isn't computationally efficient, but it's useful for debugging
        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        for this_evtdf in evtdf_div_topo[1:]:
            var, wgt = get_clipped_evts(this_evtdf, var_config.var_evt_reco_col, bins)
            univ_wgt = this_evtdf[syst_name][univ_col].copy()
            univ_wgt[np.isnan(univ_wgt)] = 1 ## IMPORTANT: make nan univ_wgt to 1. to ignore them
            background_cv, _   = np.histogram(var, bins=bins, weights=wgt)
            background_univ, _ = np.histogram(var, bins=bins, weights=wgt*univ_wgt)

            if bkgd_subtract:
                signal_univ += (background_univ - background_cv)
            else:
                signal_univ += background_univ

        signal_univ *= scale_factor
        univ_events.append(signal_univ)

    univ_events = np.array(univ_events)

    if bkgd_subtract:
        cv_events = ret["nevts_sel_reco"]
        cv_events *= scale_factor
    else:
        cv_events = ret["nevts_allsel_reco"]
        cv_events *= scale_factor 

    return univ_events, cv_events


def get_response_matrix(reco_vs_true, 
                        eff):
    denom = reco_vs_true.T.sum(axis=0)
    num = reco_vs_true.T
    response = np.divide(
        num * eff, denom,
        out=np.zeros_like(num, dtype=float),  # fill with 0 where invalid
        where=denom != 0
    )
    return response



# ====== plotting functions ======

# ==== plot additions ====
def get_textloc_x(values, bins, textloc=[0.05, 0.55]):
    textloc_x, _ = textloc
    n_firsthalf = np.sum(values[:len(bins)//2])
    n_secondhalf = np.sum(values[len(bins)//2:])
    if n_firsthalf < n_secondhalf:
        textloc_x, textloc_ha = textloc_x, 'left'
    else:
        textloc_x, textloc_ha = 1-textloc_x, 'right'
    return textloc_x, textloc_ha


def add_approval_text(approval, textloc_x, textloc_y, textloc_ha):
    if approval == "internal":
        approval_text = r"$\mathbf{SBND}$ Internal"
        textcolor = 'rosybrown'

    elif approval == "preliminary":
        approval_text = r"$\mathbf{SBND}$ Preliminary"
        textcolor = 'gray'

    else:
        return # don't add anything

    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(
        textloc_x, textloc_y, 
        approval_text, 
        transform=ax.transAxes, 
        ha=textloc_ha, va='top',
        fontsize=20, color=textcolor
    )

def add_chi2_text(chi2_val, p_val, ndof, textloc_x, textloc_y, textloc_ha):
    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(textloc_x, textloc_y, 
            f"$\chi^2$/ndof = {chi2_val:.2f}/{ndof} (p-value = {p_val:.2f})",
            transform=ax.transAxes, 
            ha=textloc_ha, va='top',
            fontsize=12, color='black')

def add_genie_version_text(textloc_x, textloc_y, textloc_ha):
    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(textloc_x, textloc_y, 
            r"GENIE v3.4.0 AR23_00i_00_000", 
            transform=ax.transAxes, 
            ha=textloc_ha, va='top',
            fontsize=12, color='gray')

def format_singlebin_plot():
    ax = plt.gcf().axes[0]
    ax.set_xticks([])


# ==== bar plot ====
def bar_plot(breakdown_type="topology", 
             mc_df=None, intime_df=None, dirt_df=None,
             show_plot=True, 
             plot_labels=["", "", ""],
             approval="internal",
             save_fig=False, save_name=None): 

    sizes = [] 
    dfs_list = []
    if mc_df is not None:
        dfs_list.append(mc_df)
    if intime_df is not None:
        dfs_list.append(intime_df)
    if dirt_df is not None:
        dfs_list.append(dirt_df)
    for df in dfs_list:
        scale = df.pot_weight.unique()[0]
        cuts, labels, colors = get_category_cuts(breakdown_type, df, ret_cuts=True)
        # this_df = df[i].groupby(level=[0,1]).head()
        this_sizes = [scale*len(df[i].groupby(level=[0,1]).head()) for i in cuts]
        sizes.append(this_sizes)
    sizes = np.array(sizes).sum(axis=0)

    ncateg = len(cuts)
    fig, ax = plt.subplots(figsize = (6, ncateg*0.6))

    bars = plt.barh(labels, sizes, align='center', color = colors)
    tot_count = np.array(sizes).sum()
    
    perc_list = []
    for bar in bars:
        width = bar.get_width()
        label_y_pos = bar.get_y() + bar.get_height() / 2
        perc = 100*(width+0.)/(tot_count+0.)
        ax.text(width+1, label_y_pos, s= ("%0.1f"%(perc) + "%"), va='center')
        perc_list.append(perc)

    plt.xlabel(plot_labels[0])
    plt.xlim(0, 1.12 * np.max(sizes))

    # === plot additions ====
    add_approval_text(approval, 0.95, 0.5, "right")

    if save_fig:
        plt.savefig(save_name, bbox_inches="tight")

    if not show_plot:
        plt.close()
        
    # # print breakdowns as a table
    # print("Breakdown ", breakdown_type, " ", generator)
    # print("-"*20)
    # for label, perc in zip(labels, perc_list):
    #     print(f"{label}: {perc:.1f}%")
    # print("-"*20)
    # print(f"Total: {np.sum(perc_list):.1f}%")

    ret_dict = {"perc_list": perc_list}
    return ret_dict


# # ==== histograms ====
# def overlay_hists(breakdown_type="topology",
#                   mc_df=None,
#                   data_df=None,
#                   intime_df=None,
#                   dirt_df=None,
#                   var_config="",
#                   plot_labels=["", "", ""],
#                   ax_ylim_ratio=1.5,
#                   ratio = False,
#                   density = False,
#                   syst = None, # fractional cov matrix
#                   textchi2 = False,
#                   vline = None,
#                   textloc=[0.05, 0.55],
#                   approval="internal",
#                   plot=True,
#                   save_fig=False, 
#                   save_name=None): 

#     # ==== prepare dfs for plotting ====

#     # MC
#     if mc_df is not None:

#         # TODO: uncomment this to append dirt_df to mc_df
#         # if dirt_df is not None:
#         #     dirt_df_ = dirt_df.copy()
#         #     # append to mc_df, bump up __ntuple index so that they are unique
#         #     ntuple_vals = mc_df.index.get_level_values(0)
#         #     ntuple_offset = ntuple_vals.max()+1
#         #     names = dirt_df_.index.names
#         #     # __ntuple should be at level 0
#         #     if "__ntuple" in names:
#         #         idx_loc = names.index("__ntuple")
#         #     else:
#         #         idx_loc = 0
#         #     new_tuples = []
#         #     for tup in dirt_df_.index:
#         #         tup = list(tup)
#         #         tup[idx_loc] = tup[idx_loc] + ntuple_offset
#         #         new_tuples.append(tuple(tup))
#         #     dirt_df_.index = pd.MultiIndex.from_tuples(new_tuples, names=names)

#         #     mc_df = pd.concat([mc_df, dirt_df])
        
#         vardf, _        = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)

#         # breakdown MC events into truth categories
#         if breakdown_type == "pdg":
#             # trk breakdown
#             labels = pdg_labels
#             colors = pdg_colors
#             cuts = get_pdg_category(mc_df, ret_cuts=True)

#         elif breakdown_type == "topology":
#             labels = topology_labels
#             colors = topology_colors
#             cuts = get_topo_category(mc_df, ret_cuts=True)

#         elif breakdown_type == "genie":
#             labels = genie_mode_labels
#             colors = genie_mode_colors
#             cuts = get_genie_category(mc_df, ret_cuts=True)

#         elif breakdown_type == "genie_sb":
#             labels = genie_sb_mode_labels
#             colors = genie_sb_mode_colors
#             cuts = get_genie_sb_category(mc_df, ret_cuts=True)
#             # hatches for marking S/B on plot
#             hatches = [None] * len(labels)
#             for i in range(5, len(labels), 2):
#                 hatches[i] = '////'
#         else:
#             raise ValueError("Invalid breakdown_type: %s, please choose between [topology, genie, or genie_sb]" % breakdown_type)
#         var_categ = [vardf[i] for i in cuts]
#         weights_categ = [list(mc_df.loc[cuts[i], 'pot_weight']) for i in range(len(cuts))] 

#         # MC stat err
#         each_mc_hist_data = []
#         each_mc_hist_err2 = []  # sum of squared weights for error
#         for v, w in zip(var_categ, weights_categ):
#             hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
#             hist_err2, _ = np.histogram(v, weights=np.square(w), bins=var_config.bins)
#             each_mc_hist_data.append(hist_vals)
#             each_mc_hist_err2.append(hist_err2)
#         total_mc = np.sum(each_mc_hist_data, axis=0)
#         total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
#         mc_stat_err = np.sqrt(total_mc_err2)

#         # if topology breakdown, add background CV
#         if breakdown_type == "topology":
#             total_mc_bkgd = total_mc - each_mc_hist_data[-1]
#         else:
#             total_mc_bkgd = None

#     else:
#         vardf = None
#         var_categ = None
#         total_mc = None
#         print("No MC data provided")

 
#     # Intime cosmics
#     if intime_df is not None:
#         vardf_intime, _ = get_clipped_evts(intime_df, var_config.var_evt_reco_col, var_config.bins)
#         total_intime, _ = np.histogram(vardf_intime, bins=var_config.bins, weights=intime_df.pot_weight)
#         # var_categ = [vardf_intime] + var_categ
#         # weights_categ = [list(intime_df.pot_weight)] + weights_categ
#         # colors = colors + ["silver"]
#         # labels = labels + ["In-time\nCosmic"]

#         # add to the cosmic item in existing list
#         var_categ[0] = pd.concat([vardf_intime, var_categ[0]])
#         weights_categ[0] = list(intime_df.pot_weight) + list(weights_categ[0])
        
#         total_mc = total_mc + total_intime

#     else:
#         vardf_intime = None
#         print("No intime cosmics provided")


#     # Dirt cosmics
#     if dirt_df is not None:
#         vardf_dirt, _ = get_clipped_evts(dirt_df, var_config.var_evt_reco_col, var_config.bins)
#         total_dirt, _ = np.histogram(vardf_dirt, bins=var_config.bins, weights=dirt_df.pot_weight)
#         var_categ = [vardf_dirt] + var_categ
#         weights_categ = [list(dirt_df.pot_weight)] + weights_categ
#         colors = colors + ["black"]
#         labels = labels + ["Low E\nDirt"]
#         total_mc = total_mc + total_dirt

#     else:
#         vardf_dirt = None
#         print("No dirt cosmics provided")


#    # Data
#     if data_df is not None:
#         vardf_data, _   = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
#         total_data, _ = np.histogram(vardf_data, bins=var_config.bins, weights=data_df.pot_weight)
#         sum_data = np.sum(total_data)
#         data_eylow, data_eyhigh = return_data_stat_err(total_data)

#         # data/MC
#         if total_mc is not None:
#             data_ratio = total_data / total_mc
#             data_ratio_eylow = data_eylow / total_mc
#             data_ratio_eyhigh = data_eyhigh / total_mc
#             data_ratio = np.nan_to_num(data_ratio, nan=-999.)
#             data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
#             data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
#     else:
#         vardf_data = None
#         total_data = None
#         print("No data data provided")



#     # if density is True, area normalize to the data
#     if mc_df is not None and data_df is not None and density == True:
#         mc_area = np.sum(total_mc)

#         if intime_df is not None:
#             intime_area = np.sum(total_intime)
#             mc_area = mc_area + intime_area

#         if dirt_df is not None:
#             dirt_area = np.sum(total_dirt)
#             mc_area = mc_area + dirt_area

#         data_area = np.sum(total_data)
#         density_factor = data_area / mc_area

#         weights_categ = [np.array(w) * density_factor for w in weights_categ]


#     # the order of cuts from get_*_category is reversed from the order of labels and colors
#     colors, labels = colors[::-1], labels[::-1]

#     # ========================================================

#     # ==== plot template ====
#     if ratio:
#         fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
#                                sharex=True, gridspec_kw={'height_ratios': [4, 1]})
#         ax, ax_r = axs[0], axs[1]
#         fig.subplots_adjust(hspace=0.1)
#         ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
#         ax_r.set_ylim(0.4, 1.6)
#         ax_r.set_xlabel(plot_labels[0])
#         ax_r.set_ylabel("Data/MC")
#         ax_r.grid(True)
#         ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)
#         ax_r.minorticks_on()

#     else:
#         fig, ax = plt.subplots(figsize=(8.5, 7))
#         ax.set_xlabel(plot_labels[0])

#     # common formatting
#     ax.set_xlim(var_config.bins[0], var_config.bins[-1])
#     ax.set_ylabel(plot_labels[1])
#     ax.set_title(plot_labels[2])

#     # ==== Plot histograms ====

#     # == rate panel ==
#     # MC
#     if var_categ is not None:
#         mc_stack, _, _ = ax.hist(var_categ,
#                                  weights=weights_categ,
#                                  bins=var_config.bins,
#                                  stacked=True,
#                                  color=colors,
#                                  label=labels,
#                                  linewidth=0,
#                                  edgecolor='none',
#                                  histtype='stepfilled')

#         breakdown_accum = [np.sum(this_mode) for this_mode in mc_stack]
#         breakdown_fractions = [breakdown_accum[0]] + [(breakdown_accum[i+1] - breakdown_accum[i]) for i in range(len(breakdown_accum) - 1)]
#         breakdown_fractions = [frac / breakdown_accum[-1] for frac in breakdown_fractions]

#         if breakdown_type == "genie_sb":
#             # hatch background portion
#             bottom = np.zeros(len(var_config.bins) - 1)
#             for i, (v, w, h) in enumerate(zip(var_categ, weights_categ, hatches)):
#                 hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
#                 ax.bar(
#                     var_config.bin_centers,
#                     hist_vals,
#                     width=np.diff(var_config.bins),
#                     bottom=bottom,
#                     color='none',
#                     hatch=h,
#                     edgecolor='white',
#                     linewidth=0.0,
#                     align='center'
#                 )
#                 bottom += hist_vals

#     if syst is not None: # list of syst uncertainties 
#         mc_stat_err_frac = mc_stat_err / total_mc
#         syst_err_frac = np.sqrt(np.diag(syst))
#         syst_err = np.sqrt(mc_stat_err_frac**2 + syst_err_frac**2) # fractional error
#         syst_err = syst_err * total_mc

#         ax.bar(
#             var_config.bin_centers,
#             2 * syst_err,
#             width=np.diff(var_config.bins),
#             bottom=total_mc - syst_err,
#             facecolor='none',             # transparent fill
#             hatch='xxx',                 # hatch pattern similar to ROOT's 3004
#             linewidth=0.0,
#             edgecolor='dimgray',            # outline color of the hatching
#             label='Syst. Unc.'
#         )

#         if data_df is not None:
#             data_stat_cov_frac = np.diag( (0.5*(data_eyhigh + data_eylow) / total_data ) ** 2)
#             mc_stat_cov_frac = np.diag(mc_stat_err_frac**2)
#             combined_cov_frac = syst + data_stat_cov_frac + mc_stat_cov_frac
#             # combined_cov_frac = data_stat_cov_frac
#             combined_cov = cov_from_fraccov(combined_cov_frac, total_mc)
#             # only calculated chi2 with nonzero bins
#             nonzero_bins = (total_data > 0)
#             chi2_val, p_val = get_chi2(total_data[nonzero_bins], total_mc[nonzero_bins], combined_cov[nonzero_bins, :][:, nonzero_bins])
 
#     else:
#         print("no syst provided")

#     # Data
#     if vardf_data is not None:
#         ax.errorbar(var_config.bin_centers, 
#                     total_data, 
#                     yerr=np.vstack((data_eylow, data_eyhigh)),
#                     color='black', 
#                     fmt='o', markersize=5, capsize=3, linewidth=1.5,
#                     label='Data')

#     # == ratio panel ==
#     if ratio:
#         # MC 
#         if syst is not None:
#             mc_content_ratio = total_mc / total_mc # dummy
#             mc_stat_err_ratio = syst_err / total_mc
#             mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
#             ax_r.bar(
#                 var_config.bin_centers,
#                 2*mc_stat_err_ratio,
#                 width=np.diff(var_config.bins),
#                 bottom=mc_content_ratio - mc_stat_err_ratio,
#                 facecolor='none',
#                 edgecolor='dimgray',
#                 hatch='xxx',
#                 linewidth=0.0,
#                 label='Syst. Unc.'
#             )
#         else:
#             pass

#         # data/MC 
#         if data_df is not None:
#             ax_r.errorbar(var_config.bin_centers, data_ratio, 
#                             yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
#                             fmt='o', color='black',
#                             markersize=5, capsize=3, linewidth=1.5)

#     # ===============================

#     # ==== Legend ====
#     # legend order: data, mc, syst
#     handles, labels_orig = ax.get_legend_handles_labels()
#     ordered_handles = []
#     ordered_labels = []

#     if data_df is not None:
#         data_handle_index = labels_orig.index('Data')
#         data_handle = handles[data_handle_index]
#         ordered_handles.extend([data_handle])
#         data_text = 'Observed ({:.0f})'.format(sum_data)
#         # if textchi2 and syst is not None:
#         #     data_text += f" \n$\chi^2$/ndof = {chi2_val:.2f}/{len(var_config.bins)-1}\np-value = {p_val:.2f}"
#         ordered_labels.extend([data_text])

#     if mc_df is not None:
#         if breakdown_type == "genie_sb":
#             # collapse S and B into a single combined legend
#             for i in range(len(genie_mode_colors)):
#                 ordered_handles.append(Patch(facecolor=genie_mode_colors[i], edgecolor='none'))

#             legend_labels = []
#             i, n = 0, len(labels)
#             while i < n:
#                 # assume paired S/B in the order of MC legend entries
#                 label_base = labels[i]
#                 if (i + 1 < n) and (labels[i + 1] == label_base): # paired S/B
#                     frac1 = breakdown_fractions[i]
#                     frac2 = breakdown_fractions[i+1]
#                     legend_label = f"{label_base} ({frac2*100:.1f}%/{frac1*100:.1f}%)"
#                     i += 2
#                 else: # single
#                     frac1 = breakdown_fractions[i]
#                     legend_label = f"{label_base} ({frac1*100:.1f}%)"
#                     i += 1
#                 legend_labels.append(legend_label)
#             ordered_labels.extend(legend_labels[::-1])

#         else:
#             if data_df is not None:
#                 mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
#             else:
#                 mc_handles = [h for i, h in enumerate(handles) if 'Unc.' not in labels_orig[i]]
#             mc_labels = [f"{label} ({frac*100:.1f}%)"
#                                 for label, frac in zip(labels, breakdown_fractions)]
#             ordered_handles.extend(mc_handles)
#             ordered_labels.extend(mc_labels[::-1]) # note the reverse order of mc_labels


#     if syst is not None:
#         unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
#         unc_label = [l for l in labels_orig if 'Unc.' in l]
#         ordered_handles.extend(unc_handle)
#         ordered_labels.extend(unc_label)

#     # adjust fontsize so that legend fits in the figure
#     fontsize = 11.3
#     ncol = 2
#     if breakdown_type == "genie_sb":
#         textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
#         fontsize = fontsize
#         ncol = 2

#         # separate legend box with S / B hatches
#         example_signal = Patch(facecolor="black", edgecolor='white', label='Signal')
#         example_background = Patch(facecolor="black", edgecolor='white', hatch='////', linewidth=0, label='Background')
#         # if textloc_x < 0.5:
#         #     box_ax = ax.inset_axes([0.17, 0.65, 0.13, 0.13], transform=ax.transAxes)
#         # else:
#         box_ax = ax.inset_axes([0.625, 0.66, 0.13, 0.13], transform=ax.transAxes)
#         box_ax.axis('off')
#         mini_legend = Legend(
#             box_ax,
#             handles=[example_signal, example_background],
#             labels=['Signal', 'Background'],
#             loc='center',
#             fontsize=fontsize,
#             frameon=False,
#             borderpad=0.7,
#             handlelength=2.1,
#             handleheight=0.9,
#             ncol=1,
#             fancybox=True,
#             framealpha=1.0
#         )
#         box_ax.add_artist(mini_legend)

#     ax.legend(
#         ordered_handles,
#         ordered_labels,
#         loc='upper left',
#         # loc='upper center',
#         fontsize=fontsize,
#         frameon=False,
#         ncol=ncol,
#         bbox_to_anchor=(0.05, 0.9, 0.8, 0.1),
#         mode='expand'
#     )

#     # ===============================

#     # ==== plot additions ====

#     # y-axis limit
#     ax.set_ylim(0., ax_ylim_ratio* np.max(total_mc))

#     # vertical lines
#     if vline is not None:
#         for v in vline:
#             ymax = ax.get_ylim()[1]
#             ax.vlines(x=v[0], ymin=0, ymax=ymax*0.75, color='red', linestyle='--')
#             # Plot arrow if v[1] is specified (0: left, 1: right)
#             if len(v) > 1:
#                 direction = v[1]
#                 arrow_params = {
#                     'y': ymax * 0.4,
#                     'dx': 0.18 * (ax.get_xlim()[1] - ax.get_xlim()[0]),  # adjustable length
#                     'width': 0.025 * (ax.get_xlim()[1] - ax.get_xlim()[0]),  # adjustable width
#                     'color': 'red',
#                     'head_width': 0.04 * (ax.get_ylim()[1] - ax.get_ylim()[0]),  # adjustable head width
#                     'head_length': 0.03 * (ax.get_xlim()[1] - ax.get_xlim()[0]),  # adjustable head length
#                     'length_includes_head': True
#                 }
#                 if direction == 0:
#                     # Left arrow
#                     ax.arrow(v[0], arrow_params['y'], -arrow_params['dx'], 0, 
#                              width=arrow_params['width'],
#                              color=arrow_params['color'],
#                              head_width=arrow_params['head_width'],
#                              head_length=arrow_params['head_length'],
#                              length_includes_head=arrow_params['length_includes_head'],
#                              clip_on=True)
#                 elif direction == 1:
#                     # Right arrow
#                     ax.arrow(v[0], arrow_params['y'], arrow_params['dx'], 0, 
#                              width=arrow_params['width'],
#                              color=arrow_params['color'],
#                              head_width=arrow_params['head_width'],
#                              head_length=arrow_params['head_length'],
#                              length_includes_head=arrow_params['length_includes_head'],
#                              clip_on=True)

#     # textboxes
#     textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
#     textloc_y = textloc[1]

#     if textchi2 and syst is not None:
#         add_chi2_text(chi2_val, p_val, len(var_config.bins)-1, textloc_x, textloc_y+0.08, textloc_ha)

#     add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

#     if breakdown_type != "pdg":
#         add_genie_version_text(textloc_x, textloc_y-0.06, textloc_ha)

#     if var_config.var_save_name == "integrated":
#         format_singlebin_plot()

#     # ===============================

#     # == save figure ==
#     if save_fig:
#         plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)

#     if plot == True:
#         plt.show()
#     else:
#         plt.close()

#     return {"cuts": cuts, 
#             "total_mc": total_mc, 
#             "total_mc_bkgd": total_mc_bkgd,
#             "total_data": total_data}


# ==== histograms ====
def overlay_hists(breakdown_type="topology",
                  mc_df=None,
                  data_df=None,
                  intime_df=None,
                  dirt_df=None,
                  var_config="",
                  plot_labels=["", "", ""],
                  ax_ylim_ratio=1.5,
                  ratio = False,
                  density = False,
                  syst = None, # fractional cov matrix
                  syst_decomp = False,
                  textchi2 = False,
                  vline = None,
                  textloc=[0.05, 0.55],
                  approval="internal",
                  plot=True,
                  save_fig=False, 
                  save_name=None): 

    # ==== prepare dfs for plotting ====

    # MC
    if mc_df is not None:

        # TODO: uncomment this to append dirt_df to mc_df
        # if dirt_df is not None:
        #     dirt_df_ = dirt_df.copy()
        #     # append to mc_df, bump up __ntuple index so that they are unique
        #     ntuple_vals = mc_df.index.get_level_values(0)
        #     ntuple_offset = ntuple_vals.max()+1
        #     names = dirt_df_.index.names
        #     # __ntuple should be at level 0
        #     if "__ntuple" in names:
        #         idx_loc = names.index("__ntuple")
        #     else:
        #         idx_loc = 0
        #     new_tuples = []
        #     for tup in dirt_df_.index:
        #         tup = list(tup)
        #         tup[idx_loc] = tup[idx_loc] + ntuple_offset
        #         new_tuples.append(tuple(tup))
        #     dirt_df_.index = pd.MultiIndex.from_tuples(new_tuples, names=names)

        #     mc_df = pd.concat([mc_df, dirt_df])
        
        vardf, _        = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)

        # breakdown MC events into truth categories
        if breakdown_type == "pdg":
            # trk breakdown
            labels = pdg_labels
            colors = pdg_colors
            cuts = get_pdg_category(mc_df, ret_cuts=True)

        elif breakdown_type == "topology":
            labels = topology_labels
            colors = topology_colors
            cuts = get_topo_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie":
            labels = genie_mode_labels
            colors = genie_mode_colors
            cuts = get_genie_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie_sb":
            labels = genie_sb_mode_labels
            colors = genie_sb_mode_colors
            cuts = get_genie_sb_category(mc_df, ret_cuts=True)
            # hatches for marking S/B on plot
            hatches = [None] * len(labels)
            for i in range(5, len(labels), 2):
                hatches[i] = '////'
        else:
            raise ValueError("Invalid breakdown_type: %s, please choose between [topology, genie, or genie_sb]" % breakdown_type)
        var_categ = [vardf[i] for i in cuts]
        weights_categ = [list(mc_df.loc[cuts[i], 'pot_weight']) for i in range(len(cuts))] 

        # MC stat err
        each_mc_hist_data = []
        each_mc_hist_err2 = []  # sum of squared weights for error
        for v, w in zip(var_categ, weights_categ):
            hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
            hist_err2, _ = np.histogram(v, weights=np.square(w), bins=var_config.bins)
            each_mc_hist_data.append(hist_vals)
            each_mc_hist_err2.append(hist_err2)
        total_mc = np.sum(each_mc_hist_data, axis=0)
        total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
        mc_stat_err = np.sqrt(total_mc_err2)

        # if topology breakdown, add background CV
        if breakdown_type == "topology":
            total_mc_bkgd = total_mc - each_mc_hist_data[-1]
        else:
            total_mc_bkgd = None

    else:
        vardf = None
        var_categ = None
        total_mc = None
        print("No MC data provided")

 
    # Intime cosmics
    if intime_df is not None:
        vardf_intime, _ = get_clipped_evts(intime_df, var_config.var_evt_reco_col, var_config.bins)
        total_intime, _ = np.histogram(vardf_intime, bins=var_config.bins, weights=intime_df.pot_weight)
        # var_categ = [vardf_intime] + var_categ
        # weights_categ = [list(intime_df.pot_weight)] + weights_categ
        # colors = colors + ["silver"]
        # labels = labels + ["In-time\nCosmic"]

        # add to the cosmic item in existing list
        var_categ[0] = pd.concat([vardf_intime, var_categ[0]])
        weights_categ[0] = list(intime_df.pot_weight) + list(weights_categ[0])
        
        total_mc = total_mc + total_intime

    else:
        vardf_intime = None
        print("No intime cosmics provided")


    # Dirt cosmics
    if dirt_df is not None:
        vardf_dirt, _ = get_clipped_evts(dirt_df, var_config.var_evt_reco_col, var_config.bins)
        total_dirt, _ = np.histogram(vardf_dirt, bins=var_config.bins, weights=dirt_df.pot_weight)
        var_categ = [vardf_dirt] + var_categ
        weights_categ = [list(dirt_df.pot_weight)] + weights_categ
        colors = colors + ["black"]
        labels = labels + ["Low E\nDirt"]
        total_mc = total_mc + total_dirt

    else:
        vardf_dirt = None
        print("No dirt cosmics provided")


   # Data
    if data_df is not None:
        vardf_data, _   = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
        total_data, _ = np.histogram(vardf_data, bins=var_config.bins, weights=data_df.pot_weight)
        sum_data = np.sum(total_data)
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        # data/MC
        if total_mc is not None:
            data_ratio = total_data / total_mc
            data_ratio_eylow = data_eylow / total_mc
            data_ratio_eyhigh = data_eyhigh / total_mc
            data_ratio = np.nan_to_num(data_ratio, nan=-999.)
            data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
            data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
    else:
        vardf_data = None
        total_data = None
        print("No data data provided")



    # if density is True, area normalize to the data
    if mc_df is not None and data_df is not None and density == True:
        mc_area = np.sum(total_mc)

        if intime_df is not None:
            intime_area = np.sum(total_intime)
            mc_area = mc_area + intime_area

        if dirt_df is not None:
            dirt_area = np.sum(total_dirt)
            mc_area = mc_area + dirt_area

        data_area = np.sum(total_data)
        density_factor = data_area / mc_area

        weights_categ = [np.array(w) * density_factor for w in weights_categ]


    # the order of cuts from get_*_category is reversed from the order of labels and colors
    colors, labels = colors[::-1], labels[::-1]

    # ========================================================

    # ==== plot template ====
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        ax, ax_r = axs[0], axs[1]
        fig.subplots_adjust(hspace=0.1)
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_r.set_ylim(0.5, 1.5)
        ax_r.set_xlabel(plot_labels[0])
        ax_r.set_ylabel("Data/MC")
        ax_r.grid(True)
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)
        ax_r.minorticks_on()

    else:
        fig, ax = plt.subplots(figsize=(8.5, 7))
        ax.set_xlabel(plot_labels[0])

    # common formatting
    ax.set_xlim(var_config.bins[0], var_config.bins[-1])
    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])

    # ==== Plot histograms ====

    # == rate panel ==
    # MC
    if var_categ is not None:
        mc_stack, _, _ = ax.hist(var_categ,
                                 weights=weights_categ,
                                 bins=var_config.bins,
                                 stacked=True,
                                 color=colors,
                                 label=labels,
                                 linewidth=0,
                                 edgecolor='none',
                                 histtype='stepfilled')

        breakdown_accum = [np.sum(this_mode) for this_mode in mc_stack]
        breakdown_fractions = [breakdown_accum[0]] + [(breakdown_accum[i+1] - breakdown_accum[i]) for i in range(len(breakdown_accum) - 1)]
        breakdown_fractions = [frac / breakdown_accum[-1] for frac in breakdown_fractions]

        if breakdown_type == "genie_sb":
            # hatch background portion
            bottom = np.zeros(len(var_config.bins) - 1)
            for i, (v, w, h) in enumerate(zip(var_categ, weights_categ, hatches)):
                hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
                ax.bar(
                    var_config.bin_centers,
                    hist_vals,
                    width=np.diff(var_config.bins),
                    bottom=bottom,
                    color='none',
                    hatch=h,
                    edgecolor='white',
                    linewidth=0.0,
                    align='center'
                )
                bottom += hist_vals

    if syst is not None: # list of syst uncertainties 

        # TODO: diff between data and mc as additional systematic
        # syst_diff = total_data - total_mc
        # syst_diff_cov = np.cov(np.array([total_data, total_mc]).T)

        # decompose into shape and norm components
        cov_norm, cov_mixed, cov_shape = Matrix_Decomp(total_mc, syst * (total_mc**2))
        syst_err_norm = np.sqrt(np.abs(np.diag(cov_norm)))
        syst_err_mixed = np.sqrt(np.abs(np.diag(cov_mixed)))
        syst_err_shape = np.sqrt(np.abs(np.diag(cov_shape)))

        mc_stat_err_frac = mc_stat_err / total_mc
        syst_err_frac = np.sqrt(np.diag(syst))
        syst_err = np.sqrt(mc_stat_err_frac**2 + syst_err_frac**2) # fractional error
        syst_err = syst_err * total_mc

        if syst_decomp == False:
            ax.bar(
                var_config.bin_centers,
                2 * syst_err,
                width=np.diff(var_config.bins),
                bottom=total_mc - syst_err,
                facecolor='none',             # transparent fill
                hatch='xxx',                 # hatch pattern similar to ROOT's 3004
                linewidth=0.0,
                edgecolor='dimgray',            # outline color of the hatching
                label='Syst. Unc.'
            )

        if syst_decomp == True:

            ax.bar(
                var_config.bin_centers,
                2*syst_err_shape,
                width=np.diff(var_config.bins),
                bottom=total_mc - syst_err_shape,
                facecolor='red',
                edgecolor='red',
                alpha=0.3,
                # hatch='xxx',
                linewidth=0.0,
                label='Syst. Unc. (Shape)'
            )

            ax.bar(
                var_config.bin_centers,
                2*syst_err_mixed,
                width=np.diff(var_config.bins),
                bottom=total_mc - syst_err_mixed,
                facecolor='none',
                edgecolor='green',
                hatch='////',
                linewidth=0.0,
                label='Syst. Unc. (Mixed)'
            )

            ax.bar(
                var_config.bin_centers,
                2*syst_err_norm,
                width=np.diff(var_config.bins),
                bottom=total_mc - syst_err_norm,
                facecolor='dimgray',
                edgecolor='dimgray',
                alpha=0.3,
                # hatch='xxx',
                linewidth=0.0,
                label='Syst. Unc. (Norm)'
            )



        if data_df is not None:
            data_stat_cov_frac = np.diag( (0.5*(data_eyhigh + data_eylow) / total_data ) ** 2)
            mc_stat_cov_frac = np.diag(mc_stat_err_frac**2)
            combined_cov_frac = syst + data_stat_cov_frac + mc_stat_cov_frac
            combined_cov = cov_from_fraccov(combined_cov_frac, total_mc)
            chi2_val, p_val = get_chi2(total_data, total_mc, combined_cov)
 
    else:
        print("no syst provided")

    # Data
    if vardf_data is not None:
        ax.errorbar(var_config.bin_centers, 
                    total_data, 
                    yerr=np.vstack((data_eylow, data_eyhigh)),
                    color='black', 
                    fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data')

    # == ratio panel ==
    if ratio:
        # MC 
        if syst is not None:
            if syst_decomp == False:
                mc_content_ratio = total_mc / total_mc # dummy
                mc_stat_err_ratio = syst_err / total_mc
                mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
                ax_r.bar(
                    var_config.bin_centers,
                    2*mc_stat_err_ratio,
                    width=np.diff(var_config.bins),
                    bottom=mc_content_ratio - mc_stat_err_ratio,
                    facecolor='none',
                    edgecolor='dimgray',
                    hatch='xxx',
                    linewidth=0.0,
                    label='Syst. Unc.'
                )

            if syst_decomp == True:
                mc_content_ratio = total_mc / total_mc # dummy
                mc_stat_err_ratio_norm = syst_err_norm / total_mc
                mc_stat_err_ratio_mixed = syst_err_mixed / total_mc
                # print("norm: ", mc_stat_err_ratio_norm)
                mc_stat_err_ratio_shape = syst_err_shape / total_mc
                # print("shape: ", mc_stat_err_ratio_shape)
                mc_stat_err_ratio_norm = np.nan_to_num(mc_stat_err_ratio_norm, nan=0.)
                mc_stat_err_ratio_mixed = np.nan_to_num(mc_stat_err_ratio_mixed, nan=0.)
                mc_stat_err_ratio_shape = np.nan_to_num(mc_stat_err_ratio_shape, nan=0.)

                ax_r.bar(
                    var_config.bin_centers,
                    2*mc_stat_err_ratio_shape,
                    width=np.diff(var_config.bins),
                    bottom=mc_content_ratio - mc_stat_err_ratio_shape,
                    facecolor='red',
                    edgecolor='red',
                    alpha=0.3,
                    # hatch='xxx',
                    linewidth=0.0,
                    label='Syst. Unc. (Shape)'
                )

                ax_r.bar(
                    var_config.bin_centers,
                    2*mc_stat_err_ratio_mixed,
                    width=np.diff(var_config.bins),
                    bottom=mc_content_ratio - mc_stat_err_ratio_mixed,
                    facecolor='none',
                    edgecolor='green',
                    hatch='////',
                    linewidth=0.0,
                    label='Syst. Unc. (Mixed)'
                )

                ax_r.bar(
                    var_config.bin_centers,
                    2*mc_stat_err_ratio_norm,
                    width=np.diff(var_config.bins),
                    bottom=mc_content_ratio - mc_stat_err_ratio_norm,
                    alpha=0.3,
                    facecolor='dimgray',
                    edgecolor='dimgray',
                    # hatch='xxx',
                    linewidth=0.0,
                    label='Syst. Unc. (Norm)'
                )


        else:
            pass

        # data/MC 
        if data_df is not None:
            ax_r.errorbar(var_config.bin_centers, data_ratio, 
                            yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                            fmt='o', color='black',
                            markersize=5, capsize=3, linewidth=1.5)

    # ===============================

    # ==== Legend ====
    # legend order: data, mc, syst
    handles, labels_orig = ax.get_legend_handles_labels()
    ordered_handles = []
    ordered_labels = []

    if data_df is not None:
        data_handle_index = labels_orig.index('Data')
        data_handle = handles[data_handle_index]
        ordered_handles.extend([data_handle])
        data_text = 'Observed ({:.0f})'.format(sum_data)
        # data_text = 'Observed'
        ordered_labels.extend([data_text])

    if mc_df is not None:
        if breakdown_type == "genie_sb":
            # collapse S and B into a single combined legend
            for i in range(len(genie_mode_colors)):
                ordered_handles.append(Patch(facecolor=genie_mode_colors[i], edgecolor='none'))

            legend_labels = []
            i, n = 0, len(labels)
            while i < n:
                # assume paired S/B in the order of MC legend entries
                label_base = labels[i]
                if (i + 1 < n) and (labels[i + 1] == label_base): # paired S/B
                    frac1 = breakdown_fractions[i]
                    frac2 = breakdown_fractions[i+1]
                    legend_label = f"{label_base} ({frac2*100:.1f}%/{frac1*100:.1f}%)"
                    i += 2
                else: # single
                    frac1 = breakdown_fractions[i]
                    legend_label = f"{label_base} ({frac1*100:.1f}%)"
                    i += 1
                legend_labels.append(legend_label)
            ordered_labels.extend(legend_labels[::-1])

        else:
            if data_df is not None:
                mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
            else:
                mc_handles = [h for i, h in enumerate(handles) if 'Unc.' not in labels_orig[i]]

            # breakdown_fractions[-1] = 0.916
            # breakdown_fractions[0] = 0.008
            mc_labels = [f"{label} ({frac*100:.1f}%)"
                                for label, frac in zip(labels, breakdown_fractions)]
            ordered_handles.extend(mc_handles)
            ordered_labels.extend(mc_labels[::-1]) # note the reverse order of mc_labels


    if syst is not None:
        unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
        unc_label = [l for l in labels_orig if 'Unc.' in l]
        ordered_handles.extend(unc_handle)
        ordered_labels.extend(unc_label)

    # adjust fontsize so that legend fits in the figure
    fontsize = 11.3
    ncol = 3
    if breakdown_type == "genie_sb":
        textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
        fontsize = fontsize
        ncol = 2

        # separate legend box with S / B hatches
        example_signal = Patch(facecolor="black", edgecolor='white', label='Signal')
        example_background = Patch(facecolor="black", edgecolor='white', hatch='////', linewidth=0, label='Background')
        # if textloc_x < 0.5:
        #     box_ax = ax.inset_axes([0.17, 0.65, 0.13, 0.13], transform=ax.transAxes)
        # else:
        box_ax = ax.inset_axes([0.625, 0.66, 0.13, 0.13], transform=ax.transAxes)
        box_ax.axis('off')
        mini_legend = Legend(
            box_ax,
            handles=[example_signal, example_background],
            labels=['Signal', 'Background'],
            loc='center',
            fontsize=fontsize,
            frameon=False,
            borderpad=0.7,
            handlelength=2.1,
            handleheight=0.9,
            ncol=1,
            fancybox=True,
            framealpha=1.0
        )
        box_ax.add_artist(mini_legend)

        ax.legend(
            ordered_handles,
            ordered_labels,
            loc='upper left',
            # loc='upper center',
            fontsize=fontsize,
            frameon=False,
            ncol=ncol,
            bbox_to_anchor=(0.05, 0.9, 0.8, 0.1),
            mode='expand'
        )

    else:
        ax.legend(
            ordered_handles,
            ordered_labels,
            loc='upper left',
            # loc='upper center',
            fontsize=fontsize,
            frameon=False,
            ncol=ncol,
        )

    # ax_r.legend(fontsize=9, ncol=2)

    # ===============================

    # ==== plot additions ====

    # y-axis limit
    ax.set_ylim(0., ax_ylim_ratio* np.max(total_mc))

    # vertical lines
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v[0], ymin=0, ymax=ymax*0.75, color='red', linestyle='--')
            # Plot arrow if v[1] is specified (0: left, 1: right)
            if len(v) > 1:
                direction = v[1]
                arrow_params = {
                    'y': ymax * 0.4,
                    'dx': 0.18 * (ax.get_xlim()[1] - ax.get_xlim()[0]),  # adjustable length
                    'width': 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0]),  # adjustable width
                    'color': 'red',
                    'head_width': 0.04 * (ax.get_ylim()[1] - ax.get_ylim()[0]),  # adjustable head width
                    'head_length': 0.03 * (ax.get_xlim()[1] - ax.get_xlim()[0]),  # adjustable head length
                    'length_includes_head': True
                }
                if direction == 0:
                    # Left arrow
                    ax.arrow(v[0], arrow_params['y'], -arrow_params['dx'], 0, 
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True)
                elif direction == 1:
                    # Right arrow
                    ax.arrow(v[0], arrow_params['y'], arrow_params['dx'], 0, 
                             width=arrow_params['width'],
                             color=arrow_params['color'],
                             head_width=arrow_params['head_width'],
                             head_length=arrow_params['head_length'],
                             length_includes_head=arrow_params['length_includes_head'],
                             clip_on=True)

    # textboxes
    textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
    textloc_y = textloc[1]

    if textchi2 and syst is not None:
        add_chi2_text(chi2_val, p_val, len(var_config.bins)-1, textloc_x, textloc_y+0.08, textloc_ha)

    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if breakdown_type != "pdg":
        add_genie_version_text(textloc_x, textloc_y-0.06, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    # ===============================

    # == save figure ==
    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()

    return {"cuts": cuts, 
            "total_mc": total_mc, 
            "total_mc_bkgd": total_mc_bkgd,
            "total_data": total_data}



def plot_efficiency(df_dict={}, 
                    stage_labels=[],
                    var_config=None, 
                    textloc=[0.05, 0.55],
                    approval="internal", 
                    legend=True,
                    plot=True,
                    save_fig=False, 
                    save_name=None):

    var_name = var_config.var_nu_col
    bins = var_config.bins
    bin_centers = var_config.bin_centers

    keys = list(df_dict.keys())
    colors = ["C" + str(i) for i in range(len(keys))]

    all_signal_df = df_dict[keys[0]][IsNuInFV_NumuCC_1p0pi(df_dict[keys[0]])]
    n_tot_integ = len(all_signal_df)
    var, wgts = get_clipped_evts(all_signal_df, var_name, bins)
    n_tot, _ = np.histogram(var, weights=wgts, bins=bins)

    fig, ax = plt.subplots()
    ax_eff = ax.twinx()

    eff_list = []
    eff_err_list = []
    for kidx, key in enumerate(keys):
        this_df = df_dict[key]
        this_signal_df = this_df[IsNuInFV_NumuCC_1p0pi(this_df)]
        var, wgts = get_clipped_evts(this_signal_df, var_name, bins)
        n, _ = np.histogram(var, weights=wgts, bins=bins)

        this_eff_integ = len(this_signal_df) / n_tot_integ
        this_label = stage_labels[kidx] + " ({:.2f}%)".format(this_eff_integ*100)
        n, bins, _ = ax.hist(var, weights=wgts, bins=bins, histtype="step", label=this_label, alpha=0.6)

        if key == keys[-1]:
            print("final purity: {:.2f}%".format(100 * len(this_signal_df) / len(this_df)))

        this_eff = n / n_tot
        stat_scale_factor = this_df.pot_weight.unique()[0]
        this_eff_err = get_eff_err(n/stat_scale_factor, n_tot/stat_scale_factor)
        ax_eff.errorbar(bin_centers, this_eff, yerr=this_eff_err, fmt="o-", color=colors[kidx], label=this_label)
        eff_list.append(this_eff)
        eff_err_list.append(this_eff_err)

    ax.set_xlabel(var_config.var_labels[0])
    ax.set_ylabel("Events")
    ax_eff.set_ylabel("Efficiency")
    ax_eff.set_ylim(0, 1.05)
    plt.xlim(bins[0], bins[-1])

    if legend:
        plt.legend(bbox_to_anchor=(1.15, 1.0))

    # ==== plot additions ====
    textloc_x, textloc_ha = get_textloc_x(n_tot, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches="tight", dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()

    return {"eff_list": eff_list, "eff_err_list": eff_err_list}


def plot_univ_hists(
                univ_events, 
                cv_events,
                syst_name, 
                var_config, 
                approval="internal",
                textloc=[0.05, 0.55],
                plot=True,
                ax_titles = ["", "", ""],
                save_fig=False, 
                save_name=None): 

    assert univ_events.shape[1] == len(cv_events) 
    n_univ = univ_events.shape[0]

    if (n_univ > 10):
        colors = ["#FDE725FF", "#1F968BFF", "#440154FF"] # viridis colors

        sorted_univs = np.sort(univ_events, axis=0)

        n_68 = int(0.68 * n_univ)
        start_68 = (n_univ - n_68) // 2
        end_68 = start_68 + n_68

        n_95 = int(0.95 * n_univ)
        start_95 = (n_univ - n_95) // 2
        end_95 = start_95 + n_95

        # Define bins & groupings of universes: [(range, color, label, skip_68)]
        segs = [
            (range(start_68, end_68), colors[0], "Universe (68%)", False),
            (range(start_95, end_95), colors[1], "Universe (95%)", True),
            ( (i for i in range(n_univ) if i not in range(start_95, end_95)), colors[2], "Universe (100%)", False)
        ]

        plotted = set()  # tracks which labels were plotted already
        for r, color, label, skip_68 in segs:
            for i in r:
                if skip_68 and i in range(start_68, end_68): 
                    continue
                show_label = label if label not in plotted else None
                plt.hist(var_config.bin_centers, bins=var_config.bins, weights=sorted_univs[i],
                        histtype="step", color=color, alpha=0.7, label=show_label)
                plotted.add(label)

    # if too few universes, just plot all of them in gray
    else:
        for i in range(n_univ):
            show_label = "Universe" if i == 0 else None
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", color="gray", label=show_label)

    # plot CV last so that it is on top of the universes
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=cv_events, histtype="step", color="k", label="Central Value")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[1])
    if ax_titles[1] != "":
        plt.ylabel(ax_titles[1])
    else:
        plt.ylabel("Events / Bin")
    plt.title(ax_titles[2])

    plt.legend(frameon=False)

    # ==== plot additions ====
    textloc_x, textloc_ha = get_textloc_x(cv_events, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()


    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


def plot_unfolded_result(unfold, 
                         measured, 
                         models,
                         var_config, 
                         chi2_list=[],
                         xsec_unit=0,
                         textloc=[0.05, 0.55],
                         approval="internal",
                         plot_labels=["", "", ""],
                         plot=True,
                         save_fig=False, 
                         save_name=None,
                         data=False,
                         closure_test=False):

    bins = var_config.bins
    bin_centers = var_config.bin_centers
    bin_widths = np.diff(bins)
    if len(var_config.bins) == 2:
        bin_widths = np.array([1.0])

    # unfolded result
    Unfolded = unfold['unfold']
    UnfoldedCov = unfold["UnfoldCov"]
    Unfolded_perwidth = Unfolded / bin_widths

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov']
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)

    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)
    UnfoldCov_syst_frac = fraccov_from_cov(UnfoldCov_syst, Unfolded)

    # --- decompose into norm and shape components
    # the first item in models dict is the nominal input model
    norm_model = list(models.keys())[0]
    SystUnfoldCov_norm, SystUnfoldCov_shape = Matrix_Decomp(models[norm_model], UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))


    # --- plot
    fig, ax = plt.subplots(figsize=(8.5, 7))
    # set err to 0 for closure test
    if closure_test:
        dummy_err = np.zeros_like(Unfolded_perwidth)
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=dummy_err, fmt='o', color='black')

    else:
        # plot shape uncertainty as error bars
        Unfold_uncert_stat_perwidth = Unfold_uncert_stat / bin_widths
        Unfold_uncert_shape_perwidth = Unfold_uncert_shape / bin_widths
        # Unfold_uncert_stat_shape_perwidth = Unfold_uncert_stat_perwidth + Unfold_uncert_shape_perwidth
        Unfold_uncert_stat_shape_perwidth = Unfold_uncert_shape_perwidth
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=Unfold_uncert_stat_shape_perwidth, fmt='o', color='black', capsize=3)

        # plot syst norm component as histogram at the bottom
        Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths
        if len(var_config.bins) != 2:
            norm_handle = plt.bar(bin_centers, Unfold_uncert_norm_perwidth, width=bin_widths, label='Syst. error (norm)', alpha=0.5, color='gray')

    if data: # get stat uncertainty for data
        # data_eylow, data_eyhigh = return_data_stat_err(measured/xsec_unit)
        # Data_frac_unc = (data_eyhigh - data_eylow) / (2 * measured/xsec_unit)
        Data_frac_unc = (1/np.sqrt(measured / xsec_unit))
        Data_stat = Unfolded_perwidth * Data_frac_unc
        # Data_stat = Unfolded * Data_frac_unc
        # factor = bin_widths.sum()/len(bin_widths)
        # Data_stat = Data_stat / factor
        Data_stat_frac_unc_smeared = (unfold['AddSmear'] @ Data_frac_unc)
        # Data_frac_unc_smeared = (1/np.sqrt(unfold['AddSmear'] @ measured / xsec_unit))
        Data_stat_smeared = Unfolded_perwidth * Data_stat_frac_unc_smeared

        Data_stat_frac_cov = np.diag(Data_frac_unc**2)
        Data_stat_cov = cov_from_fraccov(Data_stat_frac_cov, Unfolded_perwidth)
        Data_stat_frac_cov_smeared = np.diag((unfold['AddSmear'] @ Data_frac_unc)**2)
        Data_stat_cov_smeared = cov_from_fraccov(Data_stat_frac_cov_smeared, Unfolded_perwidth)

        if len(var_config.bins) == 2:
            tot_err = np.sqrt(Data_stat**2 + Unfold_uncert_norm_perwidth**2)
        else:
            tot_err = np.sqrt(Data_stat**2 + Unfold_uncert_stat_shape_perwidth**2)
        Data_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=tot_err, fmt='o', color='black', capsize=3)
        handles = [bar_handle, Data_handle]
        labels = ['SBND Development Data', 'Measured Signal']
        UnfoldCov_syst = cov_from_fraccov(UnfoldCov_syst_frac, Unfolded_perwidth)

        UnfoldCov_syst = UnfoldCov_syst + Data_stat_cov
        UnfoldCov_syst_smeared = UnfoldCov_syst + Data_stat_cov_smeared

    # divide measured & model by bin width
    measured_perwidth = measured / bin_widths
    if data == False:
        reco_handle, = plt.step(bins, np.append(measured_perwidth, measured_perwidth[-1]), where='post', label='Measured Signal (Input)')

    # --- get chi2 values for each model to compare
    if len(chi2_list) == 0:
        chi2_vals = []
        p_values = []
    else:
        chi2_vals = chi2_list
    model_handles = []
    model_labels = []
    for midx, mkey in enumerate(models.keys()):
        model_smeared = unfold['AddSmear'] @ models[mkey]
        # if "SBN" in mkey:
        model_smeared_perwidth = model_smeared / bin_widths

        # else:
        #     model_smeared_perwidth = model_smeared 

        if len(chi2_list) == 0:
            # remove bins with <= 0 events
            # Fix chi2 mask logic: mask just once, store, reuse, improve clarity
            mask = (Unfolded_perwidth > 0) & (model_smeared_perwidth > 0)
            Unfolded_perwidth_safe = Unfolded_perwidth[mask]
            model_smeared_perwidth_safe = model_smeared_perwidth[mask]
            UnfoldCov_syst_safe = UnfoldCov_syst[np.ix_(mask, mask)]
            chi2_val, p_val = get_chi2(Unfolded_perwidth_safe, model_smeared_perwidth_safe, UnfoldCov_syst_safe)
            # chi2_val, p_val = get_chi2(Unfolded_perwidth, model_smeared_perwidth, UnfoldCov_syst_smeared)
            chi2_vals.append(chi2_val)
            p_values.append(p_val)

        print("Unfolded perwidth: ", Unfolded_perwidth)
        print("Model smeared perwidth: ", model_smeared_perwidth)

        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post')
        model_handles.append(model_handle)
        # model_labels.append(f'$A_c \\otimes$ {mkey} ($\chi^2$ = {chi2_vals[midx]:.2f}/{len(bins)-1}), p-value = {p_values[midx]:.3f}')
        # model_labels.append(f'$A_c \\otimes$ {mkey} ($\chi^2$ = {chi2_vals[midx]:.2f}/{len(bins)-1})')
        model_labels.append(f'$A_c \\otimes$ {mkey}')

    # legend
    if closure_test:
        handles = [bar_handle, reco_handle] + model_handles
        labels = ['Unfolded Asimov Data', 'Measured Signal'] + model_labels
    elif data:
        if len(var_config.bins) == 2:
            handles = [bar_handle] + model_handles
            labels = ['Data (Syst. Unc. + Stat. Unc.)'] + model_labels
        else:
            handles = [bar_handle, norm_handle] + model_handles
            labels = ['Data (Shape Syst. Unc. + Stat. Unc.)', 'Norm. Syst. Unc.'] + model_labels
    else:
        handles = [bar_handle, norm_handle, reco_handle] + model_handles
        labels = ['SBND Development Data', 'Norm. Syst. Unc.', 'Measured Signal'] + model_labels
    plt.legend(handles, labels, 
               loc='upper left', fontsize=12, frameon=False, ncol=1, bbox_to_anchor=(0.02, 0.98))

    plt.xlabel(var_config.var_labels[0], fontsize=20)
    plt.ylabel(var_config.xsec_label, fontsize=20)
    plt.title(plot_labels[2])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.7)

    # ==== plot additions
    textloc_x, textloc_ha = get_textloc_x(Unfolded_perwidth, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    add_genie_version_text(textloc_x, textloc_y-0.1, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


def variation_hists(evtdfs=None, var_name=None, breakdown_type=None,
                    nevts_list=None,
                    datadf=None,
                    bins=None,
                    var_colors=None, var_labels=None,
                    plot_labels=["", "", ""],
                    vline = None,
                    textloc=[0.05, 0.55],
                    approval="internal",
                    plot=True,
                    save_fig=False, save_name=None): 

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    if evtdfs is not None:
        n_vars = len(evtdfs)
    elif nevts_list is not None:
        n_vars = len(nevts_list)
    else:
        raise ValueError("Either evtdfs or nevts_list must be provided")

    # get distribution from dfs
    if evtdfs is not None:
        vardfs, wgtdfs = [], []
        nevts_list = []
        mc_stat_err_list = []
        for df in evtdfs:
            vardf, wgtdf = get_clipped_evts(df, var_name, bins)
            vardfs.append(vardf)
            wgtdfs.append(wgtdf)

            nevts, _ = np.histogram(vardf, bins=bins, weights=wgtdf)
            total_mc_err2, _ = np.histogram(vardf, bins=bins, weights=wgtdf**2)
            mc_stat_err = np.sqrt(total_mc_err2)
            nevts_list.append(nevts)
            mc_stat_err_list.append(mc_stat_err)

        if breakdown_type is not None:
            # plot breakdown for the first variation (CV)
            if breakdown_type == "topology":
                labels = topology_labels[::-1]
                colors = topology_colors[::-1]
                cuts = get_topo_category(evtdfs[0], ret_cuts=True)

            elif breakdown_type == "genie":
                labels = genie_mode_labels
                colors = genie_mode_colors
                cuts = get_genie_category(evtdfs[0], ret_cuts=True)


    if datadf is not None:
        vardf_data, wgtdf_data = get_clipped_evts(datadf, var_name, bins)
        total_data, _ = np.histogram(vardf_data, bins=bins, weights=wgtdf_data)
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

    # ===== plot =====
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7), 
                            sharex=True, gridspec_kw={'height_ratios': [4, 1]})

    fig.subplots_adjust(hspace=0.05)
    ax = axs[0]
    ax_r = axs[1]

    for sidx in range(n_vars):
        if sidx == 0:
            if breakdown_type is not None:
                vardf, wgtdf = get_clipped_evts(evtdfs[0], var_name, bins)
                vardf_categ = [vardf[i] for i in cuts]
                weights_categ = [list(evtdfs[0].loc[cuts[i], 'pot_weight']) for i in range(len(cuts))]
                mc_stack, _, _ = ax.hist(vardf_categ,
                                        weights=weights_categ,
                                        bins=bins,
                                        stacked=True,
                                        color=colors,
                                        label=labels,
                                        linewidth=0,
                                        edgecolor='none',
                                        histtype='stepfilled')
                continue

        ax.hist(bin_centers,
                weights=nevts_list[sidx],
                bins=bins, 
                histtype="step" , 
                color=var_colors[sidx], 
                label=var_labels[sidx])

    if datadf is not None:
        ax.errorbar(bin_centers, 
                    total_data, 
                    yerr=np.vstack((data_eylow, data_eyhigh)),
                    color='black', 
                    fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data')


    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])
    ax.set_xlim(bins[0], bins[-1])

    # ==== var/CV ratio panel
    for sidx in range(n_vars):
        if sidx == 0:
            continue
        # Avoid division by zero: ignore bins where denominator is 0
        ratio = np.full_like(nevts_list[0], np.nan, dtype=float)
        nonzero_mask = nevts_list[0] != 0
        ratio[nonzero_mask] = nevts_list[sidx][nonzero_mask] / nevts_list[0][nonzero_mask]
        # this_err = np.sqrt(
        #     (mc_stat_err_list[0] / nevts_list[0])**2 + 
        #     (mc_stat_err_list[sidx] / nevts_list[sidx])**2
        # )
        # Only plot nonzero, non-nan elements in the ratio
        valid_mask = (~np.isnan(ratio)) & (ratio != 0)
        ax_r.hist(bin_centers[valid_mask], bins=bins, weights=ratio[valid_mask], linewidth=1, histtype="step", color=var_colors[sidx])

    ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
    ax_r.set_xlim(bins[0], bins[-1])
    ax_r.set_ylim(0.9, 1.1)
    ax_r.set_xlabel(plot_labels[0])
    ax_r.set_ylabel("Variation / CV")
    ax_r.grid(True)
    ax_r.minorticks_on()
    ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)

    # ==== plot additions
    # vertical lines on main panel
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    # approval rank
    textloc_x, textloc_ha = get_textloc_x(nevts_list[0], bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)
    ax.legend(loc="best")

    # == save figure ==
    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()

    return nevts_list

def signal_cut(df, detector=DETECTOR):
    print("DETECTOR: ", detector)
    signal_cut =  (df.mc.nmu_220MeVc == 1) & (df.mc.np_300MeVc == 1) & (df.mc.npi_70MeVc == 0) & (df.mc.npi0 == 0) &\
                    (np.sqrt(df.mc.mu.genp.x**2 + df.mc.mu.genp.y**2 + df.mc.mu.genp.z**2) < 1) &\
                    (np.sqrt(df.mc.p.genp.x**2 + df.mc.p.genp.y**2 + df.mc.p.genp.z**2) < 1) &\
                        InFV(df.mc.mu.start, det=detector) & InFV(df.mc.p.start, det=detector) 
    return df[signal_cut]

def signal_hists(evtdf=None,  # df with selected & reco'ed events
                 nudf=None,   # df with all MC truth
                 var_config=None,
                 return_data=False,
                 plot=True,
                 textloc=[0.05, 0.55],
                 approval="internal",
                 mode="reco",
                 save_fig=False, 
                 save_name=None):
    """
    generate / selected / reco'ed signal events

    topo_categ == 1 : signal
    topology_list has "1" as the first item, corresponding to the signal topology

    items with "sel" tags are selected events that are signal in truth 
    item with "allsel" tags are all selected events, signal + background
    """

    bins = var_config.bins
    bin_centers = var_config.bin_centers
    reco_col = var_config.var_evt_reco_col
    truth_col = var_config.var_evt_truth_col
    nu_col = var_config.var_nu_col

    # ===== all selected events =====
    # reco'ed
    var_allsel_reco, wgt_allsel_reco = get_clipped_evts(evtdf, reco_col, bins)
    var_allsel_truth, wgt_allsel_truth = get_clipped_evts(evtdf, truth_col, bins)

    # ===== true signal events =====
    evtdf_signal = evtdf[evtdf.topo_categ == 1]
    # selected, reco'ed 
    var_sel_reco, wgt_sel_reco = get_clipped_evts(evtdf_signal, reco_col, bins)
    # selected, truth (for response matrix)
    var_sel_truth, wgt_sel_truth = get_clipped_evts(evtdf_signal, truth_col, bins)

    nevts_allsel_truth, _ = np.histogram(var_allsel_truth, weights=wgt_allsel_truth, bins=bins)
    nevts_allsel_reco, _  = np.histogram(var_allsel_reco,  weights=wgt_allsel_reco,  bins=bins)
    nevts_sel_truth, _    = np.histogram(var_sel_truth,    weights=wgt_sel_truth,    bins=bins)
    nevts_sel_reco, _     = np.histogram(var_sel_reco,     weights=wgt_sel_reco,     bins=bins)

    if nudf is not None:
        # all MC, truth (for efficiency vector)
        nudf_signal = nudf[nudf.topo_categ == 1]
        var_allmc, wgt_allmc = get_clipped_evts(nudf_signal, nu_col, bins)
        nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=bins)
        if mode == "unfold":
            # don't consider containment

            nudf_signal = signal_cut(nudf)
            var_allmc, wgt_allmc = get_clipped_evts(nudf_signal, nu_col, bins)
            nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=bins)
    else:
        var_allmc = None
        wgt_allmc = None
        nevts_allmc = None

    if plot:
        plt.hist(bin_centers, weights=nevts_allsel_truth, bins=bins, histtype="step", label="Selected Events, True", color="C0")
        plt.hist(bin_centers, weights=nevts_allsel_reco,  bins=bins, histtype="step", label="Selected Events, Reco", color="C0", linestyle="--")
        plt.hist(bin_centers, weights=nevts_sel_truth,    bins=bins, histtype="step", label="Selected Signal, True", color="C1")
        plt.hist(bin_centers, weights=nevts_sel_reco,     bins=bins, histtype="step", label="Selected Signal, Reco", color="C1", linestyle="--")

        if nevts_allmc is not None:
            nevts_allmc, _, _ = plt.hist(var_allmc,     weights=wgt_allmc,     bins=bins, histtype="step", label="All Signal in MC", color="black")
        else:
            nevts_allmc = None

        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.xlim(bins[0], bins[-1])
        plt.legend()

        # ==== plot additions ====
        textloc_x, textloc_ha = get_textloc_x(nevts_allsel_truth, bins, textloc)
        textloc_y = textloc[1]
        add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

        if var_config.var_save_name == "integrated":
            format_singlebin_plot()

        if save_fig:
            plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

        if plot == True:
            plt.show()
        else:
            plt.close()

    if return_data:
        return {
            "var_allmc": var_allmc,
            "wgt_allmc": wgt_allmc,
            "nevts_allmc": nevts_allmc,

            "var_sel_truth": var_sel_truth,
            "wgt_sel_truth": wgt_sel_truth,
            "nevts_sel_truth": nevts_sel_truth,

            "var_sel_reco": var_sel_reco,
            "wgt_sel_reco": wgt_sel_reco,
            "nevts_sel_reco": nevts_sel_reco,

            "var_allsel_truth": var_allsel_truth,
            "wgt_allsel_truth": wgt_allsel_truth,
            "nevts_allsel_truth": nevts_allsel_truth,

            "var_allsel_reco": var_allsel_reco,
            "wgt_allsel_reco": wgt_allsel_reco,
            "nevts_allsel_reco": nevts_allsel_reco,
        }


# ==== fractional uncertainty plot ====
def plot_frac_unc(frac_unc_list, 
                  var_config, 
                  plot_labels=["", "", ""],
                  legends = None,
                  textloc=[0.05, 0.55],
                  approval="internal",
                  plot=True,
                  save_fig=False, 
                  save_name=None):

    for fidx, frac_unc in enumerate(frac_unc_list):
        color = "C{}".format(fidx)
        if len(frac_unc_list) == 1:
            color = "black"
        plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_unc, histtype="step", color=color)

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel("Fractional Uncertainty")
    plt.title(plot_labels[2])
    plt.grid(True)

    if legends is not None:
        plt.legend(legends)

    textloc_x, textloc_ha = get_textloc_x(frac_unc, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if var_config.var_save_name == "integrated":
        format_singlebin_plot()

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()


# ==== 2D plots ====

def get_text_color(value):
    rgba = cmap(norm(value))
    # Compute luminance (perceived brightness)
    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if luminance > 0.5 else "white"


def bin_range_labels(edges):
    return [f"{edges[i]:.2f}–{edges[i+1]:.2f}" for i in range(len(edges)-1)]


def plot_heatmap(matrix, 
                 bins,
                 plot_labels=["", "", ""],
                 approval="internal",
                 verbose=False,
                 plot=True,
                 cmap="bwr",
                 save_fig=False, 
                 save_name=None):

    nbins = len(bins)
    assert nbins-1 == matrix.shape[0] == matrix.shape[1]
    unif_bin = np.linspace(0., float(nbins - 1), nbins)
    extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

    x_edges, y_edges = np.array(bins), np.array(bins)
    x_tick_positions, y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2, (unif_bin[:-1] + unif_bin[1:]) / 2
    x_labels, y_labels = bin_range_labels(x_edges), bin_range_labels(y_edges)

    fig, ax = plt.subplots(figsize=(10, 10))
    if cmap == "bwr":
        plt.imshow(matrix, extent=extent, origin="lower", vmin=-1, vmax=1, cmap=cmap)
    else:
        plt.imshow(matrix, extent=extent, origin="lower", cmap=cmap)

    # Find the power-of-10 exponent from one of the (non-NaN) values
    exponent = 0
    flat_matrix = matrix[~np.isnan(matrix)]
    if flat_matrix.size > 0 and np.any(flat_matrix != 0):
        example_value = flat_matrix[0]
        exponent = np.floor(np.log10(abs(example_value)))
        # if exponent is infinite (matrix contains only zeros), set to 0
        if np.isinf(exponent):
            exponent = 0
        exponent = int(exponent)

        # # Round to nearest multiple of 3 for best tick label presentation
        # exponent = 3 * int(np.floor(exponent / 3))

        formatter = mpl.ticker.FuncFormatter(lambda x, _: f"{x/10**exponent:.2f}")
        cbar = plt.colorbar(shrink=0.7)
        cbar.set_label(f"{plot_labels[2]} [10$^{{{exponent}}}$]", fontsize=16)
        cbar.ax.yaxis.set_major_formatter(formatter)
    else:
        plt.colorbar(shrink=0.7, label=plot_labels[2])

    for i in range(nbins-1):      # rows (y)
        for j in range(nbins-1):  # columns (x)
            value = matrix[i, j]
            if not np.isnan(value):  # skip NaNs
                significand = value / 10**exponent
                plt.text(
                    j + 0.5, i + 0.5,
                    f"{significand:.2f}",
                    ha="center", va="center",
                    color=get_text_color(value),
                    fontsize=10
                )

    plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
    plt.yticks(y_tick_positions, y_labels)
    plt.xlabel(plot_labels[0], fontsize=20)
    plt.ylabel(plot_labels[1], fontsize=20)
    # plt.title(plot_labels[2], fontsize=20)

    if verbose:
        n_diag = np.sum(np.diag(matrix))
        diagonal_ratio = n_diag / np.sum(matrix)
        print(f"Diagonal ratio: {diagonal_ratio:.2f}")
        print(f"True ratio: {np.diag(matrix) / np.sum(matrix, axis=0)}")

        # print

    # ===== plot additions =====
    add_approval_text(approval, 0.95, 1.05, "right")

    if save_fig:
        plt.savefig(save_name+fig_ext, bbox_inches='tight', dpi=dpi)

    if plot == True:
        plt.show()
    else:
        plt.close()




####
# Exposure Accounting

def get_integrated_flux(fluxfile, plot=False):
    # flux file, units: /m^2/10^6 POT 
    # 50 MeV bins
    flux = uproot.open(fluxfile)
    numu_flux = flux["flux_sbnd_numu"].to_numpy()
    bin_edges = numu_flux[1]
    flux_vals = numu_flux[0]

    if plot:
        fig, ax = plt.subplots()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2, color="C0")
        plt.xlim(0, 3)
        plt.xlabel("Neutrino Energy [GeV]")
        plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
        plt.title("SBND $\\nu_\\mu$ Flux")
        plt.savefig("sbnd-flux.pdf", bbox_inches='tight')

    integrated_flux = flux_vals.sum() / (1e4  * 1e6) # to cm2 # to POT
    print("Integrated flux: %.3e" % integrated_flux)
    return integrated_flux


def get_active_volume(detector="SBND"):
    if detector == "SBND":
        V_SBND = 380 * 380 * 440 # cm3, the active volume of the detector 

    elif detector == "SBND_nohighyz":
        V_SBND = 380 * 380 * 440 - 380* (190 - 100) *(450-250) 

    elif detector == "SBND_face":
        V_SBND = 380 * 380 * 50

    elif detector == "SBND_face_yzcut":
        V_SBND = 380 * 380 * 50 - 380* (190 - 100) * 50 

    elif detector == "SBND_end":
        V_SBND = 380 * 380 * 50

    return V_SBND


def get_xsec_unit(tot_pot, 
                  fluxfile="/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root", 
                  detector="SBND",
                  volume=None):

    tot_flux = get_integrated_flux(fluxfile, plot=False)
    tot_flux *= tot_pot
    print("integrated flux: ", tot_flux)

    V_SBND = get_active_volume(detector)
    if volume is not None:
        print("using custom volume: ", volume)
        V_SBND = volume

    NTARGETS = RHO * V_SBND * (N_A / M_AR) #/ 40 # divide by 40 to make this per-argon nucleus
    print("# of targets: ", NTARGETS)

    xsec_unit = 1 / (tot_flux * NTARGETS)
    # # TODO: fix scalar overflow error in python v3.10+
    # if xsec_unit == 0:
    #     print("XSEC_UNIT is 0, setting to 1e-38")
    #     xsec_unit = 1e-38
    print("xsec unit: ", xsec_unit)
    return xsec_unit

    