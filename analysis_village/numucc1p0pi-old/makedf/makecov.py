import pandas as pd
import numpy as np
import sys
# local imports
sys.path.append('/exp/sbnd/app/users/munjung/xsec/wienersvd/cafpyana')
from analysis_village.numucc1p0pi.selection_definitions import *
from analysis_village.numucc1p0pi.makedf.makedf import *
from makedf.geniesyst import regen_systematics_sbnd_multisigma, regen_systematics_sbnd_morph
from pyanalib.variable_calculator import get_cc1p0pi_tki
from analysis_village.numucc1p0pi.variable_configs import VariableConfig


def make_geniewgts(df, wgt_types=["bnb","genie"]):

    var_configs = [VariableConfig.muon_momentum(),
                   VariableConfig.muon_direction(),
                   VariableConfig.proton_momentum(),
                   VariableConfig.proton_direction(),
                   VariableConfig.tki_del_p(),
                   VariableConfig.tki_del_Tp(),
                   VariableConfig.tki_del_Tp_x(),
                   VariableConfig.tki_del_Tp_y(),
                   VariableConfig.tki_del_alpha(),
                   VariableConfig.tki_del_phi()]

    df.loc[:,'nuint_categ'] = get_int_category(df)
    df.loc[:,'genie_categ'] = get_genie_category(df)

    # add tki deltas to df
    mc_mudf = df.mu
    mc_pdf = df.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)

    # make TKI if it doesn't exist
    if "mc_del_alpha" not in df.columns:
        df = multicol_add(df, tki_mc["del_alpha"].rename("mc_del_alpha"))
        df = multicol_add(df, tki_mc["del_phi"].rename("mc_del_phi"))
        df = multicol_add(df, tki_mc["del_Tp"].rename("mc_del_Tp"))
        df = multicol_add(df, tki_mc["del_p"].rename("mc_del_p"))
        df = multicol_add(df, tki_mc["del_Tp_x"].rename("mc_del_Tp_x"))
        df = multicol_add(df, tki_mc["del_Tp_y"].rename("mc_del_Tp_y"))


    dict_univ_counts = {}
    for var_config in var_configs:
        for topo_mode in mode_list:
            this_df = df[df.nuint_categ == topo_mode]

            # if df list is empty, all entries are 0
            skip = False
            if len(this_df) == 0:
                skip = True

            ## CV
            n_cv, bins = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins)
            cv_col_name = ("nuint_categ_"+str(topo_mode), 
                           "_".join([str(s) for s in var_config.var_nu_col if s]),
                           "cv",
                           "")
            dict_univ_counts[cv_col_name] = n_cv

            ## GENIE
            if "genie" in wgt_types:

                # multisim
                syst_name = "GENIE"
                n_syst_univ = 100

                n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
                for i in range(n_syst_univ):
                    weights = this_df[syst_name]["univ_{}".format(i)]
                    weights[np.isnan(weights)] = 1
                    if not skip:
                        n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
                    col_name = ("nuint_categ_"+str(topo_mode), 
                                "_".join([str(s) for s in var_config.var_nu_col if s]),
                                syst_name+"_multisim", 
                                "univ_{}".format(i))
                    dict_univ_counts[col_name] = n_univ[i, :]

                # multisigma
                n_syst_univ = 6
                for syst_name in regen_systematics_sbnd_multisigma:
                    n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
                    for i in range(n_syst_univ):
                        if (i < 3):
                            pm = "ms"
                        else:
                            pm = "ps"
                        sig = i%3 + 1
                        weights = this_df[syst_name][pm+str(sig)]
                        weights[np.isnan(weights)] = 1
                        if not skip:
                            n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
                        col_name = ("nuint_categ_"+str(topo_mode), 
                                    "_".join([str(s) for s in var_config.var_nu_col if s]),
                                    syst_name,
                                    "univ_{}".format(i))
                        dict_univ_counts[col_name] = n_univ[i, :]

                # morph
                n_syst_univ = 1
                for syst_name in regen_systematics_sbnd_morph:
                    n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
                    for i in range(n_syst_univ):
                        weights = this_df[syst_name]["morph"]
                        weights[np.isnan(weights)] = 1
                        if not skip:
                            n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
                        col_name = ("nuint_categ_"+str(topo_mode), 
                                    "_".join([str(s) for s in var_config.var_nu_col if s]),
                                    syst_name,
                                    "univ_{}".format(i))
                        dict_univ_counts[col_name] = n_univ[i, :]

            if "bnb" in wgt_types:

                ## Flux
                syst_name = "Flux"
                n_syst_univ = 500
                n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
                for i in range(n_syst_univ):
                    weights = this_df[syst_name]["univ_{}".format(i)]
                    weights[np.isnan(weights)] = 1
                    if not skip:
                        n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
                    col_name = ("nuint_categ_"+str(topo_mode), 
                                "_".join([str(s) for s in var_config.var_nu_col if s]),
                                syst_name, 
                                "univ_{}".format(i))
                    dict_univ_counts[col_name] = n_univ[i, :]


            # if "mcstat" in wgt_types:
            #     ## MCStat
            #     n_univ_mcstat = 500
            #     mc_evt_df, MCstat_univ_events = mcstat(mc_evt_df, mc_hdr_df, n_universes=n_univ_mcstat)

    df_univ_counts = pd.DataFrame(dict_univ_counts)
    if not isinstance(df_univ_counts.columns, pd.MultiIndex):
        df_univ_counts.columns = pd.MultiIndex.from_tuples(df_univ_counts.columns)

    return df_univ_counts


def make_numucc1p0pi_evtdf_final_wgts(f, mode="df", include_weights=True, multisim_nuniv=500, wgt_types=["bnb","genie"], slim=True, 
                           trkScoreCut=False, **trkArgs):
    if mode == "df":
        return make_numucc1p0pi_evtdf(f, sel_level="final", include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                                    trkScoreCut=trkScoreCut, **trkArgs)
    elif mode == "cov":
        # all weights on mcnudf
        # print("all mc nus")
        # mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)

        # print("final selection")
        mcdf = make_numucc1p0pi_evtdf(f, sel_level="final", include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                                    trkScoreCut=trkScoreCut, **trkArgs)

        mcdf_wgts = make_geniewgts(mcdf)

        hdr = loadbranches(f["recTree"], mchdrbranches).rec.hdr
        pot = hdr['pot'].sum()

        return mcdf_wgts, pot


def make_numucc1p0pi_binned_systs(filename):
    
    events = uproot.open(filename)
    try:
        ret = make_numucc1p0pi_evtdf_final_wgts(events, mode="cov")

    except Exception as e:
        print(f"Error in make_numucc1p0pi_binned_systs for file {filename}: {e}")
        return None

    return ret
