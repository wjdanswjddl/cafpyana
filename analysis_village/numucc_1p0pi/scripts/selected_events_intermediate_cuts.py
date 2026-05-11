#!/usr/bin/env python3
"""Overlay plots for intermediate cut variables using loose MC/data dfs.

Loads ``get_ana_dfs(option=\"event_selection\")`` (same sample bookkeeping as the
event-selection notebook: MC ``-sel_all-wgts``-style CAFs, ``_Fixed_all`` data).

Run uncertainties first, e.g. ``get_systematics_mcstat_flux_g4.py --mc-df-stage sel_all
--var-set intermediate``, or ``get_systematics_multisim.py`` with MC loading matched to the
same ``sel_all`` sample (outputs include ``*_syst_dict.npz`` files consumed below).

Then pass ``--syst-results-dir`` to this script pointing at the **syst disk root**
(``MCstat/``, ``Flux/``, ``G4/``, ``Cosmics/``, ``GENIE/cov_mat_dict.pkl`` — see
``syst_disk_layout``). Variables absent from the GENIE covariance pickle use a zero matrix.
"""

from __future__ import annotations

import argparse
import pickle
from datetime import datetime
from functools import partial
from os import makedirs, path

import numpy as np
import pandas as pd

import sys

sys.path.append("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")

from analysis_village.numucc_1p0pi.files_config import (
    get_ana_dfs,
    save_fig_base_dir,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import syst_disk_paths
from analysis_village.numucc_1p0pi.utils import (
    get_pot_str,
    overlay_hists,
)
import matplotlib.pyplot as plt

plt.style.use("presentation.mplstyle")


def _zero_cov(var_config):
    n = len(var_config.bin_centers)
    return np.zeros((n, n))


def load_syst_cov_bundle(var_config, syst_root: str | None, genie_pkl: str | None):
    def _cov_frac_from_npz(npz_obj, var_sn, inner_key):
        z = dict(npz_obj)
        if var_sn not in z:
            return _zero_cov(var_config)
        ret = z[var_sn].item()[inner_key]
        return ret["cov_frac"]

    vn = var_config.var_save_name
    if syst_root is None:
        raise SystemExit("Pass --syst-results-dir (syst disk root; see syst_disk_layout).")

    paths = syst_disk_paths(syst_root)
    for role in ("mcstat", "flux", "g4", "cosmics"):
        pth = paths[role]
        if not path.isfile(pth):
            raise SystemExit("Missing required %s uncertainty file: %s" % (role.upper(), pth))

    mcstat_npz = np.load(paths["mcstat"], allow_pickle=True)
    g4_npz = np.load(paths["g4"], allow_pickle=True)
    flux_npz = np.load(paths["flux"], allow_pickle=True)
    cosmics_npz = np.load(paths["cosmics"], allow_pickle=True)

    mcstat_syst = _cov_frac_from_npz(mcstat_npz, vn, "MCstat")
    g4_syst = _cov_frac_from_npz(g4_npz, vn, "G4")
    flux_syst = _cov_frac_from_npz(flux_npz, vn, "flux")
    cosmics_syst = _cov_frac_from_npz(cosmics_npz, vn, "Cosmics")

    genie_path = genie_pkl or paths["genie"]
    if not path.isfile(genie_path):
        raise SystemExit("Missing GENIE covariance pickle: %s" % genie_path)
    genie_blob = pickle.load(open(genie_path, "rb"))
    try:
        genie_syst = genie_blob[var_config.var_save_name]["genie"]
    except KeyError:
        genie_syst = _zero_cov(var_config)

    pot_frac_unc = 0.02
    ntargets_frac_unc = 0.01

    frac_uncert_total = np.zeros(len(var_config.bin_centers))
    cov_total = np.zeros((len(var_config.bin_centers), len(var_config.bin_centers)))
    systs = [mcstat_syst, genie_syst, flux_syst, g4_syst, cosmics_syst]
    for syst in systs:
        cov_total += syst
        frac_uncert_total += np.sqrt(np.diag(syst)) ** 2

    for syst in (pot_frac_unc, ntargets_frac_unc):
        cov_total += np.diag(syst * np.ones(len(var_config.bin_centers)) ** 2)
        frac_uncert_total += (syst * np.ones(len(var_config.bin_centers))) ** 2

    frac_uncert_total = np.sqrt(frac_uncert_total)
    return cov_total, frac_uncert_total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunk_idx", type=int, default=0)
    ap.add_argument("--n_time_splits", type=int, default=15)
    ap.add_argument("--syst-results-dir", type=str, required=True)
    ap.add_argument("--genie-cov-pkl", type=str, default=None)
    ap.add_argument("--save-fig", action="store_true")
    args = ap.parse_args()

    today_str = datetime.now().strftime("%Y%m%d")
    save_fig_dir = path.join(
        save_fig_base_dir,
        "selected_events-intermediate-{}/chunk{}".format(today_str, args.chunk_idx),
    )
    if args.save_fig:
        makedirs(save_fig_dir, exist_ok=True)
        print("Saving plots under", save_fig_dir)

    dfs = get_ana_dfs(option="event_selection")
    mc_evt_df = dfs["mc"].copy()
    data_evt_df_ = dfs["data"]
    intime_evt_df = dfs["intime"].copy()
    data_hdr_df_ = dfs["data_hdr"]
    intime_hdr_df = dfs["intime_hdr"]
    mc_hdr_df = dfs["mc_hdr"]

    mc_evt_df.loc[mc_evt_df.mc.iscc.isna(), ("mc", "iscc")] = 999
    data_evt_df_["mc", "iscc"] = 999

    _sorted = data_hdr_df_.sort_values(["run", "evt"], kind="mergesort")
    data_hdr_df_splits = [
        _sorted.iloc[idx]
        for idx in np.array_split(np.arange(len(_sorted)), args.n_time_splits)
    ]
    data_hdr_df = data_hdr_df_splits[args.chunk_idx]

    common_idxs = data_evt_df_.reset_index(level=[2]).index.intersection(data_hdr_df.index)
    data_evt_df = (
        data_evt_df_.reset_index(level=[2])
        .loc[common_idxs]
        .reset_index()
        .set_index(["__ntuple", "entry", "rec.slc..index"])
    )

    data_tot_pot = data_hdr_df["pot"].sum()
    data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
    pot_str = get_pot_str(data_tot_pot)
    pot_label_chunk = f"Events / Bin (POT={pot_str})"
    data_gates = data_hdr_df.nbnbinfo.sum()

    mc_tot_pot = mc_hdr_df["pot"].sum()
    mc_pot_scale = data_tot_pot / mc_tot_pot
    mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))

    intime_gates = intime_hdr_df[intime_hdr_df["first_in_subrun"] == 1]["noffbeambnb"].sum()
    f_beam = 0.0753
    scale_intime_to_lightdata = (1 - f_beam) * data_gates / intime_gates
    intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))

    plotter = partial(
        overlay_hists,
        mc_df=mc_evt_df,
        data_df=data_evt_df,
        intime_df=intime_evt_df,
        dirt_df=None,
        ax_ylim_ratio=1.5,
        ratio=True,
        textloc=[0.05, 0.55],
        approval="internal",
        plot=False,
        save_fig=args.save_fig,
        cosmic_estimate="intime",
    )

    var_configs = list(INTERMEDIATE_CUT_SYST_VARIABLE_CONFIGS)

    for var_config in var_configs:
        cov, _frac = load_syst_cov_bundle(var_config, args.syst_results_dir, args.genie_cov_pkl)
        plot_labels_hist = [var_config.var_labels[1], pot_label_chunk, ""]
        sn = path.join(save_fig_dir, "{}_topology".format(var_config.var_save_name))
        plotter(
            breakdown_type="topology",
            var_config=var_config,
            plot_labels=plot_labels_hist,
            syst=cov,
            textchi2=True,
            save_name=sn,
        )


if __name__ == "__main__":
    main()
