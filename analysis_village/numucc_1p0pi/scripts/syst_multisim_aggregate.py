#!/usr/bin/env python
"""Reduce phase for chunked Flux/G4/MCstat covariances (+ optional cosmics).

See module docstring in previous revision — summary:

* Input: ``nu__*.pkl`` from ``syst_multisim_chunk.py``.
* Sum ``univ_events`` / ``cv_events``, build covariances, optional cosmics via
  ``files_config.get_ana_dfs`` (same monolithic loads as the legacy script).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
from os import path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.covariance import get_covariance_matrix
from analysis_village.numucc_1p0pi.constants import dpi
from analysis_village.numucc_1p0pi.files_config import get_ana_dfs
from analysis_village.numucc_1p0pi.syst_disk_layout import SUB_COSMICS
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    NEUTRINO_SYST_ORDER,
    build_var_configs,
    drop_bad_g4_weights,
    save_cosmics_legacy_npz,
    save_neutrino_multisim_npzs,
)
from analysis_village.numucc_1p0pi.utils import (
    generate_tags,
    get_univ_rates,
    plot_heatmap,
    plot_univ_hists,
)


def collect_nu_chunks(chunks_dir: str) -> list[str]:
    return sorted(glob.glob(path.join(chunks_dir, "nu__*.pkl")))


def _write_covariance_manifest(
    out_dir: str,
    syst_dict: dict,
    merged_meta: list,
    mc_df_stage: str,
    var_set: str,
    cosmics_npz_dir: str | None = None,
) -> None:
    """Small JSON sidecar listing variables and NPZ outputs for downstream tools."""
    neutrino_npz = [
        "mcstat_syst_dict.npz",
        "flux_syst_dict.npz",
        "g4_syst_dict.npz",
    ]
    var_names: set[str] = set()
    for _cat, block in (syst_dict or {}).items():
        if isinstance(block, dict):
            var_names.update(block.keys())
    cosmics_block = (syst_dict or {}).get("cosmics")
    cosmics_entry = None
    if cosmics_block:
        cosmics_entry = {
            "file": "cosmics_syst_dict.npz",
            "directory": cosmics_npz_dir,
            "note": "Cosmic unisim — not neutrino multisim; separate from neutrino multisim NPZs.",
        }
    manifest = {
        "schema": "numucc_multisim_covariance_v1",
        "description": "Fractional and absolute covariance packs per category; NPZs use numpy.savez_compressed",
        "mc_df_stage": mc_df_stage,
        "var_set": var_set,
        "categories_present": sorted([k for k, v in (syst_dict or {}).items() if v]),
        "variables": sorted(var_names),
        "neutrino_multisim_npz": [
            {"file": n, "role": "per-variable dict → inner syst matrices"} for n in neutrino_npz
        ],
        "cosmics_npz": cosmics_entry,
        "map_shard_metadata": merged_meta,
    }
    outp = path.join(out_dir, "covariance_manifest.json")
    with open(outp, "w") as f:
        json.dump(manifest, f, indent=2)
    print("[multisim-agg] wrote", outp)


def merge_nu_chunks(paths: list[str]) -> dict:
    merged = None
    for fp in tqdm(paths, desc="merge nu chunks"):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        if merged is None:
            merged = {"syst": {sn: {} for sn in NEUTRINO_SYST_ORDER}, "meta": []}
        merged["meta"].append(
            {
                "df_file": d.get("df_file"),
                "splits": d.get("splits_processed"),
                "syst_names_computed": d.get("syst_names_computed"),
            }
        )
        raw_syst = d.get("syst") or {}
        for sn in NEUTRINO_SYST_ORDER:
            block = raw_syst.get(sn, {})
            for vsn, pack in block.items():
                u = np.asarray(pack["univ_events"], dtype=float)
                c = np.asarray(pack["cv_events"], dtype=float)
                if vsn not in merged["syst"][sn]:
                    merged["syst"][sn][vsn] = {"univ_events": u.copy(), "cv_events": c.copy()}
                else:
                    mu = merged["syst"][sn][vsn]
                    if mu["univ_events"].shape != u.shape:
                        raise ValueError(
                            "shape mismatch %s %s: %s vs %s"
                            % (sn, vsn, mu["univ_events"].shape, u.shape)
                        )
                    mu["univ_events"] += u
                    mu["cv_events"] += c
    if merged is None:
        raise RuntimeError("no chunks merged")
    return merged


def _syst_plot_key(sn: str):
    return ("mc", sn) if sn in ("Flux", "G4") else sn


def covariance_dict_from_merged(
    merged: dict,
    var_configs: list,
    syst_disk_root: str,
    save_fig: bool,
) -> dict:
    syst_dict = {sn: {} for sn in NEUTRINO_SYST_ORDER}
    vc_by = {v.var_save_name: v for v in var_configs}
    for sn in NEUTRINO_SYST_ORDER:
        cat_dir = path.join(syst_disk_root, sn)
        os.makedirs(cat_dir, exist_ok=True)
        sk = _syst_plot_key(sn)
        tag = sn
        for vsn, pack in tqdm(list(merged["syst"][sn].items()), desc="cov %s" % sn):
            vc = vc_by.get(vsn)
            if vc is None:
                continue
            univ = np.asarray(pack["univ_events"], dtype=float)
            cv = np.asarray(pack["cv_events"], dtype=float)
            ret = get_covariance_matrix(univ, cv)
            syst_dict[sn][vsn] = ret
            if save_fig:
                plot_univ_hists(
                    univ,
                    cv,
                    sk,
                    vc,
                    save_fig=True,
                    save_name=path.join(cat_dir, "{}-{}-universes".format(vsn, tag)),
                )
                for matrix_type in ["cov", "cov_frac", "corr"]:
                    plot_heatmap(
                        ret[matrix_type],
                        vc.bins,
                        plot_labels=[
                            vc.var_labels[1],
                            vc.var_labels[1],
                            matrix_type.capitalize(),
                        ],
                        save_fig=True,
                        save_name=path.join(
                            cat_dir,
                            "{}-{}-{}.pdf".format(vsn, tag, matrix_type),
                        ),
                    )
    return syst_dict


def process_cosmics(
    mc_evt_df: pd.DataFrame,
    offbeam_data_df: pd.DataFrame,
    intime_mc_df: pd.DataFrame,
    var_config,
    syst_dict: dict,
    save_fig_dir: str,
    save_fig: bool,
) -> dict:
    univ_events, cv_events = get_univ_rates(
        evtdf=mc_evt_df,
        var_config=var_config,
        n_univ=1,
        bkgd_subtract=True,
        syst_name=("mc", "Flux"),
    )

    data_events, _ = np.histogram(offbeam_data_df[var_config.var_evt_reco_col], bins=var_config.bins)
    mc_events, _ = np.histogram(
        intime_mc_df[var_config.var_evt_reco_col],
        weights=intime_mc_df["pot_scale"],
        bins=var_config.bins,
    )

    if var_config.var_save_name == "integrated":
        data_events = np.array([len(offbeam_data_df)])
        mc_events = np.array([len(intime_mc_df)])

    mc_events = mc_events.copy()
    data_events = data_events.copy()
    mc_events[mc_events == 0] = 1
    data_events[data_events == 0] = 1

    if save_fig:
        plt.subplots()
        plt.hist(
            var_config.bin_centers,
            weights=data_events,
            bins=var_config.bins,
            histtype="step",
            color="red",
            label="Data",
        )
        plt.hist(
            var_config.bin_centers,
            weights=mc_events,
            bins=var_config.bins,
            histtype="step",
            color="black",
            label="MC",
        )
        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.legend()
        plt.savefig(
            path.join(save_fig_dir, "{}-cosmics-universes.pdf".format(var_config.var_save_name)),
            bbox_inches="tight",
            dpi=dpi,
        )
        plt.close()

    univ_events = np.array([cv_events + (data_events - mc_events)])
    ret = get_covariance_matrix(univ_events, cv_events)
    if save_fig:
        plot_heatmap(
            ret["cov"],
            var_config.bins,
            plot_labels=[var_config.var_labels[1], var_config.var_labels[1], "Covariance"],
            save_fig=True,
            save_name=path.join(
                save_fig_dir,
                "{}-cosmics-cov.pdf".format(var_config.var_save_name),
            ),
        )
    syst_dict["cosmics"][var_config.var_save_name] = ret
    return syst_dict


def load_mc_evt_for_cosmics(mc_df_stage: str) -> pd.DataFrame:
    if mc_df_stage == "final":
        ret = get_ana_dfs(option="systs", systs_mc_df_tag="", systs_chunk_tags=None)
    else:
        ret = get_ana_dfs(
            option="systs",
            systs_mc_df_tag="-sel_all-wgts",
            systs_chunk_tags=generate_tags("ah")[1:],
        )
    return drop_bad_g4_weights(ret["evt"])


def run_syst_multisim_aggregate(
    chunks_dir: str,
    syst_disk_root: str,
    mc_df_stage: str = "final",
    var_set: str = "final",
    skip_cosmics: bool = False,
    no_plots: bool = False,
    no_legacy_npz: bool = False,
) -> None:
    os.makedirs(syst_disk_root, exist_ok=True)
    save_fig = not no_plots

    var_configs = build_var_configs(var_set)
    ck = collect_nu_chunks(chunks_dir)
    if not ck:
        raise RuntimeError("[multisim-agg] no nu__*.pkl under %s" % chunks_dir)
    print("[multisim-agg] merging %d nu chunk(s)" % len(ck))
    merged = merge_nu_chunks(ck)
    syst_dict = covariance_dict_from_merged(merged, var_configs, syst_disk_root, save_fig)

    cosmics_write_dir = None
    if not skip_cosmics:
        cosmics_write_dir = path.join(syst_disk_root, SUB_COSMICS)
        os.makedirs(cosmics_write_dir, exist_ok=True)
        print(
            "[multisim-agg] cosmics via files_config (monolithic MC + offbeam/intime) → %s"
            % cosmics_write_dir
        )
        mc_evt_full = load_mc_evt_for_cosmics(mc_df_stage)
        dfs = get_ana_dfs(option="cosmics_systs")
        syst_dict["cosmics"] = {}
        for vc in tqdm(var_configs, desc="cosmics"):
            process_cosmics(
                mc_evt_full,
                dfs["data"],
                dfs["mc"],
                vc,
                syst_dict,
                cosmics_write_dir,
                save_fig,
            )

    if not no_legacy_npz:
        save_neutrino_multisim_npzs(syst_dict, syst_disk_root)
        print("[multisim-agg] neutrino multisim NPZs →", syst_disk_root)
        if not skip_cosmics and cosmics_write_dir is not None:
            save_cosmics_legacy_npz(syst_dict, syst_disk_root)
            print("[multisim-agg] cosmics NPZ →", cosmics_write_dir)

    _write_covariance_manifest(
        syst_disk_root,
        syst_dict,
        merged["meta"],
        mc_df_stage,
        var_set,
        cosmics_npz_dir=cosmics_write_dir,
    )

    summ = path.join(syst_disk_root, "syst_multisim_aggregate_summary.txt")
    with open(summ, "w") as f:
        f.write("# merged nu chunk metadata\n")
        for m in merged["meta"]:
            f.write("%s\n" % m)
    print("[multisim-agg] wrote", summ)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chunks_dir", required=True)
    p.add_argument(
        "--syst-disk-root",
        "--out_dir",
        "--out-dir",
        dest="syst_disk_root",
        required=True,
        help="Syst disk root (see analysis_village.numucc_1p0pi.syst_disk_layout): writes "
        "MCstat/, Flux/, G4/, and Cosmics/ subfolders.",
    )
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--var-set", choices=("final", "intermediate", "both"), default="final")
    p.add_argument("--skip-cosmics", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--no-legacy-npz", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    run_syst_multisim_aggregate(
        chunks_dir=args.chunks_dir,
        syst_disk_root=args.syst_disk_root,
        mc_df_stage=args.mc_df_stage,
        var_set=args.var_set,
        skip_cosmics=args.skip_cosmics,
        no_plots=args.no_plots,
        no_legacy_npz=args.no_legacy_npz,
    )


if __name__ == "__main__":
    main()
