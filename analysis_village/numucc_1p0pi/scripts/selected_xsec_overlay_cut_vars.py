#!/usr/bin/env python3
"""Data–MC overlays of the cut variables varied in the sel_mup cut campaign.

Reuses the same ``PLOT_SETS`` / loaders as ``selected_xsec_overlay.py`` and writes
PNGs into the same ``FixedDev/selected_xsec_overlay_cuts/<tag>/`` directories so
the applied thresholds can be checked visually (distribution edges vs cut lines).
"""

from __future__ import annotations

import pickle
import sys
from datetime import datetime
from os import makedirs, path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_REPO = "/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana"
_SCRIPTS = path.join(_REPO, "analysis_village/numucc_1p0pi/scripts")
sys.path.append(_REPO)
sys.path.append(_SCRIPTS)

from analysis_village.numucc_1p0pi.selected_xsec_overlay_hist import (  # noqa: E402
    build_overlay_histdata_map,
    overlay_hists_from_counts,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
import selected_xsec_overlay as sxo  # noqa: E402

CUTVAR_PKL_NAME = "cutvar_overlay_histdata.pkl"
FORCE_REBUILD_COUNTS = False
BREAKDOWN_TYPES = ("topology", "genie_sb")

# Per-campaign applied thresholds (nominal elsewhere).
# nu_score: keep > th; chi2_avg_mu: keep < th; |mcs_range_diff|: keep < th;
# vertex_z: exclude [z0, z1] when set.
CUT_APPLIED: Dict[str, dict] = {
    "nu_score0": {
        "nu_score_th": 0.0,
        "chi2mu_th": 30.0,
        "qual_th": 0.2,
        "vz_exclude": None,
    },
    "chi2mu15": {
        "nu_score_th": 0.45,
        "chi2mu_th": 15.0,
        "qual_th": 0.2,
        "vz_exclude": None,
    },
    "chi2mu45": {
        "nu_score_th": 0.45,
        "chi2mu_th": 45.0,
        "qual_th": 0.2,
        "vz_exclude": None,
    },
    "mcs_range_diff1p0": {
        "nu_score_th": 0.45,
        "chi2mu_th": 30.0,
        "qual_th": 1.0,
        "vz_exclude": None,
    },
    "vz_exclude_200_300": {
        "nu_score_th": 0.45,
        "chi2mu_th": 30.0,
        "qual_th": 0.2,
        "vz_exclude": (200.0, 300.0),
    },
}


def _cut_var_configs() -> List[VariableConfig]:
    """Configs pointing at final-selected (mu / slc) columns on sel_mup evt."""
    return [
        VariableConfig(
            var_save_name="nu_score",
            var_plot_name=r"$\nu_{\mathrm{score}}$",
            var_labels=[
                r"Neutrino-like Score",
                r"Neutrino-like Score",
                "",
            ],
            bins=np.linspace(0.0, 1.0, 51),
            var_evt_reco_col=("slc", "nu_score", "", "", "", "", ""),
            var_evt_truth_col=("", "", "", "", "", "", ""),
            var_nu_col=("", "", ""),
            xsec_label="",
        ),
        VariableConfig(
            var_save_name="chi2_avg_mu",
            var_plot_name=r"$\chi^2_{\mu,\mathrm{avg}}$",
            var_labels=[
                r"$\mathrm{\chi^{2}_{\mu,\,avg}}$",
                r"$\mathrm{\chi^{2}_{\mu,\,avg,\,\mathrm{reco.}}}$",
                "",
            ],
            bins=np.linspace(0.0, 60.0, 61),
            var_evt_reco_col=("mu", "pfp", "trk", "chi2pid", "avg", "chi2_muon", ""),
            var_evt_truth_col=("", "", "", "", "", "", ""),
            var_nu_col=("", "", ""),
            xsec_label="",
        ),
        VariableConfig(
            var_save_name="mcs_range_diff",
            var_plot_name="MCS Range Difference",
            var_labels=[
                r"$\mathrm{(P_{Range} - P_{MCS}) \, / \, P_{Range}}$",
                r"$\mathrm{(P_{Range} - P_{MCS}) / P_{Range}}$",
                "",
            ],
            # Wide enough for QUAL_TH=1.0 selected sample.
            bins=np.linspace(-1.2, 1.2, 49),
            var_evt_reco_col=("mu", "pfp", "trk", "mcs_range_diff", "", "", ""),
            var_evt_truth_col=("", "", "", "", "", "", ""),
            var_nu_col=("", "", ""),
            xsec_label="",
        ),
        VariableConfig(
            var_save_name="vertex_z",
            var_plot_name="Neutrino Vertex Z [cm]",
            var_labels=[
                "Neutrino Vertex Z [cm]",
                "Slice Vertex Z [cm]",
                "",
            ],
            bins=np.linspace(0.0, 500.0, 51),
            var_evt_reco_col=("slc", "vertex", "z", "", "", "", ""),
            var_evt_truth_col=("", "", "", "", "", "", ""),
            var_nu_col=("", "", ""),
            xsec_label="",
        ),
    ]


def _vlines_for_tag(tag: str, var_save_name: str) -> Optional[list]:
    """Cut markers for ``overlay_hists`` ``vline`` (direction 0=left, 1=right)."""
    th = CUT_APPLIED[tag]
    if var_save_name == "nu_score":
        return [[th["nu_score_th"], 1]]
    if var_save_name == "chi2_avg_mu":
        return [[th["chi2mu_th"], 0]]
    if var_save_name == "mcs_range_diff":
        q = th["qual_th"]
        return [[-q, 1], [q, 0]]
    if var_save_name == "vertex_z":
        vz = th["vz_exclude"]
        if vz is None:
            return None
        z0, z1 = vz
        # Arrows point into the kept regions (outside the excluded window).
        return [[z0, 0], [z1, 1]]
    return None


def cutvar_pkl_path(out_dir: str) -> str:
    return path.join(out_dir, CUTVAR_PKL_NAME)


def load_cutvar_counts(out_dir: str) -> Optional[dict]:
    pkl = cutvar_pkl_path(out_dir)
    if not path.isfile(pkl):
        return None
    with open(pkl, "rb") as f:
        payload = pickle.load(f)
    if not isinstance(payload, dict) or "histdata" not in payload:
        raise ValueError(f"unrecognized cutvar counts payload: {pkl}")
    return payload


def save_cutvar_counts(
    out_dir: str,
    histdata_map: dict,
    *,
    pot_label: str,
    plot_set: dict,
    var_save_names: Sequence[str],
) -> str:
    makedirs(out_dir, exist_ok=True)
    pkl = cutvar_pkl_path(out_dir)
    payload = {
        "format": "selected_xsec_overlay_cutvar_histdata_v1",
        "pot_label": pot_label,
        "plot_set": dict(plot_set),
        "var_save_names": list(var_save_names),
        "breakdown_types": list(BREAKDOWN_TYPES),
        "histdata": dict(histdata_map),
    }
    with open(pkl, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    return pkl


def plot_cutvar_map(
    histdata_map: dict,
    var_configs: Sequence[VariableConfig],
    *,
    tag: str,
    pot_label: str,
    out_dir: str,
) -> None:
    makedirs(out_dir, exist_ok=True)
    for var_config in var_configs:
        vline = _vlines_for_tag(tag, var_config.var_save_name)
        for breakdown_type in BREAKDOWN_TYPES:
            key = (var_config.var_save_name, breakdown_type)
            hd = histdata_map.get(key)
            if hd is None:
                raise KeyError(f"missing histdata for {key}")
            save_name = path.join(out_dir, f"{key[0]}_{key[1]}")
            print(f"  plot {key[0]} ({key[1]})  vline={vline}", flush=True)
            overlay_hists_from_counts(
                hd,
                var_config=var_config,
                plot_labels=[var_config.var_labels[1], pot_label, ""],
                ax_ylim_ratio=sxo.AX_YLIM_RATIO,
                ratio=sxo.RATIO,
                textloc=sxo.TEXTLOC,
                approval=sxo.APPROVAL,
                save_fig=sxo.SAVE_FIG,
                plot=sxo.PLOT,
                textchi2=sxo.TEXTCHI2,
                save_name=save_name,
                vline=vline,
                syst=None,
            )


def run_plot_set(plot_set: dict, var_configs: Sequence[VariableConfig]) -> None:
    tag = plot_set["tag"]
    out_dir = plot_set["output_dir"]
    if tag not in CUT_APPLIED:
        raise KeyError(f"no CUT_APPLIED entry for tag={tag!r}")

    print(
        f"\n{'=' * 72}\nCut-var overlay: {tag}\n"
        f"  MC:   {plot_set['mc_dir']}\n"
        f"  data: {plot_set['data_dir']}\n"
        f"  cuts: {CUT_APPLIED[tag]}\n"
        f"  saving -> {out_dir}",
        flush=True,
    )

    payload = None if FORCE_REBUILD_COUNTS else load_cutvar_counts(out_dir)
    if payload is not None:
        print(f"  replot from counts: {cutvar_pkl_path(out_dir)}", flush=True)
        pot_label = payload["pot_label"]
        histdata_map = payload["histdata"]
    else:
        print("  filling cut-var counts from dataframes...", flush=True)
        mc_evt, mc_hdr, n_mc_loaded, n_mc_total = sxo.load_mc_sample(
            plot_set["mc_dir"], plot_set["mc_filename_str"]
        )
        if n_mc_loaded < n_mc_total:
            print(
                f"  MC subsample: {n_mc_loaded}/{n_mc_total} files",
                flush=True,
            )
        data_evt, data_hdr = sxo.load_data_sample(
            plot_set["data_dir"], plot_set["data_filename_str"]
        )
        pot_label = sxo.setup_pot_weights(mc_evt, mc_hdr, data_evt, data_hdr)
        histdata_map = build_overlay_histdata_map(
            var_configs,
            BREAKDOWN_TYPES,
            mc_df=mc_evt,
            data_df=data_evt,
        )
        pkl = save_cutvar_counts(
            out_dir,
            histdata_map,
            pot_label=pot_label,
            plot_set=plot_set,
            var_save_names=[vc.var_save_name for vc in var_configs],
        )
        print(f"  wrote counts -> {pkl}", flush=True)

    plot_cutvar_map(
        histdata_map,
        var_configs,
        tag=tag,
        pot_label=pot_label,
        out_dir=out_dir,
    )


def main() -> None:
    t0 = datetime.now()
    var_configs = _cut_var_configs()
    print(f"selected_xsec_overlay_cut_vars start {t0.isoformat()}", flush=True)
    print(
        f"vars={[vc.var_save_name for vc in var_configs]}  "
        f"mem_now={sxo.get_memory_used_frac() * 100:.1f}%",
        flush=True,
    )
    for plot_set in sxo.PLOT_SETS:
        run_plot_set(plot_set, var_configs)
    print(f"\nDone in {datetime.now() - t0}", flush=True)


if __name__ == "__main__":
    main()
