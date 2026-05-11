#!/usr/bin/env python
"""Cosmic background systematic from offbeam data vs intime MC (unisim by default).

Loads **only** offbeam and intime samples — no neutrino MC multisim weights.

For each kinematic variable (default ``--cv-mode offbeam``):

  * **CV**: histogram of **offbeam data** (counts per bin).
  * **Single unisim variation**: histogram of **intime MC** with ``pot_scale`` applied (same
    gates scaling as ``files_config.get_ana_dfs(option=\"cosmics_systs\")``).

Covariance uses ``pyanalib.covariance.get_covariance_matrix`` with one alternate universe
row (intime) relative to the CV (offbeam). Optional modes:

  * ``mean``: two universe rows (offbeam, intime) with CV = bin-wise average (legacy bracket).
  * ``intime``: CV = intime, single variation = offbeam.

Outputs ``cosmics_syst_dict.npz`` keyed by ``var_save_name``, inner key ``Cosmics`` (compatible
with ``utils.get_syst_unc``): top-level ``cov`` / ``cov_frac`` / ``corr`` plus ``univ_offbeam``
and ``univ_intime`` histograms for bookkeeping.

Example::

    export NUMUCC_SYST_DISK_ROOT=/path/to/syst_disk
    python get_systematics_cosmics.py

Or pass ``--syst-disk-root`` explicitly; plots and ``Cosmics/cosmics_syst_dict.npz`` always go under
``<root>/Cosmics/`` (see ``syst_disk_layout``).

"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from os import makedirs, path
from typing import Any, Dict, List, Mapping, MutableMapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from pyanalib.covariance import get_covariance_matrix  # noqa: E402
from pyanalib.split_df_helpers import generate_tags, load_and_concat_mc_dfs  # noqa: E402

from analysis_village.numucc_1p0pi.files_config import (  # noqa: E402
    file_dir as FILES_DEFAULT_FILE_DIR,
    get_ana_dfs,
    n_max_concat as FILES_DEFAULT_N_MAX_CONCAT,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import SUB_COSMICS, SYST_DISK_ENV  # noqa: E402
from analysis_village.numucc_1p0pi.utils import dpi, fig_ext, plot_heatmap, plot_univ_hists  # noqa: E402
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)

try:
    plt.style.use(path.join(path.dirname(__file__), "presentation.mplstyle"))
except Exception:
    pass


def _sanitize_matrix_pack(pack: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    cov = np.nan_to_num(np.asarray(pack["cov"], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    cov_frac = np.nan_to_num(np.asarray(pack["cov_frac"], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    corr = np.asarray(pack["corr"], dtype=float)
    d = np.sqrt(np.maximum(np.diag(cov), 0.0))
    outer = np.outer(d, d)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(outer > 0, cov / outer, 0.0)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return {"cov": cov, "cov_frac": cov_frac, "corr": corr}


def build_variable_configs(arg_vars: Sequence[str] | None) -> List[Any]:
    registry = {
        "integrated": VariableConfig.all_events,
        "vertex_x": VariableConfig.vertex_x,
        "vertex_y": VariableConfig.vertex_y,
        "vertex_z": VariableConfig.vertex_z,
        "muon-p": VariableConfig.muon_momentum,
        "muon-dir_z": VariableConfig.muon_direction,
        "muon-dir_x": VariableConfig.muon_direction_x,
        "muon-dir_y": VariableConfig.muon_direction_y,
        "muon-dir_phi": VariableConfig.muon_direction_phi,
        "proton-p": VariableConfig.proton_momentum,
        "proton-dir_z": VariableConfig.proton_direction,
        "proton-dir_x": VariableConfig.proton_direction_x,
        "proton-dir_y": VariableConfig.proton_direction_y,
        "proton-dir_phi": VariableConfig.proton_direction_phi,
        "muon-end_x": VariableConfig.muon_end_x,
        "muon-end_y": VariableConfig.muon_end_y,
        "muon-end_z": VariableConfig.muon_end_z,
        "opening_angle": VariableConfig.opening_angle,
        "tki-del_alpha": VariableConfig.tki_del_alpha,
        "tki-del_phi": VariableConfig.tki_del_phi,
        "tki-del_Tp": VariableConfig.tki_del_Tp,
        "tki-del_p": VariableConfig.tki_del_p,
        "tki-del_Tp_x": VariableConfig.tki_del_Tp_x,
        "tki-del_Tp_y": VariableConfig.tki_del_Tp_y,
    }
    if arg_vars:
        out = []
        for name in arg_vars:
            key = name.strip()
            if key not in registry:
                raise ValueError(f"Unknown variable key '{key}'. Choices: {sorted(registry)}")
            out.append(registry[key]())
        return out
    return with_final_selected_evt_variables(
        list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
        + [
            VariableConfig.vertex_x(),
            VariableConfig.vertex_y(),
            VariableConfig.vertex_z(),
            VariableConfig.muon_direction_x(),
            VariableConfig.muon_direction_y(),
            VariableConfig.proton_direction_x(),
            VariableConfig.proton_direction_y(),
            VariableConfig.opening_angle(),
        ]
    )


def cosmic_histograms(
    offbeam_df: pd.DataFrame,
    intime_df: pd.DataFrame,
    var_config: Any,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (h_offbeam, h_intime_scaled) with same binning as ``var_config``."""
    col = var_config.var_evt_reco_col
    bins = var_config.bins

    h_off, _ = np.histogram(offbeam_df[col], bins=bins)
    h_in, _ = np.histogram(intime_df[col], weights=intime_df["pot_scale"], bins=bins)

    if var_config.var_save_name == "integrated":
        h_off = np.array([float(len(offbeam_df))])
        h_in = np.array([float(np.sum(intime_df["pot_scale"]))])

    return h_off.astype(float), h_in.astype(float)


def univ_stack_and_cv(h_off: np.ndarray, h_in: np.ndarray, mode: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return (univ_events, cv_events) for ``get_covariance_matrix``."""
    mode = mode.strip().lower()
    h_off = np.asarray(h_off, dtype=float)
    h_in = np.asarray(h_in, dtype=float)
    if mode == "offbeam":
        # CV = data-driven cosmic estimate; intime MC is the single alternate shape
        return np.stack([h_in], axis=0), h_off.copy()
    if mode == "intime":
        return np.stack([h_off], axis=0), h_in.copy()
    if mode == "mean":
        return np.stack([h_off, h_in], axis=0), 0.5 * (h_off + h_in)
    raise ValueError(f"Unknown --cv-mode {mode!r}; use mean|intime|offbeam")


def _plot_matrices(
    ret: Mapping[str, np.ndarray],
    var_config: Any,
    save_fig_dir: str,
    save_plots: bool,
    prefix: str,
) -> None:
    labels = {
        "cov": "Covariance",
        "cov_frac": "Fractional covariance",
        "corr": "Correlation",
    }
    for matrix_type in ("cov", "cov_frac", "corr"):
        save_fig_name = path.join(save_fig_dir, f"{prefix}-{matrix_type}")
        plot_heatmap(
            ret[matrix_type],
            var_config.bins,
            plot_labels=[var_config.var_labels[1], var_config.var_labels[1], labels[matrix_type]],
            plot=False,
            save_fig=save_plots,
            save_name=save_fig_name,
        )


def process_variable_cosmics(
    offbeam_df: pd.DataFrame,
    intime_df: pd.DataFrame,
    var_config: Any,
    cv_mode: str,
    save_fig_dir: str,
    save_plots: bool,
) -> Dict[str, Any]:
    h_off, h_in = cosmic_histograms(offbeam_df, intime_df, var_config)
    univ_events, cv_events = univ_stack_and_cv(h_off, h_in, cv_mode)

    if save_plots:
        plt.figure(figsize=(8, 6))
        mode_l = cv_mode.strip().lower()
        if mode_l == "offbeam":
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                linewidth=2.0,
                label="CV (offbeam data)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                label="Unisim variation (intime MC)",
            )
        elif mode_l == "intime":
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                linewidth=2.0,
                label="CV (intime MC)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                label="Unisim variation (offbeam data)",
            )
        else:
            plt.hist(
                var_config.bin_centers,
                weights=h_off,
                bins=var_config.bins,
                histtype="step",
                color="crimson",
                label="Offbeam data",
            )
            plt.hist(
                var_config.bin_centers,
                weights=h_in,
                bins=var_config.bins,
                histtype="step",
                color="black",
                label="Intime MC (scaled)",
            )
            plt.hist(
                var_config.bin_centers,
                weights=cv_events,
                bins=var_config.bins,
                histtype="step",
                color="tab:blue",
                linestyle="--",
                label="CV (mean)",
            )
        plt.xlim(var_config.bins[0], var_config.bins[-1])
        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.legend(frameon=False)
        plt.savefig(
            path.join(save_fig_dir, f"{var_config.var_save_name}-cosmics-offbeam-vs-intime{fig_ext}"),
            bbox_inches="tight",
            dpi=dpi,
        )
        plt.close()

    ret = _sanitize_matrix_pack(get_covariance_matrix(univ_events, cv_events))

    if save_plots:
        plot_univ_hists(
            univ_events,
            cv_events,
            "Cosmics",
            var_config,
            plot=False,
            ax_titles=[
                "",
                "Events / Bin",
                f"Cosmics ({cv_mode}: CV vs variation, n_univ={univ_events.shape[0]})",
            ],
            save_fig=True,
            save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-cosmics_univ"),
        )
        _plot_matrices(
            ret,
            var_config,
            save_fig_dir,
            True,
            f"{var_config.var_save_name}-cosmics",
        )

    return {
        "cov": ret["cov"],
        "cov_frac": ret["cov_frac"],
        "corr": ret["corr"],
        "rate": ret,
        "univ_offbeam": h_off,
        "univ_intime": h_in,
        "cv_mode": cv_mode,
        "cv_histogram": cv_events,
        "n_univ": int(univ_events.shape[0]),
    }


def save_cosmics_npz(syst_dict_by_var: Mapping[str, Any], out_path: str) -> None:
    np.savez_compressed(out_path, **syst_dict_by_var)
    logging.info("Wrote %s", out_path)


def load_cosmic_samples_ana() -> Tuple[pd.DataFrame, pd.DataFrame]:
    dfs = get_ana_dfs(option="cosmics_systs")
    return dfs["data"], dfs["mc"]


def load_cosmic_samples_concat(
    file_dir: str,
    offbeam_chunk_tags: List[str],
    intime_chunk_tags: List[str],
    n_max_concat: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    data_dfs = load_and_concat_mc_dfs(
        file_dir=file_dir,
        chunk_tags=offbeam_chunk_tags,
        df_tag="",
        keys2load=["evt", "hdr"],
        n_max_concat=n_max_concat,
        sub_dir="data",
        sample_dir="OffBeam",
    )
    data_hdr = data_dfs["hdr"]
    data_evt = data_dfs["evt"]
    data_gates = data_hdr[data_hdr["first_in_subrun"] == 1]["noffbeambnb"].sum()

    mc_dfs = load_and_concat_mc_dfs(
        file_dir=file_dir,
        chunk_tags=intime_chunk_tags,
        df_tag="",
        keys2load=["evt", "hdr"],
        n_max_concat=n_max_concat,
        sub_dir="MC",
        sample_dir="intime",
    )
    mc_hdr = mc_dfs["hdr"]
    mc_evt = mc_dfs["evt"]
    mc_gates = mc_hdr[mc_hdr["first_in_subrun"] == 1]["ngenevt"].sum()

    scale = float(data_gates / mc_gates) if mc_gates else 1.0
    logging.info("Cosmic gates scale data/mc = %.6f (offbeam %.4e / intime %.4e)", scale, data_gates, mc_gates)
    mc_evt = mc_evt.copy()
    mc_evt["pot_scale"] = scale
    data_evt = data_evt.copy()
    data_evt["pot_scale"] = 1.0
    return data_evt, mc_evt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out-tag",
        default=None,
        help="Optional label for logs only (default: today); output path is always <syst-disk>/Cosmics/.",
    )
    p.add_argument(
        "--syst-disk-root",
        default=None,
        help=(
            "Syst disk layout root; writes plots + NPZ under <root>/%s/. "
            "If omitted, uses environment variable %s." % (SUB_COSMICS, SYST_DISK_ENV)
        ),
    )
    p.add_argument("--error-log", default=None, help="Per-variable failure log path")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--no-save-npz", action="store_true")
    p.add_argument(
        "--cv-mode",
        choices=("mean", "intime", "offbeam"),
        default="offbeam",
        help="offbeam: CV=offbeam, one variation=intime (default); intime: swapped; mean: 2 rows + CV=avg",
    )
    p.add_argument("--data-source", choices=("ana", "concat"), default="ana")
    p.add_argument("--file-dir", default=FILES_DEFAULT_FILE_DIR)
    p.add_argument("--n-max-concat", type=int, default=FILES_DEFAULT_N_MAX_CONCAT)
    p.add_argument(
        "--offbeam-chunk-tags",
        default=None,
        help="Comma-separated tags for OffBeam concat (default: generate_tags('ad'))",
    )
    p.add_argument(
        "--intime-chunk-tags",
        default=None,
        help="Comma-separated tags for intime concat (default: generate_tags('au'))",
    )
    p.add_argument("--vars", nargs="*", default=None)
    args = p.parse_args()
    root = args.syst_disk_root or os.environ.get(SYST_DISK_ENV)
    if not root:
        p.error(
            "Pass --syst-disk-root or set %s (outputs always go to <root>/%s/)."
            % (SYST_DISK_ENV, SUB_COSMICS)
        )
    args.syst_disk_root = path.abspath(path.expanduser(root.rstrip("/")))
    return args


def main() -> None:
    args = parse_args()
    tag = args.out_tag or datetime.now().strftime("%Y%m%d")
    save_fig_dir = path.join(args.syst_disk_root, SUB_COSMICS)
    makedirs(save_fig_dir, exist_ok=True)

    log_path = args.error_log or path.join(save_fig_dir, "cosmics_failures.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    logger = logging.getLogger("cosmics_syst")
    logger.info("Run tag=%s syst_disk_root=%s -> %s", tag, args.syst_disk_root, save_fig_dir)

    var_configs = build_variable_configs(args.vars)
    save_plots = not args.no_plots

    if args.data_source == "ana":
        offbeam_df, intime_df = load_cosmic_samples_ana()
    else:
        ot = (
            [t.strip() for t in args.offbeam_chunk_tags.split(",") if t.strip()]
            if args.offbeam_chunk_tags
            else generate_tags("ad")
        )
        it = (
            [t.strip() for t in args.intime_chunk_tags.split(",") if t.strip()]
            if args.intime_chunk_tags
            else generate_tags("au")
        )
        offbeam_df, intime_df = load_cosmic_samples_concat(
            args.file_dir, ot, it, args.n_max_concat
        )

    logger.info(
        "Loaded offbeam evt=%d intime evt=%d (cv-mode=%s)",
        len(offbeam_df),
        len(intime_df),
        args.cv_mode,
    )

    syst_dict: Dict[str, MutableMapping[str, Any]] = {}

    for var_config in tqdm(var_configs, desc="cosmics"):
        try:
            pay = process_variable_cosmics(
                offbeam_df,
                intime_df,
                var_config,
                args.cv_mode,
                save_fig_dir,
                save_plots,
            )
            syst_dict.setdefault(var_config.var_save_name, {})
            syst_dict[var_config.var_save_name]["Cosmics"] = pay
        except Exception:
            logger.error(
                "FAILED variable=%s\n%s",
                var_config.var_save_name,
                traceback.format_exc(),
            )

    if not args.no_save_npz and syst_dict:
        # Filename fixed by syst_disk_layout.Cosmics
        save_cosmics_npz(syst_dict, path.join(save_fig_dir, "cosmics_syst_dict.npz"))

    logger.info("Done -> %s (log %s)", save_fig_dir, log_path)


if __name__ == "__main__":
    main()
