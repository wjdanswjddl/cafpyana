#!/usr/bin/env python
"""Cosmic background systematic from offbeam data vs intime MC (unisim by default).

**Chunked workflow** (preferred; matches ``syst_multisim_chunk`` / ``syst_detvar_chunk``):

1. Map: ``syst_cosmics_chunk.py`` per ``.df`` from ``dataset_locations.iter_cosmics_chunk_df_paths``
   (offbeam and intime separately) → ``cosmics__<sample>__<stem>.pkl``.
2. Reduce: ``syst_cosmics_aggregate.py --chunks_dir …`` merges pickles, applies global gate
   scaling ``sum(offbeam gates)/sum(intime gates)``, builds covariances, writes
   ``<syst-disk>/Cosmics/cosmics_syst_dict.npz``.

**In-memory run** (``run-ana``): single ``get_ana_dfs(option='cosmics_systs')`` load for
notebooks — not for large-scale production I/O.

For each kinematic variable (default ``--cv-mode offbeam``):

  * **CV**: histogram of **offbeam data** (counts per bin).
  * **Single unisim variation**: histogram of **intime MC** with gate scaling applied.

Covariance uses ``pyanalib.covariance.get_covariance_matrix``. Optional ``--cv-mode``:
``mean`` | ``intime`` | ``offbeam``.

Examples::

    python syst_cosmics_chunk.py --sample offbeam --df_file PATH.df --out_dir CHUNKS
    python syst_cosmics_chunk.py --sample intime --df_file PATH.df --out_dir CHUNKS
    export NUMUCC_SYST_DISK_ROOT=/path/to/syst_disk
    python syst_cosmics_aggregate.py --chunks_dir CHUNKS --syst-disk-root $NUMUCC_SYST_DISK_ROOT

    # one-shot in-memory (same scaling as files_config cosmics_systs)
    python get_systematics_cosmics.py run-ana --syst-disk-root $NUMUCC_SYST_DISK_ROOT

"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from os import makedirs, path
from typing import Any, Dict, Mapping, MutableMapping, Tuple

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

from analysis_village.numucc_1p0pi.files_config import get_ana_dfs  # noqa: E402
from analysis_village.numucc_1p0pi.syst_disk_layout import SUB_COSMICS, SYST_DISK_ENV  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cosmics_common import (  # noqa: E402
    apply_flat_cosmic_uncertainty,
    build_variable_configs,
)
from analysis_village.numucc_1p0pi.utils import dpi, fig_ext, plot_heatmap, plot_univ_hists  # noqa: E402
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


def process_variable_cosmics_from_histograms(
    h_off: np.ndarray,
    h_in: np.ndarray,
    var_config: Any,
    cv_mode: str,
    save_fig_dir: str,
    save_plots: bool,
    flat_uncertainty: bool = True,
    blow_up_frac_unc_threshold: float = 1.0,
) -> Dict[str, Any]:
    """Covariance + optional plots from precomputed offbeam / intime histograms."""
    h_off = np.asarray(h_off, dtype=float)
    h_in = np.asarray(h_in, dtype=float)
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

    pay = {
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
    if flat_uncertainty:
        apply_flat_cosmic_uncertainty(
            pay, blow_up_frac_unc_threshold=blow_up_frac_unc_threshold
        )
    return pay


def process_variable_cosmics(
    offbeam_df: pd.DataFrame,
    intime_df: pd.DataFrame,
    var_config: Any,
    cv_mode: str,
    save_fig_dir: str,
    save_plots: bool,
    flat_uncertainty: bool = True,
    blow_up_frac_unc_threshold: float = 1.0,
) -> Dict[str, Any]:
    h_off, h_in = cosmic_histograms(offbeam_df, intime_df, var_config)
    return process_variable_cosmics_from_histograms(
        h_off,
        h_in,
        var_config,
        cv_mode,
        save_fig_dir,
        save_plots,
        flat_uncertainty=flat_uncertainty,
        blow_up_frac_unc_threshold=blow_up_frac_unc_threshold,
    )


def save_cosmics_npz(syst_dict_by_var: Mapping[str, Any], out_path: str) -> None:
    np.savez_compressed(out_path, **syst_dict_by_var)
    logging.info("Wrote %s", out_path)


def load_cosmic_samples_ana() -> Tuple[pd.DataFrame, pd.DataFrame]:
    dfs = get_ana_dfs(option="cosmics_systs")
    return dfs["data"], dfs["mc"]


def _parse_syst_disk_root(p: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    root = args.syst_disk_root or os.environ.get(SYST_DISK_ENV)
    if not root:
        p.error(
            "Pass --syst-disk-root or set %s (outputs always go to <root>/%s/)."
            % (SYST_DISK_ENV, SUB_COSMICS)
        )
    args.syst_disk_root = path.abspath(path.expanduser(root.rstrip("/")))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("run-ana", help="Single get_ana_dfs('cosmics_systs') load → Cosmics NPZ")
    pr.add_argument(
        "--out-tag",
        default=None,
        help="Optional label for logs only (default: today).",
    )
    pr.add_argument("--syst-disk-root", default=None)
    pr.add_argument("--error-log", default=None)
    pr.add_argument("--no-plots", action="store_true")
    pr.add_argument("--no-save-npz", action="store_true")
    pr.add_argument(
        "--cv-mode",
        choices=("mean", "intime", "offbeam"),
        default="offbeam",
    )
    pr.add_argument("--vars", nargs="*", default=None)

    pa = sub.add_parser(
        "aggregate",
        help="Merge chunk pickles (same as syst_cosmics_aggregate.py).",
    )
    pa.add_argument("--chunks-dir", required=True, dest="chunks_dir")
    pa.add_argument("--syst-disk-root", default=None)
    pa.add_argument("--out-tag", default=None)
    pa.add_argument("--error-log", default=None)
    pa.add_argument("--no-plots", action="store_true")
    pa.add_argument("--no-save-npz", action="store_true")
    pa.add_argument("--cv-mode", choices=("mean", "intime", "offbeam"), default="offbeam")
    pa.add_argument("--vars", nargs="*", default=None)

    args = p.parse_args()
    _parse_syst_disk_root(p, args)
    return args


def main_run_ana(args: argparse.Namespace) -> None:
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

    offbeam_df, intime_df = load_cosmic_samples_ana()
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
        save_cosmics_npz(syst_dict, path.join(save_fig_dir, "cosmics_syst_dict.npz"))

    logger.info("Done -> %s (log %s)", save_fig_dir, log_path)


def main() -> None:
    args = parse_args()
    if args.cmd == "run-ana":
        main_run_ana(args)
    elif args.cmd == "aggregate":
        from analysis_village.numucc_1p0pi.scripts.syst_cosmics_aggregate import run_aggregate

        run_aggregate(args)
    else:
        raise SystemExit("unknown cmd %r" % (args.cmd,))


if __name__ == "__main__":
    _av = sys.argv[1:]
    if not _av:
        sys.argv.append("run-ana")
    elif _av[0] not in ("run-ana", "aggregate") and _av[0].startswith("-"):
        sys.argv.insert(1, "run-ana")
    main()
