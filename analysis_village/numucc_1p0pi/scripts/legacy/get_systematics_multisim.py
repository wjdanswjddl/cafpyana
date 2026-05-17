#!/usr/bin/env python
"""Unified multisim systematic extractor for MCstat, Flux, and G4 weights.

Supersedes the removed standalone ``get_systematics-flux.py`` and
``get_systematics-G4.py`` scripts (obsolete monolithic flux/G4 drivers).

For MCstat/Flux/G4 **without** cosmics, prefer this script over legacy helpers.
Cosmic unisim remains in ``get_systematics_cosmics.py``. Chunked neutrino
multisim uses ``syst_multisim_chunk.py`` → ``syst_multisim_aggregate.py``;
``get_systematics_mcstat_flux_g4.py`` remains a thin orchestrator for that path
and legacy single-metric monolithic runs.

For each kinematic variable and systematic source this script:

  1. Builds universe histograms for **total selected event rate** (signal + all
     backgrounds in-bin, no subtraction) and **background-subtracted signal rate**
     using :func:`analysis_village.numucc_1p0pi.utils.get_univ_rates`. Background
     fluctuations enter the subtracted metric via the existing implementation:
     each universe adds ``(background_univ - background_cv)`` bin-wise for every
     non-signal topology (see ``utils.get_univ_rates``).

  2. Computes covariance, fractional covariance, and correlation matrices per
     metric via ``pyanalib.covariance.get_covariance_matrix``.

  3. Plots universe distributions with :func:`utils.plot_univ_hists` (viridis CL
     styling) and heatmaps with :func:`utils.plot_heatmap`.

  4. Writes NPZ files keyed by ``var_save_name`` (compatible with
     ``event_selection_aggregate.py`` falling back to ``utils.get_syst_unc`` when
     chunks lack universe histograms). Each systematic entry exposes top-level
     ``cov``, ``cov_frac``, and ``corr`` for the **background-subtracted** metric
     (notebook / aggregate default), plus ``total_rate`` / ``bkgd_subtracted``
     dicts holding full matrices for both metrics.

Fault tolerance: failures for individual variables are logged and skipped without
aborting the full run.

Examples
--------
Bundled Flux/G4 + MCstat (same weight columns as ``files_config.get_ana_dfs``)::

    python get_systematics_multisim.py --out-tag 20260510

Per flux knob from ``makedf.bnbsyst.bnb_systematics_beam``::

    python get_systematics_multisim.py --flux-mode knobs --flux-knob-group beam

Cosmic background systematics (offbeam vs intime) are handled separately by
``get_systematics_cosmics.py``.

"""
from __future__ import annotations

import argparse
import logging
import sys
import traceback
from datetime import datetime
from os import makedirs, path
from typing import Any, Dict, List, Mapping, MutableMapping, Sequence, Tuple, Union

import matplotlib

matplotlib.use("Agg")
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
from pyanalib.split_df_helpers import load_and_concat_mc_dfs  # noqa: E402
from makedf.bnbsyst import (  # noqa: E402
    BNB_FLUX_GROUPS,
    bnb_systematics_beam,
    bnb_systematics_hadron,
    bnb_systematics_xsec,
)
from analysis_village.numucc_1p0pi.files_config import (  # noqa: E402
    file_dir as FILES_DEFAULT_FILE_DIR,
    get_ana_dfs,
    n_max_concat as FILES_DEFAULT_N_MAX_CONCAT,
    save_fig_base_dir,
)
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    generate_tags,
    get_univ_rates,
    plot_heatmap,
    plot_univ_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (  # noqa: E402
    FILE_FLUX,
    FILE_G4,
    FILE_MCSTAT,
    SUB_FLUX,
    SUB_G4,
    SUB_MCSTAT,
    normalized_root,
)

try:
    import matplotlib.pyplot as plt

    plt.style.use(path.join(path.dirname(__file__), "presentation.mplstyle"))
except Exception:
    pass


SystKey = Union[str, Tuple[str, ...]]

G4_SYSTEMATICS_DEFAULT = [
    "reinteractions_neutron_Geant4",
    "reinteractions_piminus_Geant4",
    "reinteractions_piplus_Geant4",
    "reinteractions_proton_Geant4",
]


def _flux_group_specs(group: str) -> Tuple[List[str], int | None]:
    """Return (knob names, multisim_nuniv hint) for a named flux bundle."""
    g = group.strip().lower()
    if g in BNB_FLUX_GROUPS:
        knobs, nuniv = BNB_FLUX_GROUPS[g]
        return list(knobs), int(nuniv)
    if g == "beam":
        return list(bnb_systematics_beam), None
    if g == "hadron":
        return list(bnb_systematics_hadron), None
    if g == "xsec":
        return list(bnb_systematics_xsec), None
    raise ValueError(f"Unknown flux knob group: {group}")


def _mcstat_univ_block(mc_evt_df: pd.DataFrame):
    """MCstat weights live under ``mc.MCstat`` after ``make_pandora_evtdf`` / ``truth_match``."""
    try:
        return mc_evt_df["mc"]["MCstat"]
    except (KeyError, TypeError, AttributeError):
        return mc_evt_df["MCstat"]


def _univ_rates_syst_key(mc_evt_df: pd.DataFrame, syst_name: SystKey) -> SystKey:
    """Key passed to :func:`utils.get_univ_rates` for column indexing."""
    if syst_name == "MCstat":
        try:
            mc_evt_df["mc"]["MCstat"]
            return ("mc", "MCstat")
        except (KeyError, TypeError, AttributeError):
            return "MCstat"
    return syst_name


def infer_n_univ(mc_evt_df: pd.DataFrame, syst_name: SystKey) -> int:
    """Infer multisim count from ``univ_*`` columns under ``syst_name``."""
    if syst_name == "MCstat":
        block = _mcstat_univ_block(mc_evt_df)
    else:
        block = mc_evt_df[syst_name]
    max_i = -1
    for c in block.columns:
        leaf = c[-1] if isinstance(c, tuple) else c
        s = str(leaf)
        if s.startswith("univ_"):
            try:
                max_i = max(max_i, int(s.split("_", 1)[1]))
            except ValueError:
                continue
    if max_i < 0:
        raise ValueError(f"No univ_* columns found under syst_name={syst_name!r}")
    return max_i + 1


def drop_events_bad_g4_weights(mc_evt_df: pd.DataFrame, thresh: float = 1e3) -> pd.DataFrame:
    """Drop rows with any G4 universe weight above ``thresh`` (legacy cleanup)."""
    # Prefer consolidated mc.G4 block when present
    try:
        g4 = mc_evt_df["mc"]["G4"]
    except (KeyError, TypeError):
        return mc_evt_df
    bad_idx = set()
    for c in g4.columns:
        leaf = c[-1] if isinstance(c, tuple) else c
        if not str(leaf).startswith("univ_"):
            continue
        s = g4[c]
        bad_idx.update(s.index[s > thresh].tolist())
    if not bad_idx:
        return mc_evt_df
    logging.warning("Dropping %d events with mc.G4 univ weight > %g", len(bad_idx), thresh)
    return mc_evt_df.drop(index=list(bad_idx))


def _sanitize_matrix_pack(pack: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Replace NaN/inf in fractional covariance; rebuild correlation safely."""
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
    """Resolve CLI ``var_save_name`` strings to ``VariableConfig`` instances."""
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
    # Default: core distributions + vertex/opening_angle + merged finals (φ, μ endpoints, …)
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


def _plot_matrices(
    ret: Mapping[str, np.ndarray],
    var_config: Any,
    metric_tag: str,
    syst_tag: str,
    save_fig_dir: str,
    save_plots: bool,
) -> None:
    for matrix_type in ("cov", "cov_frac", "corr"):
        title = {"cov": "Covariance", "cov_frac": "Fractional covariance", "corr": "Correlation"}[matrix_type]
        save_fig_name = path.join(
            save_fig_dir,
            f"{var_config.var_save_name}-{syst_tag}-{metric_tag}-{matrix_type}",
        )
        plot_heatmap(
            ret[matrix_type],
            var_config.bins,
            plot_labels=[var_config.var_labels[1], var_config.var_labels[1], title],
            plot=False,
            save_fig=save_plots,
            save_name=save_fig_name,
        )


def process_variable_two_metrics(
    mc_evt_df: pd.DataFrame,
    var_config: Any,
    syst_name: SystKey,
    syst_tag: str,
    save_fig_dir: str,
    save_plots: bool,
) -> Dict[str, Any]:
    """Return payload with total-rate and background-subtracted covariance bundles."""
    n_univ = infer_n_univ(mc_evt_df, syst_name)
    rates_sk = _univ_rates_syst_key(mc_evt_df, syst_name)

    univ_total, cv_total = get_univ_rates(
        cov_type="rate",
        evtdf=mc_evt_df,
        nudf=None,
        var_config=var_config,
        syst_name=rates_sk,
        n_univ=n_univ,
        bkgd_subtract=False,
        plot=False,
    )
    univ_sub, cv_sub = get_univ_rates(
        cov_type="rate",
        evtdf=mc_evt_df,
        nudf=None,
        var_config=var_config,
        syst_name=rates_sk,
        n_univ=n_univ,
        bkgd_subtract=True,
        plot=False,
    )

    # Universe distribution plots (utils.plot_univ_hists CL colors)
    if save_plots:
        base = path.join(save_fig_dir, f"{var_config.var_save_name}-{syst_tag}")
        plot_univ_hists(
            univ_total,
            cv_total,
            syst_name,
            var_config,
            plot=False,
            ax_titles=["", "Events / Bin", f"{syst_tag}: total event rate"],
            save_fig=True,
            save_name=f"{base}_total-rate_univ",
        )
        plot_univ_hists(
            univ_sub,
            cv_sub,
            syst_name,
            var_config,
            plot=False,
            ax_titles=["", "Events / Bin", f"{syst_tag}: background-subtracted signal"],
            save_fig=True,
            save_name=f"{base}_bkgd-sub_univ",
        )

    ret_total = _sanitize_matrix_pack(get_covariance_matrix(univ_total, cv_total))
    ret_sub = _sanitize_matrix_pack(get_covariance_matrix(univ_sub, cv_sub))

    if save_plots:
        _plot_matrices(ret_total, var_config, "total_rate", syst_tag, save_fig_dir, True)
        _plot_matrices(ret_sub, var_config, "bkgd_subtracted", syst_tag, save_fig_dir, True)

    # Notebook / aggregate default: background-subtracted fractional covariance
    payload: Dict[str, Any] = {
        "cov": ret_sub["cov"],
        "cov_frac": ret_sub["cov_frac"],
        "corr": ret_sub["corr"],
        # Explicit metric blocks (also include copy under legacy 'rate' nesting)
        "rate": ret_sub,
        "total_rate": ret_total,
        "bkgd_subtracted": ret_sub,
    }
    return payload


def _merge_nested_syst_dict(dst: MutableMapping[str, Any], src: Mapping[str, Any]) -> None:
    for var_key, inner in src.items():
        if var_key not in dst:
            dst[var_key] = {}
        dst[var_key].update(inner)


def save_per_category_npz(syst_dict_by_var: Mapping[str, Any], out_path: str) -> None:
    """Write ``var_save_name``-keyed NPZ (same layout as legacy flux/G4 scripts)."""
    np.savez_compressed(out_path, **syst_dict_by_var)
    logging.info("Wrote %s", out_path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-tag", default=None, help="Date/tag subdirectory under plots base (default: today %%Y%%m%%d)")
    p.add_argument("--out-dir", default=None, help="Override output directory for plots + NPZ (default: save_fig_base_dir/systematics-multisim-<tag>)")
    p.add_argument("--error-log", default=None, help="Path for per-variable failure log")
    p.add_argument("--no-plots", action="store_true", help="Skip figure output; still write NPZ")
    p.add_argument("--no-save-npz", action="store_true", help="Plots only")
    p.add_argument("--data-source", choices=("ana", "concat"), default="ana",
                   help="Load MC via get_ana_dfs('systs') or load_and_concat_mc_dfs")
    p.add_argument("--file-dir", default=FILES_DEFAULT_FILE_DIR)
    p.add_argument("--sample-dir", default="BNB_cosmics", help="Under MC/ when using --data-source concat")
    p.add_argument("--flux-subdir", default="", help="Optional trailing subdir under sample-dir (e.g. flux_wgts-beam)")
    p.add_argument("--chunk-tags", default=None, help="Comma-separated chunk tags; default from generate_tags('bl')[-5:] for concat")
    p.add_argument("--n-max-concat", type=int, default=FILES_DEFAULT_N_MAX_CONCAT)
    p.add_argument("--drop-g4-outliers", action="store_true", help="Drop events with any mc.G4 univ weight > 1e3")
    p.add_argument("--flux-mode", choices=("bundled", "knobs"), default="bundled",
                   help="Use consolidated mc.Flux / mc.G4 columns or individual flux knobs")
    p.add_argument("--flux-knob-group", default="beam",
                   help="When --flux-mode knobs: beam | hadron | xsec | beam,hadron,... composite")
    p.add_argument("--g4-mode", choices=("bundled", "knobs"), default="bundled")
    p.add_argument("--vars", nargs="*", default=None, help="Subset of var_save_name keys (see --help in code registry)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    tag = args.out_tag or datetime.now().strftime("%Y%m%d")
    save_fig_dir = args.out_dir or path.join(save_fig_base_dir, f"systematics-multisim-{tag}")
    makedirs(save_fig_dir, exist_ok=True)

    log_path = args.error_log or path.join(save_fig_dir, "multisim_failures.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    logger = logging.getLogger("multisim")

    var_configs = build_variable_configs(args.vars)
    save_plots = not args.no_plots

    # ----- load MC -----
    if args.data_source == "ana":
        ret = get_ana_dfs(option="systs")
        mc_evt_df = ret["evt"]
        mc_hdr_df = ret["hdr"]
    else:
        chunk_tags = (
            [t.strip() for t in args.chunk_tags.split(",") if t.strip()]
            if args.chunk_tags
            else generate_tags("bl")[-5:]
        )
        sample_dir = args.sample_dir.strip("/")
        if args.flux_subdir:
            sample_dir = path.join(sample_dir, args.flux_subdir.strip("/"))
        concat_dfs = load_and_concat_mc_dfs(
            file_dir=args.file_dir,
            chunk_tags=chunk_tags,
            df_tag="",
            keys2load=["hdr", "evt"],
            n_max_concat=args.n_max_concat,
            sub_dir="MC",
            sample_dir=sample_dir,
        )
        mc_hdr_df = concat_dfs["hdr"]
        mc_evt_df = concat_dfs["evt"]

    mc_evt_df = mc_evt_df.copy()
    mc_evt_df["pot_weight"] = np.ones(len(mc_evt_df), dtype=float)

    if args.drop_g4_outliers:
        mc_evt_df = drop_events_bad_g4_weights(mc_evt_df)

    logger.info("Loaded MC evt rows=%d POT(sum hdr)=%.4e", len(mc_evt_df), float(mc_hdr_df["pot"].sum()))

    # Accumulators: separate NPZ per category for aggregate/get_syst_unc-style consumption
    dict_mcstat: Dict[str, Any] = {}
    dict_flux: Dict[str, Any] = {}
    dict_g4: Dict[str, Any] = {}

    # (columns path, plot filename tag, key inside per-variable dict for NPZ)
    flux_jobs: List[Tuple[SystKey, str, str]] = []
    if args.flux_mode == "bundled":
        flux_jobs.append((("mc", "Flux"), "Flux", "flux"))
    else:
        groups = [g.strip() for g in args.flux_knob_group.split(",") if g.strip()]
        knobs_flat: List[str] = []
        for g in groups:
            klist, _ = _flux_group_specs(g)
            knobs_flat.extend(klist)
        for knob in knobs_flat:
            flux_jobs.append((("mc", knob), knob, knob))

    g4_jobs: List[Tuple[SystKey, str, str]] = []
    if args.g4_mode == "bundled":
        g4_jobs.append((("mc", "G4"), "G4", "G4"))
    else:
        for knob in G4_SYSTEMATICS_DEFAULT:
            g4_jobs.append((("mc", knob), knob, knob))

    def run_job(
        syst_name: SystKey,
        plot_tag: str,
        npz_inner_key: str,
        target: MutableMapping[str, Any],
    ) -> None:
        for var_config in tqdm(var_configs, desc=f"{plot_tag}"):
            try:
                payload = process_variable_two_metrics(
                    mc_evt_df,
                    var_config,
                    syst_name,
                    plot_tag,
                    save_fig_dir,
                    save_plots,
                )
                target.setdefault(var_config.var_save_name, {})
                target[var_config.var_save_name][npz_inner_key] = payload
            except Exception:
                logger.error(
                    "FAILED variable=%s syst_plot=%s npz_key=%s\n%s",
                    var_config.var_save_name,
                    plot_tag,
                    npz_inner_key,
                    traceback.format_exc(),
                )

    try:
        infer_n_univ(mc_evt_df, "MCstat")
        run_job("MCstat", "MCstat", "MCstat", dict_mcstat)
    except Exception:
        logger.error("MCstat columns missing or unreadable; skipping MCstat.\n%s", traceback.format_exc())

    for sk, ptag, nk in flux_jobs:
        try:
            infer_n_univ(mc_evt_df, sk)
            run_job(sk, ptag, nk, dict_flux)
        except Exception:
            logger.error("Flux systematic unavailable %s\n%s", ptag, traceback.format_exc())

    for sk, ptag, nk in g4_jobs:
        try:
            infer_n_univ(mc_evt_df, sk)
            run_job(sk, ptag, nk, dict_g4)
        except Exception:
            logger.error("G4 systematic unavailable %s\n%s", ptag, traceback.format_exc())

    if not args.no_save_npz:
        root = normalized_root(save_fig_dir)
        if dict_mcstat:
            d = path.join(root, SUB_MCSTAT)
            makedirs(d, exist_ok=True)
            save_per_category_npz(dict_mcstat, path.join(d, FILE_MCSTAT))
        if dict_flux:
            d = path.join(root, SUB_FLUX)
            makedirs(d, exist_ok=True)
            save_per_category_npz(dict_flux, path.join(d, FILE_FLUX))
        if dict_g4:
            d = path.join(root, SUB_G4)
            makedirs(d, exist_ok=True)
            save_per_category_npz(dict_g4, path.join(d, FILE_G4))

        # Combined bundle for convenience (same per-var nesting)
        combined: Dict[str, Any] = {}
        _merge_nested_syst_dict(combined, dict_mcstat)
        _merge_nested_syst_dict(combined, dict_flux)
        _merge_nested_syst_dict(combined, dict_g4)
        if combined:
            save_per_category_npz(combined, path.join(root, "multisim_syst_dict_all.npz"))

    logger.info("Done. Outputs under %s (failures logged to %s)", save_fig_dir, log_path)


if __name__ == "__main__":
    main()
