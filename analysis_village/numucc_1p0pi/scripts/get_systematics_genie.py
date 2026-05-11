#!/usr/bin/env python
"""GENIE multisim systematics: event-rate and cross-section covariances (chunked only).

Input ``.df`` paths and per-group directory globs live in
``analysis_village.numucc_1p0pi.dataset_locations`` (``GENIE_GROUP_GLOBS`` /
:func:`iter_genie_chunk_map_tasks`). Each **knob group** (CCQE, MEC, …) uses its own sample
directory; within a group, all knobs share the same files.

**Map** (HDF splits inside each ``.df`` file):

* ``chunk-map`` — sequential reads of ``evt_{i}`` / ``mcnu_{i}``; accumulates *additive*
  histograms for the **rate** path and xsec tensors (see module docstring in git history
  for the merge math).

* ``chunk-merge`` — sums ``genie__<GROUP>__*.pkl`` from ``chunk-map`` and builds covariances.

Example::

    python get_systematics_genie.py chunk-map --df-file in.df --out-dir chunks_genie \\
        --genie-group CCQE --max-splits 0

    python get_systematics_genie.py chunk-merge --chunks-dir chunks_genie \\
        --genie-group CCQE --out-dir plots_genie --xsec-unit 1.0

Batch inputs::

    python3 -c "from analysis_village.numucc_1p0pi.dataset_locations import iter_genie_chunk_map_tasks
    for g,p in iter_genie_chunk_map_tasks(): print(g,p)"
"""
from __future__ import annotations

import argparse
import gc
import glob
import logging
import os
import pickle
import sys
from os import path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x, **kwargs):
        return x

import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

sys.path.append(path.join(path.dirname(__file__), "..", "..", ".."))

from pyanalib.covariance import get_covariance_matrix  # noqa: E402
from pyanalib.split_df_helpers import get_n_split  # noqa: E402

from analysis_village.numucc_1p0pi.categories import (  # noqa: E402
    get_genie_category,
    get_topo_category,
    topology_list,
)
from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_ORDER, GENIE_GROUP_KNOBS  # noqa: E402
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)
from analysis_village.numucc_1p0pi.utils import (  # noqa: E402
    get_clipped_evts,
    get_response_matrix,
    get_univ_rates,
    plot_heatmap,
    plot_univ_hists,
    signal_hists,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402
logger = logging.getLogger(__name__)

SystName = Tuple[str, str]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_genie_dataframes(
    dfs: Mapping[str, Any],
    *,
    context: str = "",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(evt, mcnu)`` after verifying both frames exist."""
    missing = [k for k in ("evt", "mcnu") if k not in dfs]
    if missing:
        raise ValueError(
            f"{context} GENIE setup requires both 'evt' and 'mcnu' dataframes; "
            f"missing keys: {missing}. Present: {sorted(dfs.keys())}"
        )
    evt = dfs["evt"]
    mcnu = dfs["mcnu"]
    if not isinstance(evt, pd.DataFrame) or not isinstance(mcnu, pd.DataFrame):
        raise TypeError(f"{context} 'evt' and 'mcnu' must be pandas DataFrames.")
    if len(evt) == 0 and len(mcnu) == 0:
        logger.warning("%s empty evt and mcnu frames", context)
    return evt, mcnu


def validate_split_pair(evt_df: pd.DataFrame, mcnu_df: pd.DataFrame, split_idx: int) -> None:
    if len(evt_df) == 0 and len(mcnu_df) == 0:
        logger.debug("split %d: both evt and mcnu empty", split_idx)


def _attach_phi_degrees(mc_evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame) -> None:
    """Fill mu/p track phi in degrees (same recipe as the legacy monolithic driver)."""
    for df in (mc_evt_df, mc_nu_df):
        if df is None or len(df) == 0:
            continue
        for pref in ("mu", "p"):
            try:
                df.loc[:, (pref, "pfp", "trk", "phi", "", "", "")] = np.degrees(
                    np.arctan2(
                        df[pref, "pfp", "trk", "dir", "x", "", ""],
                        df[pref, "pfp", "trk", "dir", "y", "", ""],
                    )
                )
            except (KeyError, TypeError, ValueError):
                pass


# ---------------------------------------------------------------------------
# Universe helpers
# ---------------------------------------------------------------------------


def _iter_leaf_strings(block_cols: pd.Index) -> List[str]:
    out: List[str] = []
    for c in block_cols:
        leaf = c[0] if isinstance(c, tuple) else c
        out.append(str(leaf))
    return out


def _infer_multisim_n_univ(block_cols: pd.Index) -> int:
    max_i = -1
    for leaf_s in _iter_leaf_strings(block_cols):
        if not leaf_s.startswith("univ_"):
            continue
        try:
            max_i = max(max_i, int(leaf_s.split("_", 1)[1]))
        except ValueError:
            continue
    return max_i + 1


def _ensure_univ0_from_leaf(df: pd.DataFrame, syst_key: Tuple[str, ...], src_leaf: str) -> None:
    """
    For non-multisim knobs, alias one leaf (e.g. 'ps1' or 'morph') into a synthetic 'univ_0'
    column under the same systematic block so downstream code can treat it like multisim.
    """
    block = df.loc[:, syst_key]
    src_rest = None
    for c in block.columns:
        leaf = c[-1] if isinstance(c, tuple) else c
        if str(leaf) == src_leaf:
            src_rest = c if isinstance(c, tuple) else (c,)
            break
    if src_rest is None:
        # Nothing to do; caller will raise a clearer error.
        return

    dst_rest = tuple(src_rest[:-1]) + ("univ_0",)
    full_src = tuple(syst_key) + tuple(src_rest)
    full_dst = tuple(syst_key) + tuple(dst_rest)

    if full_dst in df.columns:
        return
    df.loc[:, full_dst] = df.loc[:, full_src]


def normalize_and_infer_n_univ(mc_evt_df: pd.DataFrame, mc_nu_df: pd.DataFrame, syst_name: SystName) -> int:
    """
    Detect whether a knob is multisim vs multisigma vs unisim and normalize to a common
    interface: multisim keeps 'univ_*', while multisigma/unisim get a synthetic 'univ_0'.

    Rules (per user request):
    - multisigma: use 'ps1' as the one-universe unisim
    - unisim: use 'morph' as the one-universe unisim
    """
    key = tuple(syst_name)
    block_cols = mc_evt_df.loc[:, key].columns
    n = _infer_multisim_n_univ(block_cols)
    if n > 0:
        return n

    leaves = set(_iter_leaf_strings(block_cols))
    if "ps1" in leaves:
        _ensure_univ0_from_leaf(mc_evt_df, key, "ps1")
        _ensure_univ0_from_leaf(mc_nu_df, key, "ps1")
        return 1
    if "morph" in leaves:
        _ensure_univ0_from_leaf(mc_evt_df, key, "morph")
        _ensure_univ0_from_leaf(mc_nu_df, key, "morph")
        return 1

    raise ValueError(
        f"No univ_* columns under syst_name={syst_name!r}, and also no 'ps1' (multisigma) "
        f"or 'morph' (unisim) leaf found. Available leaves: {sorted(leaves)}"
    )


def _syst_plot_key(syst_name: SystName) -> SystName:
    return syst_name


# ---------------------------------------------------------------------------
# XSEC-path accumulation (chunk-additive)
# ---------------------------------------------------------------------------


def _empty_xsec_tensor_acc(n_univ: int, nb: int) -> Dict[str, np.ndarray]:
    return {
        "nevts_allmc": np.zeros(nb, dtype=np.float64),
        "cv_sel_reco": np.zeros(nb, dtype=np.float64),
        "cv_allsel_reco": np.zeros(nb, dtype=np.float64),
        "reco_vs_true": np.zeros((n_univ, nb, nb), dtype=np.float64),
        "signal_allmc": np.zeros((n_univ, nb), dtype=np.float64),
        "signal_sel_truth": np.zeros((n_univ, nb), dtype=np.float64),
        "bg_cv": np.zeros(nb, dtype=np.float64),
        "bg_univ": np.zeros((n_univ, nb), dtype=np.float64),
    }


def accumulate_xsec_path_chunk(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    n_univ: int,
    acc: MutableMapping[str, np.ndarray],
) -> None:
    """Add contributions from one evt/mcnu chunk into ``acc`` (in-place)."""
    bins = var_config.bins
    nb = len(bins) - 1

    evtdf_signal = mc_evt_df[mc_evt_df.topo_categ == 1]
    nudf_signal = mc_nu_df[mc_nu_df.topo_categ == 1]
    evtdf_div_topo = [mc_evt_df[mc_evt_df.topo_categ == mode] for mode in topology_list]

    ret = signal_hists(mc_evt_df, mc_nu_df, var_config, return_data=True, plot=False)
    nevts_allmc = ret["nevts_allmc"]
    if nevts_allmc is None:
        raise ValueError(
            f"xsec path requires mcnu for var={var_config.var_save_name}; nevts_allmc is None"
        )

    acc["nevts_allmc"] += np.asarray(nevts_allmc, dtype=np.float64)
    acc["cv_sel_reco"] += np.asarray(ret["nevts_sel_reco"], dtype=np.float64)
    acc["cv_allsel_reco"] += np.asarray(ret["nevts_allsel_reco"], dtype=np.float64)

    for uidx in range(n_univ):
        univ_col = f"univ_{uidx}"
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true, _, _ = np.histogram2d(
                ret["var_sel_truth"],
                ret["var_sel_reco"],
                weights=ret["wgt_sel_truth"] * evtdf_signal[syst_name][univ_col],
                bins=bins,
            )
        acc["reco_vs_true"][uidx] += reco_vs_true

        sam, _ = np.histogram(
            ret["var_allmc"],
            weights=ret["wgt_allmc"] * nudf_signal[syst_name][univ_col],
            bins=bins,
        )
        sst, _ = np.histogram(
            ret["var_sel_truth"],
            weights=ret["wgt_sel_truth"] * evtdf_signal[syst_name][univ_col],
            bins=bins,
        )
        acc["signal_allmc"][uidx] += sam
        acc["signal_sel_truth"][uidx] += sst

    for this_evtdf in evtdf_div_topo[1:]:
        var, wgt = get_clipped_evts(this_evtdf, var_config.var_evt_reco_col, bins)
        acc["bg_cv"] += np.histogram(var, bins=bins, weights=wgt)[0].astype(np.float64)
        for uidx in range(n_univ):
            uw = this_evtdf[syst_name][f"univ_{uidx}"].copy()
            uw[np.isnan(uw)] = 1.0
            acc["bg_univ"][uidx] += np.histogram(var, bins=bins, weights=wgt * uw)[0].astype(
                np.float64
            )


def finalize_xsec_univ_events(
    acc: Mapping[str, np.ndarray],
    xsec_unit: float,
) -> np.ndarray:
    """Return ``univ_events`` of shape ``(n_univ, nbins)`` for covariance (xsec path)."""
    nevts_allmc = np.asarray(acc["nevts_allmc"], dtype=np.float64)
    nb = int(nevts_allmc.shape[0])
    n_univ = int(acc["reco_vs_true"].shape[0])
    scale = float(xsec_unit)

    rows: List[np.ndarray] = []
    for uidx in range(n_univ):
        if nb == 1:
            reco_vs_true = np.array([[1.0]], dtype=np.float64)
        else:
            reco_vs_true = acc["reco_vs_true"][uidx]

        # Match utils.get_univ_rates xsec branch: ratio of summed weighted histograms.
        eff = acc["signal_sel_truth"][uidx] / acc["signal_allmc"][uidx]
        response = get_response_matrix(reco_vs_true, eff)
        signal_univ = response @ nevts_allmc
        signal_univ = signal_univ + (acc["bg_univ"][uidx] - acc["bg_cv"])
        signal_univ *= scale
        rows.append(signal_univ)
    return np.asarray(rows, dtype=np.float64)


def finalize_cv_sel_reco_xsec(
    acc: Mapping[str, np.ndarray], xsec_unit: float, *, bkgd_subtract: bool = True
) -> np.ndarray:
    scale = float(xsec_unit)
    base = acc["cv_sel_reco"] if bkgd_subtract else acc["cv_allsel_reco"]
    return np.asarray(base, dtype=np.float64) * scale


# ---------------------------------------------------------------------------
# Pickle blob layout for chunked GENIE
# ---------------------------------------------------------------------------

RATE_ACC_KEY = "rate_univ_cv"
XSEC_ACC_KEY = "xsec_accumulators"


def accumulate_chunk_into_blob_root(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    blob_root: MutableMapping[str, Any],
    var_configs: Sequence[VariableConfig],
    syst_names: Sequence[SystName],
    *,
    bkgd_subtract: bool = True,
) -> None:
    validate_split_pair(mc_evt_df, mc_nu_df, -1)
    rate_blk = blob_root.setdefault(RATE_ACC_KEY, {})

    for syst_name in syst_names:
        knob = syst_name[1]
        n_univ = normalize_and_infer_n_univ(mc_evt_df, mc_nu_df, syst_name)

        for var_config in var_configs:
            slug = var_config.var_save_name

            # ---- rate path (identical to summing get_univ_rates per chunk)
            univ_r, cv_r = get_univ_rates(
                cov_type="rate",
                syst_type="GENIE",
                evtdf=mc_evt_df,
                nudf=mc_nu_df,
                var_config=var_config,
                syst_name=syst_name,
                n_univ=n_univ,
                bkgd_subtract=bkgd_subtract,
                plot=False,
            )
            slot_r = rate_blk.setdefault(knob, {}).setdefault(
                slug,
                {"univ": np.zeros_like(univ_r, dtype=np.float64), "cv": np.zeros_like(cv_r, dtype=np.float64)},
            )
            slot_r["univ"] += np.asarray(univ_r, dtype=np.float64)
            slot_r["cv"] += np.asarray(cv_r, dtype=np.float64)

            # ---- xsec tensors
            xsec_blk = blob_root.setdefault(XSEC_ACC_KEY, {})
            knob_blk = xsec_blk.setdefault(knob, {})
            if slug not in knob_blk:
                nb = len(var_config.bins) - 1
                knob_blk[slug] = _empty_xsec_tensor_acc(n_univ, nb)
            accumulate_xsec_path_chunk(mc_evt_df, mc_nu_df, var_config, syst_name, n_univ, knob_blk[slug])


def merge_genie_chunk_pickles(paths: List[str]) -> Dict[str, Any]:
    merged: Optional[Dict[str, Any]] = None
    for fp in tqdm(paths, desc="merge GENIE chunks"):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        if merged is None:
            merged = {
                "kind": "genie_syst_merged",
                RATE_ACC_KEY: d[RATE_ACC_KEY],
                XSEC_ACC_KEY: d[XSEC_ACC_KEY],
                "meta": [d.get("meta", {})],
            }
            continue
        merged["meta"].append(d.get("meta", {}))

        # rate
        for knob, vars_d in d[RATE_ACC_KEY].items():
            for slug, pack in vars_d.items():
                mp = merged[RATE_ACC_KEY].setdefault(knob, {}).setdefault(
                    slug,
                    {"univ": pack["univ"].copy(), "cv": pack["cv"].copy()},
                )
                if mp["univ"].shape != pack["univ"].shape:
                    raise ValueError(f"rate shape mismatch {knob} {slug}")
                mp["univ"] += pack["univ"]
                mp["cv"] += pack["cv"]

        # xsec tensors
        for knob, vars_d in d[XSEC_ACC_KEY].items():
            for slug, acc_new in vars_d.items():
                mp_acc = merged[XSEC_ACC_KEY].setdefault(knob, {}).setdefault(slug, {})
                if not mp_acc:
                    for k, arr in acc_new.items():
                        mp_acc[k] = np.array(arr, dtype=np.float64, copy=True)
                    mp_acc.setdefault(
                        "cv_allsel_reco", np.zeros_like(mp_acc["cv_sel_reco"], dtype=np.float64)
                    )
                else:
                    for k in acc_new:
                        mp_acc[k] += np.asarray(acc_new[k], dtype=np.float64)
    if merged is None:
        raise RuntimeError("no GENIE chunk pickles merged")
    return merged


# ---------------------------------------------------------------------------
# Plotting / covariance packaging (matches legacy get_systematics)
# ---------------------------------------------------------------------------


def covariance_bundle_univ_events(
    univ_events: np.ndarray,
    cv_events: np.ndarray,
    syst_name: SystName,
    var_config: VariableConfig,
    cov_tag: str,
    *,
    plot: bool = False,
    save_fig: bool = False,
    save_fig_dir: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    ret = get_covariance_matrix(univ_events, cv_events)
    if save_fig and save_fig_dir:
        os.makedirs(save_fig_dir, exist_ok=True)
        sk = _syst_plot_key(syst_name)
        plot_univ_hists(
            univ_events,
            cv_events,
            sk,
            var_config,
            plot=plot,
            save_fig=True,
            save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-{syst_name[1]}_{cov_tag}-univ_hists"),
        )
        for matrix_type in ("cov", "cov_frac", "corr"):
            lab = matrix_type.replace("_", " ").title()
            plot_heatmap(
                ret[matrix_type],
                var_config.bins,
                plot_labels=[var_config.var_labels[1], var_config.var_labels[1], lab],
                plot=plot,
                save_fig=True,
                save_name=path.join(save_fig_dir, f"{var_config.var_save_name}-{syst_name[1]}_{cov_tag}-{matrix_type}"),
            )
    return ret


def get_systematics(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    var_config: VariableConfig,
    syst_name: SystName,
    syst_type: str = "GENIE",
    plot: bool = False,
    save_fig: bool = False,
    save_fig_dir: Optional[str] = None,
    *,
    xsec_unit: float = 0.0,
):
    """Legacy helper — delegates to :func:`utils.get_univ_rates` for ``rate`` and ``xsec``."""
    matrices: Dict[str, Any] = {}
    validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context="get_systematics: ")
    n_univ = normalize_and_infer_n_univ(mc_evt_df, mc_nu_df, syst_name)
    for cov_type in ["xsec", "rate"]:
        univ_events, cv_events = get_univ_rates(
            cov_type=cov_type,
            syst_type=syst_type,
            evtdf=mc_evt_df,
            nudf=mc_nu_df,
            var_config=var_config,
            syst_name=syst_name,
            n_univ=n_univ,
            xsec_unit=xsec_unit,
            plot=False,
        )
        ret = covariance_bundle_univ_events(
            univ_events,
            cv_events,
            syst_name,
            var_config,
            cov_type,
            plot=plot,
            save_fig=save_fig,
            save_fig_dir=save_fig_dir,
        )
        matrices[cov_type] = ret
    return matrices


# ---------------------------------------------------------------------------
# CLI chunk-map / chunk-merge
# ---------------------------------------------------------------------------


def run_chunk_map(
    df_file: str,
    out_dir: str,
    genie_group: str,
    var_configs: List[VariableConfig],
    syst_names: Sequence[SystName],
    n_univ_cap: int = 0,
    max_splits: int = 0,
    bkgd_subtract: bool = True,
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[genie-chunk-map] no HDF splits")

    blob_root: Dict[str, Any] = {
        "kind": "genie_syst_chunk",
        "meta": {
            "df_file": df_file,
            "splits_processed": n_use,
            "genie_group": genie_group,
        },
        RATE_ACC_KEY: {},
        XSEC_ACC_KEY: {},
    }

    for i in tqdm(range(n_use), desc="HDF splits"):
        mc_evt_df = pd.read_hdf(df_file, key=f"evt_{i}")
        mc_nu_df = pd.read_hdf(df_file, key=f"mcnu_{i}")
        validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context=f"split {i}:")
        mc_evt_df = mc_evt_df.copy()
        mc_nu_df = mc_nu_df.copy()

        # Match the df-update notebook: ensure mcnu columns have a leading "mc" level so
        # category helpers that expect df.mc.* work on both evt and mcnu frames.
        if isinstance(mc_nu_df.columns, pd.MultiIndex):
            try:
                first_level = mc_nu_df.columns.get_level_values(0)
                need_prefix = not np.all(first_level == "mc")
            except Exception:
                need_prefix = True
            if need_prefix:
                mc_nu_df.columns = pd.MultiIndex.from_tuples(
                    [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
                )

        # Event categories (required for selections/hists downstream).
        mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
        mc_nu_df.loc[:, "topo_categ"] = get_topo_category(mc_nu_df)

        # Keep parity with the notebook even if not used everywhere yet.
        mc_evt_df.loc[:, "genie_categ"] = get_genie_category(mc_evt_df)
        mc_nu_df.loc[:, "genie_categ"] = get_genie_category(mc_nu_df)

        _attach_phi_degrees(mc_evt_df, mc_nu_df)
        if n_univ_cap > 0:
            # Optionally truncate universe columns — uncommon; omitted unless weights exist
            pass

        accumulate_chunk_into_blob_root(
            mc_evt_df,
            mc_nu_df,
            blob_root,
            var_configs,
            syst_names,
            bkgd_subtract=bkgd_subtract,
        )
        del mc_evt_df, mc_nu_df
        gc.collect()

    stem = path.splitext(path.basename(df_file))[0]
    out_path = path.join(out_dir, "genie__%s__%s.pkl" % (genie_group, stem))
    with open(out_path, "wb") as f:
        pickle.dump(blob_root, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("wrote %s", out_path)
    return out_path


def run_chunk_merge(
    chunks_dir: str,
    out_dir: str,
    genie_group: str,
    var_configs: List[VariableConfig],
    xsec_unit: float,
    *,
    bkgd_subtract: bool = True,
    save_figs: bool = False,
    npz_path: Optional[str] = None,
) -> Dict[str, Any]:
    paths = sorted(glob.glob(path.join(chunks_dir, "genie__%s__*.pkl" % genie_group)))
    if not paths:
        raise SystemExit(
            "[genie-chunk-merge] no genie__%s__*.pkl under %s" % (genie_group, chunks_dir)
        )
    vc_by = {v.var_save_name: v for v in var_configs}

    merged = merge_genie_chunk_pickles(paths)
    os.makedirs(out_dir, exist_ok=True)

    syst_dict_out: Dict[str, Dict[str, Dict[str, Any]]] = {}

    for knob in tqdm(sorted(merged[RATE_ACC_KEY].keys()), desc="GENIE knobs"):
        syst_dict_out[knob] = {}
        syst_tuple: SystName = ("mc", knob)

        for slug, rate_pack in merged[RATE_ACC_KEY][knob].items():
            vc = vc_by.get(slug)
            if vc is None:
                continue
            univ_rate = np.asarray(rate_pack["univ"], dtype=np.float64)
            cv_rate = np.asarray(rate_pack["cv"], dtype=np.float64)
            rate_cov = covariance_bundle_univ_events(
                univ_rate,
                cv_rate,
                syst_tuple,
                vc,
                "rate",
                plot=False,
                save_fig=save_figs,
                save_fig_dir=out_dir if save_figs else None,
            )

            acc_x = merged[XSEC_ACC_KEY][knob][slug]
            univ_xsec = finalize_xsec_univ_events(acc_x, xsec_unit=xsec_unit)
            cv_xsec = finalize_cv_sel_reco_xsec(acc_x, xsec_unit=xsec_unit, bkgd_subtract=bkgd_subtract)
            xsec_cov = covariance_bundle_univ_events(
                univ_xsec,
                cv_xsec,
                syst_tuple,
                vc,
                "xsec",
                plot=False,
                save_fig=save_figs,
                save_fig_dir=out_dir if save_figs else None,
            )

            syst_dict_out[knob][slug] = {"rate": rate_cov, "xsec": xsec_cov}

    if npz_path:
        # Nested dict of covariance arrays: store as a 0-d object array (NumPy pickles contents).
        np.savez_compressed(npz_path, syst=np.array(syst_dict_out, dtype=object))
        logger.info("wrote %s (object array key 'syst')", npz_path)

    summ = path.join(out_dir, "genie_chunk_merge_summary_%s.txt" % genie_group)
    with open(summ, "w") as f:
        f.write("# merged chunk pickles\n")
        for p in paths:
            f.write("%s\n" % p)
    return syst_dict_out


def parse_chunk_cli(argv: Optional[Sequence[str]] = None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pm = sub.add_parser(
        "chunk-map",
        help="One .df file → genie__<GROUP>__stem.pkl accumulators for that group's knobs",
    )
    pm.add_argument("--df-file", required=True)
    pm.add_argument("--out-dir", required=True)
    pm.add_argument(
        "--genie-group",
        required=True,
        choices=list(GENIE_GROUP_ORDER),
        help="Knob group / sample layout (must match dataset_locations GENIE_GROUP_GLOBS).",
    )
    pm.add_argument("--max-splits", type=int, default=0, help="0 = all HDF splits")
    pm.add_argument(
        "--knobs",
        default=None,
        help="Comma-separated mc.* knob names (default: all knobs for this group from makedf.geniesyst)",
    )

    rg = sub.add_parser(
        "chunk-merge",
        help="Merge genie__<GROUP>__*.pkl for one group → covariance dict / NPZ",
    )
    rg.add_argument("--chunks-dir", required=True)
    rg.add_argument("--out-dir", required=True)
    rg.add_argument(
        "--genie-group",
        required=True,
        choices=list(GENIE_GROUP_ORDER),
        help="Must match the chunk-map --genie-group / pickle prefix.",
    )
    rg.add_argument("--xsec-unit", type=float, default=1.0)
    rg.add_argument("--save-figs", action="store_true")
    rg.add_argument(
        "--no-bkgd-subtract",
        action="store_true",
        help="Match get_univ_rates(..., bkgd_subtract=False) for rate+xsec CV baseline",
    )
    rg.add_argument("--npz", default=None)

    return p.parse_args(argv)


def main_cli_chunk(argv: Optional[Sequence[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_chunk_cli(argv)
    var_configs = with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))

    if args.cmd == "chunk-map":
        group = args.genie_group
        if args.knobs:
            knobs = [x.strip() for x in args.knobs.split(",") if x.strip()]
        else:
            knobs = list(GENIE_GROUP_KNOBS[group])
        if not knobs:
            raise SystemExit("[chunk-map] no knobs for group %s" % group)
        syst_names = [("mc", k) for k in knobs]
        run_chunk_map(
            args.df_file,
            args.out_dir,
            group,
            var_configs,
            syst_names,
            max_splits=args.max_splits,
        )

    elif args.cmd == "chunk-merge":
        probe_paths = sorted(glob.glob(path.join(args.chunks_dir, "genie__%s__*.pkl" % args.genie_group)))
        if probe_paths:
            with open(probe_paths[0], "rb") as f:
                probe = pickle.load(f)
            meta_g = probe.get("meta", {}).get("genie_group")
            if meta_g and meta_g != args.genie_group:
                logger.warning(
                    "chunk pickle meta genie_group=%r differs from CLI %r",
                    meta_g,
                    args.genie_group,
                )
        run_chunk_merge(
            args.chunks_dir,
            args.out_dir,
            args.genie_group,
            var_configs,
            xsec_unit=args.xsec_unit,
            bkgd_subtract=not args.no_bkgd_subtract,
            save_figs=args.save_figs,
            npz_path=args.npz,
        )


if __name__ == "__main__":
    main_cli_chunk()
