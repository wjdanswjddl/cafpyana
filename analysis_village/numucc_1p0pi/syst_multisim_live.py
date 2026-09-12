"""Live / interactive multisim systematics accumulation for notebooks.

Walks ``sel_mup`` (or similar final-selection) ``.df`` files that carry
per-knob universe weights, accumulates CV + universe rate histograms, and
refreshes a multi-panel figure (distributions, cov / corr, fractional
uncertainty) as statistics grow — same interactive pattern as
``event_selection_live.py``.

Accumulated histogram counts are written to a pickle so a later notebook
section can reload and recompute covariance products without re-reading CAFs.
"""
from __future__ import annotations

import gc
import io
import pickle
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from pyanalib.covariance import get_covariance_matrix
from pyanalib.split_df_helpers import get_n_split, load_dfs

from analysis_village.numucc_1p0pi.categories import get_topo_category
from analysis_village.numucc_1p0pi.evt_derived_kinematics import ensure_derived_trk_kinematics_cols
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    build_var_configs,
    combine_indep_knob_cov_packs,
    drop_bad_flux_knob_weights,
    drop_bad_g4_knob_weights,
    flux_mc_knob_names,
    g4_mc_knob_names,
    syst_key_for_name,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates

# Universe overlay alpha for live / light step histograms (user: "alpha=5" → 5%).
UNIV_ALPHA = 0.05
TOTAL_KNOB_KEY = "__total__"


@dataclass
class LiveMultisimConfig:
    """Notebook knobs for the live multisim accumulation run."""

    df_dir: Path | str = ""
    """Directory of ``*.df`` files (e.g. ``…/2026_09_03_022028__sel_mup-wgts_flux-corrected``)."""

    syst_name: str = "Flux"
    """``Flux``, ``G4``, or ``MCstat`` — selects knob list / weight layout."""

    flux_knob_groups: str = "all"
    """Passed to :func:`flux_mc_knob_names` when ``syst_name=='Flux'``."""

    n_universe: int = 100
    """Cap on universes read from each knob (files may store more)."""

    max_files: int | None = None
    """Cap on ``.df`` files. ``None`` = all."""

    files_per_batch: int = 200
    """Load + concat this many ``.df`` files via :func:`pyanalib.split_df_helpers.load_dfs`
    before one accumulate pass (amortizes ``get_univ_rates`` overhead)."""

    update_every_n_files: int = 200
    """Refresh the live panel whenever cumulative files processed reaches a multiple
    of this (also always on the last batch). Default matches ``files_per_batch``."""

    var_set: str = "final"
    """Variable registry passed to :func:`build_var_configs`."""

    plot_var_names: Sequence[str] | None = None
    """Subset of ``var_save_name`` shown on the live panel. ``None`` → core kinematics."""

    accumulate_var_names: Sequence[str] | None = None
    """Vars for which histogram counts are filled. ``None`` → same as ``plot_var_names``
    (keeps interactive runs tractable). Set to all final names for a full dump."""

    accumulate_total: bool = True
    """Always accumulate product-of-knobs ``__total__`` for Flux/G4 (default on)."""

    show_in_notebook: bool = True
    save_every_update: bool = True
    """Write ``accumulated_histcounts.pkl`` (+ analysis packs) on each panel refresh."""

    save_panel_png: bool = True
    """Write ``live_panel_latest.png`` (and a per-update snapshot) under ``plots_dir``."""

    work_dir: Path | str | None = None
    plots_dir: Path | str | None = None
    figsize_per_row: Tuple[float, float] = (16.0, 3.6)
    panel_dpi: int = 110
    univ_alpha: float = UNIV_ALPHA
    bkgd_subtract: bool = True
    trace: bool = False


@dataclass
class LiveMultisimResult:
    work_dir: Path
    plots_dir: Path
    payload_path: Path
    n_files_done: int
    n_files_total: int
    knobs: Tuple[str, ...]
    var_configs: List[Any]
    acc: Dict[str, Any]
    failed: List[Tuple[str, str]]


def default_work_dir(syst_name: str, tag: str | None = None) -> Path:
    today = datetime.now().strftime("%Y%m%d")
    suffix = tag or today
    return Path(
        f"/exp/sbnd/data/users/munjung/plots/numucc1p0pi/"
        f"systematics-multisim-live-{syst_name}-{suffix}"
    )


def default_plot_var_names() -> Tuple[str, ...]:
    """Compact live-panel set (xsec measurement variables)."""
    return (
        "integrated",
        "muon-p",
        "muon-dir_z",
        "proton-p",
        "proton-dir_z",
        "tki-del_Tp",
        "tki-del_p",
        "tki-del_alpha",
        "tki-del_phi",
    )


def discover_df_files(df_dir: Path | str, max_files: int | None = None) -> List[Path]:
    root = Path(df_dir).expanduser()
    files = sorted(root.glob("*.df"))
    if max_files is not None:
        files = files[: int(max_files)]
    return files


def knobs_for_syst(syst_name: str, flux_knob_groups: str = "all") -> Tuple[str, ...]:
    if syst_name == "Flux":
        return flux_mc_knob_names(flux_knob_groups)
    if syst_name == "G4":
        return g4_mc_knob_names()
    if syst_name == "MCstat":
        return ()  # bundled ``mc.MCstat`` / ``MCstat`` block — no per-knob loop
    raise ValueError(f"unsupported syst_name={syst_name!r}")


def _count_univ_columns(evt_df: pd.DataFrame, syst_col_key: Any, cap: int = 2048) -> int:
    n = 0
    for i in range(cap):
        if isinstance(syst_col_key, tuple):
            probe = syst_col_key + (f"univ_{i}",)
        else:
            probe = (syst_col_key, f"univ_{i}")
        if multicol_resolve_column_key(evt_df, probe) is None:
            break
        n += 1
    return n


def _sanitize_matrix_pack(ret: dict) -> dict:
    out = {}
    for k in ("cov", "cov_frac", "corr"):
        out[k] = np.nan_to_num(np.asarray(ret[k], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    return out


def frac_unc_from_pack(pack: dict) -> np.ndarray:
    return np.sqrt(np.maximum(np.diag(pack["cov_frac"]), 0.0))


def inject_multiplied_mc_knob_weights(
    df: pd.DataFrame,
    knobs: Sequence[str],
    bundled_tag: str,
    n_univ: int,
) -> pd.DataFrame:
    """Write ``(mc, bundled_tag, univ_i)`` as the product of per-knob weights."""
    df = df.copy()
    for uidx in range(int(n_univ)):
        w = np.ones(len(df), dtype=float)
        missing = False
        for knob in knobs:
            key = multicol_resolve_column_key(df, ("mc", knob, f"univ_{uidx}"))
            if key is None:
                missing = True
                break
            wi = np.asarray(df.loc[:, key], dtype=float)
            wi = np.nan_to_num(wi, nan=1.0, posinf=1.0, neginf=1.0)
            w *= wi
        if missing:
            continue
        col = ("mc", bundled_tag, f"univ_{uidx}", "", "", "", "")
        df.loc[:, col] = w
    return df


def _seed_or_add(bucket: dict, var_save_name: str, univ: np.ndarray, cv: np.ndarray) -> None:
    univ = np.asarray(univ, dtype=float)
    cv = np.asarray(cv, dtype=float)
    if var_save_name not in bucket:
        bucket[var_save_name] = {
            "univ_events": np.array(univ, dtype=float, copy=True),
            "cv_events": np.array(cv, dtype=float, copy=True),
        }
    else:
        entry = bucket[var_save_name]
        # Align n_univ if a later file has fewer universes.
        n_u = min(entry["univ_events"].shape[0], univ.shape[0])
        entry["univ_events"] = entry["univ_events"][:n_u] + univ[:n_u]
        entry["cv_events"] += cv


def prepare_evt_df(mc_evt_df: pd.DataFrame) -> pd.DataFrame:
    if "topo_categ" not in mc_evt_df.columns:
        mc_evt_df = mc_evt_df.copy()
        mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
    return ensure_derived_trk_kinematics_cols(mc_evt_df)


def _ntuple_level(df: pd.DataFrame) -> int:
    if isinstance(df.index, pd.MultiIndex):
        names = df.index.names or []
        if "__ntuple" in names:
            return int(names.index("__ntuple"))
        return 0
    return 0


def _remap_ntuple_index(df: pd.DataFrame, ntuple_remap: dict) -> None:
    """Rewrite ``__ntuple`` in-place so concatenated files keep unique indices."""
    if isinstance(df.index, pd.MultiIndex):
        names = list(df.index.names) if df.index.names is not None else []
        idx_loc = _ntuple_level(df)
        new_tuples = []
        for tup in df.index:
            tup = list(tup)
            old = tup[idx_loc]
            tup[idx_loc] = ntuple_remap[old]
            new_tuples.append(tuple(tup))
        df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
    elif df.index.name == "__ntuple":
        df.index = df.index.map(ntuple_remap)
    elif len(df) and isinstance(df.index[0], tuple):
        new_tuples = []
        for tup in df.index:
            tup = list(tup)
            tup[0] = ntuple_remap[tup[0]]
            new_tuples.append(tuple(tup))
        df.index = pd.Index(new_tuples)


def load_evt_batch(file_paths: Sequence[Path | str]) -> pd.DataFrame:
    """Load + concat ``evt`` tables from many ``.df`` files.

    Uses :func:`pyanalib.split_df_helpers.load_dfs` per file (all HDF splits), then
    remaps ``__ntuple`` and concatenates — same idea as ``dfs_from_dir`` /
    ``load_and_concat_mc_dfs``, but for an explicit file list.
    """
    if not file_paths:
        raise ValueError("load_evt_batch: empty file_paths")
    frames: List[pd.DataFrame] = []
    ntuple_offset = 0
    for fpath in file_paths:
        n_split = int(get_n_split(str(fpath)))
        dfs = load_dfs(str(fpath), keys2load=["evt"], n_max_concat=max(n_split, 1))
        df = dfs["evt"]
        if len(df) == 0:
            continue
        if isinstance(df.index, pd.MultiIndex):
            lvl = _ntuple_level(df)
            old_ids = np.unique(np.asarray(df.index.get_level_values(lvl)))
        elif df.index.name == "__ntuple":
            old_ids = np.unique(np.asarray(df.index))
        else:
            old_ids = np.unique(np.asarray([t[0] for t in df.index])) if len(df) else np.array([])
        remap = {old: ntuple_offset + i for i, old in enumerate(sorted(old_ids))}
        _remap_ntuple_index(df, remap)
        ntuple_offset += len(remap)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    if len(frames) == 1:
        return frames[0]
    return pd.concat(frames, axis=0, sort=False)


def accumulate_evt_df(
    mc_evt_df: pd.DataFrame,
    *,
    syst_name: str,
    knobs: Sequence[str],
    var_configs: Sequence[Any],
    acc: Dict[str, Any],
    n_universe: int,
    bkgd_subtract: bool = True,
    accumulate_total: bool = True,
    total_var_configs: Sequence[Any] | None = None,
    log_tag: str = "batch",
) -> None:
    """Add one (possibly concatenated) evt frame into ``acc`` (knob-nested)."""
    if mc_evt_df is None or len(mc_evt_df) == 0:
        return
    tot_vcs = list(total_var_configs) if total_var_configs is not None else list(var_configs)

    mc_evt_df = prepare_evt_df(mc_evt_df)
    if syst_name == "Flux" and knobs:
        mc_evt_df = drop_bad_flux_knob_weights(mc_evt_df, knobs=knobs, n_univ=n_universe)
    elif syst_name == "G4" and knobs:
        mc_evt_df = drop_bad_g4_knob_weights(mc_evt_df, knobs=knobs, n_univ=n_universe)
    if len(mc_evt_df) == 0:
        return

    if knobs:
        for knob in knobs:
            sk = ("mc", knob)
            n_u = min(int(n_universe), _count_univ_columns(mc_evt_df, sk))
            if n_u <= 0:
                continue
            knob_acc = acc.setdefault(knob, {})
            for vc in var_configs:
                try:
                    univ, cv = get_univ_rates(
                        cov_type="rate",
                        evtdf=mc_evt_df,
                        var_config=vc,
                        syst_name=sk,
                        n_univ=n_u,
                        bkgd_subtract=bkgd_subtract,
                    )
                except Exception as ex:
                    print(f"[multisim-live] skip {log_tag} knob={knob} var={vc.var_save_name}: {ex}")
                    continue
                _seed_or_add(knob_acc, vc.var_save_name, univ, cv)

        if accumulate_total and knobs and tot_vcs:
            per_knob_n = [_count_univ_columns(mc_evt_df, ("mc", k)) for k in knobs]
            n_u_tot = min(int(n_universe), min(per_knob_n) if per_knob_n else 0)
            if n_u_tot > 0:
                df_tot = inject_multiplied_mc_knob_weights(
                    mc_evt_df, knobs, bundled_tag=syst_name, n_univ=n_u_tot
                )
                tot_acc = acc.setdefault(TOTAL_KNOB_KEY, {})
                sk_tot = ("mc", syst_name)
                for vc in tot_vcs:
                    try:
                        univ, cv = get_univ_rates(
                            cov_type="rate",
                            evtdf=df_tot,
                            var_config=vc,
                            syst_name=sk_tot,
                            n_univ=n_u_tot,
                            bkgd_subtract=bkgd_subtract,
                        )
                    except Exception as ex:
                        print(f"[multisim-live] skip {log_tag} total var={vc.var_save_name}: {ex}")
                        continue
                    _seed_or_add(tot_acc, vc.var_save_name, univ, cv)
                del df_tot
    else:
        sk = syst_key_for_name(syst_name)
        if syst_name == "MCstat":
            if _count_univ_columns(mc_evt_df, ("mc", "MCstat")) > 0:
                sk = ("mc", "MCstat")
            elif _count_univ_columns(mc_evt_df, "MCstat") > 0:
                sk = "MCstat"
        n_u = min(int(n_universe), _count_univ_columns(mc_evt_df, sk))
        if n_u > 0:
            bucket = acc.setdefault(TOTAL_KNOB_KEY, {})
            for vc in var_configs:
                try:
                    univ, cv = get_univ_rates(
                        cov_type="rate",
                        evtdf=mc_evt_df,
                        var_config=vc,
                        syst_name=sk,
                        n_univ=n_u,
                        bkgd_subtract=bkgd_subtract,
                    )
                except Exception as ex:
                    print(f"[multisim-live] skip {log_tag} var={vc.var_save_name}: {ex}")
                    continue
                _seed_or_add(bucket, vc.var_save_name, univ, cv)


def accumulate_file(
    df_path: Path | str,
    *,
    syst_name: str,
    knobs: Sequence[str],
    var_configs: Sequence[Any],
    acc: Dict[str, Any],
    n_universe: int,
    bkgd_subtract: bool = True,
    accumulate_total: bool = True,
    total_var_configs: Sequence[Any] | None = None,
) -> None:
    """Backward-compatible single-file entry point (loads via :func:`load_dfs`)."""
    evt = load_evt_batch([df_path])
    accumulate_evt_df(
        evt,
        syst_name=syst_name,
        knobs=knobs,
        var_configs=var_configs,
        acc=acc,
        n_universe=n_universe,
        bkgd_subtract=bkgd_subtract,
        accumulate_total=accumulate_total,
        total_var_configs=total_var_configs,
        log_tag=Path(df_path).name,
    )
    del evt
    gc.collect()


def _iter_file_batches(files: Sequence[Path], files_per_batch: int) -> List[List[Path]]:
    n = max(1, int(files_per_batch))
    return [list(files[i : i + n]) for i in range(0, len(files), n)]


def cov_pack_for_entry(univ: np.ndarray, cv: np.ndarray) -> dict:
    return _sanitize_matrix_pack(get_covariance_matrix(np.asarray(univ), np.asarray(cv)))


def knob_display_label(knob: str) -> str:
    """Short legend label for a Flux/G4 knob name."""
    for suffix in ("_Flux", "_Geant4", "_G4"):
        if knob.endswith(suffix):
            return knob[: -len(suffix)]
    if knob == TOTAL_KNOB_KEY:
        return "product"
    return knob


def analysis_for_var(acc: dict, var_save_name: str) -> dict | None:
    """Per-variable packs: each knob, independent fractional sum, and product-of-knobs.

    Returns ``None`` if nothing is accumulated for ``var_save_name``.

    Keys
    ----
    knob_packs : dict[str, pack]
        Per-knob ``{cov, cov_frac, corr}``.
    knob_frac : dict[str, ndarray]
        Diagonal fractional uncertainty per knob.
    indep_sum_pack / indep_sum_frac
        Sum of independent fractional covariances (standard Flux/G4 combination).
    product_pack / product_frac / product_univ / product_cv
        From multiplied-knob ``__total__`` universes (when present).
    univ_events / cv_events
        Universes used for the distribution overlay — **product** when available,
        else the first knob.
    matrix_pack
        Covariance/correlation shown on the panel — product when available, else
        indep sum.
    """
    knob_packs: Dict[str, dict] = {}
    knob_order: List[str] = []
    cv_ref = None
    for knob, block in acc.items():
        if knob == TOTAL_KNOB_KEY:
            continue
        if var_save_name not in block:
            continue
        e = block[var_save_name]
        pack = cov_pack_for_entry(e["univ_events"], e["cv_events"])
        knob_packs[knob] = pack
        knob_order.append(knob)
        if cv_ref is None:
            cv_ref = np.asarray(e["cv_events"], dtype=float)

    product_pack = None
    product_univ = None
    product_cv = None
    if TOTAL_KNOB_KEY in acc and var_save_name in acc[TOTAL_KNOB_KEY]:
        e = acc[TOTAL_KNOB_KEY][var_save_name]
        product_univ = np.asarray(e["univ_events"], dtype=float)
        product_cv = np.asarray(e["cv_events"], dtype=float)
        product_pack = cov_pack_for_entry(product_univ, product_cv)
        if cv_ref is None:
            cv_ref = product_cv.copy()

    if not knob_packs and product_pack is None:
        return None

    indep_sum_pack = None
    if knob_packs and cv_ref is not None:
        packs_list = [knob_packs[k] for k in knob_order]
        indep_sum_pack = (
            combine_indep_knob_cov_packs(packs_list, cv_ref)
            if len(packs_list) > 1
            else packs_list[0]
        )

    if product_univ is not None and product_cv is not None:
        univ_events, cv_events = product_univ, product_cv
    else:
        # Fallback: first knob's universes
        first = knob_order[0]
        e = acc[first][var_save_name]
        univ_events = np.asarray(e["univ_events"], dtype=float)
        cv_events = np.asarray(e["cv_events"], dtype=float)

    matrix_pack = product_pack if product_pack is not None else indep_sum_pack
    assert matrix_pack is not None

    knob_frac = {k: frac_unc_from_pack(p) for k, p in knob_packs.items()}
    out = {
        "knob_order": knob_order,
        "knob_packs": knob_packs,
        "knob_frac": knob_frac,
        "indep_sum_pack": indep_sum_pack,
        "indep_sum_frac": frac_unc_from_pack(indep_sum_pack) if indep_sum_pack is not None else None,
        "product_pack": product_pack,
        "product_frac": frac_unc_from_pack(product_pack) if product_pack is not None else None,
        "product_univ": product_univ,
        "product_cv": product_cv,
        "univ_events": univ_events,
        "cv_events": cv_events,
        "matrix_pack": matrix_pack,
    }
    return out


def combined_pack_for_var(acc: dict, var_save_name: str) -> Tuple[dict, np.ndarray, np.ndarray] | None:
    """Backward-compatible: ``(matrix_pack, product_univ, product_cv)``."""
    a = analysis_for_var(acc, var_save_name)
    if a is None:
        return None
    return a["matrix_pack"], a["univ_events"], a["cv_events"]


def build_analysis_dict(acc: dict, var_save_names: Sequence[str]) -> Dict[str, dict]:
    """Serializable per-var analysis (packs + frac arrays) for the saved payload."""
    out: Dict[str, dict] = {}
    for vsn in var_save_names:
        a = analysis_for_var(acc, vsn)
        if a is None:
            continue
        # Drop large redundant univ arrays from nested packs? Keep histcounts in acc;
        # store packs + fracs + product univ for replot without recomputing cov.
        out[vsn] = {
            "knob_order": list(a["knob_order"]),
            "knob_packs": a["knob_packs"],
            "knob_frac": a["knob_frac"],
            "indep_sum_pack": a["indep_sum_pack"],
            "indep_sum_frac": a["indep_sum_frac"],
            "product_pack": a["product_pack"],
            "product_frac": a["product_frac"],
            "product_univ": a["product_univ"],
            "product_cv": a["product_cv"],
            "univ_events": a["univ_events"],
            "cv_events": a["cv_events"],
            "matrix_pack": a["matrix_pack"],
        }
    return out


# ---------------------------------------------------------------------------
# Drawing (axes-level; mirrors utils.plot_univ_hists / plot_frac_unc / heatmap)
# ---------------------------------------------------------------------------
def draw_univ_hists_on_ax(
    ax,
    univ_events: np.ndarray,
    cv_events: np.ndarray,
    var_config,
    *,
    univ_alpha: float = UNIV_ALPHA,
    title: str = "",
    ylabel: str = "Events",
) -> None:
    """CV in bold black; each universe as a light step histogram (``univ_alpha``)."""
    univ_events = np.asarray(univ_events, dtype=float)
    cv_events = np.asarray(cv_events, dtype=float)
    n_univ = int(univ_events.shape[0])
    bc = var_config.bin_centers
    bins = var_config.bins

    for i in range(n_univ):
        label = "Universe" if i == 0 else None
        ax.hist(
            bc,
            bins=bins,
            weights=univ_events[i],
            histtype="step",
            color="0.55",
            alpha=float(univ_alpha),
            linewidth=0.8,
            label=label,
        )

    ax.hist(
        bc,
        bins=bins,
        weights=cv_events,
        histtype="step",
        color="k",
        linewidth=2.2,
        label="Central Value",
    )
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(var_config.var_labels[0] if var_config.var_labels else "")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best")
    if getattr(var_config, "var_save_name", "") == "integrated":
        ax.set_xticks([])


def draw_frac_unc_on_ax(
    ax,
    frac_series: Sequence[Tuple[str, np.ndarray]],
    var_config,
    *,
    title: str = "",
    emphasis: Sequence[str] = ("indep sum", "product"),
) -> None:
    """Overlay fractional uncertainties. ``frac_series`` is ``[(label, frac_unc), ...]``.

    Labels in ``emphasis`` (case-insensitive match on stripped label) are drawn
    thicker / black or C0 so totals stand out over per-knob curves.
    """
    emph = {s.lower() for s in emphasis}
    for idx, (label, frac_unc) in enumerate(frac_series):
        frac_unc = np.asarray(frac_unc, dtype=float)
        is_emph = label.lower() in emph
        if is_emph and "product" in label.lower():
            color, lw = "k", 2.4
        elif is_emph and "indep" in label.lower():
            color, lw = "C0", 2.2
        else:
            color, lw = f"C{(idx % 8) + 1}", 1.0
        ax.hist(
            var_config.bin_centers,
            bins=var_config.bins,
            weights=frac_unc,
            histtype="step",
            color=color,
            linewidth=lw,
            label=label,
        )
    ax.set_xlim(var_config.bins[0], var_config.bins[-1])
    ax.set_xlabel(var_config.var_labels[0] if var_config.var_labels else "")
    ax.set_ylabel("Fractional Uncertainty")
    ax.grid(True, alpha=0.4)
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=6, loc="best", ncol=2)
    if getattr(var_config, "var_save_name", "") == "integrated":
        ax.set_xticks([])


def frac_series_from_analysis(a: dict) -> List[Tuple[str, np.ndarray]]:
    """Ordered ``(label, frac)`` for uncertainty overlays."""
    series: List[Tuple[str, np.ndarray]] = []
    for knob in a.get("knob_order") or []:
        series.append((knob_display_label(knob), a["knob_frac"][knob]))
    if a.get("indep_sum_frac") is not None:
        series.append(("indep sum", a["indep_sum_frac"]))
    if a.get("product_frac") is not None:
        series.append(("product", a["product_frac"]))
    return series

def _bin_range_labels(edges) -> List[str]:
    edges = np.asarray(edges, dtype=float)
    return [f"{edges[i]:.2f}–{edges[i+1]:.2f}" for i in range(len(edges) - 1)]


def draw_heatmap_on_ax(
    ax,
    matrix: np.ndarray,
    bins,
    *,
    title: str = "",
    cmap: str = "bwr",
    vmin=None,
    vmax=None,
    colorbar: bool = True,
) -> None:
    """Compact heatmap inspired by :func:`utils.plot_heatmap`."""
    matrix = np.asarray(matrix, dtype=float)
    nbins = len(bins)
    assert nbins - 1 == matrix.shape[0] == matrix.shape[1]
    unif = np.linspace(0.0, float(nbins - 1), nbins)
    extent = [unif[0], unif[-1], unif[0], unif[-1]]
    if cmap == "bwr" and vmin is None and vmax is None:
        vmin, vmax = -1.0, 1.0
    im = ax.imshow(matrix, extent=extent, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    tick_pos = (unif[:-1] + unif[1:]) / 2
    labels = _bin_range_labels(bins)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(labels, fontsize=7)
    if title:
        ax.set_title(title, fontsize=10)
    if colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def make_accumulation_panel(
    acc: dict,
    var_configs: Sequence[Any],
    *,
    status: str = "",
    plot_var_names: Sequence[str] | None = None,
    univ_alpha: float = UNIV_ALPHA,
    figsize_per_row: Tuple[float, float] = (16.0, 3.6),
    syst_label: str = "Flux",
) -> plt.Figure:
    """One row per variable: [product universes | cov | corr | frac unc breakdown]."""
    by_name = {vc.var_save_name: vc for vc in var_configs}
    names = list(plot_var_names) if plot_var_names is not None else list(default_plot_var_names())
    names = [n for n in names if n in by_name]
    usable = []
    analyses = {}
    for n in names:
        a = analysis_for_var(acc, n)
        if a is not None:
            usable.append(n)
            analyses[n] = a
    if not usable:
        fig, ax = plt.subplots(figsize=(8, 2.5))
        ax.text(0.5, 0.5, "No accumulated histograms yet", ha="center", va="center")
        ax.set_axis_off()
        fig.suptitle(status or "waiting…")
        return fig

    nrows = len(usable)
    fig_w, row_h = figsize_per_row
    fig, axes = plt.subplots(nrows, 4, figsize=(fig_w, max(row_h, 4.0) * nrows), squeeze=False)

    for r, vsn in enumerate(usable):
        vc = by_name[vsn]
        a = analyses[vsn]
        pack = a["matrix_pack"]
        label = getattr(vc, "var_plot_name", vsn)
        has_product = a.get("product_univ") is not None

        draw_univ_hists_on_ax(
            axes[r, 0],
            a["univ_events"],
            a["cv_events"],
            vc,
            univ_alpha=univ_alpha,
            title=(
                f"{label} — {syst_label} product universes"
                if has_product
                else f"{label} — {syst_label} universes"
            ),
        )
        draw_heatmap_on_ax(
            axes[r, 1],
            pack["cov"],
            vc.bins,
            title="Covariance (product)" if has_product else "Covariance",
            cmap="viridis",
            vmin=None,
            vmax=None,
        )
        draw_heatmap_on_ax(
            axes[r, 2],
            pack["corr"],
            vc.bins,
            title="Correlation (product)" if has_product else "Correlation",
            cmap="bwr",
            vmin=-1,
            vmax=1,
        )
        draw_frac_unc_on_ax(
            axes[r, 3],
            frac_series_from_analysis(a),
            vc,
            title="Fractional uncertainty",
        )

    fig.suptitle(status or f"{syst_label} multisim accumulation", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    return fig

def build_payload(
    acc: dict,
    *,
    cfg: LiveMultisimConfig,
    knobs: Sequence[str],
    var_configs: Sequence[Any],
    n_files_done: int,
    n_files_total: int,
    df_files: Sequence[str],
) -> dict:
    var_save_names = [vc.var_save_name for vc in var_configs]
    analysis = build_analysis_dict(acc, var_save_names)
    return {
        "schema": 2,
        "syst_name": cfg.syst_name,
        "flux_knob_groups": cfg.flux_knob_groups,
        "n_universe": cfg.n_universe,
        "var_set": cfg.var_set,
        "knobs": list(knobs),
        "var_save_names": var_save_names,
        "var_bins": {vc.var_save_name: np.asarray(vc.bins, dtype=float) for vc in var_configs},
        "var_labels": {vc.var_save_name: list(vc.var_labels) for vc in var_configs},
        "var_plot_names": {
            vc.var_save_name: getattr(vc, "var_plot_name", vc.var_save_name) for vc in var_configs
        },
        "acc": acc,
        # Derived products saved alongside raw histcounts (knob / indep sum / product).
        "analysis": analysis,
        "n_files_done": n_files_done,
        "n_files_total": n_files_total,
        "df_dir": str(cfg.df_dir),
        "df_files": list(df_files),
        "accumulate_total": True if knobs else bool(cfg.accumulate_total),
        "bkgd_subtract": cfg.bkgd_subtract,
        "created": datetime.now().isoformat(timespec="seconds"),
    }


def save_payload(payload: dict, path: Path | str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)
    return path


def load_payload(path: Path | str) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)


class LiveNotebookDisplay:
    """Update one inline image in the notebook cell (no file write)."""

    def __init__(self, dpi: int = 110):
        self.dpi = dpi
        self._display_id = None

    def __call__(self, fig: plt.Figure, status: str = "", *_args) -> None:
        from IPython.display import Image, display, update_display

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=self.dpi, bbox_inches="tight")
        buf.seek(0)
        plt.close(fig)
        img = Image(data=buf.getvalue())
        if self._display_id is None:
            handle = display(img, display_id=True)
            self._display_id = handle.display_id
        else:
            update_display(img, display_id=self._display_id)
        if status:
            print(status, flush=True)


def run_live_accumulate(
    cfg: LiveMultisimConfig | None = None,
    *,
    on_update: Optional[Callable[[plt.Figure, str, dict], None]] = None,
) -> LiveMultisimResult:
    """Load ``.df`` files in batches, accumulate, and refresh the live panel."""
    cfg = cfg or LiveMultisimConfig()
    if not cfg.df_dir:
        raise ValueError("LiveMultisimConfig.df_dir is required")

    work = Path(cfg.work_dir or default_work_dir(cfg.syst_name)).expanduser()
    plots = Path(cfg.plots_dir or (work / "plots")).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    payload_path = work / "accumulated_histcounts.pkl"

    var_configs_all = build_var_configs(cfg.var_set)
    knobs = knobs_for_syst(cfg.syst_name, cfg.flux_knob_groups)
    files = discover_df_files(cfg.df_dir, cfg.max_files)
    if not files:
        raise FileNotFoundError(f"no *.df under {cfg.df_dir}")

    plot_names = (
        list(cfg.plot_var_names)
        if cfg.plot_var_names is not None
        else list(default_plot_var_names())
    )
    accum_names = (
        list(cfg.accumulate_var_names)
        if cfg.accumulate_var_names is not None
        else list(plot_names)
    )
    by_name = {vc.var_save_name: vc for vc in var_configs_all}
    missing = [n for n in accum_names if n not in by_name]
    if missing:
        raise KeyError(f"unknown accumulate_var_names (not in var_set={cfg.var_set!r}): {missing}")
    var_configs = [by_name[n] for n in accum_names]
    do_total = bool(knobs) and bool(cfg.accumulate_total)
    total_vcs = list(var_configs) if do_total else []

    batches = _iter_file_batches(files, cfg.files_per_batch)
    acc: Dict[str, Any] = {}
    failed: List[Tuple[str, str]] = []
    display_fn = on_update
    if display_fn is None and cfg.show_in_notebook:
        display_fn = LiveNotebookDisplay(dpi=cfg.panel_dpi)

    print(
        f"[multisim-live] files={len(files)}  batches={len(batches)}  "
        f"files_per_batch={cfg.files_per_batch}  knobs={len(knobs) or 1}  "
        f"n_univ≤{cfg.n_universe}  accumulate_vars={len(var_configs)}  "
        f"plot_vars={len(plot_names)}  product_total={do_total}",
        flush=True,
    )

    def _should_update(*, done: int, is_last: bool) -> bool:
        if is_last:
            return True
        n = max(1, int(cfg.update_every_n_files))
        return done > 0 and (done % n == 0)

    def _save_panel_png(fig: plt.Figure, done: int) -> None:
        if not cfg.save_panel_png:
            return
        latest = plots / "live_panel_latest.png"
        snap = plots / f"live_panel_files{done:04d}.png"
        fig.savefig(latest, dpi=cfg.panel_dpi, bbox_inches="tight")
        fig.savefig(snap, dpi=cfg.panel_dpi, bbox_inches="tight")

    def _refresh(done: int) -> None:
        status = (
            f"{cfg.syst_name} live  files={done}/{len(files)}  "
            f"batch_size≤{cfg.files_per_batch}  knobs={len(knobs) or 1}  "
            f"n_univ≤{cfg.n_universe}  dir={Path(cfg.df_dir).name}"
        )
        fig = make_accumulation_panel(
            acc,
            var_configs,
            status=status,
            plot_var_names=plot_names,
            univ_alpha=cfg.univ_alpha,
            figsize_per_row=cfg.figsize_per_row,
            syst_label=cfg.syst_name,
        )
        payload = build_payload(
            acc,
            cfg=cfg,
            knobs=knobs,
            var_configs=var_configs,
            n_files_done=done,
            n_files_total=len(files),
            df_files=[str(p) for p in files[:done]],
        )
        payload["files_per_batch"] = int(cfg.files_per_batch)
        if cfg.save_every_update:
            save_payload(payload, payload_path)
        _save_panel_png(fig, done)
        if display_fn is not None:
            display_fn(fig, status, payload)
        else:
            plt.close(fig)

    files_done = 0
    for bi, batch in enumerate(tqdm(batches, desc=f"{cfg.syst_name} batches"), start=1):
        tag = f"batch{bi}/{len(batches)}({len(batch)} files)"
        try:
            t0 = datetime.now()
            evt = load_evt_batch(batch)
            t_load = (datetime.now() - t0).total_seconds()
            print(
                f"[multisim-live] {tag} loaded evt rows={len(evt):,} in {t_load:.1f}s — accumulating…",
                flush=True,
            )
            t1 = datetime.now()
            accumulate_evt_df(
                evt,
                syst_name=cfg.syst_name,
                knobs=knobs,
                var_configs=var_configs,
                acc=acc,
                n_universe=cfg.n_universe,
                bkgd_subtract=cfg.bkgd_subtract,
                accumulate_total=do_total,
                total_var_configs=total_vcs,
                log_tag=tag,
            )
            t_acc = (datetime.now() - t1).total_seconds()
            print(f"[multisim-live] {tag} accumulate done in {t_acc:.1f}s", flush=True)
            del evt
            gc.collect()
        except Exception as ex:
            for fpath in batch:
                failed.append((str(fpath), str(ex)))
            print(f"[multisim-live] FAILED {tag}: {ex}", flush=True)
        files_done += len(batch)
        if _should_update(done=files_done, is_last=(bi == len(batches))):
            _refresh(files_done)

    # Final save even if save_every_update is False.
    payload = build_payload(
        acc,
        cfg=cfg,
        knobs=knobs,
        var_configs=var_configs,
        n_files_done=len(files),
        n_files_total=len(files),
        df_files=[str(p) for p in files],
    )
    payload["files_per_batch"] = int(cfg.files_per_batch)
    payload["failed"] = failed
    save_payload(payload, payload_path)

    return LiveMultisimResult(
        work_dir=work,
        plots_dir=plots,
        payload_path=payload_path,
        n_files_done=len(files),
        n_files_total=len(files),
        knobs=tuple(knobs),
        var_configs=list(var_configs),
        acc=acc,
        failed=failed,
    )


def packs_from_payload(payload: dict) -> Dict[str, dict]:
    """``var_save_name → analysis dict`` (knob / indep sum / product + universes).

    Prefers the saved ``analysis`` block (schema ≥ 2); otherwise recomputes from ``acc``.
    """
    if payload.get("analysis"):
        return dict(payload["analysis"])
    acc = payload["acc"]
    names = payload.get("var_save_names") or []
    if not names:
        seen = set()
        for block in acc.values():
            seen.update(block.keys())
        names = sorted(seen)
    return build_analysis_dict(acc, names)


def analysis_from_payload(payload: dict, var_save_name: str) -> dict | None:
    packs = packs_from_payload(payload)
    return packs.get(var_save_name)

def var_configs_from_payload(payload: dict):
    """Rebuild VariableConfigs from the live ``var_set``, falling back to payload bins."""
    var_set = payload.get("var_set", "final")
    live = {vc.var_save_name: vc for vc in build_var_configs(var_set)}
    out = []
    for vsn in payload.get("var_save_names", []):
        if vsn in live:
            out.append(live[vsn])
            continue
        bins = np.asarray(payload["var_bins"][vsn], dtype=float)
        labels = list(payload.get("var_labels", {}).get(vsn, [vsn, vsn]))
        plot_name = payload.get("var_plot_names", {}).get(vsn, vsn)

        class _Stub:
            pass

        vc = _Stub()
        vc.var_save_name = vsn
        vc.var_plot_name = plot_name
        vc.var_labels = labels
        vc.bins = bins
        vc.bin_centers = (bins[:-1] + bins[1:]) / 2.0
        out.append(vc)
    return out
