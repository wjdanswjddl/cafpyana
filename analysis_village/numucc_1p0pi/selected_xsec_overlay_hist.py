"""Fill, cache, and plot selected-event MC/data overlays from histogram counts.

First pass: fill per-category MC + data histograms from event dataframes and
pickle them. Later passes: reload counts and replot without loading DFs.

Plotting uses :func:`overlay_hists_from_counts`, which renders the same stacked
overlay as :func:`analysis_village.numucc_1p0pi.utils.overlay_hists` without
modifying that DF-based path.
"""

from __future__ import annotations

import pickle
from os import makedirs, path
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from analysis_village.numucc_1p0pi.selection_framework import (
    BREAKDOWN_REGISTRY,
    OverlayHistData,
)
from analysis_village.numucc_1p0pi.utils import overlay_hists_from_histdata

HISTDATA_PKL_NAME = "overlay_histdata.pkl"
PAYLOAD_FORMAT = "selected_xsec_overlay_histdata_v1"

HistKey = Tuple[str, str]  # (var_save_name, breakdown_type)

# ``get_genie_sb_category(ret_cuts=True)`` is cosmic-first with S/B pairs after NC.
# ``get_genie_category`` is the same order with each mode's S+B summed.
_GENIE_SB_TO_GENIE_GROUPS: Tuple[Tuple[int, ...], ...] = (
    (0,),
    (1,),
    (2,),
    (3, 4),
    (5, 6),
    (7, 8),
    (9, 10),
)


def collapse_genie_sb_to_genie(hd_sb: OverlayHistData) -> OverlayHistData:
    """GENIE-mode stack: sum each genie_sb S/B pair (no per-mode topology split)."""
    n_sb = BREAKDOWN_REGISTRY["genie_sb"][0]
    n_g = BREAKDOWN_REGISTRY["genie"][0]
    if hd_sb.breakdown_type != "genie_sb":
        raise ValueError(f"expected genie_sb OverlayHistData, got {hd_sb.breakdown_type!r}")
    if hd_sb.mc_hist is None or hd_sb.mc_hist.shape[0] != n_sb:
        raise ValueError(
            f"genie_sb mc_hist has {None if hd_sb.mc_hist is None else hd_sb.mc_hist.shape[0]} "
            f"layers; expected {n_sb}"
        )
    if n_g != len(_GENIE_SB_TO_GENIE_GROUPS):
        raise RuntimeError("GENIE collapse groups do not match genie category count")

    hd = OverlayHistData(
        var_save_name=hd_sb.var_save_name,
        breakdown_type="genie",
        bins=np.asarray(hd_sb.bins, dtype=float),
    )
    for i, idxs in enumerate(_GENIE_SB_TO_GENIE_GROUPS):
        hd.mc_hist[i] = sum(np.asarray(hd_sb.mc_hist[j], dtype=float) for j in idxs)
        hd.mc_err2[i] = sum(np.asarray(hd_sb.mc_err2[j], dtype=float) for j in idxs)

    def _copy1(src, dst_attr):
        val = getattr(hd_sb, src)
        if val is not None:
            setattr(hd, dst_attr, np.asarray(val, dtype=float).copy())

    for attr in (
        "data_hist",
        "data_err2",
        "intime_hist",
        "intime_err2",
        "offbeam_hist",
        "offbeam_err2",
        "dirt_hist",
        "dirt_err2",
    ):
        _copy1(attr, attr)

    if hd_sb.mc_univ_hist:
        collapsed: Dict[str, np.ndarray] = {}
        for tag, univ in hd_sb.mc_univ_hist.items():
            arr = np.asarray(univ, dtype=float)
            if arr.ndim != 3 or arr.shape[1] != n_sb:
                continue
            out = np.zeros((arr.shape[0], n_g, arr.shape[2]), dtype=float)
            for i, idxs in enumerate(_GENIE_SB_TO_GENIE_GROUPS):
                out[:, i, :] = sum(arr[:, j, :] for j in idxs)
            collapsed[tag] = out
        hd.mc_univ_hist = collapsed or None

    hd.has_mc = hd_sb.has_mc
    hd.has_intime = hd_sb.has_intime
    hd.has_offbeam = getattr(hd_sb, "has_offbeam", False)
    hd.has_dirt = hd_sb.has_dirt
    hd.has_data = hd_sb.has_data
    return hd


def ensure_genie_breakdown(
    histdata_map: MutableMapping[HistKey, OverlayHistData],
) -> int:
    """Add ``genie`` keys by collapsing ``genie_sb`` when the mode-only stack is missing."""
    added = 0
    sb_keys = [k for k in list(histdata_map) if k[1] == "genie_sb"]
    for vsn, _bt in sb_keys:
        gkey = (vsn, "genie")
        if gkey in histdata_map:
            continue
        histdata_map[gkey] = collapse_genie_sb_to_genie(histdata_map[(vsn, "genie_sb")])
        added += 1
    return added


def histdata_pkl_path(out_dir: str) -> str:
    return path.join(out_dir, HISTDATA_PKL_NAME)


def fill_overlay_histdata(
    var_config,
    breakdown_type: str,
    mc_df: Optional[pd.DataFrame] = None,
    data_df: Optional[pd.DataFrame] = None,
    intime_df: Optional[pd.DataFrame] = None,
    offbeam_df: Optional[pd.DataFrame] = None,
    dirt_df: Optional[pd.DataFrame] = None,
) -> OverlayHistData:
    """Build stacked-ready counts for one (variable, breakdown).

    MC is filled per mode (topology / genie / genie_sb category order matches
    ``overlay_hists`` / ``BREAKDOWN_REGISTRY``).
    """
    from analysis_village.numucc_1p0pi.evt_derived_kinematics import (
        ensure_derived_trk_kinematics_cols,
    )

    hd = OverlayHistData(
        var_save_name=var_config.var_save_name,
        breakdown_type=breakdown_type,
        bins=var_config.bins,
    )
    col = var_config.var_evt_reco_col
    # Ensure φ exists when this var needs it (CAFs often omit precomputed phi).
    need_phi = str(getattr(var_config, "var_save_name", "")).endswith("dir_phi") or (
        isinstance(col, tuple) and len(col) >= 4 and col[3] == "phi"
    )
    frames = []
    for label, df in (
        ("mc", mc_df),
        ("data", data_df),
        ("intime", intime_df),
        ("offbeam", offbeam_df),
        ("dirt", dirt_df),
    ):
        if df is None:
            frames.append((label, None))
            continue
        if need_phi:
            df = ensure_derived_trk_kinematics_cols(df)
        frames.append((label, df))
    by = {k: v for k, v in frames}
    if by.get("mc") is not None:
        hd.fill_from_df(by["mc"], col, "mc")
    if by.get("data") is not None:
        hd.fill_from_df(by["data"], col, "data")
    if by.get("intime") is not None:
        hd.fill_from_df(by["intime"], col, "intime")
    if by.get("offbeam") is not None:
        hd.fill_from_df(by["offbeam"], col, "offbeam")
    if by.get("dirt") is not None:
        hd.fill_from_df(by["dirt"], col, "dirt")
    return hd


def build_overlay_histdata_map(
    var_configs: Sequence,
    breakdown_types: Sequence[str],
    mc_df: Optional[pd.DataFrame] = None,
    data_df: Optional[pd.DataFrame] = None,
    intime_df: Optional[pd.DataFrame] = None,
    offbeam_df: Optional[pd.DataFrame] = None,
    dirt_df: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> Dict[HistKey, OverlayHistData]:
    """Fill counts for every (var_config, breakdown_type)."""
    out: Dict[HistKey, OverlayHistData] = {}
    for var_config in var_configs:
        for breakdown_type in breakdown_types:
            key = (var_config.var_save_name, breakdown_type)
            if verbose:
                print(f"  fill counts {key[0]} ({key[1]})")
            out[key] = fill_overlay_histdata(
                var_config,
                breakdown_type,
                mc_df=mc_df,
                data_df=data_df,
                intime_df=intime_df,
                offbeam_df=offbeam_df,
                dirt_df=dirt_df,
            )
    return out


def save_overlay_counts(
    out_dir: str,
    histdata_map: Mapping[HistKey, OverlayHistData],
    *,
    pot_label: str,
    plot_set: Optional[Mapping[str, Any]] = None,
    var_save_names: Optional[Iterable[str]] = None,
    breakdown_types: Optional[Iterable[str]] = None,
) -> str:
    """Pickle counts next to the figure directory. Returns the pickle path."""
    makedirs(out_dir, exist_ok=True)
    pkl = histdata_pkl_path(out_dir)
    payload = {
        "format": PAYLOAD_FORMAT,
        "pot_label": pot_label,
        "plot_set": dict(plot_set) if plot_set is not None else None,
        "var_save_names": list(var_save_names)
        if var_save_names is not None
        else sorted({k[0] for k in histdata_map}),
        "breakdown_types": list(breakdown_types)
        if breakdown_types is not None
        else sorted({k[1] for k in histdata_map}),
        "histdata": dict(histdata_map),
    }
    with open(pkl, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    return pkl


def load_overlay_counts(out_dir: str) -> Optional[dict]:
    """Load a previously saved counts payload, or ``None`` if missing."""
    pkl = histdata_pkl_path(out_dir)
    if not path.isfile(pkl):
        return None
    with open(pkl, "rb") as f:
        payload = pickle.load(f)
    if not isinstance(payload, dict) or "histdata" not in payload:
        raise ValueError(f"unrecognized overlay counts payload: {pkl}")
    return payload


def overlay_hists_from_counts(
    histdata: OverlayHistData,
    var_config=None,
    **kwargs,
):
    """Plot a stacked MC + data overlay from precomputed histogram counts.

    Same visual result as ``overlay_hists(..., mc_df=..., data_df=...)`` for
    equivalent inputs. Does not touch the dataframe-based plotter path.
    """
    return overlay_hists_from_histdata(
        histdata,
        var_config=var_config,
        **kwargs,
    )


def _final_sel_vlines(var_save_name: str):
    """Cut arrows for selected-sample overlays.

    Product B overlays are the *already selected* sample (inside the pμ 0.22–1
    and pp 0.3–1 GeV analysis window). The arrows were copied from cut-stage
    plots that mark those kinematic thresholds on pre-cut distributions; they
    do not belong on the PRL data–MC overlays. Cut-stage plots still pass
    ``vline`` via ``selected_xsec_overlay_cut_vars._vlines_for_tag``.
    """
    return None


def plot_overlay_counts_map(
    histdata_map: Mapping[HistKey, OverlayHistData],
    var_configs: Sequence,
    breakdown_types: Sequence[str],
    *,
    pot_label: str,
    out_dir: str,
    get_syst=None,
    ax_ylim_ratio: float = 1.9,
    ratio: bool = True,
    textloc=None,
    approval: str = "internal",
    pot_text=None,
    save_fig: bool = True,
    plot: bool = True,
    textchi2: bool = True,
    verbose: bool = True,
    cosmic_estimate: str = "offbeam",
) -> None:
    """Render every cached (var, breakdown) overlay to ``out_dir``.

    ``cosmic_estimate`` is forwarded to ``overlay_hists_from_histdata``
    (``\"offbeam\"`` or ``\"intime\"``). Track-PDG stacks a separate
    ``Intime Cosmics`` layer; topology/genie/genie_sb fold into one ``Cosmics``.
    """
    if textloc is None:
        textloc = [0.03, 0.55]
    if save_fig:
        makedirs(out_dir, exist_ok=True)

    if isinstance(histdata_map, dict):
        n_added = ensure_genie_breakdown(histdata_map)
        if n_added and verbose:
            print(f"  synthesized {n_added} genie stacks from genie_sb", flush=True)

    for var_config in var_configs:
        syst = get_syst(var_config) if get_syst is not None else None
        vline = _final_sel_vlines(var_config.var_save_name)
        for breakdown_type in breakdown_types:
            key = (var_config.var_save_name, breakdown_type)
            hd = histdata_map.get(key)
            if hd is None:
                print(f"  skip missing histdata for {key}", flush=True)
                continue
            save_name = path.join(out_dir, f"{key[0]}_{key[1]}")
            plot_labels = [var_config.var_labels[1], pot_label, ""]
            if verbose:
                print(f"  plot from counts {key[0]} ({key[1]})")
            overlay_hists_from_counts(
                hd,
                var_config=var_config,
                plot_labels=plot_labels,
                ax_ylim_ratio=ax_ylim_ratio,
                ratio=ratio,
                textloc=textloc,
                approval=approval,
                pot_text=pot_text if pot_text is not None else pot_label,
                save_fig=save_fig,
                plot=plot,
                textchi2=textchi2,
                syst=syst,
                vline=vline,
                save_name=save_name,
                cosmic_estimate=cosmic_estimate,
            )
