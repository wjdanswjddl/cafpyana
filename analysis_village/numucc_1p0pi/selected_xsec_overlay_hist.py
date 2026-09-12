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
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import pandas as pd

from analysis_village.numucc_1p0pi.selection_framework import OverlayHistData
from analysis_village.numucc_1p0pi.utils import overlay_hists_from_histdata

HISTDATA_PKL_NAME = "overlay_histdata.pkl"
PAYLOAD_FORMAT = "selected_xsec_overlay_histdata_v1"

HistKey = Tuple[str, str]  # (var_save_name, breakdown_type)


def histdata_pkl_path(out_dir: str) -> str:
    return path.join(out_dir, HISTDATA_PKL_NAME)


def fill_overlay_histdata(
    var_config,
    breakdown_type: str,
    mc_df: Optional[pd.DataFrame] = None,
    data_df: Optional[pd.DataFrame] = None,
    intime_df: Optional[pd.DataFrame] = None,
    dirt_df: Optional[pd.DataFrame] = None,
) -> OverlayHistData:
    """Build stacked-ready counts for one (variable, breakdown).

    MC is filled per mode (topology / genie_sb category order matches
    ``overlay_hists`` / ``BREAKDOWN_REGISTRY``).
    """
    hd = OverlayHistData(
        var_save_name=var_config.var_save_name,
        breakdown_type=breakdown_type,
        bins=var_config.bins,
    )
    col = var_config.var_evt_reco_col
    if mc_df is not None:
        hd.fill_from_df(mc_df, col, "mc")
    if data_df is not None:
        hd.fill_from_df(data_df, col, "data")
    if intime_df is not None:
        hd.fill_from_df(intime_df, col, "intime")
    if dirt_df is not None:
        hd.fill_from_df(dirt_df, col, "dirt")
    return hd


def build_overlay_histdata_map(
    var_configs: Sequence,
    breakdown_types: Sequence[str],
    mc_df: Optional[pd.DataFrame] = None,
    data_df: Optional[pd.DataFrame] = None,
    intime_df: Optional[pd.DataFrame] = None,
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
    save_fig: bool = True,
    plot: bool = True,
    textchi2: bool = True,
    verbose: bool = True,
) -> None:
    """Render every cached (var, breakdown) overlay to ``out_dir``."""
    if textloc is None:
        textloc = [0.03, 0.55]
    if save_fig:
        makedirs(out_dir, exist_ok=True)

    for var_config in var_configs:
        syst = get_syst(var_config) if get_syst is not None else None
        for breakdown_type in breakdown_types:
            key = (var_config.var_save_name, breakdown_type)
            hd = histdata_map.get(key)
            if hd is None:
                raise KeyError(f"missing histdata for {key} in counts map")
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
                save_fig=save_fig,
                plot=plot,
                textchi2=textchi2,
                syst=syst,
                save_name=save_name,
            )
