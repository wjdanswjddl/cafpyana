"""Notebook helpers for MCstat / Flux / G4 Product-B multisim covariances.

Shared load → FV/topo prep → bad-weight filter → universe rates → cov packs,
plus Flux asymmetry side-by-side plots. Chunked production remains in
``scripts/syst_multisim_*.py``; live walk remains in ``syst_multisim_live.py``.
"""
from __future__ import annotations

import json
import time
from datetime import datetime
from os import makedirs, path
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from pyanalib.covariance import get_covariance_matrix
from pyanalib.pandas_helpers import pad_column_name
from pyanalib.split_df_helpers_new import dfs_from_dir
from pyanalib.variable_calculator import get_cc1p0pi_tki

from analysis_village.numucc_1p0pi.categories import get_topo_category
from analysis_village.numucc_1p0pi.dataset_locations import MULTISIM_SYST_GLOBS_FINAL
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_FLUX,
    FILE_G4,
    FILE_MCSTAT,
    SUB_FLUX,
    SUB_G4,
    SUB_MCSTAT,
    category_out_dir,
    normalized_root,
)
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    drop_bad_flux_knob_weights,
    drop_bad_g4_knob_weights,
    flux_mc_knob_names,
    g4_mc_knob_names,
    save_neutrino_multisim_npzs,
    syst_key_for_name,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

SYST_DISK_SUBFILE = {
    "MCstat": (SUB_MCSTAT, FILE_MCSTAT),
    "Flux": (SUB_FLUX, FILE_FLUX),
    "G4": (SUB_G4, FILE_G4),
}

_FLUX_KNOB_DISPLAY = {
    "expskin": "Exp.\nskin",
    "horncurrent": "Horn\ncurrent",
    "kminus": r"$K^-$",
    "kplus": r"$K^+$",
    "kzero": r"$K^0$",
    "piminus": r"$\pi^-$",
    "piplus": r"$\pi^+$",
    "pioninexsec": r"$\pi$ inel. $\sigma$",
    "pionqexsec": r"$\pi$ QE $\sigma$",
    "piontotxsec": r"$\pi$ total $\sigma$",
    "nucleoninexsec": r"N inel. $\sigma$",
    "nucleonqexsec": r"N QE $\sigma$",
    "nucleontotxsec": r"N total $\sigma$",
}


def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Load / prep
# ---------------------------------------------------------------------------

def multisim_search_dir(syst_name: str) -> str:
    """Parent directory of ``MULTISIM_SYST_GLOBS_FINAL[syst_name]``."""
    if syst_name not in MULTISIM_SYST_GLOBS_FINAL:
        raise KeyError(f"no glob for {syst_name!r}; edit MULTISIM_SYST_GLOBS_FINAL")
    return str(Path(MULTISIM_SYST_GLOBS_FINAL[syst_name]).parent)


def load_sel_mup_evt(syst_name: str, *, n_max_concat: int = 999) -> pd.DataFrame:
    dirname = multisim_search_dir(syst_name)
    log(f"loading {syst_name} MC from {dirname} ...")
    t0 = time.time()
    dfs = dfs_from_dir(
        search_dir=dirname,
        filename_str="sel_mup",
        keys2load=["hdr", "evt"],
        n_max_concat=n_max_concat,
    )
    if "evt" not in dfs:
        raise RuntimeError(f"no evt table for {syst_name} under {dirname}")
    log(f"  {syst_name}: {len(dfs['evt']):,} rows in {time.time() - t0:.1f}s")
    return dfs["evt"]


def evt_df_fixed(df: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    """TKI derived cols + FV ``|x|>10`` + ``topo_categ`` (notebook Product B prep)."""
    slc_mudf = df.mu.pfp.trk
    slc_pdf = df.p.pfp.trk
    tki_reco = get_cc1p0pi_tki(
        slc_mudf,
        slc_pdf,
        pad_column_name(("P", "p_muon"), slc_mudf),
        pad_column_name(("P", "p_proton"), slc_pdf),
    )
    df["del_Tp_x"] = tki_reco["del_Tp_x"]
    df["del_Tp_y"] = tki_reco["del_Tp_y"]

    mc_mudf = df.mu.pfp.trk.truth.p
    mc_pdf = df.p.pfp.trk.truth.p
    tki_mc = get_cc1p0pi_tki(
        mc_mudf,
        mc_pdf,
        pad_column_name(("totp",), mc_mudf),
        pad_column_name(("totp",), mc_pdf),
    )
    df["mc_del_Tp_x"] = tki_mc["del_Tp_x"]
    df["mc_del_Tp_y"] = tki_mc["del_Tp_y"]
    df[("mc", "del_Tp_x")] = tki_mc["del_Tp_x"]
    df[("mc", "del_Tp_y")] = tki_mc["del_Tp_y"]

    n_before = len(df)
    df = df[np.abs(df.slc.vertex.x) > 10]
    if "topo_categ" not in df.columns:
        df = df.copy()
        df.loc[:, "topo_categ"] = get_topo_category(df)
    return df, n_before


def load_and_prepare_evt(syst_name: str, *, n_max_concat: int = 999) -> pd.DataFrame:
    raw = load_sel_mup_evt(syst_name, n_max_concat=n_max_concat)
    log(f"running evt_df_fixed for {syst_name} ...")
    t0 = time.time()
    evtdf, n_before = evt_df_fixed(raw)
    log(
        f"  {syst_name}: {len(evtdf):,} events after FV cut ({n_before:,} before) "
        f"in {time.time() - t0:.1f}s"
    )
    if "topo_categ" in evtdf.columns:
        log("  topo_categ counts: " + repr(evtdf.topo_categ.value_counts().to_dict()))
    return evtdf


# ---------------------------------------------------------------------------
# Knobs / cleaning
# ---------------------------------------------------------------------------

def syst_knob_names(syst_name: str, *, flux_groups: str = "all") -> Tuple[str, ...]:
    """Per-knob weight columns for Flux/G4; empty for bundled MCstat."""
    if syst_name == "Flux":
        return flux_mc_knob_names(flux_groups)
    if syst_name == "G4":
        return g4_mc_knob_names()
    return ()


def drop_bad_weights(
    evtdf: pd.DataFrame,
    syst_name: str,
    knobs: Sequence[str],
    n_univ: int,
) -> pd.DataFrame:
    if syst_name == "Flux" and knobs:
        return drop_bad_flux_knob_weights(evtdf, knobs=knobs, max_wgt=1e3, n_univ=n_univ)
    if syst_name == "G4" and knobs:
        return drop_bad_g4_knob_weights(evtdf, knobs=knobs, max_wgt=1e3, n_univ=n_univ)
    return evtdf


def resolve_bundled_syst_name(evtdf: pd.DataFrame, syst_name: str):
    """Weight column root for ``get_univ_rates`` (MCstat layout variants)."""
    if syst_name == "MCstat":
        for sk in (("mc", "MCstat"), "MCstat"):
            probe = sk + ("univ_0",) if isinstance(sk, tuple) else (sk, "univ_0")
            if multicol_resolve_column_key(evtdf, probe) is not None:
                return sk
        return ("mc", "MCstat")
    return syst_key_for_name(syst_name)


def probe_weight_columns(
    evtdf: pd.DataFrame,
    syst_name: str,
    *,
    flux_groups: str = "all",
) -> None:
    knobs = syst_knob_names(syst_name, flux_groups=flux_groups)
    if knobs:
        log(f"{syst_name}: {len(knobs)} knobs; first={knobs[0]!r} ... last={knobs[-1]!r}")
        k0 = knobs[0]
        key0 = multicol_resolve_column_key(evtdf, ("mc", k0, "univ_0"))
        log(f"  weight probe knob={k0!r} univ_0 -> {key0!r}")
    else:
        sk = resolve_bundled_syst_name(evtdf, syst_name)
        probe = sk + ("univ_0",) if isinstance(sk, tuple) else (sk, "univ_0")
        key0 = multicol_resolve_column_key(evtdf, probe)
        log(f"{syst_name}: bundled weights; syst_name={sk!r} univ_0 -> {key0!r}")


# ---------------------------------------------------------------------------
# Covariance packs
# ---------------------------------------------------------------------------

def sanitize_matrix_pack(ret: Mapping[str, Any]) -> dict:
    out = {}
    for k in ("cov", "cov_frac", "corr"):
        out[k] = np.nan_to_num(np.asarray(ret[k], dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    return out


def frac_unc_from_pack(pack: Mapping[str, Any]) -> np.ndarray:
    return np.sqrt(np.maximum(np.diag(pack["cov_frac"]), 0.0))


def inject_multiplied_mc_knob_weights(
    df: pd.DataFrame,
    knobs: Sequence[str],
    bundled_tag: str = "Flux",
    n_univ: int = 100,
    verbose: bool = False,
) -> pd.DataFrame:
    """Write ``(mc, bundled_tag, univ_i)`` as the product of per-knob weights."""
    if verbose:
        log(
            f"inject multiplied weights ({bundled_tag}): "
            f"{len(knobs)} knobs x {n_univ} universes on {len(df):,} events"
        )
    t0 = time.time()
    df = df.copy()
    n_written = 0
    for uidx in tqdm(range(int(n_univ)), desc=f"inject {bundled_tag} product", disable=not verbose):
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
        n_written += 1
    if verbose:
        log(f"  wrote {n_written}/{n_univ} multiplied-universe columns in {time.time() - t0:.1f}s")
    return df


def cov_pack_for_knob(
    evtdf,
    var_config,
    knob: str,
    n_univ: int,
    cov_type: str,
    bkgd_subtract: bool,
    verbose: bool = False,
):
    univ, cv = get_univ_rates(
        cov_type,
        evtdf=evtdf,
        nudf=None,
        bkgd_subtract=bkgd_subtract,
        var_config=var_config,
        syst_name=("mc", knob),
        n_univ=n_univ,
        verbose=verbose,
    )
    pack = sanitize_matrix_pack(get_covariance_matrix(univ, cv))
    return pack, univ, cv


def cov_pack_bundled(
    evtdf,
    var_config,
    syst_name: str,
    n_univ: int,
    cov_type: str,
    bkgd_subtract: bool,
    verbose: bool = False,
):
    sk = resolve_bundled_syst_name(evtdf, syst_name)
    univ, cv = get_univ_rates(
        cov_type,
        evtdf=evtdf,
        nudf=None,
        bkgd_subtract=bkgd_subtract,
        var_config=var_config,
        syst_name=sk,
        n_univ=n_univ,
        verbose=verbose,
    )
    pack = sanitize_matrix_pack(get_covariance_matrix(univ, cv))
    return pack, univ, cv


def cov_pack_multiplied(
    evtdf,
    var_config,
    knobs: Sequence[str],
    n_univ: int,
    cov_type: str,
    bkgd_subtract: bool,
    bundled_tag: str = "Flux",
    verbose: bool = False,
):
    df_w = inject_multiplied_mc_knob_weights(
        evtdf, knobs, bundled_tag=bundled_tag, n_univ=n_univ, verbose=verbose
    )
    univ, cv = get_univ_rates(
        cov_type,
        evtdf=df_w,
        nudf=None,
        bkgd_subtract=bkgd_subtract,
        var_config=var_config,
        syst_name=("mc", bundled_tag),
        n_univ=n_univ,
        verbose=verbose,
    )
    pack = sanitize_matrix_pack(get_covariance_matrix(univ, cv))
    return pack, univ, cv


# ---------------------------------------------------------------------------
# Paths / save
# ---------------------------------------------------------------------------

def default_syst_disk_root(save_fig_base_dir: str, syst_name: str, tag: Optional[str] = None) -> str:
    today = tag or datetime.now().strftime("%Y%m%d")
    root = path.join(save_fig_base_dir, f"systematics-notebook-{syst_name.lower()}-{today}")
    makedirs(root, exist_ok=True)
    return root


def univ_hist_save_name(syst_disk_root: str, sn: str, vsn: str, tag: str) -> str:
    cat_dir = category_out_dir(syst_disk_root, sn)
    makedirs(cat_dir, exist_ok=True)
    return path.join(cat_dir, f"{vsn}-{tag}-universes")


def matrix_save_name(syst_disk_root: str, sn: str, vsn: str, tag: str) -> str:
    cat_dir = category_out_dir(syst_disk_root, sn)
    makedirs(cat_dir, exist_ok=True)
    return path.join(cat_dir, f"{vsn}-{tag}-matrix")


def write_multisim_npzs_and_manifest(
    syst_dict: dict,
    syst_disk_root: str,
    *,
    syst_name: str,
    n_univ: int,
    bkgd_subtract: bool,
    cov_type: str,
    knobs: Sequence[str],
) -> str:
    """Save NPZ(s) + ``covariance_manifest.json``; return manifest path."""
    log("writing NPZ(s) ...")
    save_neutrino_multisim_npzs(syst_dict, syst_disk_root)
    if syst_name in SYST_DISK_SUBFILE and syst_dict.get(syst_name):
        sub, fname = SYST_DISK_SUBFILE[syst_name]
        npz_path = path.join(normalized_root(syst_disk_root), sub, fname)
        log(f"  wrote {npz_path}")

    categories_present = [syst_name] if syst_dict.get(syst_name) else []
    knob_sidecars = [f"{syst_name}_by_knob"] if syst_dict.get(f"{syst_name}_by_knob") else []
    var_names = set(syst_dict.get(syst_name, {}).keys())
    manifest = {
        "schema": "numucc_multisim_covariance_v1",
        "description": (
            f"Notebook {syst_name} multisim: "
            + (
                "bundled universes only."
                if syst_name == "MCstat"
                else "per-knob packs under *_by_knob; combined entry uses multiplied per-universe weights."
            )
        ),
        "mc_df_stage": "final",
        "var_set": "notebook",
        "syst_types_run": [syst_name],
        "categories_present": categories_present,
        "knob_breakdown_sidecars": knob_sidecars,
        "variables": sorted(var_names),
        "knob_total_method": (
            {syst_name: "multiply_univ_weights"} if syst_name in ("Flux", "G4") and categories_present else {}
        ),
        "n_univ": n_univ,
        "bkgd_subtract": bkgd_subtract,
        "cov_type": cov_type,
        "knobs": {syst_name: list(knobs)} if knobs else {},
    }
    manifest_path = path.join(normalized_root(syst_disk_root), "covariance_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    log("wrote " + manifest_path)
    return manifest_path


# ---------------------------------------------------------------------------
# Flux asymmetry (integrated side-by-side)
# ---------------------------------------------------------------------------

def flux_knob_display_label(knob_name: str) -> str:
    """Compact panel titles for BNB flux multisim knobs."""
    kn = str(knob_name)
    base = kn[:-5] if kn.endswith("_Flux") else kn
    return _FLUX_KNOB_DISPLAY.get(base, base.replace("_", " "))


def plot_univ_hists_on_ax(
    ax,
    univ_events,
    cv_events,
    var_config,
    title: str = "",
    show_legend: bool = False,
):
    """Universe histograms on one axes + mean (red dotted) and CV (black)."""
    univ_events = np.asarray(univ_events, dtype=float)
    cv_events = np.asarray(cv_events, dtype=float)
    n_u = univ_events.shape[0]
    mean_events = np.mean(univ_events, axis=0)

    if n_u > 10:
        colors = ["#FDE725FF", "#1F968BFF", "#440154FF"]
        sorted_univs = np.sort(univ_events, axis=0)
        n_68 = int(0.68 * n_u)
        start_68 = (n_u - n_68) // 2
        end_68 = start_68 + n_68
        n_95 = int(0.95 * n_u)
        start_95 = (n_u - n_95) // 2
        end_95 = start_95 + n_95
        segs = [
            (range(start_68, end_68), colors[0], "Universe (68%)", False),
            (range(start_95, end_95), colors[1], "Universe (95%)", True),
            (
                (i for i in range(n_u) if i not in range(start_95, end_95)),
                colors[2],
                "Universe (100%)",
                False,
            ),
        ]
        plotted: set = set()
        for r, color, label, skip_68 in segs:
            for i in r:
                if skip_68 and i in range(start_68, end_68):
                    continue
                show_label = label if label not in plotted else None
                ax.hist(
                    var_config.bin_centers,
                    bins=var_config.bins,
                    weights=sorted_univs[i],
                    histtype="step",
                    color=color,
                    alpha=0.7,
                    label=show_label,
                )
                plotted.add(label)
    else:
        for i in range(n_u):
            ax.hist(
                var_config.bin_centers,
                bins=var_config.bins,
                weights=univ_events[i],
                histtype="step",
                color="gray",
                label="Universe" if i == 0 else None,
            )

    ax.hist(
        var_config.bin_centers,
        bins=var_config.bins,
        weights=mean_events,
        histtype="step",
        color="red",
        linestyle=":",
        linewidth=2.0,
        label="Mean Value",
    )
    ax.hist(
        var_config.bin_centers,
        bins=var_config.bins,
        weights=cv_events,
        histtype="step",
        color="k",
        label="Central Value",
    )
    ax.set_xlim(var_config.bins[0], var_config.bins[-1])
    ax.set_title(title, fontsize=9)
    if getattr(var_config, "var_save_name", "") == "integrated":
        ax.set_xticks([])
    if show_legend:
        ax.legend(frameon=False, fontsize=7, loc="upper right")


def collect_flux_asymmetry_panels(
    mc_evt_df: pd.DataFrame,
    knobs: Sequence[str],
    *,
    var_config=None,
    n_univ: int = 100,
    cov_type: str = "rate",
    bkgd_subtract: bool = True,
    verbose: bool = False,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    """Per-knob + multiplied Flux-total universe rates for the asymmetry figure."""
    if var_config is None:
        var_config = VariableConfig.all_events()
    panel_specs: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for iknob, knob in enumerate(tqdm(knobs, desc="Flux knobs (integrated)")):
        log(f"knob [{iknob + 1}/{len(knobs)}] {knob}")
        univ, cv = get_univ_rates(
            cov_type,
            evtdf=mc_evt_df,
            nudf=None,
            bkgd_subtract=bkgd_subtract,
            var_config=var_config,
            syst_name=("mc", knob),
            n_univ=n_univ,
            verbose=verbose,
        )
        panel_specs.append(
            (flux_knob_display_label(knob), np.asarray(univ, float), np.asarray(cv, float))
        )

    log("combined Flux (multiplied per-universe weights) ...")
    df_tot = inject_multiplied_mc_knob_weights(
        mc_evt_df, knobs, bundled_tag="Flux", n_univ=n_univ, verbose=verbose
    )
    univ_tot, cv_tot = get_univ_rates(
        cov_type,
        evtdf=df_tot,
        nudf=None,
        bkgd_subtract=bkgd_subtract,
        var_config=var_config,
        syst_name=("mc", "Flux"),
        n_univ=n_univ,
        verbose=verbose,
    )
    panel_specs.append(
        ("Flux (total)", np.asarray(univ_tot, float), np.asarray(cv_tot, float))
    )
    log(f"computed {len(panel_specs)} panels")
    return panel_specs


def draw_flux_asymmetry_figure(
    panel_specs: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    var_config,
    *,
    save_path: Optional[str] = None,
    panel_w: float = 0.675,
    fig_h: float = 3.6,
    dpi: int = 140,
    fig_ext: str = ".png",
):
    """Side-by-side integrated universe hists; y-lim from Flux (total) panel."""
    n_panels = len(panel_specs)
    _, univ_tot_plot, cv_tot_plot = panel_specs[-1]
    y_vals = np.concatenate([
        univ_tot_plot.ravel(),
        cv_tot_plot.ravel(),
        np.mean(univ_tot_plot, axis=0).ravel(),
    ])
    y_lo, y_hi = float(np.min(y_vals)), float(np.max(y_vals))
    y_pad = 0.06 * (y_hi - y_lo) if y_hi > y_lo else 0.02 * max(abs(y_hi), 1.0)
    ylim = (y_lo - y_pad, y_hi + y_pad)

    fig, axes = plt.subplots(
        1, n_panels, figsize=(panel_w * n_panels, fig_h),
        sharey=True, constrained_layout=True,
    )
    if n_panels == 1:
        axes = [axes]

    legend_handles = legend_labels = None
    for i, (ax, (label, univ, cv)) in enumerate(zip(axes, panel_specs)):
        plot_univ_hists_on_ax(
            ax, univ, cv, var_config, title=label, show_legend=(i == n_panels - 1),
        )
        ax.set_ylim(ylim)
        if i == n_panels - 1:
            legend_handles, legend_labels = ax.get_legend_handles_labels()
            if ax.get_legend() is not None:
                ax.get_legend().remove()

    axes[0].set_ylabel("Events")
    offset = axes[0].yaxis.get_offset_text()
    offset.set_x(-0.22)
    offset.set_ha("right")
    fig.supxlabel("Integrated")
    fig.legend(
        legend_handles, legend_labels,
        loc="upper center", bbox_to_anchor=(0.5, 1.08),
        ncol=5, frameon=False, fontsize=8,
    )
    if save_path:
        out = save_path if str(save_path).endswith(fig_ext) else save_path + fig_ext
        makedirs(path.dirname(out) or ".", exist_ok=True)
        fig.savefig(out, bbox_inches="tight", dpi=dpi)
        log("saved " + out)
    return fig, axes
