#!/usr/bin/env python3
"""
Xsec category-totals summary plots with DENT unisim overlaid.

Mirrors ``systematics-summary.ipynb`` ``plot_category_totals_single_panel`` with
``genie_kind="xsec"`` (filtered combined GENIE **xsec** knobs — not ``genie_rate``),
and adds a **DENT** curve from matched CV/DENT ``sel_mup`` histograms:

    unc_DENT[%] = 100 * |DENT - CV| / CV

Total syst. folds DENT in as a fully correlated unisim
(``cov_frac += outer(frac, frac)``).

Example:
    python plot_xsec_summary_with_dent.py
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_KNOBS, GENIE_GROUP_ORDER
from analysis_village.numucc_1p0pi.scripts.dent_compare import load_hists
from analysis_village.numucc_1p0pi.syst_category_summary import (
    CAT_COSMICS,
    CAT_DETECTOR,
    CAT_FLUX,
    CAT_G4,
    CAT_GENIE_RATE,
    CAT_GENIE_XSEC,
    CAT_MCSTAT,
    CAT_NTARGETS,
    CAT_POT,
    NTARGETS_FRAC_UNC_PCT,
    POT_FRAC_UNC_PCT,
    frac_unc_pct,
    frac_weights_for_plot,
    load_category_syst_summary,
    sum_cov_frac_matrices,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import category_summary_npz_path
from analysis_village.numucc_1p0pi.utils import dpi as DEFAULT_DPI
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

SYST_DISK_ROOT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final")
DENT_CACHE = SYST_DISK_ROOT / "DENT" / "cache" / "dent_sel_mup_hists.pkl"
GENIE_PKL = SYST_DISK_ROOT / "GENIE" / "cov_mat_dict.pkl"
OUT_DIR = SYST_DISK_ROOT / "DENT-on-xsec-summary" / "plots"

VARS: Sequence[VariableConfig] = (
    VariableConfig.tki_del_Tp(),
    VariableConfig.tki_del_alpha(),
    VariableConfig.tki_del_phi(),
)

CATEGORY_PLOT_ORDER = (
    "Flux",
    "GENIE",
    "G4",
    "Detector",
    "DENT",
    "Exposure",
    "Targets",
    "Cosmics",
    "MC stat.",
)

# Non-GENIE categories from the exported category summary NPZ.
CATEGORY_KEYS = {
    "Flux": CAT_FLUX,
    "G4": CAT_G4,
    "Detector": CAT_DETECTOR,
    "Cosmics": CAT_COSMICS,
    "MC stat.": CAT_MCSTAT,
    "Exposure": CAT_POT,
    "Targets": CAT_NTARGETS,
}

BREAKDOWN_FIGSIZE = (8.0, 6.0)
BREAKDOWN_FIG_DPI = 100
BREAKDOWN_SAVE_DPI = int(DEFAULT_DPI)
LEGEND_NCOL = 3


# ---------------------------------------------------------------------------
# GENIE combined xsec (same conventions as systematics-summary.ipynb)
# ---------------------------------------------------------------------------
def _genie_knob_excluded(kn: str, mode: str) -> bool:
    """Mirror ``genie_knob_excluded_from_combined_breakdown`` in the summary notebook."""
    if mode == "Other" and (kn.endswith("_pi") or kn.endswith("_N")):
        return True
    if mode == "Ar23p":
        if "D_ZExp" in kn:
            return True
        if "q0bin5" in kn:
            return True
        if "EDepFSI_DecayAngMEC" in kn:
            return True
        if "EDepFSI_NormCCMEC" in kn:
            return True
    return False


def _genie_group_key(knob: str) -> str:
    """Mirror ``genie_combined_breakdown_group_key``."""
    import re

    kn = str(knob)
    if re.search(r"_b\d+$", kn) and "ZExp" in kn:
        return "ZExp"
    m = re.search(r"_dial_\d+$", kn)
    if m:
        return kn[: m.start()].rsplit("_", 1)[-1]
    m = re.search(r"_q0bin\d+$", kn)
    if m:
        prefix = kn[: m.start()]
        if "Martini" in prefix:
            return "MEC Martini"
        if "Valenica" in prefix:
            return "MEC Valencia"
        return prefix.rsplit("_", 1)[-1]
    m = re.search(r"bin\d+$", kn)
    if m:
        prefix = kn[: m.start()]
        if "Martini" in prefix:
            return "MEC Martini"
        if "Valenica" in prefix:
            return "MEC Valencia"
        return prefix.rsplit("_", 1)[-1] if prefix else kn
    return kn


def genie_combined_xsec_cov(genie_blob: Mapping, var_name: str) -> np.ndarray:
    """Filtered sum of GENIE **xsec** (non-``_rate``) knobs — notebook xsec total."""
    gd = genie_blob.get(var_name)
    if not isinstance(gd, dict):
        raise KeyError(f"GENIE pickle missing variable {var_name!r}")

    xsec_parts: Dict[str, np.ndarray] = {}
    for mode in GENIE_GROUP_ORDER:
        for kn in GENIE_GROUP_KNOBS.get(mode, []):
            if _genie_knob_excluded(kn, mode):
                continue
            if kn not in gd:
                continue
            mat = gd[kn]
            if not isinstance(mat, np.ndarray) or np.asarray(mat).ndim != 2:
                continue
            gkey = _genie_group_key(kn)
            arr = np.asarray(mat, dtype=np.float64)
            xsec_parts[gkey] = arr if gkey not in xsec_parts else xsec_parts[gkey] + arr

    if not xsec_parts:
        if "genie" in gd:
            return np.asarray(gd["genie"], dtype=np.float64)
        raise KeyError(f"No GENIE xsec knobs for {var_name!r}")

    total = sum_cov_frac_matrices(xsec_parts.values())
    if total is None:
        raise KeyError(f"Empty GENIE xsec sum for {var_name!r}")
    return np.asarray(total, dtype=np.float64)


def _step_colors(n: int) -> List:
    tab10 = list(plt.cm.tab10.colors)
    set2 = list(plt.cm.Set2.colors)
    pool = tab10 + set2
    if n <= len(pool):
        return pool[:n]
    out = list(pool)
    k = 0
    while len(out) < n:
        out.append(pool[k % len(pool)])
        k += 1
    return out[:n]


def dent_frac_shift(cv: np.ndarray, dent: np.ndarray) -> np.ndarray:
    cv = np.asarray(cv, dtype=np.float64)
    dent = np.asarray(dent, dtype=np.float64)
    out = np.zeros_like(cv, dtype=np.float64)
    m = np.abs(cv) > 0
    out[m] = (dent[m] - cv[m]) / cv[m]
    return out


def dent_unisim_cov_frac(cv: np.ndarray, dent: np.ndarray) -> np.ndarray:
    f = dent_frac_shift(cv, dent)
    return np.outer(f, f)


def style_uncertainty_axis(ax, var_config, pct_series_list) -> None:
    flat = [
        np.asarray(s, dtype=np.float64).ravel()
        for s in pct_series_list
        if s is not None and len(s)
    ]
    ymax = max(
        (float(np.nanmax(np.where(np.isfinite(x), x, 0.0))) for x in flat),
        default=0.0,
    )
    ax.set_xlim(float(var_config.bins[0]), float(var_config.bins[-1]))
    ax.set_ylim(bottom=0.0, top=max(ymax * 1.25, 1.0))
    xlab = (
        var_config.var_labels[1]
        if getattr(var_config, "var_labels", None)
        else var_config.var_save_name
    )
    ax.set_xlabel(xlab, fontsize=22)
    ax.set_ylabel("Uncertainty [%]", fontsize=22)
    ax.grid(which="major", linestyle="-", linewidth=0.7, alpha=0.7)
    ax.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.minorticks_on()
    ax.legend(
        loc="upper center",
        ncol=LEGEND_NCOL,
        fontsize=10,
        frameon=True,
        bbox_to_anchor=(0.5, 0.98),
    )


def plot_xsec_summary_with_dent(
    var_config,
    summary,
    dent_hists: dict,
    genie_blob: Mapping,
    *,
    out_dir: Path,
) -> Path:
    vsn = var_config.var_save_name
    bins = np.asarray(var_config.bins, dtype=float)
    bc = np.asarray(var_config.bin_centers, dtype=float)
    pack = summary["by_var"][vsn]

    category_pct: Dict[str, np.ndarray] = {}
    covs: Dict[str, np.ndarray] = {}

    for lab, key in CATEGORY_KEYS.items():
        if key not in pack["categories"]:
            continue
        cf = np.asarray(pack["categories"][key]["cov_frac"], dtype=np.float64)
        covs[lab] = cf
        category_pct[lab] = frac_weights_for_plot(cf, var_config)

    n = len(bc)
    covs.setdefault(
        "Exposure",
        np.diag(np.full(n, (POT_FRAC_UNC_PCT / 100.0) ** 2, dtype=np.float64)),
    )
    covs.setdefault(
        "Targets",
        np.diag(np.full(n, (NTARGETS_FRAC_UNC_PCT / 100.0) ** 2, dtype=np.float64)),
    )
    category_pct.setdefault("Exposure", np.full(n, POT_FRAC_UNC_PCT, dtype=float))
    category_pct.setdefault("Targets", np.full(n, NTARGETS_FRAC_UNC_PCT, dtype=float))

    # GENIE **xsec** (filtered combined), never genie_rate.
    genie_xsec = genie_combined_xsec_cov(genie_blob, vsn)
    genie_rate = np.asarray(pack["categories"][CAT_GENIE_RATE]["cov_frac"], dtype=np.float64)
    genie_xsec_pct = frac_weights_for_plot(genie_xsec, var_config)
    genie_rate_pct = frac_weights_for_plot(genie_rate, var_config)
    if float(np.nanmean(genie_xsec_pct)) > float(np.nanmean(genie_rate_pct)):
        raise RuntimeError(
            f"{vsn}: GENIE xsec mean ({genie_xsec_pct.mean():.2f}%) > rate "
            f"({genie_rate_pct.mean():.2f}%) — refusing to plot (likely swapped)."
        )
    covs["GENIE"] = genie_xsec
    category_pct["GENIE"] = genie_xsec_pct

    cv = np.asarray(dent_hists["cv"][vsn], dtype=float)
    dent = np.asarray(dent_hists["dent"][vsn], dtype=float)
    if len(cv) != n:
        raise ValueError(f"{vsn}: hist len {len(cv)} != nbins {n}")
    dent_cov = dent_unisim_cov_frac(cv, dent)
    covs["DENT"] = dent_cov
    category_pct["DENT"] = frac_weights_for_plot(dent_cov, var_config)

    total_cov = sum_cov_frac_matrices(covs.values())
    if total_cov is None:
        raise RuntimeError(f"{vsn}: empty total covariance")
    total_pct = frac_weights_for_plot(total_cov, var_config)

    fig, ax = plt.subplots(figsize=BREAKDOWN_FIGSIZE, dpi=BREAKDOWN_FIG_DPI)
    plot_labels = [lab for lab in CATEGORY_PLOT_ORDER if lab in category_pct]
    colors = _step_colors(len(plot_labels))
    pct_list: List[np.ndarray] = []
    for lab, color in zip(plot_labels, colors):
        w = category_pct[lab]
        pct_list.append(w)
        ax.hist(
            bc, bins=bins, weights=w, histtype="step", linewidth=2, color=color, label=lab,
        )
    pct_list.append(total_pct)
    ax.hist(
        bc, bins=bins, weights=total_pct, histtype="step", linewidth=2, color="k",
        label="Total syst.",
    )
    style_uncertainty_axis(ax, var_config, pct_list)
    fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"syst_break_totals_xsec_with_DENT__{vsn}.png"
    fig.savefig(out_path, dpi=BREAKDOWN_SAVE_DPI, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(
        f"{vsn}: GENIE xsec max={genie_xsec_pct.max():.2f}% "
        f"(rate would be {genie_rate_pct.max():.2f}%)  "
        f"DENT max={category_pct['DENT'].max():.2f}%  "
        f"Total max={total_pct.max():.2f}%  → {out_path}",
        flush=True,
    )
    return out_path


def main() -> int:
    summary_path = category_summary_npz_path(str(SYST_DISK_ROOT))
    for p in (summary_path, DENT_CACHE, GENIE_PKL):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

    summary = load_category_syst_summary(summary_path)
    dent_payload = load_hists(str(DENT_CACHE))
    dent_hists = dent_payload["hists"]
    with open(GENIE_PKL, "rb") as fh:
        genie_blob = pickle.load(fh)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for vc in VARS:
        vsn = vc.var_save_name
        if vsn not in summary["by_var"]:
            print(f"skip {vsn}: not in category summary", flush=True)
            continue
        if vsn not in dent_hists.get("cv", {}) or vsn not in dent_hists.get("dent", {}):
            print(f"skip {vsn}: missing in DENT hist cache", flush=True)
            continue
        plot_xsec_summary_with_dent(
            vc, summary, dent_hists, genie_blob, out_dir=OUT_DIR,
        )

    print(f"Done → {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
