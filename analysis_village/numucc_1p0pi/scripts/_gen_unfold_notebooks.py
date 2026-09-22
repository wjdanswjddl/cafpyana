#!/usr/bin/env python3
"""One-shot generator for PRL Product B unfold notebooks. Not part of the runtime path."""

from __future__ import annotations

import json
from pathlib import Path

NB_DIR = Path(__file__).resolve().parents[1] / "notebooks"


def md(src: str) -> dict:
    lines = src.split("\n")
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [ln + "\n" for ln in lines],
    }


def code(src: str) -> dict:
    lines = src.split("\n")
    out = []
    for i, line in enumerate(lines):
        if i < len(lines) - 1:
            out.append(line + "\n")
        elif line or not lines:
            out.append(line if line else "\n")
    if not out:
        out = [""]
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": out,
    }


def notebook(cells: list) -> dict:
    return {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "venv_py310_cafpyana",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "cells": cells,
    }


def write_nb(name: str, cells: list) -> Path:
    path = NB_DIR / name
    path.write_text(json.dumps(notebook(cells), indent=1))
    print("wrote", path)
    return path


# =============================================================================
# unfolding-prepare.ipynb
# =============================================================================

PREP = []

PREP.append(
    md(
        """# Unfolding prepare (PRL Product B) — DFs, overlays, response matrices

**Step 1** of the Product B data-release unfold:

1. Load beam-quality data + `sel_mup` MC (`evt`, `mcnu`, `hdr`) — Sep-1 `fvfix` campaign
2. Recompute data–MC overlays and **assert** against `PRL/data_mc_overlays/productB_sel_mup/counts_report.npz`
3. Build signal efficiency + response matrices from `mcnu` / `evt`; plot and save under `PRL/response_matrices/`

**Next:** [`unfolding.ipynb`](unfolding.ipynb) (closure + data Wiener-SVD).

**Legacy Gen1 (May recovered cov):** [`unfolding-legacy-gen1.ipynb`](unfolding-legacy-gen1.ipynb)

## Physics notes (do not silently "fix")

- Product B Sep-1 sample ≠ May Gen1 cache — unfolded χ² vs generators will differ from 34.5/12.
- Overlay bkg = topology non-signal layers (`mc_background` in counts report). Cosmics also enter via CategorySummary; Product B overlay has `mc_other=0` (no intime).
- Migration uses `wgt_sel_truth` (canonical). Do **not** copy the archive fake-data `wgt_sel_reco` path.
- Response formula unchanged: `R = get_response_matrix(reco_vs_true, eff)` with column sums = efficiency."""
    )
)

PREP.append(code("%load_ext autoreload\n%autoreload 2"))

PREP.append(
    code(
        r'''import sys
import json
import gc
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path
from os import makedirs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(SCRIPTS))

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from pyanalib.split_df_helpers_new import dfs_from_dir
from analysis_village.numucc_1p0pi.final_selected_evt_vars import CORE_SELECTED_EVT_VARIABLE_CONFIGS
from analysis_village.numucc_1p0pi.utils import (
    get_response_matrix,
    get_topo_category,
    plot_heatmap,
    signal_hists,
    strip_pot_from_ylabel,
    format_pot_corner_text,
    fig_ext,
    dpi,
    get_category_summary_syst_unc,
)
from analysis_village.numucc_1p0pi.selected_xsec_overlay_hist import (
    build_overlay_histdata_map,
    plot_overlay_counts_map,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import SYST_DISK_ENV

import selected_xsec_overlay as sxo
import data_mc_overlay_products as dmo
'''
    )
)

PREP.append(
    code(
        r'''# ---- Product B paths (same as data_mc_overlay_products.py) ----
PRL_ROOT = Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL")
OVERLAY_DIR = PRL_ROOT / "data_mc_overlays/productB_sel_mup"
COUNTS_NPZ = OVERLAY_DIR / "counts_report.npz"
RESP_DIR = PRL_ROOT / "response_matrices"
FIG_DIR = RESP_DIR / "plots"

DATA_DIR = dmo.PRODUCT_B_DATA_DIR
DATA_FN = dmo.PRODUCT_B_DATA_FN
MC_DIR = dmo.PRODUCT_B_MC_DIR
MC_FN = dmo.PRODUCT_B_MC_FN
EXPECTED_DATA_N = dmo.EXPECTED_DATA_N
SYST_ROOT = dmo.PRODUCT_B_SYST_ROOT

# Overlay vars match Product B counts report (no integrated)
OVERLAY_VAR_CONFIGS = [
    vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS if vc.var_save_name != "integrated"
]
# Unfold / response vars include integrated
UNFOLD_VAR_CONFIGS = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
BREAKDOWN_TYPES = ("topology", "genie_sb")

makedirs(RESP_DIR, exist_ok=True)
makedirs(FIG_DIR, exist_ok=True)
makedirs(OVERLAY_DIR, exist_ok=True)
os.environ[SYST_DISK_ENV] = SYST_ROOT

print("data dir :", DATA_DIR)
print("mc dir   :", MC_DIR)
print("counts   :", COUNTS_NPZ, "exists=", COUNTS_NPZ.is_file())
print("resp out :", RESP_DIR)
print("n overlay vars:", len(OVERLAY_VAR_CONFIGS), " n unfold vars:", len(UNFOLD_VAR_CONFIGS))
'''
    )
)

PREP.append(md("## 1. Load data (beam quality) and MC (`evt`, `mcnu`, `hdr`)"))

PREP.append(
    code(
        r'''# Align overlay module knobs with Product B
sxo.APPLY_BEAM_QUALITY = True
sxo.LOAD_SYST = True
sxo.SYST_DISK_ROOT = SYST_ROOT
sxo.SAVE_FIG = True
sxo.PLOT = True
sxo.APPROVAL = ""
sxo.TEXTCHI2 = True

print("Loading data (live beam-quality cuts)...")
data_evt, data_hdr = sxo.load_data_sample(DATA_DIR, DATA_FN)
n_evt_good = int(len(data_evt))
assert n_evt_good == EXPECTED_DATA_N, (
    f"n_evt_good={n_evt_good} != expected {EXPECTED_DATA_N}; check DF dir {DATA_DIR}"
)
print(f"beam-quality OK: n_evt_good={n_evt_good}")

print("Loading MC evt+hdr+mcnu...")
mc_dfs = dfs_from_dir(
    MC_DIR,
    filename_str=MC_FN,
    keys2load=["hdr", "evt", "mcnu"],
    n_max_concat=sxo.N_MAX_CONCAT,
)
mc_evt = mc_dfs["evt"]
mc_hdr = mc_dfs["hdr"]
mc_nu = mc_dfs["mcnu"]
if "mc" in mc_evt.columns.get_level_values(0):
    mc_evt.loc[mc_evt.mc.iscc.isna(), ("mc", "iscc")] = 999
if "mc" in mc_nu.columns.get_level_values(0):
    mc_nu.loc[mc_nu.mc.iscc.isna(), ("mc", "iscc")] = 999

print(f"MC: evt={len(mc_evt):,}  hdr={len(mc_hdr):,}  mcnu={len(mc_nu):,}")

pot_label_raw = sxo.setup_pot_weights(mc_evt, mc_hdr, data_evt, data_hdr)
mc_scale = float(mc_evt["pot_weight"].iloc[0])
mc_nu["pot_weight"] = mc_scale * np.ones(len(mc_nu))

if "topo_categ" not in mc_evt.columns:
    mc_evt = mc_evt.copy()
    mc_evt.loc[:, "topo_categ"] = get_topo_category(mc_evt)
if "topo_categ" not in mc_nu.columns:
    mc_nu = mc_nu.copy()
    mc_nu.loc[:, "topo_categ"] = get_topo_category(mc_nu)

data_tot_pot = float(data_hdr["pot"].sum())
mc_tot_pot = float(mc_hdr["pot"].sum())
print(f"data_tot_pot={data_tot_pot:.6e}  mc_tot_pot={mc_tot_pot:.6e}  mc_scale={mc_scale:.6e}")
print("pot_label:", pot_label_raw)
'''
    )
)

PREP.append(md("## 2. Recompute data–MC overlays and assert vs `counts_report.npz`"))

PREP.append(
    code(
        r'''def _product_b_syst_cov(var_config):
    try:
        _unc, cov = get_category_summary_syst_unc(
            var_config, syst_kind="rate", syst_disk_root=SYST_ROOT,
        )
        if cov is None or not np.any(np.isfinite(cov)):
            return None
        return cov
    except Exception as ex:
        print(f"  [syst] skip {var_config.var_save_name}: {ex}")
        return None


# Snapshot reference counts BEFORE export overwrites the NPZ
assert COUNTS_NPZ.is_file(), f"missing reference counts: {COUNTS_NPZ}"
ref_before = {k: np.asarray(v) for k, v in np.load(COUNTS_NPZ).items()}
print(f"loaded reference counts: {len(ref_before)} arrays from {COUNTS_NPZ}")

print("Building overlay histdata from live DFs...")
histdata_map = build_overlay_histdata_map(
    OVERLAY_VAR_CONFIGS,
    BREAKDOWN_TYPES,
    mc_df=mc_evt,
    data_df=data_evt,
)

pot_text = format_pot_corner_text(pot_label_raw)
pot_label = strip_pot_from_ylabel(pot_label_raw) or "Events / Bin"

plot_overlay_counts_map(
    histdata_map,
    OVERLAY_VAR_CONFIGS,
    BREAKDOWN_TYPES,
    pot_label=pot_label,
    out_dir=str(OVERLAY_DIR),
    get_syst=_product_b_syst_cov,
    ax_ylim_ratio=sxo.AX_YLIM_RATIO,
    ratio=sxo.RATIO,
    textloc=sxo.TEXTLOC,
    approval="",
    pot_text=pot_text,
    save_fig=True,
    plot=True,
    textchi2=True,
)

items = []
for vc in OVERLAY_VAR_CONFIGS:
    key = (vc.var_save_name, "topology")
    if key in histdata_map:
        items.append((vc.var_save_name, histdata_map[key]))

# Build live count arrays without overwriting the reference yet
from analysis_village.numucc_1p0pi.categories import topology_labels

overlay_counts = {}
mismatches = []
for vc in OVERLAY_VAR_CONFIGS:
    slug = vc.var_save_name
    hd = histdata_map[(slug, "topology")]
    data = np.asarray(hd.data_hist, dtype=float).ravel()
    bins = np.asarray(hd.bins, dtype=float)
    mh = np.asarray(hd.mc_hist, dtype=float)
    sig = mh[-1].copy()
    bkg = mh.sum(axis=0) - sig
    mc_tot = mh.sum(axis=0)
    overlay_counts[slug] = {
        "bins": bins,
        "data": data,
        "mc_signal": sig,
        "mc_background": bkg,
        "mc_total": mc_tot,
    }
    for suffix, arr in (
        ("bins", bins),
        ("data", data),
        ("mc_signal", sig),
        ("mc_background", bkg),
        ("mc_total", mc_tot),
    ):
        key = f"{slug}__{suffix}"
        if key not in ref_before:
            mismatches.append(f"missing in ref: {key}")
            continue
        a = np.asarray(ref_before[key], dtype=float)
        b = np.asarray(arr, dtype=float)
        if a.shape != b.shape:
            mismatches.append(f"{key}: shape {a.shape} vs {b.shape}")
            continue
        if not np.allclose(a, b, rtol=1e-5, atol=1e-4):
            maxdiff = float(np.max(np.abs(a - b)))
            mismatches.append(f"{key}: max|Δ|={maxdiff:.3e}")

if mismatches:
    raise AssertionError(
        "Overlay counts disagree with counts_report.npz:\n  - "
        + "\n  - ".join(mismatches[:40])
        + (f"\n  ... ({len(mismatches)} total)" if len(mismatches) > 40 else "")
    )

print(f"ASSERT OK: {len(OVERLAY_VAR_CONFIGS)} vars match {COUNTS_NPZ.name}")

# Refresh on-disk report (same content; keeps provenance timestamp)
live_counts_npz = dmo.export_counts_report(
    str(OVERLAY_DIR),
    product="B",
    histdata_items=items,
    extra_meta={
        "data_dir": DATA_DIR,
        "mc_dir": MC_DIR,
        "syst_disk_root": SYST_ROOT,
        "syst_kind": "rate",
        "n_evt_good": n_evt_good,
        "expected_n_evt_good": EXPECTED_DATA_N,
        "source": "unfolding-prepare.ipynb",
    },
)
print("refreshed counts report:", live_counts_npz)
'''
    )
)

PREP.append(
    md(
        """## 3. Signal efficiency + response matrices

Uses `signal_hists(evt, mcnu)` → migration `histogram2d(truth, reco; wgt_sel_truth)` →
`eff = nevts_sel_truth / nevts_allmc` → `get_response_matrix`.

Efficiency curve: final-selection signal efficiency vs truth (denom = `mcnu` signal),
styled like the historical `event_selection.ipynb` / `plot_efficiency` twin-axis idea."""
    )
)

PREP.append(
    code(
        r'''def build_response(sig, bins):
    """Same recipe as Gen1 / archive unfolding-data (wgt_sel_truth)."""
    if len(bins) == 2:
        reco_vs_true = np.array([[float(np.sum(sig["nevts_sel_truth"]))]], dtype=float)
    else:
        reco_vs_true, _, _ = np.histogram2d(
            sig["var_sel_truth"],
            sig["var_sel_reco"],
            weights=sig["wgt_sel_truth"],
            bins=[bins, bins],
        )
    eff = np.asarray(sig["nevts_sel_truth"], dtype=float) / np.asarray(
        sig["nevts_allmc"], dtype=float
    )
    eff = np.where(np.isfinite(eff), eff, 0.0)
    response = get_response_matrix(reco_vs_true, eff)
    return reco_vs_true, eff, response


def plot_signal_efficiency(var_config, eff, nevts_sel, nevts_allmc, save_name):
    """Final-selection efficiency vs truth (response denominator)."""
    bins = var_config.bins
    centers = var_config.bin_centers
    fig, ax = plt.subplots()
    ax_eff = ax.twinx()
    ax.hist(
        centers, bins=bins, weights=nevts_allmc, histtype="step",
        color="C0", label="All signal (mcnu)", linewidth=1.5,
    )
    ax.hist(
        centers, bins=bins, weights=nevts_sel, histtype="step",
        color="C1", label="Selected signal", linewidth=1.5,
    )
    ax_eff.errorbar(centers, eff, fmt="ko-", markersize=4, label="Efficiency")
    ax.set_xlabel(var_config.var_labels[0])
    ax.set_ylabel("Events / Bin (POT-scaled)")
    ax_eff.set_ylabel("Efficiency")
    ax_eff.set_ylim(0, min(1.05, max(0.2, float(np.nanmax(eff)) * 1.3 if np.any(eff) else 0.2)))
    ax.set_xlim(bins[0], bins[-1])
    ax.set_title(f"Signal efficiency — {var_config.var_save_name}")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_eff.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=10)
    fig.savefig(save_name + fig_ext, bbox_inches="tight", dpi=dpi)
    plt.show()
    plt.close(fig)


response_pack = {
    "meta": {
        "schema": "prl_productB_response_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "data_dir": DATA_DIR,
        "mc_dir": MC_DIR,
        "n_evt_good": n_evt_good,
        "data_tot_pot": data_tot_pot,
        "mc_tot_pot": mc_tot_pot,
        "mc_pot_scale": mc_scale,
        "counts_report": str(COUNTS_NPZ),
        "syst_disk_root": SYST_ROOT,
    },
    "variables": {},
}

for vc in UNFOLD_VAR_CONFIGS:
    vsn = vc.var_save_name
    print("=" * 60)
    print("response:", vsn)
    sig = signal_hists(
        evtdf=mc_evt,
        nudf=mc_nu,
        var_config=vc,
        return_data=True,
        plot=False,
        mode="reco",
    )
    reco_vs_true, eff, response = build_response(sig, vc.bins)

    if vsn in overlay_counts:
        n_data = overlay_counts[vsn]["data"]
        n_bkg = overlay_counts[vsn]["mc_background"]
    else:
        reco_col = vc.var_evt_reco_col
        n_data, _ = np.histogram(
            data_evt[reco_col], bins=vc.bins, weights=data_evt["pot_weight"]
        )
        mc_bkg = mc_evt[mc_evt.topo_categ != 1]
        n_bkg, _ = np.histogram(
            mc_bkg[reco_col], bins=vc.bins, weights=mc_bkg["pot_weight"]
        )
        n_data = np.asarray(n_data, dtype=float)
        n_bkg = np.asarray(n_bkg, dtype=float)

    n_sel_data = n_data - n_bkg

    plot_signal_efficiency(
        vc,
        eff,
        np.asarray(sig["nevts_sel_truth"], dtype=float),
        np.asarray(sig["nevts_allmc"], dtype=float),
        str(FIG_DIR / f"{vsn}-efficiency"),
    )
    if reco_vs_true.shape[0] > 1:
        plot_heatmap(
            reco_vs_true,
            vc.bins,
            plot_labels=["True", "Reco", f"{vsn} migration"],
            plot=True,
            save_fig=True,
            save_name=str(FIG_DIR / f"{vsn}-reco_vs_true"),
        )
        plot_heatmap(
            response,
            vc.bins,
            plot_labels=["True", "Reco", f"{vsn} response"],
            plot=True,
            save_fig=True,
            save_name=str(FIG_DIR / f"{vsn}-response"),
        )
    else:
        print(f"  integrated: eff={eff}, response={response}")

    response_pack["variables"][vsn] = {
        "bins": np.asarray(vc.bins, dtype=float),
        "bin_centers": np.asarray(vc.bin_centers, dtype=float),
        "reco_vs_true": np.asarray(reco_vs_true, dtype=float),
        "eff": np.asarray(eff, dtype=float),
        "response": np.asarray(response, dtype=float),
        "nevts_allmc": np.asarray(sig["nevts_allmc"], dtype=float),
        "nevts_sel_truth": np.asarray(sig["nevts_sel_truth"], dtype=float),
        "nevts_sel_reco": np.asarray(sig["nevts_sel_reco"], dtype=float),
        "nevts_allsel_reco": np.asarray(sig["nevts_allsel_reco"], dtype=float),
        "n_data": np.asarray(n_data, dtype=float),
        "n_mc_bkg": np.asarray(n_bkg, dtype=float),
        "n_sel_data": np.asarray(n_sel_data, dtype=float),
        "var_labels": list(vc.var_labels),
    }
    print(
        f"  eff mean={float(np.nanmean(eff)):.4f}  "
        f"n_sel_data sum={float(np.sum(n_sel_data)):.1f}  "
        f"nevts_allmc sum={float(np.sum(sig['nevts_allmc'])):.1f}"
    )

print("built responses for", list(response_pack["variables"]))
'''
    )
)

PREP.append(
    code(
        r'''# Save lean NPZ + manifest for unfolding.ipynb
savez_kw = {"meta_json": np.array([json.dumps(response_pack["meta"])], dtype=object)}
for vsn, pack in response_pack["variables"].items():
    for key, arr in pack.items():
        if key == "var_labels":
            savez_kw[f"{vsn}::var_labels"] = np.array(arr, dtype=object)
        else:
            savez_kw[f"{vsn}::{key}"] = np.asarray(arr)

out_npz = RESP_DIR / "response_matrices.npz"
np.savez_compressed(out_npz, **savez_kw)

manifest = {
    "schema": "prl_productB_response_v1",
    "created_utc": response_pack["meta"]["created_utc"],
    "npz_path": str(out_npz),
    "variables": sorted(response_pack["variables"]),
    "meta": response_pack["meta"],
    "key_format": "{var_save_name}::{field}",
    "fields": [
        "bins", "bin_centers", "reco_vs_true", "eff", "response",
        "nevts_allmc", "nevts_sel_truth", "nevts_sel_reco", "nevts_allsel_reco",
        "n_data", "n_mc_bkg", "n_sel_data", "var_labels",
    ],
}
man_path = RESP_DIR / "response_matrices_manifest.json"
with open(man_path, "w") as f:
    json.dump(manifest, f, indent=2)

print("wrote", out_npz)
print("wrote", man_path)
print("Done prepare.")
'''
    )
)

write_nb("unfolding-prepare.ipynb", PREP)


# =============================================================================
# unfolding.ipynb (PRL Product B step 2)
# =============================================================================

UNF = []

UNF.append(
    md(
        """# Unfolding (PRL Product B Wiener-SVD)

**Step 2** of the Product B data-release unfold.

1. Load response matrices from [`unfolding-prepare.ipynb`](unfolding-prepare.ipynb)
2. Load Product B CategorySummary **`total_xsec`** covariance
3. Compute `xsec_unit` from Gen1 ray-traced flux × `FV_split_truncY` × data POT
4. **Closure test** (Asimov / MC): unfold selected MC reco; compare to `A_c @ model`
5. **Data unfold**; save under `PRL/unfolded/`

Wiener-SVD settings (unchanged): `C_type=2`, `Norm_type=0.0`, `stat_scaling=xsec_unit`.

**Legacy Gen1 recovered-cov path:** [`unfolding-legacy-gen1.ipynb`](unfolding-legacy-gen1.ipynb)

## Physics notes

- Cov recipe (archive `unfolding-data`, **not** Gen1 `CovRotation` recovery):
  `Covariance = cov_from_fraccov(total_xsec_frac, nevts_sel_reco) * xsec_unit**2`
- `measured = (n_data - n_mc_bkg) * xsec_unit`, `model = nevts_allmc * xsec_unit`
- Product B `total_xsec` categories: flux, g4, mcstat, detector, cosmics, **genie_xsec**, pot, ntargets
- Do **not** use Gen1 `mean(model/nevts_allmc)` for `xsec_unit`
- Sep-1 Product B ≠ May Gen1 → χ² vs generators will differ from 34.5/12"""
    )
)

UNF.append(code("%load_ext autoreload\n%autoreload 2"))

UNF.append(
    code(
        r'''import sys
import json
import shutil
import warnings
from datetime import datetime, timezone
from pathlib import Path
from os import makedirs

import numpy as np
import matplotlib.pyplot as plt

REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "analysis_village/numucc_1p0pi/scripts"))

warnings.filterwarnings("ignore", category=FutureWarning)

from analysis_village.numucc_1p0pi.constants import M_AR, N_A, RHO
from analysis_village.numucc_1p0pi.final_selected_evt_vars import CORE_SELECTED_EVT_VARIABLE_CONFIGS
from analysis_village.numucc_1p0pi.syst_category_summary import (
    load_category_syst_summary,
    total_cov_frac,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import category_summary_npz_path
from analysis_village.numucc_1p0pi.utils import (
    cov_from_fraccov,
    get_chi2,
    get_integrated_flux,
    plot_heatmap,
    plot_unfolded_result,
    fig_ext,
    dpi,
)
from analysis_village.unfolding.wienersvd import WienerSVD
from analysis_village.flux.raytrace_volume_defs import (
    FV_SPLIT_TRUNCY_BOXES,
    RAYTRACE_VOLUME_LABEL,
)
from unfolding_data import pack_unfold_results
'''
    )
)

UNF.append(
    code(
        r'''PRL_ROOT = Path("/exp/sbnd/data/users/munjung/xsec/numucc_1p0pi/PRL")
RESP_NPZ = PRL_ROOT / "response_matrices/response_matrices.npz"
SYST_ROOT = PRL_ROOT / "systematics/productB_sel_mup"
CAT_NPZ = Path(category_summary_npz_path(str(SYST_ROOT)))
OUT_DIR = PRL_ROOT / "unfolded"
FIG_DIR = OUT_DIR / "plots"
FLUX_FILE = Path("/exp/sbnd/data/users/munjung/flux/SBND_gsimple_raytrace/Gen1.root")

C_TYPE = 2
NORM_TYPE = 0.0

UNFOLD_VAR_CONFIGS = list(CORE_SELECTED_EVT_VARIABLE_CONFIGS)
vc_by = {vc.var_save_name: vc for vc in UNFOLD_VAR_CONFIGS}

makedirs(OUT_DIR, exist_ok=True)
makedirs(FIG_DIR, exist_ok=True)

print("response :", RESP_NPZ, "exists=", RESP_NPZ.is_file())
print("category :", CAT_NPZ, "exists=", CAT_NPZ.is_file())
print("flux     :", FLUX_FILE, "exists=", FLUX_FILE.is_file())
print("out      :", OUT_DIR)
'''
    )
)

UNF.append(md("## 1. Load response pack + Product B `total_xsec`"))

UNF.append(
    code(
        r'''blob = np.load(RESP_NPZ, allow_pickle=True)
meta = json.loads(str(blob["meta_json"][0]))
data_tot_pot = float(meta["data_tot_pot"])

variables = {}
# discover var names from keys
vsns = sorted({k.split("::", 1)[0] for k in blob.files if "::" in k})
for vsn in vsns:
    pack = {}
    for key in blob.files:
        if key.startswith(vsn + "::"):
            field = key.split("::", 1)[1]
            pack[field] = blob[key]
    variables[vsn] = pack

print("meta:", json.dumps({k: meta[k] for k in ("schema", "n_evt_good", "data_tot_pot", "mc_pot_scale")}, indent=2))
print("variables:", vsns)

summary = load_category_syst_summary(str(CAT_NPZ))
print("category summary variables:", sorted(summary["by_var"].keys())[:8], "...")
'''
    )
)

UNF.append(md("## 2. Cross-section unit (flux × FV_split_truncY × data POT)"))

UNF.append(
    code(
        r'''print("Fiducial volume (FV_split_truncY):", RAYTRACE_VOLUME_LABEL["FV_split_truncY"])
v_sbnd = 0.0
for i, box in enumerate(FV_SPLIT_TRUNCY_BOXES):
    dx = box["x_range"][1] - box["x_range"][0]
    dy = box["y_range"][1] - box["y_range"][0]
    dz = box["z_range"][1] - box["z_range"][0]
    v_box = dx * dy * dz
    v_sbnd += v_box
    print(f"  slab {i+1}: -> {v_box:.4e} cm3")
print(f"V_SBND = {v_sbnd:.6e} cm3")

integrated_flux_per_pot = get_integrated_flux(str(FLUX_FILE), plot=False)
integrated_flux = integrated_flux_per_pot * data_tot_pot
n_targets = (RHO * v_sbnd / M_AR) * N_A
XSEC_UNIT = 1.0 / (integrated_flux * n_targets)

flux_info = {
    "flux_file": str(FLUX_FILE),
    "flux_fv": "FV_split_truncY",
    "flux_fv_label": RAYTRACE_VOLUME_LABEL["FV_split_truncY"],
    "volume_cm3": float(v_sbnd),
    "integrated_flux_per_pot": float(integrated_flux_per_pot),
    "data_tot_pot": float(data_tot_pot),
    "integrated_flux": float(integrated_flux),
    "n_targets": float(n_targets),
    "xsec_unit": float(XSEC_UNIT),
}
print(f"integrated flux x POT = {integrated_flux:.6e} /cm2")
print(f"N_targets = {n_targets:.4e}")
print(f"xsec_unit = {XSEC_UNIT:.6e} cm2/nucleon")
'''
    )
)

UNF.append(
    md(
        """## 3. Closure test (Asimov / MC)

Unfold the selected signal reco spectrum as if it were data:
`measured = nevts_sel_reco * xsec_unit`, with MC-stat-only diagonal covariance
(plus a tiny floor). Compare unfolded result to `A_c @ model`."""
    )
)

UNF.append(
    code(
        r'''def asimov_covariance(nevts_sel_reco, xsec_unit):
    """MC-stat-only cov for closure (diagonal in event counts, scaled to xsec)."""
    n = np.asarray(nevts_sel_reco, dtype=float)
    # fractional MC-stat ~ 1/n with floor
    frac = np.where(n > 0, 1.0 / n, 0.0)
    frac_cov = np.diag(frac)
    return cov_from_fraccov(frac_cov, n) * (xsec_unit ** 2)


closure_results = {}
for vsn in vsns:
    vc = vc_by.get(vsn)
    if vc is None:
        print("skip (no VariableConfig):", vsn)
        continue
    pack = variables[vsn]
    response = np.asarray(pack["response"], dtype=float)
    nevts_allmc = np.asarray(pack["nevts_allmc"], dtype=float)
    nevts_sel_reco = np.asarray(pack["nevts_sel_reco"], dtype=float)
    model = nevts_allmc * XSEC_UNIT
    measured = nevts_sel_reco * XSEC_UNIT
    cov = asimov_covariance(nevts_sel_reco, XSEC_UNIT)

    unfold = WienerSVD(
        response, model, measured, cov,
        C_TYPE, NORM_TYPE, stat_scaling=XSEC_UNIT,
    )
    ac_model = np.asarray(unfold["AddSmear"], dtype=float) @ model
    u = np.asarray(unfold["unfold"], dtype=float)
    ucov = np.asarray(unfold["UnfoldCov"], dtype=float)

    # chi2 vs smeared truth
    mask = (u > 0) & (ac_model > 0)
    ndof = int(np.count_nonzero(mask))
    if ndof >= 1:
        chi2, p_val = get_chi2(u[mask], ac_model[mask], ucov[np.ix_(mask, mask)])
    else:
        chi2, p_val = float("nan"), float("nan")

    print(f"[closure] {vsn}: chi2/ndof = {chi2:.3f}/{ndof}  p={p_val:.3g}")
    plot_unfolded_result(
        unfold,
        measured,
        {"GENIE (smeared truth)": ac_model, "GENIE (truth)": model},
        vc,
        plot=True,
        save_fig=True,
        save_name=str(FIG_DIR / f"{vsn}-closure"),
        closure_test=True,
    )
    closure_results[vsn] = {
        "chi2": float(chi2),
        "ndof": int(ndof),
        "unfold": u,
        "AddSmear": np.asarray(unfold["AddSmear"], dtype=float),
    }

print("closure done")
'''
    )
)

UNF.append(
    md(
        """## 4. Data unfold

`measured = (n_data - n_mc_bkg) * xsec_unit`  
`Covariance = cov_from_fraccov(total_xsec_frac, nevts_sel_reco) * xsec_unit**2`"""
    )
)

UNF.append(
    code(
        r'''data_results = {}
ingredients = {}

for vsn in vsns:
    vc = vc_by.get(vsn)
    if vc is None:
        continue
    pack = variables[vsn]
    response = np.asarray(pack["response"], dtype=float)
    nevts_allmc = np.asarray(pack["nevts_allmc"], dtype=float)
    nevts_sel_reco = np.asarray(pack["nevts_sel_reco"], dtype=float)
    n_sel_data = np.asarray(pack["n_sel_data"], dtype=float)

    model = nevts_allmc * XSEC_UNIT
    measured = n_sel_data * XSEC_UNIT

    try:
        frac_cov = total_cov_frac(summary, vsn, kind="xsec")
    except KeyError as ex:
        print(f"[data] skip {vsn}: no total_xsec ({ex})")
        continue

    Covariance = cov_from_fraccov(frac_cov, nevts_sel_reco) * (XSEC_UNIT ** 2)

    unfold = WienerSVD(
        response, model, measured, Covariance,
        C_TYPE, NORM_TYPE, stat_scaling=XSEC_UNIT,
    )

    ac_model = np.asarray(unfold["AddSmear"], dtype=float) @ model
    u = np.asarray(unfold["unfold"], dtype=float)
    ucov = np.asarray(unfold["UnfoldCov"], dtype=float)
    mask = (u > 0) & (ac_model > 0)
    ndof = int(np.count_nonzero(mask))
    if ndof >= 1:
        chi2, p_val = get_chi2(u[mask], ac_model[mask], ucov[np.ix_(mask, mask)])
    else:
        chi2, p_val = float("nan"), float("nan")

    print(f"[data] {vsn}: chi2/ndof = {chi2:.3f}/{ndof}  p={p_val:.3g}")
    plot_unfolded_result(
        unfold,
        measured,
        {"GENIE AR23": ac_model},
        vc,
        plot=True,
        save_fig=True,
        save_name=str(FIG_DIR / f"{vsn}-data_unfold"),
    )
    if response.shape[0] > 1:
        plot_heatmap(
            np.asarray(unfold["AddSmear"], dtype=float),
            vc.bins,
            plot_labels=["True", "True", f"{vsn} AddSmear (A_c)"],
            plot=True,
            save_fig=True,
            save_name=str(FIG_DIR / f"{vsn}-AddSmear"),
        )

    packed = pack_unfold_results(unfold, vc)
    packed["chi2_vs_genie"] = float(chi2)
    packed["ndof_vs_genie"] = int(ndof)
    packed["measured"] = measured
    packed["model"] = model
    packed["response"] = response
    packed["Covariance_input"] = Covariance
    packed["frac_cov_total_xsec"] = np.asarray(frac_cov, dtype=float)
    packed["n_sel_data"] = n_sel_data
    packed["nevts_sel_reco"] = nevts_sel_reco
    packed["nevts_allmc"] = nevts_allmc
    packed["xsec_unit"] = float(XSEC_UNIT)

    data_results[vsn] = packed
    ingredients[vsn] = {
        "bins": np.asarray(pack["bins"], dtype=float),
        "response": response,
        "eff": np.asarray(pack["eff"], dtype=float),
        "reco_vs_true": np.asarray(pack["reco_vs_true"], dtype=float),
        "n_data": np.asarray(pack["n_data"], dtype=float),
        "n_mc_bkg": np.asarray(pack["n_mc_bkg"], dtype=float),
        "n_sel_data": n_sel_data,
        "nevts_allmc": nevts_allmc,
        "nevts_sel_reco": nevts_sel_reco,
        "measured": measured,
        "model": model,
        "Covariance_input": Covariance,
        "frac_cov_total_xsec": np.asarray(frac_cov, dtype=float),
    }

print("data unfold done for", list(data_results))
'''
    )
)

UNF.append(md("## 5. Save flux + unfolded products for data release"))

UNF.append(
    code(
        r'''# Copy flux file into release dir + write flux_info.json
flux_dest = OUT_DIR / "Gen1_flux.root"
if FLUX_FILE.is_file():
    if not flux_dest.exists() or flux_dest.stat().st_size != FLUX_FILE.stat().st_size:
        shutil.copy2(FLUX_FILE, flux_dest)
        print("copied flux ->", flux_dest)
    else:
        print("flux already present:", flux_dest)
else:
    print("WARN: flux file missing:", FLUX_FILE)

flux_json = OUT_DIR / "flux_info.json"
with open(flux_json, "w") as f:
    json.dump(flux_info, f, indent=2)
print("wrote", flux_json)

# Combined pickle + flat npz
import pickle

release = {
    "meta": {
        "schema": "prl_productB_unfolded_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "response_npz": str(RESP_NPZ),
        "category_summary_npz": str(CAT_NPZ),
        "c_type": C_TYPE,
        "norm_type": NORM_TYPE,
        "xsec_unit": float(XSEC_UNIT),
        "data_tot_pot": float(data_tot_pot),
        "flux_info": flux_info,
        "response_meta": meta,
        "closure": {k: {"chi2": v["chi2"], "ndof": v["ndof"]} for k, v in closure_results.items()},
    },
    "ingredients": ingredients,
    "results": data_results,
}

pkl_path = OUT_DIR / "unfolding_ingredients_and_results.pkl"
with open(pkl_path, "wb") as f:
    pickle.dump(release, f, protocol=pickle.HIGHEST_PROTOCOL)
print("wrote", pkl_path)

# Flat NPZ for easy loading without pickle
savez = {
    "meta_json": np.array([json.dumps(release["meta"], default=str)], dtype=object),
    "xsec_unit": np.float64(XSEC_UNIT),
    "data_tot_pot": np.float64(data_tot_pot),
}
for vsn, res in data_results.items():
    for key in (
        "bins", "bin_centers", "bin_widths",
        "unfold", "unfold_per_bin_width",
        "stat_err", "syst_err", "total_err",
        "stat_err_per_bin_width", "syst_err_per_bin_width", "total_err_per_bin_width",
        "AddSmear", "UnfoldCov", "StatUnfoldCov", "SystUnfoldCov",
        "measured", "model", "response", "Covariance_input",
        "frac_cov_total_xsec", "n_sel_data", "nevts_sel_reco", "nevts_allmc",
    ):
        if key in res:
            savez[f"{vsn}::{key}"] = np.asarray(res[key])
    savez[f"{vsn}::chi2_vs_genie"] = np.float64(res.get("chi2_vs_genie", np.nan))
    savez[f"{vsn}::ndof_vs_genie"] = np.int32(res.get("ndof_vs_genie", -1))

npz_path = OUT_DIR / "unfolding_ingredients_and_results.npz"
np.savez_compressed(npz_path, **savez)
print("wrote", npz_path)

# Per-variable slim files
per_var_dir = OUT_DIR / "by_variable"
makedirs(per_var_dir, exist_ok=True)
for vsn, res in data_results.items():
    p = per_var_dir / f"{vsn}.npz"
    np.savez_compressed(
        p,
        bins=res["bins"],
        bin_centers=res["bin_centers"],
        unfold=res["unfold"],
        UnfoldCov=res["UnfoldCov"],
        StatUnfoldCov=res["StatUnfoldCov"],
        SystUnfoldCov=res["SystUnfoldCov"],
        AddSmear=res["AddSmear"],
        measured=res["measured"],
        model=res["model"],
        response=res["response"],
        xsec_unit=np.float64(XSEC_UNIT),
        chi2_vs_genie=np.float64(res.get("chi2_vs_genie", np.nan)),
        ndof_vs_genie=np.int32(res.get("ndof_vs_genie", -1)),
    )
    print("wrote", p)

manifest = {
    "schema": "prl_productB_unfolded_v1",
    "created_utc": release["meta"]["created_utc"],
    "pkl": str(pkl_path),
    "npz": str(npz_path),
    "flux_info": str(flux_json),
    "flux_file_copy": str(flux_dest),
    "variables": sorted(data_results),
    "closure_chi2": {k: {"chi2": v["chi2"], "ndof": v["ndof"]} for k, v in closure_results.items()},
}
with open(OUT_DIR / "unfolded_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)
print("wrote", OUT_DIR / "unfolded_manifest.json")
print("Done.")
'''
    )
)

write_nb("unfolding.ipynb", UNF)
print("all notebooks written")
