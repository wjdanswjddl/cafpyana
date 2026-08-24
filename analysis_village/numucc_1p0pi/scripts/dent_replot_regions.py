#!/usr/bin/env python3
"""
Remake DENT region plots from a saved ``dent_sel_*_regions.pkl`` cache
without re-reading matched dataframe files.

Example:
    python dent_replot_regions.py \\
        --cache .../DENT-regions/cache/dent_sel_mup_regions.pkl \\
        --fig-dir .../DENT-regions/plots/sel_mup \\
        --regions all half_E half_W
"""

from __future__ import annotations

import argparse
import sys
from os import makedirs, path
from typing import Dict, List, Tuple

import numpy as np

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc


def _ordered_region_keys(payload: dict, kind: str) -> List[str]:
    """Prefer cached half/octant meta order; fall back to sorted region keys."""
    meta_key = {"half": "half_meta", "octant": "octant_meta"}.get(kind)
    if meta_key and meta_key in payload:
        return [m["slug"] for m in payload[meta_key]]
    keys = []
    for key, blob in payload.get("regions", {}).items():
        if blob.get("meta", {}).get("kind") == kind:
            keys.append(key)
    return sorted(keys)


def _panel_title(meta: dict, key: str) -> str:
    """Panel title; strip legacy East/West suffixes from half titles."""
    title = meta.get("title") or meta.get("name") or key
    for suffix in (" (East)", " (West)"):
        if title.endswith(suffix):
            title = title[: -len(suffix)]
    return title


def _panel_entries(
    payload: dict, keys: List[str],
) -> List[Tuple[str, Dict[str, Dict[str, np.ndarray]]]]:
    out = []
    for key in keys:
        blob = payload["regions"].get(key)
        if blob is None:
            continue
        meta = blob.get("meta", {})
        out.append((_panel_title(meta, key), {"cv": blob["hists"]["cv"], "dent": blob["hists"]["dent"]}))
    return out


# Plane-I2 chi2 vars superseded by plane-averaged chi2 in DENT diagnostics.
_SKIP_PANEL_VARS = frozenset({"chi2_mu", "chi2_p"})


def plot_region_panels(payload: dict, fig_dir: str, var_defs: dict) -> None:
    """Wide slide-friendly panels: 1×2 halves and 2×4 octants."""
    half_keys = _ordered_region_keys(payload, "half")
    oct_keys = _ordered_region_keys(payload, "octant")
    half_entries = _panel_entries(payload, half_keys)
    oct_entries = _panel_entries(payload, oct_keys)

    for var_name, cfg in var_defs.items():
        if var_name in _SKIP_PANEL_VARS:
            continue
        plot_cfg = {
            "label": cfg.get("label", var_name),
            "bins": np.asarray(cfg["bins"]),
        }
        if half_entries and all(var_name in e[1]["cv"] for e in half_entries):
            odir = path.join(fig_dir, "halves")
            makedirs(odir, exist_ok=True)
            print(f"panel halves → {odir}/panel_diagnostic_{var_name}", flush=True)
            dc.plot_var_comparison_panel(
                var_name, plot_cfg, half_entries,
                fig_dir=odir,
                fig_name=f"panel_diagnostic_{var_name}",
                nrows=1, ncols=2, figsize=(14.0, 4.8),
            )
        if oct_entries and all(var_name in e[1]["cv"] for e in oct_entries):
            odir = path.join(fig_dir, "octants")
            makedirs(odir, exist_ok=True)
            print(f"panel octants → {odir}/panel_diagnostic_{var_name}", flush=True)
            dc.plot_var_comparison_panel(
                var_name, plot_cfg, oct_entries,
                fig_dir=odir,
                fig_name=f"panel_diagnostic_{var_name}",
                nrows=2, ncols=4, figsize=(16.0, 7.2),
            )


def main() -> int:
    p = argparse.ArgumentParser(description="Replot DENT regions from cache")
    p.add_argument("--cache", required=True, help="Path to dent_sel_*_regions.pkl")
    p.add_argument("--fig-dir", required=True, help="Output plots directory for this stage")
    p.add_argument(
        "--regions",
        nargs="*",
        default=None,
        help="Region keys to replot (default: all regions in cache)",
    )
    p.add_argument("--skip-overlays", action="store_true")
    p.add_argument("--skip-frac", action="store_true")
    p.add_argument("--skip-ratio", action="store_true")
    p.add_argument("--skip-panels", action="store_true")
    p.add_argument("--panels-only", action="store_true",
                   help="Only make half/octant panel plots")
    args = p.parse_args()

    payload = dc.load_hists(args.cache)
    var_defs = payload.get("var_defs", {})
    frac_bins = np.asarray(payload.get("frac_bins", dc.FRAC_DIFF_BINS))
    ratio_bins = np.asarray(payload.get("ratio_bins", dc.RATIO_BINS))
    regions = args.regions or list(payload.get("regions", {}).keys())

    if not args.panels_only:
        for key in regions:
            if key not in payload["regions"]:
                print(f"skip unknown region {key}", flush=True)
                continue
            blob = payload["regions"][key]
            meta = blob.get("meta", {"name": key, "slug": key, "kind": "all"})
            kind = meta.get("kind", "all")
            slug = meta.get("slug", key)
            if kind == "all":
                odir = path.join(args.fig_dir, "all")
            elif kind == "octant":
                odir = path.join(args.fig_dir, "octants", slug)
            elif kind == "half":
                odir = path.join(args.fig_dir, "halves", slug)
            else:
                odir = path.join(args.fig_dir, slug)

            print(f"replot {key} → {odir}", flush=True)
            if not args.skip_overlays:
                all_hists = {"cv": blob["hists"]["cv"], "dent": blob["hists"]["dent"]}
                for var_name, cfg in var_defs.items():
                    if var_name in _SKIP_PANEL_VARS:
                        continue
                    if var_name not in all_hists["cv"]:
                        continue
                    plot_cfg = {
                        "label": cfg.get("label", var_name),
                        "bins": np.asarray(cfg["bins"]),
                    }
                    dc.plot_var_comparison(
                        var_name, plot_cfg, all_hists, fig_dir=odir, pot_scales=None,
                    )

            paired = blob.get("paired", {})
            tag = meta.get("name", key)
            for var_name, pdata in paired.items():
                label = var_defs.get(var_name, {}).get("label", var_name)
                cfg = {"label": f"{label}  [{tag}]"}
                if not args.skip_frac:
                    makedirs(path.join(odir, "fracdiff"), exist_ok=True)
                    dc.plot_frac_diff(
                        var_name, cfg, np.asarray(pdata["frac"]),
                        fig_dir=path.join(odir, "fracdiff"), bins=frac_bins,
                    )
                if not args.skip_ratio:
                    makedirs(path.join(odir, "ratio"), exist_ok=True)
                    dc.plot_ratio_dist(
                        var_name, cfg, np.asarray(pdata["ratio"]),
                        fig_dir=path.join(odir, "ratio"), bins=ratio_bins,
                    )

    if not args.skip_panels:
        plot_region_panels(payload, args.fig_dir, var_defs)

    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
