#!/usr/bin/env python
"""Turn stream-sum ``rate_covs.pkl`` files into Product A/B Flux/G4 NPZs + plots.

Product **A** = cut-stage slugs (``CUT_STAGE_RATE_ONLY_SLUGS``).
Product **B** = final-stage measurement variables (everything else).

Expected pickle (from ``syst_histcounts_stream_sum.py --build-covs``)::

    {"family": "Flux"|"G4", "by_var": {slug: {knob: pack, family.lower(): combined}}, ...}

Writes::

    <syst-disk-root>/productA/Flux/flux_syst_dict.npz
    <syst-disk-root>/productA/G4/g4_syst_dict.npz
    <syst-disk-root>/productB/Flux/flux_syst_dict.npz
    <syst-disk-root>/productB/G4/g4_syst_dict.npz
    <syst-disk-root>/plots/product{A,B}/{Flux,G4}/...
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from os import path
from typing import Any, Dict, Mapping, Optional, Sequence

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import numpy as np

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.syst_histcounts import histcounts_var_configs
from analysis_village.numucc_1p0pi.syst_multisim_common import save_neutrino_multisim_npzs
from analysis_village.numucc_1p0pi.syst_pipeline_walker import CUT_STAGE_RATE_ONLY_SLUGS
from analysis_village.numucc_1p0pi.utils import plot_frac_unc, plot_heatmap

SLIM_KEYS = frozenset(
    {
        "slim",
        "slim_multisim",
        "Flux_slim",
        "Flux_slim_multisim",
        "G4_slim",
        "G4_slim_multisim",
        "flux_total",
        "g4_total",
        "Flux",
        "G4",
        "GENIE",
    }
)


def parse_cli(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cov-pkl",
        action="append",
        required=True,
        help="rate_covs.pkl (repeat for Flux and G4)",
    )
    p.add_argument(
        "--syst-disk-root",
        required=True,
        help="Root that will contain productA/ and productB/",
    )
    p.add_argument(
        "--no-plots",
        action="store_true",
        help="Write NPZs only",
    )
    return p.parse_args(argv)


def _family_norm(raw: str) -> str:
    s = str(raw).strip()
    if s.lower() == "flux":
        return "Flux"
    if s.lower() == "g4":
        return "G4"
    return s


def _combined_key(family: str) -> str:
    return "flux" if family == "Flux" else "G4"


def _load_cov(pkl_path: str) -> tuple[str, Dict[str, Dict[str, Any]]]:
    with open(pkl_path, "rb") as fh:
        blob = pickle.load(fh)
    family = _family_norm(blob.get("family", ""))
    by_var = blob.get("by_var")
    if family not in ("Flux", "G4") or not isinstance(by_var, dict):
        raise SystemExit("bad cov pickle %s (need family Flux|G4 and by_var)" % pkl_path)
    return family, by_var


def _split_a_b(
    by_var: Mapping[str, Mapping[str, Any]],
) -> tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    a: Dict[str, Dict[str, Any]] = {}
    b: Dict[str, Dict[str, Any]] = {}
    for slug, cell in by_var.items():
        if slug in CUT_STAGE_RATE_ONLY_SLUGS:
            a[slug] = dict(cell)
        else:
            b[slug] = dict(cell)
    return a, b


def _to_syst_dict(family: str, by_var: Mapping[str, Mapping[str, Any]]) -> dict:
    comb = _combined_key(family)
    block: Dict[str, Any] = {}
    by_knob: Dict[str, Dict[str, Any]] = {}
    for slug, cell in by_var.items():
        combined = cell.get(comb) or cell.get(family.lower()) or cell.get(family)
        if combined is None:
            continue
        block[slug] = combined
        knobs = {
            k: v
            for k, v in cell.items()
            if k not in SLIM_KEYS
            and k not in (comb, family, family.lower())
            and isinstance(v, dict)
            and "cov_frac" in v
        }
        if knobs:
            by_knob[slug] = knobs
    out = {family: block}
    if by_knob:
        out["%s_by_knob" % family] = by_knob
    return out


def _var_config_map() -> Dict[str, Any]:
    return {vc.var_save_name: vc for vc in histcounts_var_configs()}


def _frac_unc(pack: Mapping[str, Any]) -> np.ndarray:
    cf = np.asarray(pack["cov_frac"], dtype=float)
    return np.sqrt(np.clip(np.diag(cf), 0.0, None))


def _plot_product(
    *,
    family: str,
    product: str,
    by_var: Mapping[str, Mapping[str, Any]],
    vc_map: Mapping[str, Any],
    out_dir: str,
) -> int:
    os.makedirs(out_dir, exist_ok=True)
    comb = _combined_key(family)
    n = 0
    for slug, cell in sorted(by_var.items()):
        pack = cell.get(comb) or cell.get(family.lower()) or cell.get(family)
        if pack is None:
            continue
        vc = vc_map.get(slug)
        bins = np.asarray(pack["cov_frac"]).shape[0]
        dummy_bins = np.arange(bins + 1, dtype=float) if vc is None else np.asarray(vc.bins)
        labels = ["", "", "%s %s %s" % (family, product, slug)]
        stem = path.join(out_dir, slug)
        plot_frac_unc(
            [_frac_unc(pack)],
            vc
            if vc is not None
            else type("VC", (), {"bins": dummy_bins, "bin_centers": 0.5 * (dummy_bins[:-1] + dummy_bins[1:]), "var_save_name": slug, "var_labels": [slug, slug]})(),
            plot_labels=labels,
            plot=False,
            save_fig=True,
            save_name=stem + "-frac_unc",
        )
        cf = np.asarray(pack["cov_frac"], dtype=float)
        corr = np.asarray(pack.get("corr", pack.get("corr_frac", np.eye(cf.shape[0]))), dtype=float)
        if cf.shape[0] >= 1:
            plot_heatmap(
                cf,
                dummy_bins,
                plot_labels=["", "", "%s cov_frac" % slug],
                plot=False,
                cmap="viridis",
                save_fig=True,
                save_name=stem + "-cov_frac",
            )
            plot_heatmap(
                corr,
                dummy_bins,
                plot_labels=["", "", "%s corr" % slug],
                plot=False,
                save_fig=True,
                save_name=stem + "-corr",
            )
        n += 1
    print("[covs-to-npz] plots %s %s  n=%d  dir=%s" % (family, product, n, out_dir), flush=True)
    return n


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_cli(argv)
    root = path.abspath(path.expanduser(args.syst_disk_root.rstrip(os.sep)))
    vc_map = _var_config_map()
    for pkl in args.cov_pkl:
        family, by_var = _load_cov(pkl)
        a, b = _split_a_b(by_var)
        print(
            "[covs-to-npz] %s  vars=%d  A=%d  B=%d  from %s"
            % (family, len(by_var), len(a), len(b), pkl),
            flush=True,
        )
        for product, block in (("productA", a), ("productB", b)):
            dest = path.join(root, product)
            save_neutrino_multisim_npzs(_to_syst_dict(family, block), dest)
            npz_name = "flux_syst_dict.npz" if family == "Flux" else "g4_syst_dict.npz"
            print(
                "[covs-to-npz] wrote %s"
                % path.join(dest, family, npz_name),
                flush=True,
            )
            if not args.no_plots:
                _plot_product(
                    family=family,
                    product=product,
                    by_var=block,
                    vc_map=vc_map,
                    out_dir=path.join(root, "plots", product, family),
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
