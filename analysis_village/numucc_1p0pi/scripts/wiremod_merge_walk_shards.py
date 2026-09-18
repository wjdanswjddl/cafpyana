#!/usr/bin/env python3
"""Merge WireMod walk shards + external CV walk → products/NPZs (no POT scale).

Envelope baseline is the matched Sep-4 CV sample Product A/B hists, not WireMod
in-file ``cv``. YZ/XTXW multi-univ walks are left unscaled — event matching
already defines the common sample.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.scripts.reprocess_wiremod import _merge_univ_products
from analysis_village.numucc_1p0pi.syst_detvar_common import (
    WIREMOD_ENVELOPE_SHIFTED,
    build_wiremod_detector_dict,
    log,
    save_detector_npz,
    wiremod_geometry_hists_for_envelope,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_DETECTOR,
    FILE_DETECTOR_SEL,
    SUB_DETECTOR,
    SUB_DETECTOR_SEL,
)


def _load_acc(path: Path) -> dict:
    with open(path, "rb") as fh:
        state = pickle.load(fh)
    if isinstance(state, dict) and "acc" in state:
        return state.get("acc") or {}
    return state


def _merge_cv_acc(acc: dict, chunk: dict) -> dict:
    if not acc:
        return {
            "hists_cut": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_cut"].items()},
            "hists_final": {k: np.asarray(v, dtype=float).copy() for k, v in chunk["hists_final"].items()},
            "pot": float(chunk["pot"]),
            "cut_var_names": list(chunk["cut_var_names"]),
            "final_var_names": list(chunk["final_var_names"]),
        }
    acc["pot"] = float(acc["pot"]) + float(chunk["pot"])
    for key in ("hists_cut", "hists_final"):
        for var, hist in chunk[key].items():
            acc[key][var] = np.asarray(acc[key].get(var, 0.0), dtype=float) + np.asarray(
                hist, dtype=float
            )
    return acc


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-base", required=True)
    p.add_argument("--yz-base-ckpt", default=None)
    p.add_argument("--yz-shard-ckpts", nargs="+", required=True)
    p.add_argument("--xtxw-shard-ckpts", nargs="+", required=True)
    p.add_argument("--cv-shard-ckpts", nargs="+", required=True)
    p.add_argument("--cv-campaign", default="2026_09_04_172912__sel_all-mc-CV")
    p.add_argument(
        "--mu-chi2mu-th",
        type=float,
        default=None,
        help="Recorded in products metadata (walks must already use this cut)",
    )
    args = p.parse_args(argv)

    out_base = Path(args.out_base)
    cache = out_base / "cache"
    cache.mkdir(parents=True, exist_ok=True)

    yz: dict = {}
    if args.yz_base_ckpt and Path(args.yz_base_ckpt).is_file():
        yz = _load_acc(Path(args.yz_base_ckpt))
        log(f"YZ base pot={yz.get('pot', 0):.3e}")
    for ck in args.yz_shard_ckpts:
        chunk = _load_acc(Path(ck))
        log(f"YZ merge {Path(ck).name} pot={chunk.get('pot', 0):.3e}")
        yz = _merge_univ_products(yz, chunk)

    xtxw: dict = {}
    for ck in args.xtxw_shard_ckpts:
        chunk = _load_acc(Path(ck))
        log(f"XTXW merge {Path(ck).name} pot={chunk.get('pot', 0):.3e}")
        xtxw = _merge_univ_products(xtxw, chunk)

    cv: dict = {}
    for ck in args.cv_shard_ckpts:
        chunk = _load_acc(Path(ck))
        log(f"CV merge {Path(ck).name} pot={chunk.get('pot', 0):.3e}")
        cv = _merge_cv_acc(cv, chunk)

    by_geom = {"YZ": yz, "XTXW": xtxw}
    pot_by = {
        "CV": float(cv["pot"]),
        "YZ": float(yz["pot"]),
        "XTXW": float(xtxw["pot"]),
    }
    log(f"POT (informational only — no scale applied): {pot_by}")

    mu_th = args.mu_chi2mu_th
    if mu_th is None:
        mu_th = (yz.get("mu_p_candidate_kwargs") or {}).get("mu_chi2mu_th")
    payload = {
        "by_geom": by_geom,
        "cv": cv,
        "pot_by_variation": pot_by,
        "pot_scales": {k: 1.0 for k in pot_by},
        "match_stage": "sel_all",
        "envelope": "calo_plus_efield_vs_external_cv",
        "cv_role": "envelope_baseline",
        "cv_campaign": args.cv_campaign,
        "batch_size": "sharded",
        "mu_chi2mu_th": mu_th,
        "mu_p_candidate_kwargs": {"mu_chi2mu_th": mu_th} if mu_th is not None else {},
    }
    product_cache = cache / "wiremod_sel_all_products.pkl"
    with open(product_cache, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"wrote {product_cache}")

    cut_names = yz["cut_var_names"]
    final_names = yz["final_var_names"]
    cv_cut = cv["hists_cut"]
    cv_final = cv["hists_final"]

    def _all_hists(product: str):
        return {
            lab: wiremod_geometry_hists_for_envelope(prod["by_universe"], product=product)
            for lab, prod in by_geom.items()
        }

    dict_a = build_wiremod_detector_dict(
        _all_hists("cut"),
        cut_names,
        wiremod_labels=("YZ", "XTXW"),
        shifted_univs=WIREMOD_ENVELOPE_SHIFTED,
        cv_hists=cv_cut,
    )
    dict_b = build_wiremod_detector_dict(
        _all_hists("final"),
        final_names,
        wiremod_labels=("YZ", "XTXW"),
        shifted_univs=WIREMOD_ENVELOPE_SHIFTED,
        cv_hists=cv_final,
    )
    npz_a = out_base / SUB_DETECTOR_SEL / FILE_DETECTOR_SEL
    npz_b = out_base / SUB_DETECTOR / FILE_DETECTOR
    save_detector_npz(
        dict_a,
        npz_a,
        manifest={
            "source": "WireMod",
            "product": "A_selection",
            "method": "maxabs_dev_unc_actual_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_a.get("detector", {})),
            "mu_chi2mu_th": mu_th,
        },
    )
    save_detector_npz(
        dict_b,
        npz_b,
        manifest={
            "source": "WireMod",
            "product": "B_measurement",
            "method": "maxabs_dev_unc_actual_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_b.get("detector", {})),
            "mu_chi2mu_th": mu_th,
        },
    )
    log(f"Product A → {npz_a}")
    log(f"Product B → {npz_b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
