#!/usr/bin/env python
"""Walk matched WireMod + DENT sel_all dfs for chi2 track-subset vars only.

Requires env ``NUMUCC_CUT_STAGE_SLUGS`` (set automatically). Writes
``<out-dir>/detector_sel_syst_dict.npz``.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path
from typing import List

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.syst_detvar_common import (
    accumulate_matched_sel_all_cv_vs_var_products,
    accumulate_matched_sel_all_products,
    accumulate_wiremod_matched_products,
    build_dent_detector_dict,
    build_wiremod_detector_dict,
    combine_wiremod_dent_detector_dict,
    log,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import CHI2_TRACK_SUBSET_SLUGS

WM_ROOT = Path(
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final-archive/WireMod/matched"
)
DENT_ROOT = Path(
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final-archive/DENT/matched"
)


def _list_dfs(root: Path, max_files: int = 0) -> List[str]:
    files = sorted(str(p) for p in root.glob("*.df"))
    if max_files > 0:
        files = files[:max_files]
    return files


def _wiremod_univ_cut_hists(acc: dict) -> dict:
    """``{univ: {var: hist}}`` from accumulate_wiremod_matched_products."""
    by_u = acc.get("by_universe") or {}
    out = {}
    for u, info in by_u.items():
        if isinstance(info, dict) and "hists_cut" in info:
            out[u] = info["hists_cut"]
        elif isinstance(info, dict):
            out[u] = info
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--workers", type=int, default=4, help="reserved (walk is sequential)")
    p.add_argument("--max-files", type=int, default=0)
    args = p.parse_args(argv)

    os.environ["NUMUCC_CUT_STAGE_SLUGS"] = ",".join(CHI2_TRACK_SUBSET_SLUGS)
    os.environ["NUMUCC_SKIP_FINAL_STAGE"] = "1"
    os.environ.setdefault("MPLBACKEND", "Agg")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    yz_files = _list_dfs(WM_ROOT / "yz", args.max_files)
    xt_files = _list_dfs(WM_ROOT / "xtxw", args.max_files)
    cv_files = _list_dfs(WM_ROOT / "cv", args.max_files)
    log(f"WireMod yz={len(yz_files)} xtxw={len(xt_files)} cv={len(cv_files)}")

    yz_acc = accumulate_wiremod_matched_products(
        yz_files, final_var_defs={}, include_cut_stage=True
    )
    xt_acc = accumulate_wiremod_matched_products(
        xt_files, final_var_defs={}, include_cut_stage=True
    )
    # Matched CV sample is CV-only (no calo univs) — use plain sel_all walk.
    cv_acc = accumulate_matched_sel_all_products(
        cv_files, final_var_defs={}, include_cut_stage=True
    )

    all_hists = {
        "YZ": _wiremod_univ_cut_hists(yz_acc),
        "XTXW": _wiremod_univ_cut_hists(xt_acc),
    }
    cv_hists = {
        s: np.asarray(cv_acc["hists_cut"][s], dtype=float)
        for s in CHI2_TRACK_SUBSET_SLUGS
        if s in cv_acc.get("hists_cut", {})
    }
    with open(out_dir / "wiremod_hists.pkl", "wb") as f:
        pickle.dump({"all_hists": all_hists, "cv_hists": cv_hists}, f)

    wm_dict = build_wiremod_detector_dict(
        all_hists,
        list(CHI2_TRACK_SUBSET_SLUGS),
        cv_hists=cv_hists or None,
    )

    dent_cv = _list_dfs(DENT_ROOT / "cv", args.max_files)
    dent_var = _list_dfs(DENT_ROOT / "dent", args.max_files)
    log(f"DENT cv={len(dent_cv)} dent={len(dent_var)}")
    dent_prod = accumulate_matched_sel_all_cv_vs_var_products(
        dent_cv,
        dent_var,
        final_var_defs={},
        include_cut_stage=True,
    )
    with open(out_dir / "dent_hists.pkl", "wb") as f:
        pickle.dump(dent_prod, f)

    h_cv = dent_prod["cv"]["hists_cut"]
    h_var = dent_prod["var"]["hists_cut"]
    dent_dict = build_dent_detector_dict(h_cv, h_var, list(CHI2_TRACK_SUBSET_SLUGS))

    combined = combine_wiremod_dent_detector_dict(wm_dict, dent_dict)
    out_npz = out_dir / "detector_sel_syst_dict.npz"
    np.savez_compressed(out_npz, **combined)
    log(f"wrote {out_npz} detector vars={sorted((combined.get('detector') or {}).keys())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
