#!/usr/bin/env python3
"""WireMod matched write + batched Product A/B walk (memory-safe).

Uses Sep-14 YZ/XTXW + Sep-4 CV common keys. Walks YZ/XTXW updatecalo universes
and the matched Sep-4 CV sample (envelope baseline). Event matching defines the
common sample — do **not** POT-scale. Envelope = max |univ − CV| per geometry.
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.dataset_locations import PLOTS_BASE, SPRING_GEN1_ROOT
from analysis_village.numucc_1p0pi.scripts import dent_compare as dc
from analysis_village.numucc_1p0pi.syst_detvar_common import (
    WIREMOD_ENVELOPE_SHIFTED,
    accumulate_matched_sel_all_products,
    accumulate_wiremod_matched_products,
    build_wiremod_detector_dict,
    glob_matched_dfs,
    log,
    run_dent_match,
    save_detector_npz,
    wiremod_geometry_hists_for_envelope,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_DETECTOR,
    FILE_DETECTOR_SEL,
    SUB_DETECTOR,
    SUB_DETECTOR_SEL,
)

DEFAULT_RAW = {
    "yz": "2026_09_14_025216__sel_all-mc-BNB_cosmics-WireModYZ",
    "xtxw": "2026_09_14_024629__sel_all-mc-BNB_cosmics-WireModXTXW",
    "cv": "2026_09_04_172912__sel_all-mc-CV",
}


def _rss_gb() -> float:
    try:
        with open(f"/proc/{os.getpid()}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1]) / (1024.0 * 1024.0)
    except Exception:
        pass
    return float("nan")


def _merge_univ_products(acc: dict, chunk: dict) -> dict:
    """Add chunk universe hists into accumulator (in place)."""
    if not acc:
        return {
            "by_universe": {
                u: {
                    "hists_cut": {k: np.asarray(v, dtype=float).copy() for k, v in p["hists_cut"].items()},
                    "hists_final": {k: np.asarray(v, dtype=float).copy() for k, v in p["hists_final"].items()},
                }
                for u, p in chunk["by_universe"].items()
            },
            "pot": float(chunk["pot"]),
            "cut_var_names": list(chunk["cut_var_names"]),
            "final_var_names": list(chunk["final_var_names"]),
            "universes": list(chunk["universes"]),
        }
    acc["pot"] = float(acc["pot"]) + float(chunk["pot"])
    for u, payload in chunk["by_universe"].items():
        dst = acc["by_universe"].setdefault(
            u,
            {
                "hists_cut": {k: np.zeros_like(v, dtype=float) for k, v in payload["hists_cut"].items()},
                "hists_final": {k: np.zeros_like(v, dtype=float) for k, v in payload["hists_final"].items()},
            },
        )
        for key in ("hists_cut", "hists_final"):
            for var, hist in payload[key].items():
                dst[key][var] = np.asarray(dst[key].get(var, 0.0), dtype=float) + np.asarray(
                    hist, dtype=float
                )
    return acc


def walk_wiremod_batched(
    files: Sequence[str],
    *,
    batch_size: int,
    checkpoint: Path,
    final_var_defs: Mapping[str, dict],
    rss_limit_gb: float,
) -> dict:
    """Walk matched WireMod files in batches; checkpoint after each batch."""
    files = list(files)
    n = len(files)
    start_batch = 0
    acc: dict = {}
    if checkpoint.is_file():
        with open(checkpoint, "rb") as fh:
            state = pickle.load(fh)
        acc = state.get("acc") or {}
        start_batch = int(state.get("next_batch", 0))
        log(f"  resume {checkpoint.name}: next_batch={start_batch} pot={acc.get('pot', 0):.3e}")

    n_batches = (n + batch_size - 1) // batch_size
    for bi in range(start_batch, n_batches):
        lo = bi * batch_size
        hi = min(n, lo + batch_size)
        batch = files[lo:hi]
        rss0 = _rss_gb()
        log(f"  batch {bi + 1}/{n_batches} files[{lo}:{hi}] rss={rss0:.2f} GiB")
        if rss0 > rss_limit_gb:
            raise RuntimeError(
                f"RSS {rss0:.2f} GiB exceeds limit {rss_limit_gb:.2f} GiB before batch {bi}"
            )
        chunk = accumulate_wiremod_matched_products(
            batch,
            final_var_defs=final_var_defs,
            include_cut_stage=True,
        )
        acc = _merge_univ_products(acc, chunk)
        del chunk
        gc.collect()
        rss1 = _rss_gb()
        log(f"    pot_acc={acc['pot']:.3e} rss={rss1:.2f} GiB")
        if rss1 > rss_limit_gb:
            raise RuntimeError(
                f"RSS {rss1:.2f} GiB exceeds limit {rss_limit_gb:.2f} GiB after batch {bi}"
            )
        with open(checkpoint, "wb") as fh:
            pickle.dump({"acc": acc, "next_batch": bi + 1, "n_files": n}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return acc


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--skip-match", action="store_true", help="Reuse existing matched dfs")
    p.add_argument("--skip-hist", action="store_true", help="Skip hist walk / NPZ")
    p.add_argument("--force-hist", action="store_true", help="Ignore product cache / batch checkpoints")
    p.add_argument("--batch-size", type=int, default=int(os.environ.get("WIREMOD_BATCH_SIZE", "50")))
    p.add_argument("--rss-limit-gb", type=float, default=float(os.environ.get("WIREMOD_RSS_LIMIT_GB", "20")))
    p.add_argument("--max-files", type=int, default=None, help="Cap matched inputs (smoke test)")
    args = p.parse_args(argv)

    dfs = Path(os.environ.get("NUMUCC_SPRING_GEN1_ROOT", SPRING_GEN1_ROOT))
    out_base = Path(os.environ.get("WIREMOD_OUT_BASE", str(Path(PLOTS_BASE) / "systematics-final" / "WireMod")))
    matched_out = out_base / "matched"
    cache = out_base / "cache"
    det_a = out_base / SUB_DETECTOR_SEL
    det_b = out_base / SUB_DETECTOR
    for d in (matched_out / "yz", matched_out / "xtxw", matched_out / "cv", cache, det_a, det_b):
        d.mkdir(parents=True, exist_ok=True)

    raw = {k: dfs / DEFAULT_RAW[k] for k in DEFAULT_RAW}
    keys_pkl = cache / "wiremod_common_keys_sel_all.pkl"
    if not keys_pkl.is_file():
        raise SystemExit(f"missing common keys: {keys_pkl} (run overlap scan first)")

    if not args.skip_match:
        log("=== WireMod matched write (phase=write) ===")
        rc = run_dent_match(
            {k: str(v) for k, v in raw.items()},
            fmt="sel_all",
            filename_str="sel_all",
            matched_out_dir=str(matched_out),
            matched_suffix="_matched",
            common_keys_pkl=str(keys_pkl),
            phase="write",
            summary_csv=str(cache / "wiremod_matched_summary-sel_all.csv"),
            n_workers=1,
            skip_existing_matched=True,
        )
        if rc != 0:
            return int(rc)

    if args.skip_hist:
        log("skip hist walk")
        return 0

    matched = {
        "YZ": glob_matched_dfs(matched_out / "yz", filename_str="sel_all"),
        "XTXW": glob_matched_dfs(matched_out / "xtxw", filename_str="sel_all"),
        "CV": glob_matched_dfs(matched_out / "cv", filename_str="sel_all"),
    }
    if args.max_files is not None:
        matched = {k: v[: args.max_files] for k, v in matched.items()}
    for lab, files in matched.items():
        log(f"{lab}: {len(files)} matched files")
        if not files:
            raise SystemExit(f"no matched files for {lab} under {matched_out}")

    product_cache = cache / "wiremod_sel_all_products.pkl"
    if product_cache.is_file() and not args.force_hist:
        log(f"reuse product cache {product_cache}")
        with open(product_cache, "rb") as fh:
            payload = pickle.load(fh)
    else:
        if args.force_hist:
            for lab in ("YZ", "XTXW"):
                ck = cache / f"wiremod_walk_checkpoint_{lab.lower()}.pkl"
                if ck.is_file():
                    ck.unlink()
                    log(f"removed {ck.name}")

        final_defs = dc.build_final_var_defs()
        by_geom: Dict[str, dict] = {}
        pot_by: Dict[str, float] = {}

        # External matched CV: full Product A/B walk (envelope baseline; no POT scale)
        log(f"[CV] matched sel_all walk n={len(matched['CV'])} (envelope baseline)")
        cv_prod = accumulate_matched_sel_all_products(
            matched["CV"], final_var_defs=final_defs, include_cut_stage=True
        )
        pot_by["CV"] = float(cv_prod["pot"])
        log(f"[CV] done POT={cv_prod['pot']:.3e} (informational; no scale applied)")

        for lab in ("YZ", "XTXW"):
            ck = cache / f"wiremod_walk_checkpoint_{lab.lower()}.pkl"
            log(f"[{lab}] batched walk batch_size={args.batch_size} rss_limit={args.rss_limit_gb} GiB")
            prod = walk_wiremod_batched(
                matched[lab],
                batch_size=max(int(args.batch_size), 1),
                checkpoint=ck,
                final_var_defs=final_defs,
                rss_limit_gb=float(args.rss_limit_gb),
            )
            by_geom[lab] = prod
            pot_by[lab] = float(prod["pot"])
            log(f"  {lab} done POT={prod['pot']:.3e} univs={prod.get('universes')}")

        payload = {
            "by_geom": by_geom,
            "cv": cv_prod,
            "pot_by_variation": pot_by,
            "pot_scales": {k: 1.0 for k in pot_by},
            "match_stage": "sel_all",
            "envelope": "calo_plus_efield_vs_external_cv",
            "cv_role": "envelope_baseline",
            "cv_campaign": DEFAULT_RAW["cv"],
            "batch_size": args.batch_size,
        }
        with open(product_cache, "wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
        log(f"wrote {product_cache}")

    by_geom = payload["by_geom"]
    if "cv" not in payload:
        raise SystemExit(
            f"{product_cache} lacks external CV products — rebuild with force-hist "
            "(envelope must use matched Sep-4 CV, not WireMod in-file cv)"
        )
    cv_prod = payload["cv"]
    cut_names = next(iter(by_geom.values()))["cut_var_names"]
    final_names = next(iter(by_geom.values()))["final_var_names"]

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
        cv_hists=cv_prod["hists_cut"],
    )
    dict_b = build_wiremod_detector_dict(
        _all_hists("final"),
        final_names,
        wiremod_labels=("YZ", "XTXW"),
        shifted_univs=WIREMOD_ENVELOPE_SHIFTED,
        cv_hists=cv_prod["hists_final"],
    )
    npz_a = det_a / FILE_DETECTOR_SEL
    npz_b = det_b / FILE_DETECTOR
    save_detector_npz(
        dict_a,
        npz_a,
        manifest={
            "source": "WireMod",
            "product": "A_selection",
            "method": "total_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_a.get("detector", {})),
        },
    )
    save_detector_npz(
        dict_b,
        npz_b,
        manifest={
            "source": "WireMod",
            "product": "B_measurement",
            "method": "total_envelope_vs_matched_cv",
            "cv_role": "envelope_baseline",
            "shifted_univs": list(WIREMOD_ENVELOPE_SHIFTED),
            "n_vars": len(dict_b.get("detector", {})),
        },
    )
    log(f"Product A → {npz_a} ({len(dict_a.get('detector', {}))} vars)")
    log(f"Product B → {npz_b} ({len(dict_b.get('detector', {}))} vars)")
    log(f"done rss={_rss_gb():.2f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
