#!/usr/bin/env python3
"""Pair-wise DENT match (CV1↔DENT1, CV2↔DENT2) then combined unisim products.

Match within each pair so batch metadata overlap cannot cross-match, then treat
the union of CV matched files as CV and DENT as DENT for Products A/B.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.dataset_locations import PLOTS_BASE, SPRING_GEN1_ROOT
from analysis_village.numucc_1p0pi.syst_detvar_common import (
    accumulate_matched_sel_all_cv_vs_var_products,
    assert_variations_matched,
    build_dent_detector_dict,
    glob_matched_dfs,
    load_or_build_cache,
    log,
    pot_scales_to_cv,
    run_dent_match,
    save_detector_npz,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (
    FILE_DETECTOR,
    FILE_DETECTOR_SEL,
    SUB_DETECTOR,
    SUB_DETECTOR_SEL,
)

DEFAULT_PAIRS = [
    {
        "name": "pair1",
        "cv": "2026_09_14_142556__sel_all-mc-CV1",
        "dent": "2026_09_14_142722__sel_all-mc-DENT1",
    },
    {
        "name": "pair2",
        "cv": "2026_09_14_153348__sel_all-mc-CV2",
        "dent": "2026_09_14_142928__sel_all-mc-DENT2",
    },
]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--skip-match", action="store_true", help="Reuse existing matched files")
    p.add_argument("--force-hist", action="store_true", help="Rewalk hists even if cache exists")
    p.add_argument(
        "--only-pair",
        action="append",
        default=None,
        help="Only match this pair name (repeatable). Default: all pairs.",
    )
    p.add_argument("--n-workers", type=int, default=int(os.environ.get("DENT_MATCH_N_WORKERS", "8")))
    p.add_argument("--file-timeout", type=float, default=float(os.environ.get("DENT_MATCH_FILE_TIMEOUT", "600")))
    p.add_argument("--meta-retries", type=int, default=int(os.environ.get("DENT_MATCH_META_RETRIES", "3")))
    p.add_argument("--assert-max-files", type=int, default=20)
    args = p.parse_args(argv)

    dfs = Path(os.environ.get("NUMUCC_SPRING_GEN1_ROOT", SPRING_GEN1_ROOT))
    out_base = Path(PLOTS_BASE) / "systematics-final" / "DENT-highstats"
    matched_out = out_base / "matched"
    cache = out_base / "cache"
    match_cache = Path(PLOTS_BASE) / "systematics-final" / "DetectorMatch" / "cache"
    for d in (matched_out / "cv", matched_out / "dent", cache, match_cache):
        d.mkdir(parents=True, exist_ok=True)

    pairs = [
        {
            "name": spec["name"],
            "cv": dfs / spec["cv"],
            "dent": dfs / spec["dent"],
        }
        for spec in DEFAULT_PAIRS
    ]
    if args.only_pair:
        wanted = set(args.only_pair)
        pairs = [pair for pair in pairs if pair["name"] in wanted]
        if not pairs:
            raise SystemExit(f"--only-pair matched nothing; known: {[s['name'] for s in DEFAULT_PAIRS]}")
    suffix = "_matched_hs"

    if not args.skip_match:
        pair_key_sets = []
        for pair in pairs:
            if not pair["cv"].is_dir() or not pair["dent"].is_dir():
                raise FileNotFoundError(f"missing dirs for {pair['name']}: {pair['cv']} / {pair['dent']}")
            keys_pkl = cache / f"dent_common_keys_sel_all_{pair['name']}.pkl"
            log(f"=== DENT match [{pair['name']}] cv={pair['cv'].name} dent={pair['dent'].name} ===")
            rc = run_dent_match(
                {"cv": str(pair["cv"]), "dent": str(pair["dent"])},
                fmt="sel_all",
                filename_str="sel_all",
                summary_csv=str(match_cache / f"dent_matched_summary-sel_all-{pair['name']}.csv"),
                matched_out_dir=str(matched_out),
                matched_suffix=suffix,
                common_keys_pkl=str(keys_pkl),
                n_workers=args.n_workers,
                file_timeout=args.file_timeout,
                meta_retries=args.meta_retries,
            )
            log(f"match [{pair['name']}] exit={rc}")
            if rc != 0:
                return int(rc)
            with open(keys_pkl, "rb") as fh:
                pair_key_sets.append(set(pickle.load(fh)))

        # Union across *all* pair key pkls on disk (so --only-pair still updates the global set).
        all_pair_sets = []
        for spec in DEFAULT_PAIRS:
            pkl = cache / f"dent_common_keys_sel_all_{spec['name']}.pkl"
            if pkl.is_file():
                with open(pkl, "rb") as fh:
                    all_pair_sets.append(set(pickle.load(fh)))
        if not all_pair_sets:
            all_pair_sets = pair_key_sets
        union_keys = set.union(*all_pair_sets) if all_pair_sets else set()
        union_pkl = cache / "dent_common_keys_sel_all.pkl"
        with open(union_pkl, "wb") as fh:
            pickle.dump(union_keys, fh, protocol=pickle.HIGHEST_PROTOCOL)
        log(
            f"union common keys {len(union_keys):,} ← "
            + " + ".join(f"{len(s):,}" for s in all_pair_sets)
        )
    else:
        log("skip match — using existing matched files / common-keys pkl")

    dent_dirs = {"CV": matched_out / "cv", "DENT": matched_out / "dent"}
    variations = {
        lab: glob_matched_dfs(path, filename_str="sel_all", matched_suffix=suffix)
        for lab, path in dent_dirs.items()
    }
    for lab, files in variations.items():
        log(f"{lab}: {len(files)} matched files")
    if not variations["CV"] or not variations["DENT"]:
        raise RuntimeError("no matched files — run without --skip-match first")

    union_pkl = cache / "dent_common_keys_sel_all.pkl"
    assert_variations_matched(
        variations,
        max_files_per_var=args.assert_max_files,
        format="sel_all",
        common_keys_path=union_pkl if union_pkl.is_file() else None,
    )

    cache_path = cache / "dent_sel_all_products.pkl"

    def _build():
        files_cv, files_dent = variations["CV"], variations["DENT"]
        log(f"walk matched sel_all: CV={len(files_cv)} DENT={len(files_dent)}")
        products = accumulate_matched_sel_all_cv_vs_var_products(files_cv, files_dent)
        pot_by = {"CV": products["cv"]["pot"], "DENT": products["var"]["pot"]}
        scales = pot_scales_to_cv(pot_by)
        sc = float(scales.get("DENT", 1.0))
        for key in ("hists_cut", "hists_final"):
            products["var"][key] = {k: (v * sc) for k, v in products["var"][key].items()}
        return {
            "products": products,
            "pot_by_variation": pot_by,
            "scales": scales,
            "match_stage": "sel_all",
            "pairs": [
                {"name": pair["name"], "cv": str(pair["cv"]), "dent": str(pair["dent"])}
                for pair in pairs
            ],
        }

    payload = load_or_build_cache(cache_path, _build, force=args.force_hist or not args.skip_match)
    # Always record the full campaign pair list in products/manifests.
    payload["pairs"] = [
        {
            "name": spec["name"],
            "cv": str(dfs / spec["cv"]),
            "dent": str(dfs / spec["dent"]),
        }
        for spec in DEFAULT_PAIRS
    ]
    products = payload["products"]
    log(f"POT scales: {payload.get('scales')}")

    dict_a = build_dent_detector_dict(
        products["cv"]["hists_cut"],
        products["var"]["hists_cut"],
        list(products["cv"]["hists_cut"].keys()),
    )
    dict_b = build_dent_detector_dict(
        products["cv"]["hists_final"],
        products["var"]["hists_final"],
        list(products["cv"]["hists_final"].keys()),
    )

    npz_a = out_base / SUB_DETECTOR_SEL / FILE_DETECTOR_SEL
    npz_b = out_base / SUB_DETECTOR / FILE_DETECTOR
    save_detector_npz(
        dict_a,
        npz_a,
        manifest={
            "source": "DENT",
            "product": "A_selection",
            "match_stage": "sel_all",
            "cache": str(cache_path),
            "n_vars": len(dict_a.get("detector-DENT", {})),
            "pairs": payload.get("pairs"),
        },
    )
    save_detector_npz(
        dict_b,
        npz_b,
        manifest={
            "source": "DENT",
            "product": "B_measurement",
            "match_stage": "sel_all",
            "cache": str(cache_path),
            "n_vars": len(dict_b.get("detector-DENT", {})),
            "pairs": payload.get("pairs"),
        },
    )
    log(f"Product A vars: {len(dict_a.get('detector-DENT', {}))}")
    log(f"Product B vars: {len(dict_b.get('detector-DENT', {}))}")
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
