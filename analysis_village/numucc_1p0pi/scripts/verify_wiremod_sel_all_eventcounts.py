#!/usr/bin/env python3
"""Verify sel_all matched unique (run,subrun,evt) counts for CV / YZ / XTXW.

Hdr-only parallel scan. For XTXW also applies drop-map to get effective walk count.
Writes JSON summary under WireMod/cache/.
"""
from __future__ import annotations

import json
import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import pandas as pd

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.syst_detvar_common import get_n_split, glob_matched_dfs, log

ArtKey = Tuple[int, int, int]
OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"


def _hdr_one(fpath: str):
    art: Set[ArtKey] = set()
    n_hdr = 0
    try:
        n_split = get_n_split(fpath)
    except Exception as ex:
        return fpath, None, str(ex)
    try:
        with pd.HDFStore(fpath, "r") as store:
            keys = {k.lstrip("/") for k in store.keys()}
            for i in range(n_split):
                key = f"hdr_{i}"
                if key not in keys:
                    continue
                hr = store[key].reset_index()
                n_hdr += len(hr)
                for r, sr, ev in zip(
                    hr["run"].astype(int), hr["subrun"].astype(int), hr["evt"].astype(int)
                ):
                    art.add((int(r), int(sr), int(ev)))
    except Exception as ex:
        return fpath, None, str(ex)
    return fpath, (n_hdr, art), None


def scan_lab(lab: str, workers: int = 8):
    files = glob_matched_dfs(OUT / "matched" / lab, filename_str="sel_all")
    log(f"[{lab}] scanning {len(files)} files workers={workers}")
    art: Set[ArtKey] = set()
    n_hdr = 0
    n_ok = 0
    n_err = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_hdr_one, f) for f in files]
        for k, fut in enumerate(as_completed(futs), 1):
            fpath, data, err = fut.result()
            if err or data is None:
                n_err += 1
                if n_err <= 3:
                    log(f"  [{lab}] SKIP {Path(fpath).name}: {err}")
            else:
                nh, a = data
                n_hdr += nh
                art |= a
                n_ok += 1
            if k % 500 == 0 or k == len(files):
                log(f"  [{lab}] {k}/{len(files)} ok={n_ok} unique={len(art)} hdr={n_hdr}")
    return {
        "lab": lab,
        "n_files": len(files),
        "n_ok": n_ok,
        "n_err": n_err,
        "n_hdr_rows": n_hdr,
        "n_unique_art": len(art),
        "art": art,
    }


def main() -> int:
    with open(CACHE / "wiremod_common_keys_sel_all.pkl", "rb") as fh:
        common = pickle.load(fh)
    # common keys are (E,run,subrun,evt); project to art
    common_art = {(int(r), int(sr), int(ev)) for (_e, r, sr, ev) in common}
    log(f"common_keys={len(common)} common_art={len(common_art)}")

    drop_payload = None
    drop_path = CACHE / "wiremod_xtxw_drop_map.pkl"
    if drop_path.is_file():
        with open(drop_path, "rb") as fh:
            drop_payload = pickle.load(fh)
        log(
            f"XTXW drop_map unique_claimed={drop_payload.get('n_unique_claimed')} "
            f"n_keep={drop_payload.get('n_keep')} n_drop={drop_payload.get('n_drop')}"
        )

    results = {}
    for lab in ("cv", "yz", "xtxw"):
        results[lab] = scan_lab(lab, workers=8)

    summary: Dict = {
        "common_keys": len(common),
        "common_art": len(common_art),
        "by_lab": {},
        "set_diffs_art": {},
        "equal_unique_art": None,
        "xtxw_effective": None,
    }
    for lab, r in results.items():
        summary["by_lab"][lab] = {
            "n_files": r["n_files"],
            "n_ok": r["n_ok"],
            "n_err": r["n_err"],
            "n_hdr_rows": r["n_hdr_rows"],
            "n_unique_art": r["n_unique_art"],
            "hdr_minus_unique": r["n_hdr_rows"] - r["n_unique_art"],
            "vs_common_art": {
                "only_lab": len(r["art"] - common_art),
                "only_common": len(common_art - r["art"]),
                "intersection": len(r["art"] & common_art),
            },
        }
        log(
            f"[{lab}] hdr={r['n_hdr_rows']} unique_art={r['n_unique_art']} "
            f"dup_rows={r['n_hdr_rows'] - r['n_unique_art']} "
            f"|lab-common|={len(r['art'] - common_art)} |common-lab|={len(common_art - r['art'])}"
        )

    for a, b in (("cv", "yz"), ("cv", "xtxw"), ("yz", "xtxw")):
        A, B = results[a]["art"], results[b]["art"]
        d = {
            "a": a,
            "b": b,
            "n_a": len(A),
            "n_b": len(B),
            "only_a": len(A - B),
            "only_b": len(B - A),
            "intersection": len(A & B),
        }
        summary["set_diffs_art"][f"{a}_vs_{b}"] = d
        log(
            f"{a} vs {b}: |A|={d['n_a']} |B|={d['n_b']} "
            f"|A-B|={d['only_a']} |B-A|={d['only_b']} |A&B|={d['intersection']}"
        )

    u_cv = results["cv"]["n_unique_art"]
    u_yz = results["yz"]["n_unique_art"]
    u_xtxw = results["xtxw"]["n_unique_art"]
    # XTXW unique set should already equal keep count; drop-map only removes duplicate rows
    xtxw_eff = drop_payload.get("n_unique_claimed") if drop_payload else u_xtxw
    summary["xtxw_effective"] = {
        "unique_art_on_disk": u_xtxw,
        "drop_map_unique_claimed": xtxw_eff,
        "drop_map_n_keep": None if not drop_payload else drop_payload.get("n_keep"),
        "drop_map_n_drop": None if not drop_payload else drop_payload.get("n_drop"),
    }
    summary["equal_unique_art"] = u_cv == u_yz == u_xtxw == len(common_art)
    summary["equal_unique_art_values"] = {
        "cv": u_cv,
        "yz": u_yz,
        "xtxw": u_xtxw,
        "common_art": len(common_art),
    }

    out_json = CACHE / "wiremod_sel_all_eventcount_verify.json"
    # strip non-serializable
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    log(f"wrote {out_json}")
    log(f"EQUAL unique art across CV/YZ/XTXW/common? {summary['equal_unique_art']}")
    print(json.dumps(summary["equal_unique_art_values"], indent=2))
    print(json.dumps(summary["set_diffs_art"], indent=2))
    return 0 if summary["equal_unique_art"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
