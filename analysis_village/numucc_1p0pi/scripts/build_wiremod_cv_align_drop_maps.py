#!/usr/bin/env python3
"""Align WireMod matched samples to the (incomplete) chi2_new CV event set.

The updatecalo-cvonly CV rematch is missing ~3.3% of the original common keys.
Patching those events is not practical (nearly every updatecalo file is short).
This script:

1. Collects keep artkeys ``(run, subrun, evt)`` from matched CV.
2. Builds walk drop-maps for YZ / XTXW that drop any hdr row whose artkey is
   not in that keep set.
3. For XTXW, also first-claim-dedupes remaining artkeys (same as before).

Walks then see the same event set as CV.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.syst_detvar_common import get_n_split, glob_matched_dfs, log

ArtKey = Tuple[int, int, int]
EntryKey = Tuple[int, int]

warnings.filterwarnings(
    "ignore",
    message="object name is not a valid Python identifier",
    module="tables",
)


def _read_file_hdrs(fpath: str):
    try:
        n_split = get_n_split(fpath)
        out = {}
        with pd.HDFStore(fpath, "r") as store:
            keys = {k.lstrip("/") for k in store.keys()}
            for i in range(n_split):
                key = f"hdr_{i}"
                if key not in keys:
                    continue
                hr = store[key].reset_index()
                rows = list(
                    zip(
                        hr["run"].astype(int),
                        hr["subrun"].astype(int),
                        hr["evt"].astype(int),
                        hr["__ntuple"].astype(int),
                        hr["entry"].astype(int),
                    )
                )
                out[i] = [((r, sr, ev), (nt, en)) for r, sr, ev, nt, en in rows]
        return fpath, out, None
    except Exception as ex:
        return fpath, None, str(ex)


def _collect_artkeys(files: List[str], *, workers: int, label: str) -> Set[ArtKey]:
    arts: Set[ArtKey] = set()
    n_err = 0
    with ProcessPoolExecutor(max_workers=max(workers, 1)) as ex:
        futs = [ex.submit(_read_file_hdrs, f) for f in files]
        for k, fut in enumerate(as_completed(futs), 1):
            fpath, data, err = fut.result()
            if err or data is None:
                n_err += 1
                if n_err <= 5:
                    log(f"  [{label}] SKIP {Path(fpath).name}: {err}")
                continue
            for rows in data.values():
                for art, _ in rows:
                    arts.add(art)
            if k % 200 == 0 or k == len(files):
                log(f"  [{label}] read {k}/{len(files)} arts={len(arts)} err={n_err}")
    return arts


def _load_hdrs(files: List[str], *, workers: int, label: str) -> Dict[str, dict]:
    hdrs: Dict[str, dict] = {}
    n_err = 0
    with ProcessPoolExecutor(max_workers=max(workers, 1)) as ex:
        futs = [ex.submit(_read_file_hdrs, f) for f in files]
        for k, fut in enumerate(as_completed(futs), 1):
            fpath, data, err = fut.result()
            if err or data is None:
                n_err += 1
                if n_err <= 5:
                    log(f"  [{label}] SKIP {Path(fpath).name}: {err}")
            else:
                hdrs[fpath] = data
            if k % 200 == 0 or k == len(files):
                log(f"  [{label}] read {k}/{len(files)} ok={len(hdrs)} err={n_err}")
    return hdrs


def _build_drop_map(
    files: List[str],
    hdrs: Dict[str, dict],
    keep_art: Set[ArtKey],
    *,
    dedupe: bool,
) -> dict:
    claimed: Set[ArtKey] = set()
    drop_map: Dict[str, Dict[int, Set[EntryKey]]] = {}
    n_drop = n_keep = n_files_with_drops = 0
    n_drop_not_in_cv = n_drop_dup = 0

    for fpath in files:
        data = hdrs.get(fpath)
        if not data:
            continue
        file_drop: Dict[int, Set[EntryKey]] = {}
        n_before = 0
        for i, rows in data.items():
            n_before += len(rows)
            drop_entries: Set[EntryKey] = set()
            for art, entry in rows:
                if art not in keep_art:
                    drop_entries.add(entry)
                    n_drop_not_in_cv += 1
                    continue
                if dedupe:
                    if art in claimed:
                        drop_entries.add(entry)
                        n_drop_dup += 1
                    else:
                        claimed.add(art)
                else:
                    claimed.add(art)
            if drop_entries:
                file_drop[i] = drop_entries
        n_file_drop = sum(len(s) for s in file_drop.values())
        n_drop += n_file_drop
        n_keep += n_before - n_file_drop
        if file_drop:
            n_files_with_drops += 1
            drop_map[os.path.abspath(fpath)] = file_drop

    return {
        "drop_map": drop_map,
        "n_files": len(files),
        "n_files_with_drops": n_files_with_drops,
        "n_drop": n_drop,
        "n_keep": n_keep,
        "n_drop_not_in_cv": n_drop_not_in_cv,
        "n_drop_dup": n_drop_dup,
        "n_unique_claimed": len(claimed),
        "kind": "cv_align" + ("+dedupe" if dedupe else ""),
    }


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-base", required=True)
    p.add_argument("--workers", type=int, default=12)
    args = p.parse_args(argv)

    out = Path(args.out_base)
    cache = out / "cache"
    cache.mkdir(parents=True, exist_ok=True)

    cv_files = glob_matched_dfs(out / "matched" / "cv", filename_str="sel_all")
    yz_files = glob_matched_dfs(out / "matched" / "yz", filename_str="sel_all")
    xtxw_files = glob_matched_dfs(out / "matched" / "xtxw", filename_str="sel_all")
    log(f"files cv={len(cv_files)} yz={len(yz_files)} xtxw={len(xtxw_files)}")

    keep_art = _collect_artkeys(cv_files, workers=args.workers, label="cv")
    keep_pkl = cache / "wiremod_keep_artkeys_cv.pkl"
    with open(keep_pkl, "wb") as fh:
        pickle.dump(keep_art, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"keep_art from CV: {len(keep_art)} -> {keep_pkl}")

    # Document restricted common keys (EventKey set) if original common exists.
    common_path = cache / "wiremod_common_keys_sel_all.pkl"
    if common_path.is_file():
        with open(common_path, "rb") as fh:
            common = pickle.load(fh)
        restricted = {k for k in common if (k[1], k[2], k[3]) in keep_art}
        rest_path = cache / "wiremod_common_keys_cv_aligned.pkl"
        with open(rest_path, "wb") as fh:
            pickle.dump(restricted, fh, protocol=pickle.HIGHEST_PROTOCOL)
        log(
            f"common={len(common)} cv_aligned_eventkeys={len(restricted)} "
            f"dropped={len(common) - len(restricted)} -> {rest_path}"
        )

    yz_hdrs = _load_hdrs(yz_files, workers=args.workers, label="yz")
    yz_payload = _build_drop_map(yz_files, yz_hdrs, keep_art, dedupe=False)
    yz_out = cache / "wiremod_yz_cv_align_drop_map.pkl"
    with open(yz_out, "wb") as fh:
        pickle.dump(yz_payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(
        f"YZ drop_map: drop={yz_payload['n_drop']} keep={yz_payload['n_keep']} "
        f"not_in_cv={yz_payload['n_drop_not_in_cv']} -> {yz_out}"
    )

    xtxw_hdrs = _load_hdrs(xtxw_files, workers=args.workers, label="xtxw")
    xtxw_payload = _build_drop_map(xtxw_files, xtxw_hdrs, keep_art, dedupe=True)
    xtxw_out = cache / "wiremod_xtxw_cv_align_drop_map.pkl"
    # Also overwrite the canonical XTXW drop-map name used by launch scripts.
    xtxw_canon = cache / "wiremod_xtxw_drop_map.pkl"
    with open(xtxw_out, "wb") as fh:
        pickle.dump(xtxw_payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    with open(xtxw_canon, "wb") as fh:
        pickle.dump(xtxw_payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(
        f"XTXW drop_map: drop={xtxw_payload['n_drop']} keep={xtxw_payload['n_keep']} "
        f"not_in_cv={xtxw_payload['n_drop_not_in_cv']} dup={xtxw_payload['n_drop_dup']} "
        f"claimed={xtxw_payload['n_unique_claimed']} -> {xtxw_out} (+ {xtxw_canon.name})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
