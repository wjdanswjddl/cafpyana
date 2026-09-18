#!/usr/bin/env python3
"""Build a cross-file artkey drop map for matched WireMod dfs (hdr-only, fast).

Phase 1: parallel hdr reads. Phase 2: sequential claim in sorted file order.
Writes a pickle used by ``wiremod_walk_shard.py --drop-map-pkl``.
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
    """Return (fpath, {split_i: [(art, nt, en), ...]}) or (fpath, None, err)."""
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


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matched-dir", required=True)
    p.add_argument("--drop-map-out", required=True)
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args(argv)

    files = glob_matched_dfs(args.matched_dir, filename_str="sel_all")
    log(f"dedupe-scan {len(files)} files workers={args.workers}")

    hdrs: Dict[str, dict] = {}
    n_err = 0
    with ProcessPoolExecutor(max_workers=max(int(args.workers), 1)) as ex:
        futs = [ex.submit(_read_file_hdrs, f) for f in files]
        for k, fut in enumerate(as_completed(futs), 1):
            fpath, data, err = fut.result()
            if err or data is None:
                n_err += 1
                if n_err <= 5:
                    log(f"  SKIP {Path(fpath).name}: {err}")
            else:
                hdrs[fpath] = data
            if k % 200 == 0 or k == len(files):
                log(f"  read {k}/{len(files)} ok={len(hdrs)} err={n_err}")

    claimed: Set[ArtKey] = set()
    drop_map: Dict[str, Dict[int, Set[EntryKey]]] = {}
    n_drop = 0
    n_keep = 0
    n_files_with_drops = 0

    for fi, fpath in enumerate(files):
        data = hdrs.get(fpath)
        if not data:
            continue
        file_drop: Dict[int, Set[EntryKey]] = {}
        n_before = 0
        for i, rows in data.items():
            n_before += len(rows)
            drop_entries: Set[EntryKey] = set()
            for art, entry in rows:
                if art in claimed:
                    drop_entries.add(entry)
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

    payload = {
        "drop_map": drop_map,
        "n_files": len(files),
        "n_files_with_drops": n_files_with_drops,
        "n_drop": n_drop,
        "n_keep": n_keep,
        "n_unique_claimed": len(claimed),
        "matched_dir": os.path.abspath(args.matched_dir),
    }
    out = Path(args.drop_map_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(
        f"done keep={n_keep} drop={n_drop} files_with_drops={n_files_with_drops} "
        f"unique_claimed={len(claimed)} wrote {out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
