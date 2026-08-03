#!/usr/bin/env python3
"""
Find events common to multiple WireMod (or other) MC variations and write per-file
``_matched`` HDF5 outputs.

Mirrors ``notebooks/wiremod.ipynb`` but supports an arbitrary number of variations
(typically three: e.g. YZ, X-ThetaXW, and CV/nominal).  Memory use stays bounded:
phase 1 scans only ``meta``; phase 3 processes one file split at a time and never
concatenates across input files.

Usage
-----
    python wiremod_match_common_events.py \\
        --use-default-variations \\
        --variation cv /pnfs/.../THIRD_SAMPLE_DIR \\
        --common-keys-pkl /pnfs/.../wiremod_common_keys.pkl

    # or pass all three explicitly:
    python wiremod_match_common_events.py \\
        --variation yz  /pnfs/.../WireModYZ/merged_perTPC \\
        --variation xtxw /pnfs/.../WireModXTXW \\
        --variation cv  /pnfs/.../THIRD_SAMPLE_DIR

    # meta scan only (write common-key pickle for a later write pass)
    python wiremod_match_common_events.py --phase meta ... --common-keys-pkl common.pkl

    # write pass using a saved key set
    python wiremod_match_common_events.py --phase write --common-keys-pkl common.pkl ...
"""

from __future__ import annotations

import argparse
import glob
import os
import pickle
import sys
import warnings
from os import path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from pyanalib.split_df_helpers_new import get_n_split, load_dfs

EventKey = Tuple[float, int, int, int]  # (E, run, subrun, evt)

DEFAULT_KEYS2LOAD = [
    "meta",
    "evt_cv",
    "evt_ccal_p",
    "evt_ccal_m",
    "evt_alpha_p",
    "evt_alpha_m",
    "evt_beta_p",
    "evt_beta_m",
    "evt_R_p",
    "evt_R_m",
]

DEFAULT_VARIATIONS = {
    "yz": "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModYZ/merged_perTPC",
    "xtxw": "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_05_16_161556__sel_2prong-mc-BNB_cosmics-WireModXTXW",
    # set the third sample directory for your CV / nominal / calovar production:
    "cv": None,
}


def matched_out_path(fpath: str, suffix: str = "_matched") -> str:
    base, ext = path.splitext(fpath)
    return f"{base}{suffix}{ext}"


def list_df_files(search_dir: str, filename_str: str) -> List[str]:
    pattern = path.join(search_dir, f"*{filename_str}*.df")
    files = sorted(glob.glob(pattern))
    # skip outputs from a previous run
    files = [f for f in files if "_matched" not in path.basename(f)]
    return files


def collect_event_keys(files: Sequence[str]) -> Set[EventKey]:
    """Load only ``meta`` per file; return unique (E, run, subrun, evt) tuples."""
    keys: Set[EventKey] = set()
    for fpath in tqdm(files, desc="meta scan"):
        try:
            meta = load_dfs(fpath, ["meta"], n_max_concat=10000)["meta"]
        except Exception as exc:
            print(f"Error loading meta from {fpath}: {exc}", flush=True)
            continue
        meta_r = meta.reset_index()
        keys.update(map(tuple, meta_r[["E", "run", "subrun", "evt"]].values))
    return keys


def intersect_event_keys(per_variation: Dict[str, Set[EventKey]]) -> Set[EventKey]:
    names = list(per_variation.keys())
    common = set(per_variation[names[0]])
    for name in names[1:]:
        common &= per_variation[name]
    return common


def _filter_meta_split(meta_df: pd.DataFrame, common_event_keys: Set[EventKey]):
    meta_r = meta_df.reset_index()
    mask = np.array(
        [
            (e, r, sr, ev) in common_event_keys
            for e, r, sr, ev in zip(meta_r["E"], meta_r["run"], meta_r["subrun"], meta_r["evt"])
        ],
        dtype=bool,
    )
    if not mask.any():
        return None, None
    match_idx = pd.MultiIndex.from_arrays(
        [meta_r.loc[mask, "__ntuple"], meta_r.loc[mask, "entry"], meta_r.loc[mask, "E"]],
        names=["__ntuple", "entry", "E"],
    )
    meta_out = meta_df.iloc[np.flatnonzero(mask)]
    return meta_out, match_idx


def _filter_evt_split(evt_df: pd.DataFrame, match_idx: pd.MultiIndex):
    evt_r = evt_df.reset_index()
    evt_r["E"] = evt_df.mc.E.values
    evt_idx = pd.MultiIndex.from_arrays(
        [evt_r["__ntuple"], evt_r["entry"], evt_r["E"]],
        names=["__ntuple", "entry", "E"],
    )
    mask = evt_idx.isin(match_idx)
    if not mask.any():
        return None
    return evt_df.iloc[np.flatnonzero(mask)]


def save_matched_files(
    files: Sequence[str],
    common_event_keys: Set[EventKey],
    keys2load: Sequence[str],
    *,
    matched_suffix: str = "_matched",
    variation_name: str = "",
) -> pd.DataFrame:
    keys_evt = [k for k in keys2load if k != "meta"]
    summary = []
    desc = f"write {variation_name}" if variation_name else "write matched"

    for fpath in tqdm(files, desc=desc):
        out_path = matched_out_path(fpath, suffix=matched_suffix)
        n_split = get_n_split(fpath)
        written_splits = []

        for i in range(n_split):
            try:
                meta = pd.read_hdf(fpath, key=f"meta_{i}")
            except Exception as exc:
                print(f"Error loading meta_{i} from {fpath}: {exc}", flush=True)
                continue

            meta_out, match_idx = _filter_meta_split(meta, common_event_keys)
            if meta_out is None:
                continue

            split_out = {"meta": meta_out}
            n_evt_cv = 0

            for key in keys_evt:
                try:
                    evt = pd.read_hdf(fpath, key=f"{key}_{i}")
                except Exception as exc:
                    print(f"Error loading {key}_{i} from {fpath}: {exc}", flush=True)
                    continue
                evt_out = _filter_evt_split(evt, match_idx)
                if evt_out is not None:
                    split_out[key] = evt_out
                    if key == "evt_cv":
                        n_evt_cv += len(evt_out)

            written_splits.append((split_out, len(meta_out), n_evt_cv))

        if not written_splits:
            continue

        if path.exists(out_path):
            os.remove(out_path)

        with pd.HDFStore(out_path, mode="w") as store:
            store.put("split", pd.DataFrame({"n_split": [len(written_splits)]}), format="fixed")
            for out_i, (split_out, _, _) in enumerate(written_splits):
                for key, df in split_out.items():
                    store.put(f"{key}_{out_i}", df, format="fixed")

        summary.append(
            {
                "variation": variation_name,
                "input": fpath,
                "output": out_path,
                "n_splits_written": len(written_splits),
                "n_meta_rows": sum(s[1] for s in written_splits),
                "n_evt_cv_rows": sum(s[2] for s in written_splits),
            }
        )

    return pd.DataFrame(summary)


def parse_variations(
    variation_args: Sequence[Sequence[str]] | None,
    *,
    use_defaults: bool,
) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if use_defaults:
        out.update({k: v for k, v in DEFAULT_VARIATIONS.items() if v})

    if variation_args:
        for pair in variation_args:
            if len(pair) != 2:
                raise ValueError("--variation requires pairs: NAME DIR")
            name, dir_path = pair
            out[name] = dir_path

    if len(out) < 2:
        raise ValueError(
            "need at least two variations: pass --variation NAME DIR twice, or "
            "--use-default-variations plus --variation cv DIR for the third sample"
        )
    return out


def run_meta_phase(
    variations: Dict[str, str],
    filename_str: str,
    max_files: int | None,
) -> Tuple[Dict[str, Set[EventKey]], Set[EventKey], Dict[str, List[str]]]:
    per_var_keys: Dict[str, Set[EventKey]] = {}
    file_lists: Dict[str, List[str]] = {}

    for name, search_dir in variations.items():
        files = list_df_files(search_dir, filename_str)
        if max_files is not None:
            files = files[:max_files]
        file_lists[name] = files
        print(f"[{name}] {len(files)} files under {search_dir}", flush=True)
        per_var_keys[name] = collect_event_keys(files)
        print(f"[{name}] unique events: {len(per_var_keys[name])}", flush=True)

    common = intersect_event_keys(per_var_keys)
    print(f"common events (all {len(variations)} variations): {len(common)}", flush=True)
    for name, keys in per_var_keys.items():
        only_here = len(keys - common)
        print(f"  {name} only: {only_here}", flush=True)

    return per_var_keys, common, file_lists


def run_write_phase(
    variations: Dict[str, str],
    file_lists: Dict[str, List[str]],
    common_keys: Set[EventKey],
    keys2load: Sequence[str],
    matched_suffix: str,
) -> pd.DataFrame:
    summaries = []
    for name, files in file_lists.items():
        if not files:
            files = list_df_files(variations[name], "")
            files = [f for f in files if "_matched" not in path.basename(f)]
        print(f"\nWriting matched files for [{name}] ...", flush=True)
        summary = save_matched_files(
            files,
            common_keys,
            keys2load,
            matched_suffix=matched_suffix,
            variation_name=name,
        )
        summaries.append(summary)
        if len(summary):
            print(
                f"  [{name}] files written: {len(summary)}, "
                f"evt_cv rows: {int(summary['n_evt_cv_rows'].sum())}",
                flush=True,
            )
        else:
            print(f"  [{name}] no matched files written", flush=True)
    if not summaries:
        return pd.DataFrame()
    return pd.concat(summaries, ignore_index=True)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Select events common to multiple MC variations and write _matched .df files.",
    )
    p.add_argument(
        "--variation",
        nargs=2,
        action="append",
        metavar=("NAME", "DIR"),
        help="Variation label and directory containing *.df files (repeat for each sample).",
    )
    p.add_argument(
        "--use-default-variations",
        action="store_true",
        help="Pre-fill yz and xtxw dirs from this script; add others with --variation.",
    )
    p.add_argument(
        "--filename-str",
        default="sel_2prong",
        help="Substring matched in input .df filenames (default: sel_2prong).",
    )
    p.add_argument(
        "--keys",
        default=",".join(DEFAULT_KEYS2LOAD),
        help="Comma-separated HDF keys to copy into matched files.",
    )
    p.add_argument(
        "--phase",
        choices=("all", "meta", "write"),
        default="all",
        help="Run meta scan only, write only, or both (default: all).",
    )
    p.add_argument(
        "--common-keys-pkl",
        default=None,
        help="Pickle path for the common (E, run, subrun, evt) set (written in meta phase).",
    )
    p.add_argument(
        "--summary-csv",
        default=None,
        help="Optional CSV path for per-file write summary.",
    )
    p.add_argument(
        "--matched-suffix",
        default="_matched",
        help="Suffix inserted before .df on output files (default: _matched).",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Cap number of input files per variation (for tests).",
    )
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    variations = parse_variations(args.variation, use_defaults=args.use_default_variations)
    if len(variations) < 2:
        raise SystemExit("need at least two --variation entries")

    keys2load = [k.strip() for k in args.keys.split(",") if k.strip()]
    if "meta" not in keys2load:
        keys2load = ["meta"] + keys2load

    common_keys: Set[EventKey] | None = None
    file_lists: Dict[str, List[str]] = {}

    if args.phase in ("all", "meta"):
        _, common_keys, file_lists = run_meta_phase(
            variations, args.filename_str, args.max_files
        )
        if args.common_keys_pkl:
            with open(args.common_keys_pkl, "wb") as fh:
                pickle.dump(common_keys, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"wrote common keys to {args.common_keys_pkl}", flush=True)

    if args.phase in ("all", "write"):
        if args.phase == "write":
            if not args.common_keys_pkl or not path.isfile(args.common_keys_pkl):
                raise SystemExit("--phase write requires an existing --common-keys-pkl")
            with open(args.common_keys_pkl, "rb") as fh:
                common_keys = pickle.load(fh)
            print(f"loaded {len(common_keys)} common keys from {args.common_keys_pkl}", flush=True)
            for name, search_dir in variations.items():
                files = list_df_files(search_dir, args.filename_str)
                if args.max_files is not None:
                    files = files[: args.max_files]
                file_lists[name] = files

        assert common_keys is not None
        summary = run_write_phase(
            variations,
            file_lists,
            common_keys,
            keys2load,
            args.matched_suffix,
        )
        if args.summary_csv and len(summary):
            summary.to_csv(args.summary_csv, index=False)
            print(f"wrote summary to {args.summary_csv}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
