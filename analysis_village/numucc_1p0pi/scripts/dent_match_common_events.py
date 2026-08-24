#!/usr/bin/env python3
"""
Find events common to CV and DENT MC variations and write per-file ``_matched`` HDF5
outputs.

DENT is a detector unisim (like WireMod). Unlike WireMod / SCE comparisons that
start from ``sel_mup`` or ``sel_2prong``, the DENT study matches events at
``sel_all`` so early selection variables can be compared.

Two input layouts are supported:

* ``sel_all`` — raw ``evt`` / ``trk`` / ``hdr`` tables (no ``meta``). Event keys
  are built from ``hdr`` (run, subrun, evt) plus generator neutrino energy ``E`` from
  ``evt.mc``.
* ``sel_mup`` — final-selection layout with ``meta`` and ``evt_cv`` (delegates to
  ``wiremod_match_common_events``).

Usage
-----
    # sel_all (primary DENT workflow)
    python dent_match_common_events.py --format sel_all \\
        --variation cv   /pnfs/.../2026_08_18_125707__sel_all-mc-CV \\
        --variation dent /pnfs/.../2026_08_18_125707__sel_all-mc-DENT \\
        --summary-csv /pnfs/.../dent_matched_summary-sel_all.csv

    # sel_mup (final-selected variables, optional)
    python dent_match_common_events.py --format sel_mup \\
        --variation cv   /pnfs/.../2026_08_18_130158__sel_mup-mc-CV \\
        --variation dent /pnfs/.../2026_08_18_130158__sel_mup-mc-DENT \\
        --filename-str sel_mup \\
        --summary-csv /pnfs/.../dent_matched_summary-sel_mup.csv
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

from pyanalib.split_df_helpers_new import get_n_split

from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series

EventKey = Tuple[float, int, int, int]  # (E, run, subrun, evt)

_DFS_ROOT = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

DEFAULT_VARIATIONS_SEL_ALL = {
    "cv": f"{_DFS_ROOT}/2026_08_19_031254__sel_all-mc-CV",
    "dent": f"{_DFS_ROOT}/2026_08_18_120607__sel_all-mc-DENT",
}

DEFAULT_VARIATIONS_SEL_MUP = {
    "cv": f"{_DFS_ROOT}/2026_08_19_031423__sel_mup-mc-CV",
    "dent": f"{_DFS_ROOT}/2026_08_18_120753__sel_mup-mc-DENT",
}

SEL_ALL_KEYS = ["evt", "trk", "hdr"]
SEL_MUP_KEYS = ["meta", "evt_cv", "evt"]


def matched_out_path(fpath: str, suffix: str = "_matched") -> str:
    base, ext = path.splitext(fpath)
    return f"{base}{suffix}{ext}"


def list_df_files(search_dir: str, filename_str: str) -> List[str]:
    pattern = path.join(search_dir, f"*{filename_str}*.df")
    files = sorted(glob.glob(pattern))
    return [f for f in files if "_matched" not in path.basename(f)]


def _event_energy_map(evt: pd.DataFrame) -> pd.Series:
    """Map (__ntuple, entry) -> generator neutrino energy from evt slices."""
    E = multicol_get_series(evt, ("mc", "E", "", "", "", ""))
    return E.groupby(level=["__ntuple", "entry"]).first().rename("E")


def hdr_event_keys(hdr: pd.DataFrame, evt: pd.DataFrame) -> Set[EventKey]:
    """Unique (E, run, subrun, evt) keys for one split."""
    if hdr is None or len(hdr) == 0 or evt is None or len(evt) == 0:
        return set()
    hdr_r = hdr.reset_index()
    E_map = _event_energy_map(evt)
    hdr_r = hdr_r.merge(
        E_map.reset_index(),
        on=["__ntuple", "entry"],
        how="inner",
    )
    hdr_r = hdr_r.dropna(subset=["E"])
    return set(
        zip(
            hdr_r["E"].astype(float),
            hdr_r["run"].astype(int),
            hdr_r["subrun"].astype(int),
            hdr_r["evt"].astype(int),
        )
    )


def collect_sel_all_event_keys(files: Sequence[str]) -> Set[EventKey]:
    keys: Set[EventKey] = set()
    for fpath in tqdm(files, desc="sel_all meta scan"):
        n_split = get_n_split(fpath)
        for i in range(n_split):
            try:
                hdr = pd.read_hdf(fpath, key=f"hdr_{i}")
                evt = pd.read_hdf(fpath, key=f"evt_{i}")
            except Exception as exc:
                print(f"Error loading split {i} from {fpath}: {exc}", flush=True)
                continue
            keys.update(hdr_event_keys(hdr, evt))
    return keys


def _matched_art_entries(hdr: pd.DataFrame, evt: pd.DataFrame, common_keys: Set[EventKey]):
    """Return boolean mask on hdr rows and matching MultiIndex for evt/trk filtering."""
    if hdr is None or len(hdr) == 0:
        return None, None
    hdr_r = hdr.reset_index()
    E_map = _event_energy_map(evt) if evt is not None and len(evt) else pd.Series(dtype=float)
    hdr_r = hdr_r.merge(
        E_map.reset_index(),
        on=["__ntuple", "entry"],
        how="left",
    )
    mask = np.array(
        [
            (e, r, sr, ev) in common_keys
            for e, r, sr, ev in zip(hdr_r["E"], hdr_r["run"], hdr_r["subrun"], hdr_r["evt"])
        ],
        dtype=bool,
    )
    if not mask.any():
        return None, None
    matched_entries = set(zip(hdr_r.loc[mask, "__ntuple"], hdr_r.loc[mask, "entry"]))
    return mask, matched_entries


def _filter_df_by_entries(df: pd.DataFrame, matched_entries: Set[Tuple[int, int]]):
    if df is None or len(df) == 0:
        return None
    df_r = df.reset_index()
    mask = [
        (nt, en) in matched_entries
        for nt, en in zip(df_r["__ntuple"], df_r["entry"])
    ]
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return None
    return df.iloc[np.flatnonzero(mask)]


def save_matched_sel_all_files(
    files: Sequence[str],
    common_event_keys: Set[EventKey],
    *,
    matched_suffix: str = "_matched",
    variation_name: str = "",
) -> pd.DataFrame:
    summary = []
    desc = f"write sel_all {variation_name}" if variation_name else "write sel_all matched"

    for fpath in tqdm(files, desc=desc):
        out_path = matched_out_path(fpath, suffix=matched_suffix)
        n_split = get_n_split(fpath)
        written_splits = []

        for i in range(n_split):
            try:
                hdr = pd.read_hdf(fpath, key=f"hdr_{i}")
                evt = pd.read_hdf(fpath, key=f"evt_{i}")
            except Exception as exc:
                print(f"Error loading split {i} from {fpath}: {exc}", flush=True)
                continue

            _, matched_entries = _matched_art_entries(hdr, evt, common_event_keys)
            if matched_entries is None:
                continue

            hdr_out = _filter_df_by_entries(hdr, matched_entries)
            evt_out = _filter_df_by_entries(evt, matched_entries)
            if hdr_out is None or evt_out is None:
                continue

            split_out = {"hdr": hdr_out, "evt": evt_out}
            try:
                trk = pd.read_hdf(fpath, key=f"trk_{i}")
                trk_out = _filter_df_by_entries(trk, matched_entries)
                if trk_out is not None:
                    split_out["trk"] = trk_out
            except Exception:
                pass

            written_splits.append((split_out, len(hdr_out), len(evt_out)))

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
                "n_hdr_rows": sum(s[1] for s in written_splits),
                "n_evt_rows": sum(s[2] for s in written_splits),
            }
        )

    return pd.DataFrame(summary)


def intersect_event_keys(per_variation: Dict[str, Set[EventKey]]) -> Set[EventKey]:
    names = list(per_variation.keys())
    common = set(per_variation[names[0]])
    for name in names[1:]:
        common &= per_variation[name]
    return common


def parse_variations(
    variation_args: Sequence[Sequence[str]] | None,
    *,
    use_defaults: bool,
    default_variations: Dict[str, str],
) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if use_defaults:
        out.update(default_variations)
    if variation_args:
        for pair in variation_args:
            if len(pair) != 2:
                raise ValueError("--variation requires pairs: NAME DIR")
            name, dir_path = pair
            out[name] = dir_path
    if len(out) < 2:
        raise ValueError("need at least two --variation entries (cv and dent)")
    return out


def run_sel_all_match(args) -> int:
    variations = parse_variations(
        args.variation,
        use_defaults=args.use_default_variations,
        default_variations=DEFAULT_VARIATIONS_SEL_ALL,
    )
    file_lists: Dict[str, List[str]] = {}
    per_var_keys: Dict[str, Set[EventKey]] = {}

    if args.phase in ("all", "meta"):
        for name, search_dir in variations.items():
            files = list_df_files(search_dir, args.filename_str)
            if args.max_files is not None:
                files = files[: args.max_files]
            file_lists[name] = files
            print(f"[{name}] {len(files)} files under {search_dir}", flush=True)
            per_var_keys[name] = collect_sel_all_event_keys(files)
            print(f"[{name}] unique events: {len(per_var_keys[name])}", flush=True)

        common_keys = intersect_event_keys(per_var_keys)
        print(f"common events (CV ∩ DENT): {len(common_keys)}", flush=True)
        if len(common_keys) == 0:
            print(
                "WARNING: zero common events — CV and DENT inputs may be from different "
                "production runs. Paired unisim samples must share the same underlying MC "
                "(e.g. 2026_08_18_125707 for both CV and DENT).",
                flush=True,
            )
        for name, keys in per_var_keys.items():
            print(f"  {name} only: {len(keys - common_keys)}", flush=True)

        if args.common_keys_pkl:
            with open(args.common_keys_pkl, "wb") as fh:
                pickle.dump(common_keys, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"wrote common keys to {args.common_keys_pkl}", flush=True)
    else:
        common_keys = None

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
        if len(common_keys) == 0:
            print("Skipping matched write (empty intersection).", flush=True)
        else:
            summaries = []
            for name, files in file_lists.items():
                print(f"\nWriting matched sel_all files for [{name}] ...", flush=True)
                summary = save_matched_sel_all_files(
                    files,
                    common_keys,
                    matched_suffix=args.matched_suffix,
                    variation_name=name,
                )
                summaries.append(summary)
                if len(summary):
                    print(
                        f"  [{name}] files written: {len(summary)}, "
                        f"evt rows: {int(summary['n_evt_rows'].sum())}",
                        flush=True,
                    )
                else:
                    print(f"  [{name}] no matched files written", flush=True)

            if summaries:
                summary = pd.concat(summaries, ignore_index=True)
                if args.summary_csv and len(summary):
                    summary.to_csv(args.summary_csv, index=False)
                    print(f"wrote summary to {args.summary_csv}", flush=True)

    return 0


def run_sel_mup_match(args) -> int:
    import wiremod_match_common_events as _core

    _core.DEFAULT_KEYS2LOAD = SEL_MUP_KEYS
    _core.DEFAULT_VARIATIONS = DEFAULT_VARIATIONS_SEL_MUP

    argv = []
    if args.use_default_variations:
        argv.append("--use-default-variations")
    if args.variation:
        for pair in args.variation:
            argv.extend(["--variation", pair[0], pair[1]])
    argv.extend(["--filename-str", args.filename_str])
    argv.extend(["--phase", args.phase])
    if args.common_keys_pkl:
        argv.extend(["--common-keys-pkl", args.common_keys_pkl])
    if args.summary_csv:
        argv.extend(["--summary-csv", args.summary_csv])
    argv.extend(["--matched-suffix", args.matched_suffix])
    if args.max_files is not None:
        argv.extend(["--max-files", str(args.max_files)])
    return _core.main(argv)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Match CV and DENT MC events and write _matched .df files.",
    )
    p.add_argument(
        "--format",
        choices=("sel_all", "sel_mup"),
        default="sel_all",
        help="Input layout: sel_all (evt/trk/hdr) or sel_mup (meta/evt_cv).",
    )
    p.add_argument(
        "--variation",
        nargs=2,
        action="append",
        metavar=("NAME", "DIR"),
        help="Variation label and directory (repeat: cv DIR dent DIR).",
    )
    p.add_argument(
        "--use-default-variations",
        action="store_true",
        help="Use bundled CV/DENT directory paths for the chosen --format.",
    )
    p.add_argument(
        "--filename-str",
        default="sel_all",
        help="Substring matched in input .df filenames (default: sel_all).",
    )
    p.add_argument(
        "--phase",
        choices=("all", "meta", "write"),
        default="all",
        help="Run meta scan only, write only, or both (default: all).",
    )
    p.add_argument("--common-keys-pkl", default=None)
    p.add_argument("--summary-csv", default=None)
    p.add_argument("--matched-suffix", default="_matched")
    p.add_argument("--max-files", type=int, default=None)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.use_default_variations and not args.variation:
        raise SystemExit("pass --use-default-variations or --variation NAME DIR (twice)")

    if args.format == "sel_all":
        return run_sel_all_match(args)
    return run_sel_mup_match(args)


if __name__ == "__main__":
    raise SystemExit(main())
