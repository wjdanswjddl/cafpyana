"""
Helper functions for working with the HDF5 df files.  

Each dataset (identified by a *key*) is stored as a set of smaller dfs (one per split),
so that very large samples can be handled in chunks.  
"""

from __future__ import annotations

import glob
from os import path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm


def get_n_split(file):
    this_split_df = pd.read_hdf(file, key="split")
    return int(this_split_df.n_split.iloc[0])


def print_keys(file):
    with pd.HDFStore(file, mode="r") as store:
        keys = store.keys()  # list of all keys in the file
        # print("Keys:", keys)


def _coerce_hdf_frame_for_concat(df: pd.DataFrame) -> pd.DataFrame:
    """Cast object columns that store only bools — pandas 2.2+ concat can raise otherwise."""
    out = df.copy()
    for col in out.columns:
        s = out[col]
        if s.dtype != object:
            continue
        non_na = s.dropna()
        if non_na.empty:
            continue
        sample = non_na.iloc[: min(len(non_na), 4096)]
        if not sample.map(lambda x: isinstance(x, (bool, np.bool_))).all():
            continue
        try:
            out[col] = s.astype(bool)
        except (TypeError, ValueError):
            pass
    return out


def _ntuple_level(df: pd.DataFrame) -> Optional[int]:
    """Index of the ``__ntuple`` level, or 0 for unnamed MultiIndex; ``None`` if not MultiIndex."""
    if not isinstance(df.index, pd.MultiIndex):
        return None
    names = df.index.names
    if names and "__ntuple" in names:
        return names.index("__ntuple")
    return 0


def _ntuple_values(df: pd.DataFrame) -> np.ndarray:
    """Per-row ``__ntuple`` scalars (for union/remap), regardless of index representation."""
    if isinstance(df.index, pd.MultiIndex):
        lvl = _ntuple_level(df)
        return np.asarray(df.index.get_level_values(lvl)).reshape(-1)
    if df.index.name == "__ntuple":
        return np.asarray(df.index).reshape(-1)
    if len(df) and isinstance(df.index[0], tuple):
        return np.fromiter((t[0] for t in df.index), dtype=np.int64, count=len(df))
    return np.asarray(df.index).reshape(-1)


def _unique_ntuple_values_across_keys(
    mc_dfs: dict[str, pd.DataFrame], keys2load: Sequence[str]
) -> np.ndarray:
    """Distinct ``__ntuple`` ids across all loaded keys in one HDF file."""
    parts: List[np.ndarray] = []
    for k in keys2load:
        df = mc_dfs.get(k)
        if df is None or len(df) == 0:
            continue
        parts.append(_ntuple_values(df))
    if not parts:
        return np.array([], dtype=np.int64)
    return np.sort(np.unique(np.concatenate(parts)))


def _remap_ntuple_index(df: pd.DataFrame, ntuple_remap: dict) -> None:
    """Rewrite ``__ntuple`` in-place (MultiIndex, named index, or tuple plain Index)."""
    if isinstance(df.index, pd.MultiIndex):
        names = list(df.index.names) if df.index.names is not None else []
        idx_loc = _ntuple_level(df)
        new_tuples = []
        for tup in df.index:
            tup = list(tup)
            old = tup[idx_loc]
            if old not in ntuple_remap:
                raise KeyError(
                    f"ntuple {old!r} missing from remap; known ntuples: {sorted(ntuple_remap)}"
                )
            tup[idx_loc] = ntuple_remap[old]
            new_tuples.append(tuple(tup))
        df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
    elif df.index.name == "__ntuple":
        df.index = df.index.map(ntuple_remap)
    elif len(df) and isinstance(df.index[0], tuple):
        new_tuples = []
        for tup in df.index:
            tup = list(tup)
            old = tup[0]
            if old not in ntuple_remap:
                raise KeyError(
                    f"ntuple {old!r} missing from remap; known ntuples: {sorted(ntuple_remap)}"
                )
            tup[0] = ntuple_remap[old]
            new_tuples.append(tuple(tup))
        df.index = pd.Index(new_tuples)


def _concat_hdf_frames(frames: Sequence[pd.DataFrame], *, label: str = "") -> pd.DataFrame:
    if not frames:
        suffix = f" ({label})" if label else ""
        raise ValueError(f"No dataframes to concatenate{suffix}")
    fixed = [_coerce_hdf_frame_for_concat(f) for f in frames]
    if len(fixed) == 1:
        return fixed[0]
    return pd.concat(fixed, axis=0, sort=False)


def load_dfs(file, keys2load, n_max_concat=100):
    out_df_dict = {}
    this_n_keys = get_n_split(file)
    n_concat = min(int(n_max_concat), int(this_n_keys))
    for key in keys2load:
        dfs: List[pd.DataFrame] = []
        for i in range(n_concat):
            dfs.append(pd.read_hdf(file, key=f"{key}_{i}"))
        out_df_dict[key] = _concat_hdf_frames(dfs, label=f"{file}:{key}")

    return out_df_dict


def dfs_from_dir(
    search_dir,
    filename_str=None,
    keys2load=("hdr", "evt"),
    n_max_concat=3,
    n_max_splits_per_file: Optional[int] = None,
):
    """
    Loop over ``*<filename_str>*.df`` files under ``search_dir``.

    * ``n_max_concat`` — maximum number of **files** to read.
    * ``n_max_splits_per_file`` — maximum HDF **splits** per file (``evt_0``, …).
      Default ``None`` loads every split in each file (recommended for merged GENIE grids).
    """
    df_lists = {k: [] for k in keys2load}
    ntuple_offset = np.int64(0)

    if filename_str is None:
        raise ValueError("No filename string provided, aborting...")

    pattern = path.join(search_dir, f"*{filename_str}*.df")
    matched_files = glob.glob(pattern)
    files_to_process = sorted(matched_files)[: int(n_max_concat)]
    print(f"Found {len(files_to_process)} files to process")
    print(f"Files to process: {files_to_process}")

    if not files_to_process:
        raise FileNotFoundError(
            f"No files matching {pattern!r} under {search_dir!r}"
        )

    load_errors = []

    for mc_file in tqdm(files_to_process):
        mc_n_split = get_n_split(mc_file)
        splits_cap = int(mc_n_split) if n_max_splits_per_file is None else int(n_max_splits_per_file)
        try:
            mc_dfs = load_dfs(mc_file, keys2load, n_max_concat=splits_cap)
        except Exception as e:
            load_errors.append((mc_file, e))
            print(f"Error loading file {mc_file}: {e}")
            continue

        unique_ntuples = _unique_ntuple_values_across_keys(mc_dfs, keys2load)
        ntuple_remap = {
            old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)
        }
        n_unique = np.int64(len(ntuple_remap))

        for df_key in keys2load:
            df = mc_dfs[df_key]
            _remap_ntuple_index(df, ntuple_remap)
            df_lists[df_key].append(df)

        ntuple_offset += n_unique

    if load_errors and not any(df_lists[k] for k in keys2load):
        detail = "\n".join(f"  {fp}: {err}" for fp, err in load_errors[:5])
        raise RuntimeError(
            f"Failed to load any HDF splits from {search_dir!r} (pattern {pattern!r}).\n{detail}"
        )

    concat_dfs = {
        k: _concat_hdf_frames(df_lists[k], label=k) for k in keys2load if df_lists[k]
    }
    if set(concat_dfs.keys()) != set(keys2load):
        missing = [k for k in keys2load if k not in concat_dfs]
        raise RuntimeError(
            f"Missing keys after concat {missing!r} from {search_dir!r}; "
            f"{len(load_errors)} file(s) failed to load"
        )

    print("REMEMBER TO RECALCULATE TKI AND CHECK FV!!")
    return concat_dfs
