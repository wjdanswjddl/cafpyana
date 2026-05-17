#!/usr/bin/env python3
"""
Merge grid-job ``.df`` (HDF5) outputs into ``<df_dir>/merged/``, following
``notebooks/combine_dfs.ipynb``:

- Concatenate per-job files with dense ``__ntuple`` remapping (same as
  ``pyanalib.split_df_helpers_new.dfs_from_dir`` / ``selected_events_cumulative``).
- Discover HDF keys from **one readable** input (grid jobs share the same schema);
  merge **all** keys found there.
- **Plan chunks from disk sizes**: read each input's ``os.path.getsize`` up front,
  then greedily group files so the sum of input sizes per group stays under a
  budget derived from ``--max-chunk-mb`` and ``--size-budget-fraction`` (merged
  output is only approximate; the fraction leaves headroom so the written HDF
  usually stays under the limit without probing every step).
- **Merge one planned chunk at a time**: load only that chunk's files, concat,
  add derived columns, write (with optional row-split if a chunk still exceeds
  ``--max-chunk-mb``), then drop the merged dict and ``gc.collect()`` — never
  retain the full sample in memory. Row-splitting advances ``iloc`` only on a
  **reference key** (``evt`` if present, else the longest table); other HDF keys
  may have different row counts and are written **in full** in every part.
- **Merge log** (JSONL): plan + each output records which input paths were merged.

Use ``--max-files`` to cap inputs for tests.
"""

from __future__ import annotations

import argparse
import gc
import glob
import json
import os
import re
import sys
import warnings
from os import path

import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
from tqdm import tqdm

warnings.filterwarnings("ignore", category=PerformanceWarning)

# Repository root (…/cafpyana)
_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.categories import get_genie_category, get_topo_category 
from pyanalib.pandas_helpers import multicol_add, pad_column_name 
from pyanalib.split_df_helpers_new import load_dfs 
from pyanalib.variable_calculator import get_cc1p0pi_tki 
from makedf.util import InFV 

_TKI_VAR_NAMES = (
    "del_alpha",
    "del_phi",
    "del_Tp",
    "del_p",
    "del_Tp_x",
    "del_Tp_y",
)


def _prefix_mcnu_columns(mc_nu_df: pd.DataFrame) -> None:
    """Leading ``mc`` level on ``mcnu`` columns (same as ``event_selection_chunk``)."""
    if not isinstance(mc_nu_df.columns, pd.MultiIndex):
        return
    try:
        first_level = mc_nu_df.columns.get_level_values(0)
        need_prefix = not np.all(first_level == "mc")
    except Exception:
        need_prefix = True
    if need_prefix:
        mc_nu_df.columns = pd.MultiIndex.from_tuples(
            [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
        )

CATHODE_INSET = 10
def perTPC_cut(df):
    in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=CATHODE_INSET) & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=CATHODE_INSET) & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=CATHODE_INSET)
    in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=CATHODE_INSET) & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=CATHODE_INSET) & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=CATHODE_INSET)
    perTPC_cut = in_TPC1_cut | in_TPC2_cut
    return perTPC_cut

def contTPC1_cut(df):
    cont_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=CATHODE_INSET) & InFV(df.mu.pfp.trk.end, det="SBND_TPC1", incathode=CATHODE_INSET) & InFV(df.p.pfp.trk.end, det="SBND_TPC1", incathode=CATHODE_INSET)
    return cont_TPC1_cut

def contTPC2_cut(df):
    cont_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=CATHODE_INSET) & InFV(df.mu.pfp.trk.end, det="SBND_TPC2", incathode=CATHODE_INSET) & InFV(df.p.pfp.trk.end, det="SBND_TPC2", incathode=CATHODE_INSET)
    return cont_TPC2_cut


def discover_dataframe_keys(hdf_path: str) -> list[str]:
    """
    List logical HDF keys (e.g. ``hdr``, ``evt``) from a ``.df`` file, excluding ``split``.

    Only includes bases that use the usual split layout (``<base>_0`` present), so
    ``load_dfs`` can read them.

    Uses ``h5py`` for a lightweight root key listing when available (faster on NFS
    than ``pd.HDFStore``); falls back to pandas if needed.
    """
    pat = re.compile(r"^(.+)_(\d+)$")
    raw_set: set[str] = set()
    try:
        import h5py

        with h5py.File(hdf_path, "r") as h5:
            raw_set = {str(k) for k in h5.keys()}
    except Exception:
        raw_set = set()

    def _bases_from_raw(names: set[str]) -> list[str]:
        bases_local: set[str] = set()
        for name in names:
            if name == "split":
                continue
            m = pat.match(name)
            if m:
                bases_local.add(m.group(1))
        return sorted(b for b in bases_local if f"{b}_0" in names)

    out = _bases_from_raw(raw_set)
    if out:
        return out

    with pd.HDFStore(hdf_path, mode="r") as store:
        raw_set = {str(k).strip("/") for k in store.keys()}
    return _bases_from_raw(raw_set)


def discover_df_keys_from_one_good_file(files: list[str]) -> tuple[list[str] | None, str | None, int]:
    """
    Read HDF dataset keys from the first input that opens cleanly.

    Assumes all grid outputs share the same keys (typical for one workflow).
    Tries ``files`` in sorted order until ``discover_dataframe_keys`` succeeds
    and returns a non-empty list.

    Returns ``(keys, probe_path, 0)`` on success, or ``(None, None, 1)`` on failure.
    """
    if not files:
        return None, None, 1
    for fp in files:
        try:
            keys = discover_dataframe_keys(fp)
        except Exception as e:
            print(f"[merge] WARNING: key probe skipped corrupt/unreadable file {fp}: {e}", flush=True)
            continue
        if not keys:
            print(f"[merge] WARNING: key probe found no logical keys in {fp}, trying next.", flush=True)
            continue
        print(
            f"[merge] HDF keys from {path.basename(fp)} "
            f"(same schema assumed for all {len(files)} input(s)): {keys}",
            flush=True,
        )
        return keys, fp, 0

    print("[merge] ERROR: could not read HDF keys from any input file.", flush=True)
    return None, None, 1


def _unique_ntuple_values_across_keys(
    mc_dfs: dict[str, pd.DataFrame], keys2load: list[str]
) -> np.ndarray:
    """All distinct ``__ntuple`` (or fallback index level) values across keys in one file."""
    parts: list[np.ndarray] = []
    for k in keys2load:
        df = mc_dfs[k]
        if df is None or len(df) == 0:
            continue
        if isinstance(df.index, pd.MultiIndex):
            names = list(df.index.names) if df.index.names is not None else []
            idx_loc = names.index("__ntuple") if "__ntuple" in names else 0
            v = df.index.get_level_values(idx_loc)
        elif df.index.name == "__ntuple":
            v = df.index
        else:
            v = df.index
        parts.append(np.asarray(v).reshape(-1))
    if not parts:
        return np.array([], dtype=np.int64)
    return np.sort(np.unique(np.concatenate(parts)))


def _remap_ntuple_indices(mc_dfs: dict[str, pd.DataFrame], keys2load: list[str], ntuple_offset: np.int64):
    """
    Dense ``__ntuple`` remapping in-place.

    Builds the remap from the **union** of ``__ntuple`` values over all keys in
    ``mc_dfs``, so a key that lists extra ntuples (e.g. ``evt`` vs ``hdr``) does
    not hit ``KeyError`` when the first key omits them.
    """
    unique_ntuples = _unique_ntuple_values_across_keys(mc_dfs, keys2load)
    ntuple_remap = {old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)}
    n_unique = np.int64(len(unique_ntuples))

    for df_key in keys2load:
        df = mc_dfs[df_key]
        if len(df) == 0:
            continue
        if isinstance(df.index, pd.MultiIndex):
            names = list(df.index.names) if df.index.names is not None else []
            idx_loc = names.index("__ntuple") if "__ntuple" in names else 0
            new_tuples = []
            for tup in df.index:
                tup = list(tup)
                tup[idx_loc] = ntuple_remap[tup[idx_loc]]
                new_tuples.append(tuple(tup))
            df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
        else:
            if df.index.name == "__ntuple":
                df.index = df.index.map(ntuple_remap)

    return n_unique


def gather_input_sizes(files: list[str]) -> list[tuple[str, int]]:
    """``(path, size_bytes)`` in the same order as ``files``; skip missing paths."""
    out: list[tuple[str, int]] = []
    for fp in files:
        if not path.isfile(fp):
            print(f"[merge] WARNING: skip missing file: {fp}", flush=True)
            continue
        out.append((fp, int(os.path.getsize(fp))))
    return out


def plan_file_groups(
    sized_files: list[tuple[str, int]],
    budget_bytes: int,
) -> list[list[str]]:
    """
    Greedily pack ``sized_files`` in order into groups whose summed input sizes
    stay ``<= budget_bytes``. Each oversized file becomes its own group.
    """
    if not sized_files:
        return []
    groups: list[list[str]] = []
    cur: list[str] = []
    cur_sum = 0
    for fp, sz in sized_files:
        if not cur:
            cur = [fp]
            cur_sum = sz
            continue
        if cur_sum + sz <= budget_bytes:
            cur.append(fp)
            cur_sum += sz
        else:
            groups.append(cur)
            cur = [fp]
            cur_sum = sz
    if cur:
        groups.append(cur)
    return groups


def concat_files_for_chunk(
    files: list[str],
    keys: list[str],
    n_max_splits: int,
    ntuple_offset: np.int64,
) -> tuple[dict[str, pd.DataFrame], np.int64]:
    """
    Load and vertically concat only ``files`` into one merged dict per key,
    applying dense ``__ntuple`` remapping with starting offset ``ntuple_offset``.
    Returns ``(merged, next_ntuple_offset)``.
    """
    df_lists: dict[str, list[pd.DataFrame]] = {k: [] for k in keys}
    off = np.int64(ntuple_offset)

    for mc_file in tqdm(files, desc="merge chunk inputs"):
        try:
            mc_dfs = load_dfs(mc_file, keys, n_max_concat=n_max_splits)
        except Exception as e:
            print(f"[merge] WARNING: failed to load {mc_file}: {e}", flush=True)
            continue

        n_unique = _remap_ntuple_indices(mc_dfs, keys, off)
        off += n_unique
        for k in keys:
            df_lists[k].append(mc_dfs[k])

    merged = {
        k: pd.concat(df_lists[k], axis=0, sort=False) if df_lists[k] else pd.DataFrame()
        for k in keys
    }
    return merged, off


def _row_split_ref_key(keys: list[str], merged: dict[str, pd.DataFrame]) -> str | None:
    """
    HDF row-splitting advances ``iloc`` on this key only.

    Prefer ``evt`` when present; otherwise the non-empty key with the most rows.
    Other keys may have different lengths and are always written in full in each part.
    """
    if "evt" in merged and len(merged["evt"]) > 0:
        return "evt"
    nonempty = [(k, len(merged[k])) for k in keys if k in merged and len(merged[k]) > 0]
    if not nonempty:
        return None
    return max(nonempty, key=lambda x: x[1])[0]


def _build_row_split_slice(
    merged: dict[str, pd.DataFrame],
    keys: list[str],
    ref_key: str,
    i0: int,
    i1: int,
) -> dict[str, pd.DataFrame]:
    """Slice ``ref_key`` with ``iloc[i0:i1]``; every other non-empty key is included in full."""
    out: dict[str, pd.DataFrame] = {}
    for k in keys:
        if k not in merged:
            continue
        df = merged[k]
        if len(df) == 0:
            continue
        if k == ref_key:
            out[k] = df.iloc[i0:i1]
        else:
            out[k] = df
    return out


def iter_row_slices_under_size(
    merged: dict[str, pd.DataFrame],
    keys: list[str],
    max_bytes: int,
    tmp_dir: str,
) -> tuple[str | None, list[tuple[int, int]]]:
    """
    Build ``(i0, i1)`` ranges on the **reference key** only so each written HDF stays
    under ``max_bytes``. Non-reference keys are included **in full** in every part
    (they may have different row counts than ``ref_key``).
    """
    ref_key = _row_split_ref_key(keys, merged)
    if ref_key is None or len(merged[ref_key]) == 0:
        return None, []

    n = len(merged[ref_key])
    tmp_measure = path.join(tmp_dir, ".merge_measure_tmp.df")
    ranges: list[tuple[int, int]] = []
    i0 = 0
    while i0 < n:
        lo, hi = i0 + 1, n
        best = i0
        while lo <= hi:
            mid = (lo + hi) // 2
            part = _build_row_split_slice(merged, keys, ref_key, i0, mid)
            try:
                sz = measure_hdf_bytes(part, tmp_measure)
            except Exception as e:
                print(f"[merge] WARNING: size probe failed at {ref_key} rows [{i0}:{mid}]: {e}", flush=True)
                sz = max_bytes + 1
            if sz <= max_bytes:
                best = mid
                lo = mid + 1
            else:
                hi = mid - 1

        if best <= i0:
            part = _build_row_split_slice(merged, keys, ref_key, i0, i0 + 1)
            try:
                sz_one = measure_hdf_bytes(part, tmp_measure)
            except Exception:
                sz_one = max_bytes + 1
            if sz_one > max_bytes:
                print(
                    f"[merge] WARNING: a single-row slice of {ref_key!r} (plus full other keys) "
                    f"exceeds max_bytes={max_bytes}; writing it anyway.",
                    flush=True,
                )
            best = i0 + 1

        ranges.append((i0, best))
        i0 = best

    if path.exists(tmp_measure):
        try:
            os.remove(tmp_measure)
        except OSError:
            pass
    return ref_key, ranges


def add_notebook_derived_columns(merged: dict[str, pd.DataFrame], fv_cut: str) -> None:
    """In-place: same additions as ``combine_dfs.ipynb`` cell 2 (best-effort)."""
    if "mcnu" in merged and len(merged["mcnu"]) > 0:
        mc_nu_df = merged["mcnu"]
        try:
            _prefix_mcnu_columns(mc_nu_df)
            mc_nu_df.loc[:, "topo_categ"] = get_topo_category(mc_nu_df)
            mc_nu_df.loc[:, "genie_categ"] = get_genie_category(mc_nu_df)
        except Exception as e:
            print(f"[merge] WARNING: mcnu topo/genie columns skipped: {e}", flush=True)

        try:
            mc_mudf = mc_nu_df["mc"]["mu"]
            mc_pdf = mc_nu_df["mc"]["p"]
            mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
            mc_P_p_col = pad_column_name(("totp",), mc_pdf)
            tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
            for var_name in _TKI_VAR_NAMES:
                merged["mcnu"] = multicol_add(merged["mcnu"], tki_mc[var_name].rename(f"{var_name}"))
        except Exception as e:
            print(f"[merge] WARNING: mcnu TKI variables skipped: {e}", flush=True)

    if "evt" in merged and len(merged["evt"]) > 0:
        mc_evt_df = merged["evt"]

        if fv_cut == "perTPC":
            mc_evt_df = mc_evt_df[perTPC_cut(mc_evt_df)]
        elif fv_cut == "contTPC1":
            mc_evt_df = mc_evt_df[contTPC1_cut(mc_evt_df)]
        elif fv_cut == "contTPC2":
            mc_evt_df = mc_evt_df[contTPC2_cut(mc_evt_df)]
        elif fv_cut == "nominal":
            pass
        else:
            print(f"[merge] WARNING: invalid FV cut: {fv_cut}", flush=True)
            return

        merged["evt"] = mc_evt_df

        try:
            mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
            mc_evt_df.loc[:, "genie_categ"] = get_genie_category(mc_evt_df)
        except Exception as e:
            print(f"[merge] WARNING: evt topo/genie columns skipped: {e}", flush=True)

        try:
            slc_mudf = mc_evt_df["mu"]["pfp"]["trk"]["truth"]["p"]
            slc_pdf = mc_evt_df["p"]["pfp"]["trk"]["truth"]["p"]
            slc_P_mu_col = pad_column_name(("totp",), slc_mudf)
            slc_P_p_col = pad_column_name(("totp",), slc_pdf)
            tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
            for var_name in _TKI_VAR_NAMES:
                merged["evt"] = multicol_add(merged["evt"], tki_reco[var_name].rename("mc_" + var_name))
        except Exception as e:
            print(f"[merge] WARNING: evt TKI variables skipped: {e}", flush=True)


def measure_hdf_bytes(slice_dfs: dict[str, pd.DataFrame], tmp_path: str) -> int:
    """Write ``slice_dfs`` to ``tmp_path`` (overwrite) and return file size in bytes."""
    if path.exists(tmp_path):
        os.remove(tmp_path)
    with pd.HDFStore(tmp_path, mode="w") as store:
        store.put("split", pd.DataFrame({"n_split": [1]}), format="fixed")
        for k, df in slice_dfs.items():
            if df is None or len(df) == 0:
                continue
            store.put(f"{k}_0", df, format="fixed")
    return os.path.getsize(tmp_path)


def write_merged_chunk(out_path: str, slice_dfs: dict[str, pd.DataFrame]) -> int:
    """Write a single merged ``.df`` (``split`` + ``<key>_0``). Returns size in bytes."""
    if path.exists(out_path):
        os.remove(out_path)
    with pd.HDFStore(out_path, mode="w") as store:
        store.put("split", pd.DataFrame({"n_split": [1]}), format="fixed")
        for k, df in slice_dfs.items():
            if df is None or len(df) == 0:
                continue
            store.put(f"{k}_0", df, format="fixed")
    return os.path.getsize(out_path)


def _merged_nonempty(merged: dict[str, pd.DataFrame], keys: list[str]) -> bool:
    return any(len(merged[k]) > 0 for k in keys)


def flush_merged_to_outputs(
    merged: dict[str, pd.DataFrame],
    input_files: list[str],
    *,
    planned_bytes: int,
    plan_chunk_index: int,
    keys: list[str],
    out_dir: str,
    stem: str,
    max_bytes: int,
    skip_derived: bool,
    log_fh,
    chunk_idx: int,
    fv_cut: str,
) -> int:
    """
    Apply optional derived columns, split by ``max_bytes`` if needed, write HDF(s),
    append JSONL log lines. Returns next ``chunk_idx``.
    """
    if not _merged_nonempty(merged, keys):
        return chunk_idx

    if not skip_derived:
        add_notebook_derived_columns(merged, fv_cut)

    ref_key, ranges = iter_row_slices_under_size(merged, keys, max_bytes=max_bytes, tmp_dir=out_dir)
    if not ranges or ref_key is None:
        return chunk_idx

    for r0, r1 in ranges:
        slice_dfs = _build_row_split_slice(merged, keys, ref_key, r0, r1)
        n_rows_ref = int(r1 - r0)
        row_iloc_log: list[str | int] = [ref_key, int(r0), int(r1)]

        out_path = path.join(out_dir, f"{stem}_merged_{chunk_idx:04d}.df")
        sz = write_merged_chunk(out_path, slice_dfs)
        entry = {
            "chunk_index": chunk_idx,
            "output_path": path.abspath(out_path),
            "plan_chunk_index": int(plan_chunk_index),
            "input_files": list(input_files),
            "planned_input_bytes": int(planned_bytes),
            "row_split_ref_key": ref_key,
            "buffer_merged_row_iloc": row_iloc_log,
            "n_rows_ref_key": n_rows_ref,
            "size_bytes": int(sz),
        }
        log_fh.write(json.dumps(entry, sort_keys=True) + "\n")
        log_fh.flush()
        print(
            f"[merge] wrote {out_path}  plan_chunk={plan_chunk_index}  "
            f"inputs={len(input_files)}  {ref_key} rows=[{r0}:{r1})  "
            f"(other keys written in full)  size={sz / (1024 ** 2):.2f} MiB",
            flush=True,
        )
        chunk_idx += 1

    return chunk_idx


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df-dir", required=True, help="Grid output directory containing ``*.df`` files.")
    p.add_argument(
        "--filename-str",
        required=True,
        help="Substring for ``glob('*<filename-str>*.df')`` (same as ``dfs_from_dir``).",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="If set, only merge the first N files after sorting (default: all matches).",
    )
    p.add_argument(
        "--n-max-splits",
        type=int,
        default=999,
        help="Max internal HDF splits read per input file (passed to ``load_dfs``).",
    )
    p.add_argument(
        "--max-chunk-mb",
        type=float,
        default=1024.0,
        help="Target upper bound for each written ``.df`` (default: 1024).",
    )
    p.add_argument(
        "--size-budget-fraction",
        type=float,
        default=0.92,
        help=(
            "When packing input files by ``os.path.getsize``, the sum of sizes per "
            "group must stay <= ``max_chunk_mb * fraction`` (default 0.92). Lower if "
            "merged HDFs often exceed the sum of inputs (e.g. many new columns)."
        ),
    )
    p.add_argument(
        "--fv",
        type=str,
        default="perTPC",
        help="FV cut to apply to the merged dataframe (default: perTPC).",
    )
    p.add_argument(
        "--output-stem",
        default=None,
        help="Stem for output files ``<stem>_merged_<NNNN>.df`` (default: basename of --df-dir).",
    )
    p.add_argument(
        "--merge-log",
        default=None,
        help="JSONL log path (default: ``<df-dir>/merged/merge_log.jsonl``).",
    )
    p.add_argument(
        "--skip-derived",
        action="store_true",
        help="Skip topo/genie/TKI columns from ``combine_dfs.ipynb``.",
    )
    args = p.parse_args()

    df_dir = path.abspath(args.df_dir)
    pattern = path.join(df_dir, f"*{args.filename_str}*.df")
    files = sorted(glob.glob(pattern))
    if args.max_files is not None:
        files = files[: int(args.max_files)]

    if not files:
        print(f"[merge] ERROR: no files matched: {pattern}", flush=True)
        return 1

    print(f"[merge] matched {len(files)} file(s) via {pattern}", flush=True)

    sized = gather_input_sizes(files)
    if not sized:
        print("[merge] ERROR: no readable input files on disk.", flush=True)
        return 1
    if len(sized) != len(files):
        print(f"[merge] WARNING: using {len(sized)} on-disk file(s) (skipped missing).", flush=True)

    sized_paths = [fp for fp, _ in sized]
    print(
        "[merge] Discovering HDF keys from a single readable input "
        "(all jobs assumed to share the same schema)…",
        flush=True,
    )
    keys, keys_probe_file, dcode = discover_df_keys_from_one_good_file(sized_paths)
    if dcode != 0 or keys is None:
        return 1

    max_bytes = int(float(args.max_chunk_mb) * 1024 * 1024)
    frac = float(args.size_budget_fraction)
    if not (0.0 < frac <= 1.0):
        print("[merge] ERROR: --size-budget-fraction must be in (0, 1].", flush=True)
        return 1
    budget_bytes = int(max_bytes * frac)

    file_groups = plan_file_groups(sized, budget_bytes)
    plan_summary = []
    for i, g in enumerate(file_groups):
        pb = sum(os.path.getsize(f) for f in g)
        plan_summary.append({"plan_chunk_index": i, "n_inputs": len(g), "planned_input_bytes": pb, "input_files": g})

    print(
        f"[merge] planned {len(file_groups)} merge chunk(s) from input sizes "
        f"(budget={budget_bytes / (1024 ** 2):.2f} MiB ≈ {args.max_chunk_mb} MiB × {frac})",
        flush=True,
    )
    for row in plan_summary:
        print(
            f"  plan_chunk {row['plan_chunk_index']}: {row['n_inputs']} file(s), "
            f"{row['planned_input_bytes'] / (1024 ** 2):.2f} MiB on disk",
            flush=True,
        )

    out_dir = path.join(df_dir, "merged" + "_" + args.fv)
    os.makedirs(out_dir, exist_ok=True)

    stem = args.output_stem or path.basename(df_dir.rstrip(path.sep))
    n_planned_outputs = len(file_groups)
    if n_planned_outputs > 0:
        print(
            f"[merge] Planned output file(s): {n_planned_outputs} "
            f"({stem}_merged_0000.df … {stem}_merged_{n_planned_outputs - 1:04d}.df if no row-splits). "
            "If a merged chunk still exceeds --max-chunk-mb on disk, it is split into extra files.",
            flush=True,
        )
    else:
        print("[merge] Planned output file(s): 0 (nothing to merge).", flush=True)

    log_path = args.merge_log or path.join(out_dir, "merge_log.jsonl")
    tmp_measure = path.join(out_dir, ".merge_measure_tmp.df")

    ntuple_offset = np.int64(0)
    chunk_idx = 0

    with open(log_path, "w", encoding="utf-8") as log_fh:
        meta = {
            "record_type": "run_header",
            "df_dir": df_dir,
            "filename_glob": pattern,
            "n_input_files": len(sized_paths),
            "keys": keys,
            "keys_discovery": "single_good_file",
            "keys_probe_file": path.abspath(keys_probe_file) if keys_probe_file else None,
            "max_chunk_mb": float(args.max_chunk_mb),
            "max_chunk_bytes": int(max_bytes),
            "size_budget_fraction": frac,
            "size_budget_bytes": int(budget_bytes),
            "skip_derived": bool(args.skip_derived),
            "planned_merge_chunks": plan_summary,
            "planned_n_output_files": int(n_planned_outputs),
            "planned_n_output_files_note": (
                "One HDF per planned size chunk unless row-splitting (on the ref key, "
                "evt if present) adds files when merged size exceeds --max-chunk-mb."
            ),
        }
        log_fh.write(json.dumps(meta, sort_keys=True) + "\n")

        for plan_i, group in enumerate(tqdm(file_groups, desc="merge plan chunks")):
            planned_bytes = sum(os.path.getsize(f) for f in group)
            # try:
            merged, ntuple_offset = concat_files_for_chunk(
                group,
                keys,
                int(args.n_max_splits),
                ntuple_offset,
            )
            # except Exception as e:
            #     print(f"[merge] ERROR: failed merging plan_chunk {plan_i} ({group}): {e}", flush=True)
            #     return 1

            chunk_idx = flush_merged_to_outputs(
                merged,
                group,
                planned_bytes=planned_bytes,
                plan_chunk_index=plan_i,
                keys=keys,
                out_dir=out_dir,
                stem=stem,
                max_bytes=max_bytes,
                skip_derived=bool(args.skip_derived),
                log_fh=log_fh,
                chunk_idx=chunk_idx,
                fv_cut=args.fv,
            )

            del merged
            gc.collect()

        if path.exists(tmp_measure):
            try:
                os.remove(tmp_measure)
            except OSError:
                pass

        summary = {
            "record_type": "run_footer",
            "n_output_chunks": int(chunk_idx),
            "merge_log_path": path.abspath(log_path),
        }
        log_fh.write(json.dumps(summary, sort_keys=True) + "\n")

    print(f"[merge] log written to {path.abspath(log_path)}", flush=True)
    print("[merge] done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
