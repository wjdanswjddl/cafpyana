#!/usr/bin/env python3
"""Add-only enrichment of split HDF5 ``.df`` / ``.h5`` tables (evt_*, mcnu_*).

Derived from ``notebooks/update_dfs.ipynb``. Existing columns and keys are never
renamed, dropped, or overwritten. New columns are skipped if already present.

Safety model (write-aside + verify + publish to sibling ``*_updated``)
----------------------------------------------------------------------
In-place ``HDFStore`` edits (even on a copy) can leave a half-updated file if a
``put`` fails mid-way. Instead:

1. **Process** — read each original file (read-only), enrich in memory, write a
   *complete* new HDF5 under a staging tmp tree (same filesystem when possible).
2. **Verify** — reopen the tmp file; check keys, row counts, and that every
   original column is still present with matching values on a sample of columns.
3. **Commit** — move verified files into a parallel sibling directory
   ``<original_directory_name>_updated/`` (same basenames). Originals are never
   replaced. Failed files leave originals untouched and keep the tmp + a log
   entry for inspection.

Usage
-----
    python update_dfs_add_columns.py --dirs /path/A /path/B
    python update_dfs_add_columns.py --dirs /path/A --dry-run
    python update_dfs_add_columns.py --dirs /path/A --process-only
    python update_dfs_add_columns.py --dirs /path/A --commit-only
    python update_dfs_add_columns.py --dirs /path/A --log-file /path/run.log
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import traceback
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from os import path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
try:
    import tables

    warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)
except Exception:
    pass

# cafpyana root (…/cafpyana)
_CAFPYANA_ROOT = path.dirname(
    path.dirname(path.dirname(path.dirname(path.abspath(__file__))))
)
if _CAFPYANA_ROOT not in sys.path:
    sys.path.insert(0, _CAFPYANA_ROOT)

from analysis_village.numucc_1p0pi.categories import (  # noqa: E402
    get_genie_category,
    get_topo_category,
)
from pyanalib.pandas_helpers import multicol_add, pad_column_name  # noqa: E402
from pyanalib.variable_calculator import get_cc1p0pi_tki  # noqa: E402

LOG = logging.getLogger("update_dfs_add_columns")

TKI_VAR_NAMES = ("del_alpha", "del_phi", "del_Tp", "del_p", "del_Tp_x", "del_Tp_y")
_KEY_BASE_RE = re.compile(r"^/?(.+)_(\d+)$")
_DEFAULT_GLOBS = ("*.df", "*.h5")
# No leading dot: pnfs/dCache often rejects hidden (.* ) directory names with EPERM.
_TMP_DIRNAME = "update_dfs_tmp"
_UPDATED_SUFFIX = "_updated"
_MANIFEST_NAME = "manifest.json"
_VERIFY_SAMPLE_COLS = 8


def updated_dir_for(src_or_dir: str) -> str:
    """Sibling publish dir: ``/path/foo`` → ``/path/foo_updated``."""
    d = path.abspath(src_or_dir)
    if path.isfile(d):
        d = path.dirname(d)
    return d.rstrip(os.sep) + _UPDATED_SUFFIX


def dest_path_for(src: str) -> str:
    """Final path under ``<srcdir>_updated/`` with the same basename."""
    return path.join(updated_dir_for(src), path.basename(src))


# ---------------------------------------------------------------------------
# Column helpers (add-only)
# ---------------------------------------------------------------------------


def _level0_names(df: pd.DataFrame) -> set:
    if not isinstance(df.columns, pd.MultiIndex):
        return set(df.columns)
    return set(df.columns.get_level_values(0))


def _has_level0(df: pd.DataFrame, name: str) -> bool:
    return name in _level0_names(df)


def _has_multicolumn(df: pd.DataFrame, parts: Tuple[str, ...]) -> bool:
    """True if a MultiIndex column matching ``parts`` (trailing '' ok) exists."""
    if not isinstance(df.columns, pd.MultiIndex):
        return len(parts) == 1 and parts[0] in df.columns
    try:
        key = pad_column_name(parts, df)
        return key in df.columns
    except Exception:
        return False


def _assign_level0(df: pd.DataFrame, name: str, values) -> pd.DataFrame:
    """Assign a top-level column if missing (add-only)."""
    if _has_level0(df, name):
        return df
    df = df.copy()
    df.loc[:, name] = values
    return df


def _key_base(store_key: str) -> Optional[str]:
    m = _KEY_BASE_RE.match(str(store_key).strip("/"))
    return m.group(1) if m else None


def _list_store_keys(hdf_path: str) -> List[str]:
    with pd.HDFStore(hdf_path, mode="r") as store:
        return [str(k).strip("/") for k in store.keys()]


# ---------------------------------------------------------------------------
# Enrichment (mirrors update_dfs.ipynb; add-only)
# ---------------------------------------------------------------------------


def add_opening_angle(df: pd.DataFrame, *, truth: bool = False, nu: bool = False) -> pd.DataFrame:
    """Opening angle in **degrees** (same as ``update_dfs.ipynb``)."""
    if nu:
        if _has_level0(df, "theta_mu_p"):
            return df
        opening = (
            df.mc.mu.dir[["x", "y", "z"]].to_numpy() * df.mc.p.dir[["x", "y", "z"]].to_numpy()
        ).sum(axis=1)
        return _assign_level0(df, "theta_mu_p", np.arccos(np.clip(opening, -1.0, 1.0)) * 180.0 / np.pi)

    if truth:
        if _has_level0(df, "mc_theta_mu_p"):
            return df
        opening = (
            df.mc.mu.dir[["x", "y", "z"]].to_numpy() * df.mc.p.dir[["x", "y", "z"]].to_numpy()
        ).sum(axis=1)
        return _assign_level0(
            df, "mc_theta_mu_p", np.arccos(np.clip(opening, -1.0, 1.0)) * 180.0 / np.pi
        )

    if _has_level0(df, "theta_mu_p"):
        return df
    opening = (
        df.mu.pfp.trk.dir[["x", "y", "z"]].to_numpy() * df.p.pfp.trk.dir[["x", "y", "z"]].to_numpy()
    ).sum(axis=1)
    return _assign_level0(df, "theta_mu_p", np.arccos(np.clip(opening, -1.0, 1.0)) * 180.0 / np.pi)


def add_track_direction(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``(mu|p).pfp.trk.truth.p.dir.{x,y,z}`` from genp / totp when missing."""
    targets = [
        (("mu", "pfp", "trk", "truth", "p", "dir", "x"), ("mu", "pfp", "trk", "truth", "p", "genp", "x")),
        (("mu", "pfp", "trk", "truth", "p", "dir", "y"), ("mu", "pfp", "trk", "truth", "p", "genp", "y")),
        (("mu", "pfp", "trk", "truth", "p", "dir", "z"), ("mu", "pfp", "trk", "truth", "p", "genp", "z")),
        (("p", "pfp", "trk", "truth", "p", "dir", "x"), ("p", "pfp", "trk", "truth", "p", "genp", "x")),
        (("p", "pfp", "trk", "truth", "p", "dir", "y"), ("p", "pfp", "trk", "truth", "p", "genp", "y")),
        (("p", "pfp", "trk", "truth", "p", "dir", "z"), ("p", "pfp", "trk", "truth", "p", "genp", "z")),
    ]
    need = [t for t, _ in targets if not _has_multicolumn(df, t)]
    if not need:
        return df

    out = df.copy()
    mu_totp = out.mu.pfp.trk.truth.p.totp
    p_totp = out.p.pfp.trk.truth.p.totp
    for dest, src in targets:
        if _has_multicolumn(out, dest):
            continue
        axis = dest[-1]
        if dest[0] == "mu":
            out.loc[:, pad_column_name(dest, out)] = out.mu.pfp.trk.truth.p.genp[axis] / mu_totp
        else:
            out.loc[:, pad_column_name(dest, out)] = out.p.pfp.trk.truth.p.genp[axis] / p_totp
    return out


def _with_mc_prefix(df: pd.DataFrame) -> pd.DataFrame:
    """Temporary view with leading ``mc`` level for category helpers (not persisted)."""
    if not isinstance(df.columns, pd.MultiIndex):
        raise TypeError("mcnu enrichment expects MultiIndex columns")
    first = df.columns.get_level_values(0)
    if np.all(first == "mc"):
        return df
    out = df.copy()
    out.columns = pd.MultiIndex.from_tuples([tuple(["mc"] + list(c)) for c in out.columns])
    return out


def enrich_evt(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """Add evt columns from ``update_dfs.ipynb``. Returns (df, notes)."""
    notes: List[str] = []
    if df is None or len(df) == 0:
        return df, notes
    if not isinstance(df.columns, pd.MultiIndex):
        notes.append("skip evt: not MultiIndex")
        return df, notes

    if not _has_level0(df, "topo_categ"):
        df = _assign_level0(df, "topo_categ", get_topo_category(df))
        notes.append("added topo_categ")
    else:
        notes.append("skip topo_categ (exists)")

    if not _has_level0(df, "genie_categ"):
        df = _assign_level0(df, "genie_categ", get_genie_category(df))
        notes.append("added genie_categ")
    else:
        notes.append("skip genie_categ (exists)")

    before = set(df.columns)
    df = add_track_direction(df)
    if set(df.columns) - before:
        notes.append("added truth track dir")
    else:
        notes.append("skip truth track dir (exists)")

    before = set(df.columns)
    df = add_opening_angle(df)
    df = add_opening_angle(df, truth=True)
    if set(df.columns) - before:
        notes.append("added opening angles")
    else:
        notes.append("skip opening angles (exist)")

    # Truth-particle TKI on selected tracks → top-level ``mc_<var>`` (notebook naming)
    missing_tki = [v for v in TKI_VAR_NAMES if not _has_level0(df, "mc_" + v)]
    if missing_tki:
        slc_mudf = df.mu.pfp.trk.truth.p
        slc_pdf = df.p.pfp.trk.truth.p
        slc_P_mu_col = pad_column_name(("totp",), slc_mudf)
        slc_P_p_col = pad_column_name(("totp",), slc_pdf)
        tki = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
        for var_name in missing_tki:
            df = multicol_add(df, tki[var_name].rename("mc_" + var_name))
        notes.append(f"added mc_tki: {missing_tki}")
    else:
        notes.append("skip mc_tki (exist)")

    return df, notes


def enrich_mcnu(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """Add mcnu columns from ``update_dfs.ipynb`` without persisting an ``mc`` prefix.

    Category helpers need ``df.mc.*``; we compute on a temporary prefixed copy and
    only write the new scalar columns onto the original layout (add-only).
    Opening angle / TKI use top-level ``mu`` / ``p`` when present, else ``mc.mu`` /
    ``mc.p`` (already-prefixed files).
    """
    notes: List[str] = []
    if df is None or len(df) == 0:
        return df, notes
    if not isinstance(df.columns, pd.MultiIndex):
        notes.append("skip mcnu: not MultiIndex")
        return df, notes

    # Categories via temporary mc-prefixed view
    need_topo = not _has_level0(df, "topo_categ")
    need_genie = not _has_level0(df, "genie_categ")
    if need_topo or need_genie:
        view = _with_mc_prefix(df)
        if need_topo:
            df = _assign_level0(df, "topo_categ", get_topo_category(view))
            notes.append("added topo_categ")
        else:
            notes.append("skip topo_categ (exists)")
        if need_genie:
            df = _assign_level0(df, "genie_categ", get_genie_category(view))
            notes.append("added genie_categ")
        else:
            notes.append("skip genie_categ (exists)")
    else:
        notes.append("skip categories (exist)")

    # Opening angle: prefer already-prefixed layout (notebook), else temp prefix
    if not _has_level0(df, "theta_mu_p"):
        levels0 = _level0_names(df)
        if "mc" in levels0 and "mu" in set(df["mc"].columns.get_level_values(0)):
            df = add_opening_angle(df, nu=True)
        else:
            view = _with_mc_prefix(df)
            opening = (
                view.mc.mu.dir[["x", "y", "z"]].to_numpy()
                * view.mc.p.dir[["x", "y", "z"]].to_numpy()
            ).sum(axis=1)
            df = _assign_level0(
                df, "theta_mu_p", np.arccos(np.clip(opening, -1.0, 1.0)) * 180.0 / np.pi
            )
        notes.append("added theta_mu_p")
    else:
        notes.append("skip theta_mu_p (exists)")

    missing_tki = [v for v in TKI_VAR_NAMES if not _has_level0(df, v)]
    if missing_tki:
        levels0 = _level0_names(df)
        if "mc" in levels0:
            mc_blk = df["mc"]
            sub0 = set(mc_blk.columns.get_level_values(0))
            if "mu" in sub0 and "p" in sub0:
                mc_mudf, mc_pdf = mc_blk["mu"], mc_blk["p"]
            elif "mu" in levels0 and "p" in levels0:
                mc_mudf, mc_pdf = df["mu"], df["p"]
            else:
                raise KeyError("mcnu missing mu/p for TKI")
        elif "mu" in levels0 and "p" in levels0:
            mc_mudf, mc_pdf = df["mu"], df["p"]
        else:
            raise KeyError("mcnu missing mu/p for TKI")

        mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
        mc_P_p_col = pad_column_name(("totp",), mc_pdf)
        tki = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
        for var_name in missing_tki:
            df = multicol_add(df, tki[var_name].rename(var_name))
        notes.append(f"added tki: {missing_tki}")
    else:
        notes.append("skip tki (exist)")

    return df, notes


def enrich_by_key(key: str, df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    base = _key_base(key)
    if base == "evt":
        return enrich_evt(df)
    if base == "mcnu":
        return enrich_mcnu(df)
    return df, [f"passthrough {key}"]


# ---------------------------------------------------------------------------
# Verify / I/O
# ---------------------------------------------------------------------------


def _sample_original_columns(df: pd.DataFrame, n: int = _VERIFY_SAMPLE_COLS) -> list:
    cols = list(df.columns)
    if len(cols) <= n:
        return cols
    # Prefer stable early columns (not newly added tails)
    step = max(1, len(cols) // n)
    return [cols[i] for i in range(0, min(len(cols), step * n), step)][:n]


def verify_add_only(
    original: dict,
    written_path: str,
) -> None:
    """Raise if tmp file lost keys/rows/columns or changed sampled original values."""
    with pd.HDFStore(written_path, mode="r") as store:
        new_keys = {str(k).strip("/") for k in store.keys()}
    old_keys = set(original.keys())
    if old_keys - new_keys:
        raise RuntimeError(f"tmp missing keys: {sorted(old_keys - new_keys)}")
    if new_keys - old_keys:
        raise RuntimeError(f"tmp has unexpected keys: {sorted(new_keys - old_keys)}")

    for key, old_df in original.items():
        new_df = pd.read_hdf(written_path, key=key)
        if len(new_df) != len(old_df):
            raise RuntimeError(f"{key}: row count {len(old_df)} → {len(new_df)}")
        if list(new_df.index.names) != list(old_df.index.names):
            raise RuntimeError(f"{key}: index names changed")

        old_cols = set(old_df.columns)
        new_cols = set(new_df.columns)
        missing = old_cols - new_cols
        if missing:
            raise RuntimeError(f"{key}: lost columns (sample): {list(missing)[:5]}")

        for col in _sample_original_columns(old_df):
            if col not in new_df.columns:
                raise RuntimeError(f"{key}: missing sampled column {col}")
            a = old_df[col]
            b = new_df[col]
            if isinstance(a, pd.DataFrame) or isinstance(b, pd.DataFrame):
                continue
            if not a.equals(b):
                # Allow pure dtype widening that still compares equal numerically
                try:
                    if np.allclose(
                        np.asarray(a, dtype=float),
                        np.asarray(b, dtype=float),
                        equal_nan=True,
                    ):
                        continue
                except Exception:
                    pass
                raise RuntimeError(f"{key}: values changed for column {col}")


def write_hdf_dict(out_path: str, frames: dict) -> None:
    os.makedirs(path.dirname(out_path) or ".", exist_ok=True)
    part = out_path + ".partial"
    if path.exists(part):
        os.remove(part)
    try:
        with pd.HDFStore(part, mode="w") as store:
            for key, df in frames.items():
                store.put(key, df, format="fixed")
        os.replace(part, out_path)
    except Exception:
        if path.exists(part):
            try:
                os.remove(part)
            except OSError:
                pass
        raise


@dataclass
class FileResult:
    src: str
    tmp: str
    status: str  # ok | failed | skipped
    message: str = ""
    notes: List[str] = field(default_factory=list)
    dst: str = ""  # path under <dir>_updated/ after commit


@dataclass
class Manifest:
    created_utc: str
    dirs: List[str]
    tmp_root: Optional[str]
    results: List[FileResult] = field(default_factory=list)

    def save(self, path_out: str) -> None:
        payload = {
            "created_utc": self.created_utc,
            "dirs": self.dirs,
            "tmp_root": self.tmp_root,
            "results": [asdict(r) for r in self.results],
        }
        os.makedirs(path.dirname(path_out) or ".", exist_ok=True)
        # Write via partial + replace: avoids leaving a stuck/unreadable pnfs
        # object if a direct open(path, "w") fails mid-way (seen as EPERM).
        partial = path_out + ".partial"
        if path.exists(partial):
            try:
                os.remove(partial)
            except OSError:
                pass
        with open(partial, "w") as f:
            json.dump(payload, f, indent=2)
            f.flush()
            os.fsync(f.fileno())
        try:
            os.replace(partial, path_out)
        except OSError:
            # pnfs can leave an unreadable target; remove and retry once
            if path.exists(path_out):
                os.remove(path_out)
            os.replace(partial, path_out)


def _default_tmp_for(src: str, tmp_root: Optional[str]) -> str:
    src_dir = path.dirname(path.abspath(src))
    name = path.basename(src)
    if tmp_root:
        # Mirror absolute path under tmp_root to avoid collisions
        rel = path.abspath(src).lstrip(os.sep)
        return path.join(path.abspath(tmp_root), rel)
    return path.join(src_dir, _TMP_DIRNAME, name)


def discover_files(dirs: Sequence[str], patterns: Sequence[str]) -> List[str]:
    import glob as _glob

    out: List[str] = []
    for d in dirs:
        d = path.abspath(d)
        if not path.isdir(d):
            raise FileNotFoundError(f"not a directory: {d}")
        base = path.basename(d.rstrip(os.sep))
        if base.endswith(_UPDATED_SUFFIX) or base == _TMP_DIRNAME:
            LOG.warning("skipping input dir that looks like an output/tmp tree: %s", d)
            continue
        for pat in patterns:
            out.extend(sorted(_glob.glob(path.join(d, pat))))
    # de-dupe preserve order
    seen = set()
    uniq = []
    for f in out:
        if f not in seen and path.isfile(f):
            parent = path.basename(path.dirname(f))
            # skip files already inside staging or publish dirs
            if parent == _TMP_DIRNAME or parent.endswith(_UPDATED_SUFFIX):
                continue
            if f"/{_TMP_DIRNAME}/" in f:
                continue
            seen.add(f)
            uniq.append(f)
    return uniq


def process_one_file(src: str, tmp: str) -> FileResult:
    notes: List[str] = []
    dst = dest_path_for(src)
    try:
        keys = _list_store_keys(src)
        if not keys:
            return FileResult(src, tmp, "skipped", "no keys", notes, dst=dst)

        originals: dict = {}
        frames: dict = {}
        for key in keys:
            df = pd.read_hdf(src, key=key)
            originals[key] = df
            try:
                new_df, key_notes = enrich_by_key(key, df)
                notes.extend([f"{key}: {n}" for n in key_notes])
                frames[key] = new_df
            except Exception as exc:
                # Keep original table for this key; record and continue
                notes.append(
                    f"{key}: ENRICH_FAILED ({type(exc).__name__}: {exc}); kept original"
                )
                frames[key] = df

        write_hdf_dict(tmp, frames)
        verify_add_only(originals, tmp)
        return FileResult(src, tmp, "ok", "verified", notes, dst=dst)
    except Exception as exc:
        tb = traceback.format_exc(limit=8)
        LOG.error("FAILED %s\n%s", src, tb)
        return FileResult(src, tmp, "failed", f"{type(exc).__name__}: {exc}", notes, dst=dst)


def commit_one(result: FileResult) -> FileResult:
    """Publish verified tmp → ``<srcdir>_updated/<basename>`` (never replace src)."""
    if result.status != "ok":
        return result
    if not path.isfile(result.tmp):
        return FileResult(
            result.src,
            result.tmp,
            "failed",
            "tmp missing at commit",
            result.notes,
            dst=result.dst,
        )
    dst = result.dst or dest_path_for(result.src)
    os.makedirs(path.dirname(dst) or ".", exist_ok=True)
    # Atomic publish onto the same filesystem as tmp when possible
    os.replace(result.tmp, dst)
    # clean empty tmp parent if possible
    tmp_parent = path.dirname(result.tmp)
    try:
        if path.isdir(tmp_parent) and not os.listdir(tmp_parent):
            os.rmdir(tmp_parent)
    except OSError:
        pass
    return FileResult(
        result.src, result.tmp, "ok", "committed", result.notes, dst=dst
    )


def load_manifest(manifest_path: str) -> Manifest:
    with open(manifest_path) as f:
        raw = json.load(f)
    results = [FileResult(**r) for r in raw.get("results", [])]
    return Manifest(
        created_utc=raw["created_utc"],
        dirs=raw["dirs"],
        tmp_root=raw.get("tmp_root"),
        results=results,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Add-only column enrichment for split HDF5 dataframe files. "
            f"Publishes to sibling <dir>{_UPDATED_SUFFIX}/; originals are never replaced."
        )
    )
    p.add_argument(
        "--dirs",
        nargs="+",
        required=True,
        help="Subdirectories containing .df / .h5 files to process",
    )
    p.add_argument(
        "--glob",
        dest="globs",
        action="append",
        default=None,
        help="Filename glob (repeatable). Default: *.df and *.h5",
    )
    p.add_argument(
        "--tmp-root",
        default=None,
        help=(
            "Optional root for staging tmp outputs (mirrors absolute paths). "
            f"Default: <each file's dir>/{_TMP_DIRNAME}/"
        ),
    )
    p.add_argument(
        "--manifest",
        default=None,
        help=(
            "Manifest JSON path "
            f"(default: <first_dir>{_UPDATED_SUFFIX}/{_MANIFEST_NAME})"
        ),
    )
    p.add_argument(
        "--process-only",
        action="store_true",
        help="Write+verify staging tmp files only; do not publish to *_updated",
    )
    p.add_argument(
        "--commit-only",
        action="store_true",
        help="Publish verified tmp files from an existing manifest into *_updated",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="List files that would be processed and their *_updated destinations",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Process at most N files (0 = all)",
    )
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument(
        "--log-file",
        default=None,
        help=(
            "Write verbose (DEBUG) detail to this file. Console stays at INFO "
            "(progress bar only); per-file paths go to the log file."
        ),
    )
    return p.parse_args(argv)


def _configure_logging(*, verbose: bool, log_file: Optional[str]) -> None:
    """Console stays quiet (INFO); verbose detail goes to ``log_file`` when set.

    With ``--log-file``, the file is always DEBUG. ``-v`` without a log file
    raises the console to DEBUG; with a log file, ``-v`` does not spam the
    terminal (progress bars are enough there).
    """
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(logging.DEBUG)

    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    if log_file:
        console.setLevel(logging.INFO)
    else:
        console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(fmt)
    root.addHandler(console)

    if log_file:
        log_path = path.abspath(log_file)
        os.makedirs(path.dirname(log_path) or ".", exist_ok=True)
        fh = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        logging.info("verbose log file: %s", log_path)


def _manifest_path(args: argparse.Namespace) -> str:
    if args.manifest:
        return path.abspath(args.manifest)
    if args.tmp_root:
        return path.join(path.abspath(args.tmp_root), _MANIFEST_NAME)
    # Live with the published outputs
    return path.join(updated_dir_for(args.dirs[0]), _MANIFEST_NAME)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    _configure_logging(verbose=args.verbose, log_file=args.log_file)
    if args.process_only and args.commit_only:
        LOG.error("choose at most one of --process-only / --commit-only")
        return 2

    patterns = tuple(args.globs) if args.globs else _DEFAULT_GLOBS
    files = discover_files(args.dirs, patterns)
    if args.limit > 0:
        files = files[: args.limit]

    LOG.info("dirs=%s  files=%d  patterns=%s", args.dirs, len(files), patterns)
    for f in files:
        LOG.debug("  %s  →  %s", f, dest_path_for(f))
    if args.dry_run:
        LOG.info("dry-run: no writes (%d files listed in log at DEBUG)", len(files))
        return 0

    manifest_path = _manifest_path(args)

    if args.commit_only:
        if not path.isfile(manifest_path):
            LOG.error("manifest not found: %s", manifest_path)
            return 1
        man = load_manifest(manifest_path)
        n_ok = n_fail = 0
        updated: List[FileResult] = []
        for r in tqdm(man.results, desc="commit", unit="file"):
            if r.status != "ok" or r.message == "committed":
                if r.message == "committed":
                    updated.append(r)
                    n_ok += 1
                    continue
                updated.append(r)
                if r.status != "ok":
                    n_fail += 1
                    LOG.warning("skip commit (%s): %s — %s", r.status, r.src, r.message)
                continue
            out = commit_one(r)
            updated.append(out)
            if out.message == "committed":
                n_ok += 1
                LOG.debug("committed %s → %s", out.src, out.dst)
            else:
                n_fail += 1
                LOG.error("commit failed %s: %s", out.src, out.message)
        man.results = updated
        man.save(manifest_path)
        LOG.info(
            "commit done: %d published to *_updated, %d not published  manifest=%s",
            n_ok,
            n_fail,
            manifest_path,
        )
        return 0 if n_fail == 0 else 1

    # process (and optionally commit / publish)
    man = Manifest(
        created_utc=datetime.now(timezone.utc).isoformat(),
        dirs=[path.abspath(d) for d in args.dirs],
        tmp_root=path.abspath(args.tmp_root) if args.tmp_root else None,
    )

    n_fail_running = 0
    pbar = tqdm(files, desc="update_dfs", unit="file")
    for src in pbar:
        tmp = _default_tmp_for(src, args.tmp_root)
        LOG.debug("process %s → staging %s", src, tmp)
        result = process_one_file(src, tmp)
        if result.status == "ok" and not args.process_only:
            result = commit_one(result)
            if result.message == "committed":
                LOG.debug("committed %s → %s", src, result.dst)
            else:
                n_fail_running += 1
                LOG.error("commit failed %s: %s", src, result.message)
        elif result.status == "ok":
            LOG.debug("verified (not committed) %s  (would publish to %s)", src, result.dst)
        else:
            n_fail_running += 1
            LOG.error("%s: %s — %s", result.status, src, result.message)
        man.results.append(result)
        pbar.set_postfix(
            file=path.basename(src)[:36],
            ok=sum(1 for r in man.results if r.status == "ok"),
            fail=n_fail_running,
            refresh=False,
        )

    man.save(manifest_path)
    n_ok = sum(1 for r in man.results if r.status == "ok")
    n_fail = sum(1 for r in man.results if r.status == "failed")
    n_skip = sum(1 for r in man.results if r.status == "skipped")
    LOG.info(
        "done: ok=%d failed=%d skipped=%d  manifest=%s",
        n_ok,
        n_fail,
        n_skip,
        manifest_path,
    )
    if n_fail:
        LOG.error(
            "Failed files leave originals intact; inspect staging tmp + manifest. "
            "Re-run with --commit-only after fixing, or delete bad tmp files."
        )
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
