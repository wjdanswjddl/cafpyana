#!/usr/bin/env python3
"""SCE detector-variation selection comparison (Sept 2026 sel_mup samples).

Find truth-level common events across CV / 0xSCE / 2xSCE, compare which pass
the mup selection (``evt``), accumulate kinematic histograms, and write
exclusive-selection event lists.

Memory-safe: never loads the full multi-million-entry common-key set into RAM.
Common keys live in SQLite (built one variation at a time). Each variation is
processed file-by-file with SQLite probes for membership checks.
"""

from __future__ import annotations

import argparse
import gc
import glob
import os
import pickle
import sqlite3
import warnings
from os import makedirs, path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

import sys

_REPO_ROOT = path.normpath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from pyanalib.split_df_helpers_new import get_n_split
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
    with_final_selected_evt_variables,
)

EventKey = Tuple[float, int, int, int]

DEFAULT_BASE = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"
DEFAULT_VARIATIONS = {
    "CV": "2026_09_01_195328__sel_mup-mc-BNB_cosmics-CV",
    "0xSCE": "2026_09_01_194801__sel_mup-mc-BNB_cosmics-0xSCE",
    "2xSCE": "2026_09_01_195043__sel_mup-mc-BNB_cosmics-2xSCE",
}
DEFAULT_OUT = (
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/SCE/"
    "sce_selection_sept2026"
)

META_COLS = ["E", "run", "subrun", "evt"]


def list_df_files(search_dir: str) -> List[str]:
    files = sorted(glob.glob(path.join(search_dir, "*.df")))
    return [f for f in files if "_matched" not in path.basename(f)]


def memory_status_gb() -> Tuple[float, float, float]:
    meminfo: Dict[str, int] = {}
    with open("/proc/meminfo", encoding="ascii") as fh:
        for line in fh:
            key, value = line.split(":", 1)
            meminfo[key] = int(value.strip().split()[0])
    total_kb = meminfo["MemTotal"]
    avail_kb = meminfo.get("MemAvailable", meminfo.get("MemFree", 0))
    used_kb = total_kb - avail_kb
    total_gb = total_kb / (1024**2)
    used_gb = used_kb / (1024**2)
    return total_gb, used_gb, used_gb / total_gb if total_gb else 0.0


def process_rss_gb() -> float:
    try:
        with open("/proc/self/status", encoding="ascii") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / (1024**2)
    except OSError:
        pass
    return 0.0


def _var_cache_path(cache_dir: str, stem: str, variation: str) -> str:
    safe = variation.replace("/", "_").replace(" ", "_")
    return path.join(cache_dir, f"{stem}_{safe}.pkl")


def assert_memory_ok(max_used_fraction: float, label: str = "") -> None:
    gc.collect()
    _, _, frac = memory_status_gb()
    if frac <= max_used_fraction:
        return
    gc.collect()
    _, used_gb, frac = memory_status_gb()
    if frac > max_used_fraction:
        raise MemoryError(
            f"System memory usage {100 * frac:.1f}% ({used_gb:.1f} GiB) exceeds limit "
            f"{100 * max_used_fraction:.1f}%"
            + (f" during {label}" if label else "")
        )


def _flatten_index_cols(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    flat = df.reset_index()
    if isinstance(flat.columns, pd.MultiIndex):
        flat = flat.copy()
        flat.columns = flat.columns.get_level_values(0)
    missing = [c for c in cols if c not in flat.columns]
    if missing:
        raise KeyError(f"Missing columns after flatten: {missing}")
    return flat[list(cols)]


def _meta_keys_from_file(fpath: str) -> Set[EventKey]:
    keys: Set[EventKey] = set()
    n_split = get_n_split(fpath)
    for i in range(n_split):
        try:
            meta = pd.read_hdf(fpath, key=f"meta_{i}", columns=META_COLS)
        except Exception:
            meta = pd.read_hdf(fpath, key=f"meta_{i}")
            meta = meta[META_COLS]
        keys.update(map(tuple, meta[META_COLS].values))
        del meta
    return keys


def _meta_lookup_df(fpath: str, split_i: int) -> pd.DataFrame:
    return _flatten_index_cols(
        pd.read_hdf(fpath, key=f"meta_{split_i}"),
        ["__ntuple", "entry", *META_COLS],
    )


def _file_tag(fpath: str) -> str:
    return path.basename(fpath)


class CommonKeyStore:
    """Disk-backed common truth-level event keys."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")

    def close(self) -> None:
        self.conn.close()

    def is_built(self) -> bool:
        row = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='common_keys'"
        ).fetchone()
        return row is not None

    def count(self) -> int:
        if not self.is_built():
            return 0
        row = self.conn.execute("SELECT COUNT(*) FROM common_keys").fetchone()
        return int(row[0]) if row else 0

    def import_from_pickle(self, pkl_path: str, batch_size: int = 100_000) -> None:
        self.reset()
        keys: Set[EventKey] = pickle.load(open(pkl_path, "rb"))
        self.conn.execute(
            """
            CREATE TABLE common_keys (
                E REAL NOT NULL,
                run INTEGER NOT NULL,
                subrun INTEGER NOT NULL,
                evt INTEGER NOT NULL,
                PRIMARY KEY (run, subrun, evt, E)
            )
            """
        )
        rows = [(float(e), int(r), int(sr), int(ev)) for e, r, sr, ev in keys]
        del keys
        for i in range(0, len(rows), batch_size):
            self.conn.executemany(
                "INSERT OR IGNORE INTO common_keys VALUES (?, ?, ?, ?)",
                rows[i : i + batch_size],
            )
        del rows
        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_common_evt "
            "ON common_keys (run, subrun, evt, E)"
        )
        self.conn.commit()
        gc.collect()

    def reset(self) -> None:
        self.conn.execute("DROP TABLE IF EXISTS common_keys")
        self.conn.commit()

    def build(
        self,
        base_dir: str,
        variations: Dict[str, str],
        *,
        max_used_fraction: float,
    ) -> None:
        self.reset()
        var_names = list(variations.keys())

        for i, name in enumerate(var_names):
            tbl = f"keys_{i}"
            self.conn.execute(
                f"""
                CREATE TABLE {tbl} (
                    E REAL NOT NULL,
                    run INTEGER NOT NULL,
                    subrun INTEGER NOT NULL,
                    evt INTEGER NOT NULL,
                    PRIMARY KEY (run, subrun, evt, E)
                )
                """
            )
            files = list_df_files(path.join(base_dir, variations[name]))
            for fpath in tqdm(files, desc=f"{name} meta->sqlite"):
                assert_memory_ok(max_used_fraction, label=f"{name} meta")
                keys = _meta_keys_from_file(fpath)
                if keys:
                    self.conn.executemany(
                        f"INSERT OR IGNORE INTO {tbl} VALUES (?, ?, ?, ?)",
                        [(float(e), int(r), int(sr), int(ev)) for e, r, sr, ev in keys],
                    )
                del keys
                gc.collect()

            if i == 0:
                self.conn.execute(
                    f"CREATE TABLE common_keys AS SELECT E, run, subrun, evt FROM {tbl}"
                )
            else:
                self.conn.execute(
                    f"""
                    CREATE TABLE common_keys_new AS
                    SELECT E, run, subrun, evt FROM common_keys
                    INTERSECT
                    SELECT E, run, subrun, evt FROM {tbl}
                    """
                )
                self.conn.execute("DROP TABLE common_keys")
                self.conn.execute(
                    "ALTER TABLE common_keys_new RENAME TO common_keys"
                )

            self.conn.execute(f"DROP TABLE {tbl}")
            self.conn.commit()
            gc.collect()

        self.conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_common_evt "
            "ON common_keys (run, subrun, evt, E)"
        )
        self.conn.commit()

    def probe(self, keys: Sequence[EventKey]) -> Set[EventKey]:
        if not keys:
            return set()
        self.conn.execute("DROP TABLE IF EXISTS _probe")
        self.conn.execute(
            "CREATE TEMP TABLE _probe (E REAL, run INT, subrun INT, evt INT)"
        )
        self.conn.executemany(
            "INSERT INTO _probe VALUES (?, ?, ?, ?)",
            [(float(e), int(r), int(sr), int(ev)) for e, r, sr, ev in keys],
        )
        rows = self.conn.execute(
            """
            SELECT p.E, p.run, p.subrun, p.evt
            FROM _probe p
            INNER JOIN common_keys c
              ON p.E = c.E AND p.run = c.run AND p.subrun = c.subrun AND p.evt = c.evt
            """
        ).fetchall()
        self.conn.execute("DROP TABLE IF EXISTS _probe")
        return {(float(e), int(r), int(sr), int(ev)) for e, r, sr, ev in rows}


class SelectedKeyStore:
    """Disk-backed store for selected event keys and per-file progress."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS selected (
                variation TEXT NOT NULL,
                E REAL NOT NULL,
                run INTEGER NOT NULL,
                subrun INTEGER NOT NULL,
                evt INTEGER NOT NULL,
                PRIMARY KEY (variation, E, run, subrun, evt)
            )
            """
        )
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS progress (
                variation TEXT NOT NULL,
                file_tag TEXT NOT NULL,
                PRIMARY KEY (variation, file_tag)
            )
            """
        )
        self.conn.commit()

    def close(self) -> None:
        self.conn.close()

    def reset(self) -> None:
        self.conn.execute("DELETE FROM selected")
        self.conn.execute("DELETE FROM progress")
        self.conn.commit()

    def completed_files(self, variation: str) -> Set[str]:
        rows = self.conn.execute(
            "SELECT file_tag FROM progress WHERE variation=?",
            (variation,),
        ).fetchall()
        return {r[0] for r in rows}

    def mark_file_done(self, variation: str, file_tag: str) -> None:
        self.conn.execute(
            "INSERT OR IGNORE INTO progress (variation, file_tag) VALUES (?, ?)",
            (variation, file_tag),
        )
        self.conn.commit()

    def insert_keys(self, variation: str, keys: Iterable[EventKey]) -> int:
        rows = [(variation, float(e), int(r), int(sr), int(ev)) for e, r, sr, ev in keys]
        if not rows:
            return 0
        self.conn.executemany(
            "INSERT OR IGNORE INTO selected (variation, E, run, subrun, evt) "
            "VALUES (?, ?, ?, ?, ?)",
            rows,
        )
        self.conn.commit()
        return len(rows)

    def count_selected(self, variation: str) -> int:
        row = self.conn.execute(
            "SELECT COUNT(*) FROM selected WHERE variation=?",
            (variation,),
        ).fetchone()
        return int(row[0]) if row else 0

    def count_all_three(self, variations: Sequence[str]) -> int:
        if len(variations) != 3:
            return 0
        row = self.conn.execute(
            """
            SELECT COUNT(*) FROM (
                SELECT E, run, subrun, evt FROM selected WHERE variation=?
                INTERSECT
                SELECT E, run, subrun, evt FROM selected WHERE variation=?
                INTERSECT
                SELECT E, run, subrun, evt FROM selected WHERE variation=?
            )
            """,
            tuple(variations),
        ).fetchone()
        return int(row[0]) if row else 0

    def exclusive_dataframe(self, variation: str, others: Sequence[str]) -> pd.DataFrame:
        if not others:
            query = (
                "SELECT E, run, subrun, evt FROM selected WHERE variation=? "
                "ORDER BY run, subrun, evt, E"
            )
            params: Sequence[str] = (variation,)
        else:
            placeholders = ",".join("?" for _ in others)
            query = f"""
            SELECT E, run, subrun, evt FROM selected WHERE variation=?
            EXCEPT
            SELECT E, run, subrun, evt FROM selected WHERE variation IN ({placeholders})
            ORDER BY run, subrun, evt, E
            """
            params = (variation, *others)
        return pd.read_sql_query(query, self.conn, params=params)


def _per_evt_col(evt_df: pd.DataFrame, col_tuple: tuple) -> np.ndarray:
    try:
        return multicol_get_series(evt_df, col_tuple).to_numpy(dtype=float)
    except Exception:
        return np.array([], dtype=float)


def build_plot_var_defs() -> Dict[str, dict]:
    var_defs: Dict[str, dict] = {}
    final_configs = with_final_selected_evt_variables(list(CORE_SELECTED_EVT_VARIABLE_CONFIGS))
    ne_cfg = VariableConfig.neutrino_energy()
    if ne_cfg.var_save_name not in {c.var_save_name for c in final_configs}:
        final_configs.append(ne_cfg)

    for vc in final_configs:
        if vc.var_save_name == "integrated":
            continue
        col = vc.var_evt_reco_col
        var_defs[vc.var_save_name] = {
            "label": vc.var_labels[0] if vc.var_labels else vc.var_plot_name,
            "bins": np.asarray(vc.bins),
            "col": col,
        }
    return var_defs


def _empty_hists(var_defs: Dict[str, dict]) -> Dict[str, np.ndarray]:
    return {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}


def _fill_hists_from_evt(
    evt_df: pd.DataFrame,
    var_defs: Dict[str, dict],
    hists: Dict[str, np.ndarray],
) -> int:
    if evt_df is None or evt_df.empty:
        return 0
    n = len(evt_df)
    for var_name, cfg in var_defs.items():
        vals = _per_evt_col(evt_df, cfg["col"])
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        bins = cfg["bins"]
        eps = (bins[-1] - bins[0]) * 1e-9
        vals = np.clip(vals, bins[0], bins[-1] - eps)
        counts, _ = np.histogram(vals, bins=bins)
        hists[var_name] += counts
    return n


def _process_evt_file(
    fpath: str,
    common_store: CommonKeyStore,
    var_defs: Dict[str, dict],
    hists: Dict[str, np.ndarray],
) -> Tuple[Set[EventKey], int]:
    selected: Set[EventKey] = set()
    n_hist_rows = 0
    n_split = get_n_split(fpath)

    for i in range(n_split):
        try:
            meta = _meta_lookup_df(fpath, i)
            evt = pd.read_hdf(fpath, key=f"evt_{i}")
        except Exception:
            continue

        evt_flat = evt.reset_index()
        if isinstance(evt_flat.columns, pd.MultiIndex):
            evt_flat.columns = evt_flat.columns.get_level_values(0)

        merged = evt_flat.merge(meta, on=["__ntuple", "entry"], how="inner")
        del meta, evt_flat
        if merged.empty:
            del merged, evt
            continue

        probe_keys = [
            (float(e), int(r), int(sr), int(ev))
            for e, r, sr, ev in merged[META_COLS].drop_duplicates().itertuples(
                index=False, name=None
            )
        ]
        common_in_split = common_store.probe(probe_keys)
        del probe_keys
        if not common_in_split:
            del merged, evt
            continue

        common_df = pd.DataFrame(list(common_in_split), columns=META_COLS)
        sel_merged = merged.merge(common_df, on=META_COLS, how="inner")
        del merged, common_df, common_in_split
        if sel_merged.empty:
            del evt, sel_merged
            continue

        for e, r, sr, ev in sel_merged[META_COLS].itertuples(index=False, name=None):
            selected.add((float(e), int(r), int(sr), int(ev)))

        pair_set = set(zip(sel_merged["__ntuple"].to_numpy(), sel_merged["entry"].to_numpy()))
        nt = evt.index.get_level_values("__ntuple")
        en = evt.index.get_level_values("entry")
        mask = np.array([(a, b) in pair_set for a, b in zip(nt, en)], dtype=bool)
        evt_sel = evt.iloc[np.flatnonzero(mask)]
        del evt, sel_merged, pair_set, mask
        n_hist_rows += _fill_hists_from_evt(evt_sel, var_defs, hists)
        del evt_sel

    return selected, n_hist_rows


def _load_partial_hists(
    cache_dir: str, variation: str, var_defs: Dict[str, dict]
) -> Tuple[Dict[str, np.ndarray], int]:
    partial = _var_cache_path(cache_dir, "variation_hists_partial", variation)
    if path.exists(partial):
        payload = pickle.load(open(partial, "rb"))
        return payload["hists"], payload["n_evts"]
    return _empty_hists(var_defs), 0


def _save_partial_hists(
    cache_dir: str,
    variation: str,
    hists: Dict[str, np.ndarray],
    n_evts: int,
) -> None:
    pickle.dump(
        {"hists": hists, "n_evts": n_evts},
        open(_var_cache_path(cache_dir, "variation_hists_partial", variation), "wb"),
    )


def process_variation(
    variation: str,
    files: Sequence[str],
    common_store: CommonKeyStore,
    selected_store: SelectedKeyStore,
    var_defs: Dict[str, dict],
    *,
    cache_dir: str,
    max_used_fraction: float,
    completed: Optional[Set[str]] = None,
) -> Tuple[Dict[str, np.ndarray], int]:
    completed = completed or set()
    hists, n_evts = _load_partial_hists(cache_dir, variation, var_defs)
    todo = [f for f in files if _file_tag(f) not in completed]

    for fpath in tqdm(todo, desc=f"{variation} evt"):
        assert_memory_ok(max_used_fraction, label=f"{variation} {_file_tag(fpath)}")
        sel_i, n_rows = _process_evt_file(fpath, common_store, var_defs, hists)
        selected_store.insert_keys(variation, sel_i)
        n_evts += n_rows
        del sel_i
        selected_store.mark_file_done(variation, _file_tag(fpath))
        _save_partial_hists(cache_dir, variation, hists, n_evts)
        gc.collect()

    return hists, n_evts


def summary_table_from_store(
    store: SelectedKeyStore,
    variations: Dict[str, str],
    n_common: int,
) -> pd.DataFrame:
    rows = []
    for name in variations:
        n_sel = store.count_selected(name)
        others = [v for v in variations if v != name]
        n_excl = len(store.exclusive_dataframe(name, others))
        rows.append(
            {
                "variation": name,
                "n_selected_common": n_sel,
                "frac_common_selected": n_sel / n_common if n_common else np.nan,
                "n_exclusive": n_excl,
                "frac_common_exclusive": n_excl / n_common if n_common else np.nan,
                "frac_selected_exclusive": n_excl / n_sel if n_sel else np.nan,
            }
        )
    n_all = store.count_all_three(list(variations))
    rows.append(
        {
            "variation": "all_three",
            "n_selected_common": n_all,
            "frac_common_selected": n_all / n_common if n_common else np.nan,
            "n_exclusive": 0,
            "frac_common_exclusive": 0.0,
            "frac_selected_exclusive": np.nan,
        }
    )
    return pd.DataFrame(rows)


def run_analysis(
    *,
    base_dir: str = DEFAULT_BASE,
    variations: Optional[Dict[str, str]] = None,
    out_dir: str = DEFAULT_OUT,
    rerun_meta: bool = False,
    rerun_selection: bool = False,
    rerun_hists: bool = False,
    max_used_fraction: float = 0.70,
) -> dict:
    variations = variations or DEFAULT_VARIATIONS
    cache_dir = path.join(out_dir, "cache")
    excl_dir = path.join(out_dir, "exclusive_events")
    makedirs(cache_dir, exist_ok=True)
    makedirs(excl_dir, exist_ok=True)
    makedirs(path.join(out_dir, "plots"), exist_ok=True)

    total_gb, used_gb, used_frac = memory_status_gb()
    print(
        f"Memory at start: {used_gb:.1f}/{total_gb:.1f} GiB used "
        f"({100 * used_frac:.1f}%), process RSS {process_rss_gb():.2f} GiB",
        flush=True,
    )
    assert_memory_ok(max_used_fraction, label="start")

    common_db = path.join(cache_dir, "common_keys.sqlite")
    common_pkl = path.join(cache_dir, "common_keys.pkl")
    common_store = CommonKeyStore(common_db)
    if rerun_meta or not common_store.is_built():
        if path.exists(common_pkl) and not rerun_meta:
            print(f"Importing common keys from {common_pkl} -> SQLite...", flush=True)
            common_store.import_from_pickle(common_pkl)
        else:
            print("Building common-key SQLite index (one variation at a time)...", flush=True)
            common_store.build(base_dir, variations, max_used_fraction=max_used_fraction)
    n_common = common_store.count()
    print(f"Common truth events: {n_common}", flush=True)

    selected_db = path.join(cache_dir, "selected_keys.sqlite")
    selected_store = SelectedKeyStore(selected_db)
    if rerun_selection or rerun_meta:
        selected_store.reset()
        for name in variations:
            partial = _var_cache_path(cache_dir, "variation_hists_partial", name)
            if path.exists(partial):
                os.remove(partial)

    var_defs = build_plot_var_defs()
    var_defs_payload = {
        k: {"label": v["label"], "bins": v["bins"]} for k, v in var_defs.items()
    }
    hists_cache = path.join(cache_dir, "variation_hists.pkl")
    all_hists: Dict[str, Dict[str, np.ndarray]] = {}
    all_nevts: Dict[str, int] = {}

    if rerun_hists and not rerun_selection:
        for name in variations:
            partial = _var_cache_path(cache_dir, "variation_hists_partial", name)
            if path.exists(partial):
                os.remove(partial)
            selected_store.conn.execute(
                "DELETE FROM progress WHERE variation=?", (name,)
            )
        selected_store.conn.commit()

    need_pass = rerun_selection or rerun_meta or rerun_hists
    for name, subdir in variations.items():
        files = list_df_files(path.join(base_dir, subdir))
        completed = selected_store.completed_files(name)
        if not need_pass and len(completed) >= len(files):
            print(f"Using cached {name} selection/hists", flush=True)
            all_hists[name], all_nevts[name] = _load_partial_hists(
                cache_dir, name, var_defs
            )
            continue

        if rerun_hists and not rerun_selection:
            hists, n_evts = _empty_hists(var_defs), 0
        else:
            hists, n_evts = _load_partial_hists(cache_dir, name, var_defs)

        hists, n_evts = process_variation(
            name,
            files,
            common_store,
            selected_store,
            var_defs,
            cache_dir=cache_dir,
            max_used_fraction=max_used_fraction,
            completed=completed if not (rerun_selection or rerun_meta) else set(),
        )
        all_hists[name] = hists
        all_nevts[name] = n_evts
        pickle.dump(
            {"hists": hists, "n_evts": n_evts},
            open(_var_cache_path(cache_dir, "variation_hists", name), "wb"),
        )

    hist_payload = {
        "var_defs": var_defs_payload,
        "hists": all_hists,
        "n_evts": all_nevts,
    }
    pickle.dump(hist_payload, open(hists_cache, "wb"))

    summary = summary_table_from_store(selected_store, variations, n_common)
    summary.to_csv(path.join(out_dir, "selection_summary.csv"), index=False)

    for name in variations:
        others = [v for v in variations if v != name]
        selected_store.exclusive_dataframe(name, others).to_csv(
            path.join(excl_dir, f"exclusive_{name}.csv"), index=False
        )

    common_store.close()
    selected_store.close()

    total_gb, used_gb, used_frac = memory_status_gb()
    print(
        f"Memory at end: {used_gb:.1f}/{total_gb:.1f} GiB used "
        f"({100 * used_frac:.1f}%), process RSS {process_rss_gb():.2f} GiB",
        flush=True,
    )

    return {
        "n_common": n_common,
        "summary": summary,
        "hist_payload": hist_payload,
        "out_dir": out_dir,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-dir", default=DEFAULT_BASE)
    p.add_argument("--out-dir", default=DEFAULT_OUT)
    p.add_argument(
        "--rerun-meta",
        action="store_true",
        help="Rebuild common-key SQLite index.",
    )
    p.add_argument(
        "--rerun-selection",
        action="store_true",
        help="Reset selected-key store and reprocess all variations.",
    )
    p.add_argument(
        "--rerun-hists",
        action="store_true",
        help="Re-accumulate histograms for all variations.",
    )
    p.add_argument(
        "--max-mem-fraction",
        type=float,
        default=0.70,
        help="Abort if system memory usage exceeds this fraction (default: 0.70).",
    )
    args = p.parse_args(argv)

    result = run_analysis(
        base_dir=args.base_dir,
        out_dir=args.out_dir,
        rerun_meta=args.rerun_meta,
        rerun_selection=args.rerun_selection,
        rerun_hists=args.rerun_hists,
        max_used_fraction=args.max_mem_fraction,
    )
    print("common truth events:", result["n_common"])
    print(result["summary"].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
