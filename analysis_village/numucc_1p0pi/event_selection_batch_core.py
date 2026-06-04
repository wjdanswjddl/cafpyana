"""Shared helpers for the batched event-selection map phase.

Loads multiple ``.df`` files into one in-memory bundle (with ``__ntuple`` remapping),
runs the notebook pipeline on that bundle, and writes histogram/count pickles compatible
with ``scripts/event_selection_aggregate.py``.
"""
from __future__ import annotations

import gc
import os
from os import path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from pyanalib.split_df_helpers_new import (
    _concat_hdf_frames,
    _remap_ntuple_index,
    _unique_ntuple_values_across_keys,
    get_n_split,
    load_dfs,
)
from pyanalib.pandas_helpers import pad_column_name

from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_runner
from analysis_village.numucc_1p0pi.evt_derived_kinematics import (
    ensure_derived_trk_kinematics_cols,
    ensure_mc_level_phi_mcnu,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key


def hdf_has_mcnu(df_file: str) -> bool:
    try:
        with pd.HDFStore(df_file, mode="r") as store:
            keys = store.keys()
        return any(str(k).startswith("/mcnu_") for k in keys)
    except Exception:
        return False


def hdr_chunk_pot(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "pot" not in hdr_df.columns:
        return 0.0
    return float(hdr_df["pot"].sum())


def hdr_data_gates_bnb(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "nbnbinfo" not in hdr_df.columns:
        return 0.0
    return float(hdr_df["nbnbinfo"].sum())


def hdr_cosmic_gates_intime(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "ngenevt" not in hdr_df.columns:
        return 0.0
    return float(hdr_df.loc[hdr_df["first_in_subrun"] == 1, "ngenevt"].sum())


def hdr_cosmic_gates_offbeam(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "noffbeambnb" not in hdr_df.columns:
        return 0.0
    return float(hdr_df.loc[hdr_df["first_in_subrun"] == 1, "noffbeambnb"].sum())


def file_hdr_meta(sample: str, hdr_df: pd.DataFrame | None, df_file: str) -> Dict[str, Any]:
    meta = {
        "path": df_file,
        "size_bytes": os.path.getsize(df_file) if os.path.isfile(df_file) else 0,
        "pot": hdr_chunk_pot(hdr_df),
    }
    if sample == "data":
        meta["gates_bnb"] = hdr_data_gates_bnb(hdr_df)
    elif sample == "intime":
        meta["cosmic_gates_intime"] = hdr_cosmic_gates_intime(hdr_df)
    elif sample == "offbeam":
        meta["cosmic_gates_offbeam"] = hdr_cosmic_gates_offbeam(hdr_df)
    return meta


def accumulate_hdr_meta(
    sample: str,
    hdr_df: pd.DataFrame | None,
    chunk_pot: List[float],
    chunk_gates_bnb: List[float],
    chunk_cosmic_gates_intime: List[float],
    chunk_cosmic_gates_offbeam: List[float],
) -> None:
    chunk_pot[0] += hdr_chunk_pot(hdr_df)
    if sample == "data":
        chunk_gates_bnb[0] += hdr_data_gates_bnb(hdr_df)
    elif sample == "intime":
        chunk_cosmic_gates_intime[0] += hdr_cosmic_gates_intime(hdr_df)
    elif sample == "offbeam":
        chunk_cosmic_gates_offbeam[0] += hdr_cosmic_gates_offbeam(hdr_df)


def intrinsic_weight_series(
    df: pd.DataFrame | None,
    sample: str,
    use_mc_genweight: bool = False,
) -> np.ndarray:
    if df is None or len(df) == 0:
        return np.ones(0, dtype=float)
    if sample == "data":
        return np.ones(len(df), dtype=float)
    if sample in ("intime", "offbeam"):
        return np.ones(len(df), dtype=float)
    if sample in ("mc", "dirt"):
        if not use_mc_genweight:
            return np.ones(len(df), dtype=float)
        try:
            gw = df["mc"]["genweight"]
            w = np.asarray(gw, dtype=float).reshape(-1)
            return np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
        except Exception:
            return np.ones(len(df), dtype=float)
    raise ValueError(sample)


def attach_intrinsic_weights(
    evt_df: pd.DataFrame | None,
    trk_df: pd.DataFrame | None,
    sample: str,
    use_mc_genweight: bool = False,
) -> None:
    if evt_df is not None and len(evt_df) > 0:
        evt_df["pot_weight"] = intrinsic_weight_series(evt_df, sample, use_mc_genweight)
    if trk_df is not None and len(trk_df) > 0:
        trk_df["pot_weight"] = intrinsic_weight_series(trk_df, sample, use_mc_genweight)


def prefix_mcnu_columns(mc_nu_df: pd.DataFrame) -> None:
    if isinstance(mc_nu_df.columns, pd.MultiIndex):
        try:
            first_level = mc_nu_df.columns.get_level_values(0)
            need_prefix = not np.all(first_level == "mc")
        except Exception:
            need_prefix = True
        if need_prefix:
            mc_nu_df.columns = pd.MultiIndex.from_tuples(
                [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
            )


def ensure_trk_phi_col(trk_df: pd.DataFrame | None) -> None:
    if trk_df is None or len(trk_df) == 0:
        return
    if not isinstance(trk_df.columns, pd.MultiIndex):
        return
    if multicol_resolve_column_key(trk_df, ("pfp", "trk", "phi", "", "", "")) is not None:
        return
    kx = multicol_resolve_column_key(trk_df, ("pfp", "trk", "dir", "x", ""))
    ky = multicol_resolve_column_key(trk_df, ("pfp", "trk", "dir", "y", ""))
    if kx is None or ky is None:
        return
    phi_col = pad_column_name(("pfp", "trk", "phi", "", "", ""), trk_df)
    trk_df.loc[:, phi_col] = np.degrees(
        np.arctan2(
            np.asarray(trk_df.loc[:, kx], dtype=float),
            np.asarray(trk_df.loc[:, ky], dtype=float),
        )
    )


def ensure_phi_and_kinematics_cols(
    evt_df: pd.DataFrame,
    trk_df: pd.DataFrame | None,
    mcnu_df: pd.DataFrame | None,
) -> Tuple[pd.DataFrame, pd.DataFrame | None]:
    evt_df = ensure_derived_trk_kinematics_cols(evt_df)
    ensure_trk_phi_col(trk_df)
    if mcnu_df is not None and len(mcnu_df) > 0:
        prefix_mcnu_columns(mcnu_df)
        mcnu_df = ensure_mc_level_phi_mcnu(mcnu_df)
    return evt_df, mcnu_df


def load_and_concat_df_files(
    df_files: Sequence[str],
    keys2load: Sequence[str],
    *,
    max_splits_per_file: int | None = None,
) -> Tuple[Dict[str, pd.DataFrame], List[Dict[str, Any]]]:
    """Load and concatenate multiple HDF ``.df`` files (notebook-style index remapping)."""
    if not df_files:
        raise ValueError("df_files is empty")
    df_lists: Dict[str, List[pd.DataFrame]] = {k: [] for k in keys2load}
    per_file_meta: List[Dict[str, Any]] = []
    ntuple_offset = np.int64(0)

    for df_file in df_files:
        n_splits = get_n_split(df_file)
        cap = n_splits if max_splits_per_file is None else min(n_splits, max_splits_per_file)
        file_dfs = load_dfs(df_file, list(keys2load), n_max_concat=cap)
        per_file_meta.append(file_hdr_meta(sample, hdr_df, df_file))

        unique_ntuples = _unique_ntuple_values_across_keys(file_dfs, keys2load)
        ntuple_remap = {
            old: np.int64(ntuple_offset + i) for i, old in enumerate(unique_ntuples)
        }
        ntuple_offset += np.int64(len(ntuple_remap))

        for k in keys2load:
            df = file_dfs[k]
            _remap_ntuple_index(df, ntuple_remap)
            df_lists[k].append(df)

    out = {k: _concat_hdf_frames(df_lists[k], label=k) for k in keys2load}
    return out, per_file_meta


def run_batch_selection(
    sample: str,
    df_files: Sequence[str],
    out_path: str,
    *,
    job_id: str = "",
    use_mc_genweight: bool = False,
    mc_univ_syst_tags: Sequence[str] = (),
    max_splits_per_file: int | None = None,
    pipeline_trace=None,
) -> Dict[str, Any]:
    """Load a file batch, run the notebook pipeline, write one pickle."""
    keys = ["evt", "trk", "hdr"]
    load_mcnu = sample == "mc" and any(hdf_has_mcnu(f) for f in df_files)
    keys_load = keys + (["mcnu"] if load_mcnu else [])

    mc_univ_tags = tuple(mc_univ_syst_tags) if sample == "mc" else ()
    runner = build_runner(sample, mc_univ_syst_tags=mc_univ_tags or None)

    chunk_pot = [0.0]
    chunk_gates_bnb = [0.0]
    chunk_cosmic_gates_intime = [0.0]
    chunk_cosmic_gates_offbeam = [0.0]
    per_file_meta: List[Dict[str, Any]] = []
    n_evt_total = 0

    for df_file in df_files:
        n_splits = get_n_split(df_file)
        cap = n_splits if max_splits_per_file is None else min(n_splits, max_splits_per_file)
        file_keys = keys_load
        file_dfs = load_dfs(df_file, file_keys, n_max_concat=cap)

        hdr_df = file_dfs.get("hdr")
        per_file_meta.append(file_hdr_meta(sample, hdr_df, df_file))
        accumulate_hdr_meta(
            sample,
            hdr_df,
            chunk_pot,
            chunk_gates_bnb,
            chunk_cosmic_gates_intime,
            chunk_cosmic_gates_offbeam,
        )

        evt_df = file_dfs["evt"]
        trk_df = file_dfs["trk"]
        mcnu_df = file_dfs.get("mcnu") if load_mcnu else None

        attach_intrinsic_weights(evt_df, trk_df, sample, use_mc_genweight)
        evt_df, mcnu_df = ensure_phi_and_kinematics_cols(evt_df, trk_df, mcnu_df)
        n_evt_total += int(len(evt_df))

        runner.run(
            {"evt": evt_df, "trk": trk_df, "hdr": hdr_df, "mcnu": mcnu_df},
            pipeline_trace=pipeline_trace,
        )

        del file_dfs, evt_df, trk_df, hdr_df, mcnu_df
        gc.collect()

    meta = {
        "weight_scheme": "intrinsic",
        "use_mc_genweight": bool(use_mc_genweight),
        "mc_univ_syst_tags": list(mc_univ_tags),
        "df_files": list(df_files),
        "per_file": per_file_meta,
        "sample": sample,
        "job_id": job_id,
        "chunk_pot": chunk_pot[0],
        "chunk_gates_bnb": chunk_gates_bnb[0],
        "chunk_cosmic_gates_intime": chunk_cosmic_gates_intime[0],
        "chunk_cosmic_gates_offbeam": chunk_cosmic_gates_offbeam[0],
        "n_evt": n_evt_total,
        "n_files": len(df_files),
        "workflow": "batched_notebook",
        "mc_efficiency_enabled": load_mcnu,
    }
    runner.save(out_path, extra_meta=meta)
    return meta
