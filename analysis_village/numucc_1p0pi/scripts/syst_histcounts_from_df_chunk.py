#!/usr/bin/env python
"""Map one sel_all weight ``.df`` → histcounts HDF (``syst_hists`` + ``var_configs``).

Starts from cached ``sel_all`` tables that already carry Flux / G4 multisim
weights (no CAF / xrootd). Re-runs the numuCC 1p0pi selection walk and
accumulates **rate** histcounts with a vectorized weight-matrix path (same
idea as ``syst_multisim_chunk --input-stage sel_all``), then packs the long
``syst_hists`` schema expected by ``systematics-histcounts.ipynb``.

Typical inputs::

    .../2026_09_04_172234__sel_all-wgts_flux-corrected_updated/*.df
    .../2026_09_04_175940__sel_all-wgts_g4-corrected_updated/*.df

Usage::

    python syst_histcounts_from_df_chunk.py \\
      --df-file PATH.df --out-dir OUTDIR --family Flux \\
      [--n-universe 1000] [--include-slim/--no-include-slim] [--sample mc]
"""
from __future__ import annotations

import argparse
import gc
import os
import shutil
import sys
import tempfile
import warnings
from os import path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from pyanalib.split_df_helpers import get_n_split

from analysis_village.numucc_1p0pi.evt_derived_kinematics import ensure_derived_trk_kinematics_cols
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from analysis_village.numucc_1p0pi.syst_histcounts import (
    attach_family_slim_products,
    empty_histcounts_df,
    histcounts_var_configs_df,
    pack_blob_to_df,
    slim_product_names,
    sum_histcounts_dfs,
)
from analysis_village.numucc_1p0pi.syst_multisim_common import (
    drop_bad_flux_knob_weights,
    drop_bad_g4_knob_weights,
    flux_mc_knob_names,
    g4_mc_knob_names,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    final_stage_var_configs,
    get_var_series,
    histogram_var,
    walk_pipeline,
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--df-file", required=True, help="Input sel_all weight HDF (.df)")
    p.add_argument("--out-dir", required=True, help="Directory for output histcounts .df")
    p.add_argument(
        "--family",
        required=True,
        choices=("Flux", "G4", "flux", "g4"),
        help="Weight family on the DF (Flux or G4)",
    )
    p.add_argument(
        "--n-universe",
        type=int,
        default=1000,
        help="Multisim universe count (default 1000; must match DF build)",
    )
    p.add_argument(
        "--include-slim",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Attach slim_multisim / slim products before fill (default: on)",
    )
    p.add_argument(
        "--sample",
        default="mc",
        choices=("mc", "dirt"),
        help="Selection sample tag for walk_pipeline (default mc)",
    )
    p.add_argument(
        "--flux-knob-groups",
        default="all",
        help="Passed to flux_mc_knob_names when --family Flux (default all)",
    )
    p.add_argument(
        "--max-splits",
        type=int,
        default=0,
        help="If >0, only process the first N HDF splits (debug)",
    )
    p.add_argument(
        "--out-name",
        default="",
        help="Output basename (default: hist_mc_<family>__<input_stem>.df)",
    )
    p.add_argument(
        "--hist-backend",
        default="vector",
        choices=("vector", "loop"),
        help="Universe histogram backend: vector (default) or loop (legacy, slow)",
    )
    return p.parse_args(argv)


def _family_norm(raw: str) -> str:
    s = str(raw).strip()
    if s.lower() == "flux":
        return "Flux"
    if s.lower() == "g4":
        return "G4"
    return s


def _knob_names(family: str, flux_knob_groups: str) -> Tuple[str, ...]:
    if family == "Flux":
        return flux_mc_knob_names(flux_knob_groups)
    if family == "G4":
        return g4_mc_knob_names()
    raise ValueError("unsupported family %r" % family)


def _drop_bad(evt: pd.DataFrame, family: str, knobs: Sequence[str], n_univ: int) -> pd.DataFrame:
    if family == "Flux":
        return drop_bad_flux_knob_weights(evt, knobs=knobs, n_univ=n_univ)
    if family == "G4":
        return drop_bad_g4_knob_weights(evt, knobs=knobs, n_univ=n_univ)
    return evt


def _out_path(args: argparse.Namespace, family: str) -> str:
    stem = path.splitext(path.basename(args.df_file))[0]
    if args.out_name:
        leaf = args.out_name if args.out_name.endswith(".df") else args.out_name + ".df"
    else:
        tag = "flux" if family == "Flux" else "g4"
        leaf = "hist_mc_%s__%s.df" % (tag, stem)
    return path.join(args.out_dir, leaf)


def _count_univ(evt_df: pd.DataFrame, knob: str, cap: int = 2048) -> int:
    n = 0
    for i in range(cap):
        if multicol_resolve_column_key(evt_df, ("mc", knob, "univ_%d" % i)) is None:
            break
        n += 1
    return n


def _weight_matrix(evt_df: pd.DataFrame, knob: str, n_univ: int) -> Optional[np.ndarray]:
    rows = []
    for i in range(n_univ):
        key = multicol_resolve_column_key(evt_df, ("mc", knob, "univ_%d" % i))
        if key is None:
            break
        w = np.asarray(evt_df[key], dtype=float).reshape(-1)
        rows.append(np.nan_to_num(w, nan=1.0, posinf=1.0, neginf=1.0))
    if not rows:
        return None
    return np.vstack(rows)


def _stage_plots(stage_key: str, final_vcs: Sequence[Any]) -> List[Tuple[Any, str]]:
    out: List[Tuple[Any, str]] = []
    for spec in CUT_STAGE_VAR_SPECS:
        if spec.stage_key == stage_key:
            out.append((spec.var_config, spec.target))
    if stage_key == FINAL_STAGE_KEY:
        for vc in final_vcs:
            out.append((vc, "evt"))
    return out


def _add_rate(
    rate: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    knob: str,
    slug: str,
    cv: np.ndarray,
    univ: np.ndarray,
) -> None:
    slot = rate.setdefault(knob, {}).setdefault(
        slug,
        {
            "cv": np.zeros_like(cv, dtype=np.float64),
            "univ": np.zeros_like(univ, dtype=np.float64),
        },
    )
    # Grow univ axis if needed (first write may have fewer rows).
    if slot["univ"].shape[0] < univ.shape[0]:
        grown = np.zeros((univ.shape[0], slot["univ"].shape[1]), dtype=np.float64)
        grown[: slot["univ"].shape[0]] = slot["univ"]
        slot["univ"] = grown
    n = min(slot["univ"].shape[0], univ.shape[0])
    slot["cv"] += cv
    slot["univ"][:n] += univ[:n]


def _prep_hist_values(
    values: np.ndarray,
    bins: np.ndarray,
    idx: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Match :func:`histogram_var` clipping; return ``(v, bi, idx_kept)``.

    ``bi[i]`` is the bin index for ``v[i]`` (same edges as ``np.histogram`` after
    the shared clip-to-``bins[-1]-eps`` treatment).
    """
    bins = np.asarray(bins, dtype=float)
    nb = int(len(bins) - 1)
    v = np.asarray(values, dtype=float).reshape(-1)
    idx_a = np.asarray(idx, dtype=np.int64).reshape(-1)
    if v.shape[0] != idx_a.shape[0]:
        raise ValueError("values/idx length mismatch: %d vs %d" % (v.shape[0], idx_a.shape[0]))
    if v.size == 0 or nb <= 0:
        return (
            np.zeros(0, dtype=np.float64),
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=np.int64),
        )
    finite = np.isfinite(v)
    if not finite.all():
        v = v[finite]
        idx_a = idx_a[finite]
    if v.size == 0:
        return v, np.zeros(0, dtype=np.int64), idx_a
    eps = (float(bins[-1]) - float(bins[0])) * 1e-9
    v = np.clip(v, bins[0], bins[-1] - eps)
    # After clipping below the upper edge, all bins are half-open [e_i, e_{i+1}).
    bi = np.searchsorted(bins, v, side="right") - 1
    bi = np.clip(bi, 0, nb - 1).astype(np.int64, copy=False)
    return v, bi, idx_a


def _univ_hist_from_bins(
    wmat: np.ndarray,
    idx_kept: np.ndarray,
    bi: np.ndarray,
    nb: int,
) -> np.ndarray:
    """Weighted fills for all universes: ``out[u, b] = sum_j wmat[u, idx[j]]`` with ``bi[j]==b``."""
    n_univ = int(wmat.shape[0])
    if bi.size == 0 or nb <= 0:
        return np.zeros((n_univ, nb), dtype=np.float64)
    w = np.asarray(wmat[:, idx_kept], dtype=np.float64)
    w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    # one-hot (n, nb) @ — equivalent to per-universe np.histogram(..., weights=...)
    n = int(bi.shape[0])
    oh = np.zeros((n, nb), dtype=np.float64)
    oh[np.arange(n, dtype=np.int64), bi] = 1.0
    return w @ oh


def _univ_hist_loop(
    values: np.ndarray,
    bins: np.ndarray,
    wmat: np.ndarray,
    idx: np.ndarray,
) -> np.ndarray:
    """Reference: per-universe :func:`histogram_var` (slow; for validation only)."""
    nb = int(len(bins) - 1)
    n_row = int(wmat.shape[0])
    univ_h = np.zeros((n_row, nb), dtype=np.float64)
    for u in range(n_row):
        univ_h[u] = histogram_var(values, bins, weights=wmat[u][idx])
    return univ_h


def _accumulate_split(
    evt: pd.DataFrame,
    trk: pd.DataFrame,
    *,
    family: str,
    knobs: Sequence[str],
    n_univ: int,
    sample: str,
    hist_backend: str = "vector",
) -> Dict[str, Any]:
    """Walk selection; return histcounts rate blob for one split.

    ``hist_backend``: ``vector`` (default, digitize once + weight matmul) or
    ``loop`` (legacy per-universe ``histogram_var`` — for numerical checks).
    """
    backend = str(hist_backend or "vector").strip().lower()
    if backend not in ("vector", "loop"):
        raise ValueError("hist_backend must be 'vector' or 'loop', got %r" % hist_backend)

    blob: Dict[str, Any] = {"family": family, "rate": {}, "xsec": {}}
    final_vcs = list(final_stage_var_configs())
    state: Dict[str, Any] = {"evt": evt, "trk": trk, "hdr": None, "mcnu": None}

    for stage_key, post_state in walk_pipeline(state, sample=sample):
        plots = _stage_plots(stage_key, final_vcs)
        if not plots:
            continue
        post_evt = post_state.get("evt")
        if post_evt is None or len(post_evt) == 0:
            continue
        post_evt = ensure_derived_trk_kinematics_cols(post_evt)
        post_state = dict(post_state)
        post_state["evt"] = post_evt

        # Resolve weight matrices once per stage (knobs present on this frame).
        wmats: Dict[str, np.ndarray] = {}
        for knob in knobs:
            n_u = min(n_univ, _count_univ(post_evt, knob))
            if n_u <= 0:
                continue
            wmat = _weight_matrix(post_evt, knob, n_u)
            if wmat is not None:
                wmats[knob] = wmat
        if not wmats:
            continue

        for vc, target in plots:
            slug = vc.var_save_name
            bins = np.asarray(vc.bins, dtype=float)
            nb = len(bins) - 1
            if slug == "integrated":
                cv_h = np.array([float(len(post_evt))], dtype=np.float64)
                for knob, wmat in wmats.items():
                    univ_h = np.sum(wmat, axis=1).reshape(-1, 1).astype(np.float64)
                    _add_rate(blob["rate"], knob, slug, cv_h, univ_h)
                continue

            got = get_var_series(post_state, vc, target)
            if got is None:
                continue
            values, idx = got
            cv_h = histogram_var(values, bins)
            if backend == "loop":
                for knob, wmat in wmats.items():
                    univ_h = _univ_hist_loop(values, bins, wmat, idx)
                    _add_rate(blob["rate"], knob, slug, cv_h, univ_h)
            else:
                _v, bi, idx_kept = _prep_hist_values(values, bins, idx)
                for knob, wmat in wmats.items():
                    univ_h = _univ_hist_from_bins(wmat, idx_kept, bi, nb)
                    _add_rate(blob["rate"], knob, slug, cv_h, univ_h)

    return blob


def process_one_file(args: argparse.Namespace) -> str:
    family = _family_norm(args.family)
    knobs = list(_knob_names(family, args.flux_knob_groups))
    n_univ = int(args.n_universe)
    sample = str(args.sample)

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if int(args.max_splits) <= 0 else min(int(args.max_splits), n_keys)
    if n_use <= 0:
        raise SystemExit("[histcounts-from-df] no splits in %s" % args.df_file)

    parts: List[pd.DataFrame] = []
    for i in range(n_use):
        evt = pd.read_hdf(args.df_file, key="evt_%d" % i)
        try:
            trk = pd.read_hdf(args.df_file, key="trk_%d" % i)
        except (KeyError, ValueError) as ex:
            print(
                "[histcounts-from-df] missing trk_%d in %s (%s); skip split"
                % (i, args.df_file, ex)
            )
            del evt
            gc.collect()
            continue

        evt = _drop_bad(evt, family, knobs, n_univ)
        if evt is None or len(evt) == 0:
            del evt, trk
            gc.collect()
            continue
        try:
            trk = trk.loc[trk.index.intersection(evt.index)]
        except Exception:
            pass

        live_knobs = list(knobs)
        if args.include_slim:
            evt, slim_names = attach_family_slim_products(
                evt, family=family, n_univ=n_univ, knob_names=list(knobs)
            )
            live_knobs = list(knobs) + list(slim_names)
            # Also accept short names if attach used GENIE-style slim labels.
            name_ms, name_full = slim_product_names(family)
            for extra in (name_ms, name_full, "slim", "slim_multisim"):
                if extra not in live_knobs:
                    live_knobs.append(extra)

        blob = _accumulate_split(
            evt,
            trk,
            family=family,
            knobs=live_knobs,
            n_univ=n_univ,
            sample=sample,
            hist_backend=getattr(args, "hist_backend", "vector"),
        )
        hist = pack_blob_to_df(blob)
        if hist is not None and len(hist) > 0:
            parts.append(hist)
        del evt, trk, blob, hist
        gc.collect()

    out = sum_histcounts_dfs(parts) if parts else empty_histcounts_df()
    os.makedirs(args.out_dir, exist_ok=True)
    dest = _out_path(args, family)

    # Write on local disk first — PyTables + dCache/pnfs often leaves 0-byte or
    # corrupt files when using mkstemp+replace directly on /pnfs.
    local_dir = tempfile.mkdtemp(prefix="histcounts_from_df_")
    tmp = path.join(local_dir, path.basename(dest))
    try:
        with pd.HDFStore(tmp, mode="w", complevel=5, complib="zlib") as store:
            store.put("syst_hists_0", out, format="table")
            store.put("var_configs_0", histcounts_var_configs_df(), format="table")
            store.put("split", pd.DataFrame({"n_split": [1]}), format="table")
        shutil.copy2(tmp, dest)
    finally:
        try:
            shutil.rmtree(local_dir, ignore_errors=True)
        except OSError:
            pass

    n_knobs = out["knob"].nunique() if len(out) else 0
    n_vars = out["var"].nunique() if len(out) else 0
    print(
        "[histcounts-from-df] wrote %s  rows=%d  family=%s  knobs=%d  vars=%d"
        % (dest, len(out), family, n_knobs, n_vars),
        flush=True,
    )
    return dest


def run_with_args(args: argparse.Namespace, *, skip_existing: bool = False):
    family = _family_norm(args.family)
    dest = _out_path(args, family)
    if skip_existing and path.isfile(dest) and path.getsize(dest) > 0:
        print("[histcounts-from-df] skip existing %s" % dest, flush=True)
        return dest, "skipped"
    out = process_one_file(args)
    return out, "ok"


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    run_with_args(args, skip_existing=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
