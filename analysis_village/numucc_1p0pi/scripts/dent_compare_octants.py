#!/usr/bin/env python3
"""
Per-octant DENT vs CV comparison (distribution overlays + per-event frac-diff).

Octants are defined from the reco slice vertex with boundaries
    x = 0,  y = 0,  z = 250
matching ``data_mc_comparison.ipynb`` / ``selected_events._add_sbnd_octant_labels``.

Saves a recombination-friendly pickle under ``cache/`` so octants can later be
merged into quadrants (or any union) without re-reading the dfs:

    hists[octant][variation][var]  — sum across octants for overlay plots
    frac_vals[octant][var]         — concatenate then re-histogram for frac-diff

Example:
    python dent_compare_octants.py --skip-match --max-files 3 \\
        --out-base /exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT-fracdiff-octants
"""

from __future__ import annotations

import argparse
import gc
import glob
import os
import sys
import warnings
from os import makedirs, path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
os.environ.setdefault("MPLBACKEND", "Agg")

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
_REPO_ROOT = path.normpath(path.join(_SCRIPT_DIR, "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.event_selection_batch_core import (
    attach_intrinsic_weights,
    ensure_phi_and_kinematics_cols,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_get_series
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    get_var_series,
    histogram_var,
    walk_pipeline,
)
from pyanalib.split_df_helpers_new import get_n_split

from analysis_village.numucc_1p0pi.scripts import dent_compare as dc

_OUT_BASE = (
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT-fracdiff-octants"
)

# Boundaries (cm): x=0, y=0, z=250
X0, Y0, Z0 = 0.0, 0.0, 250.0

# Bit packing: id = (x>=0)<<2 | (y>=0)<<1 | (z>=250)<<0  (same as data_mc_comparison)
OCTANT_META: List[dict] = []
for _i in range(8):
    _x_ge = bool((_i >> 2) & 1)
    _y_ge = bool((_i >> 1) & 1)
    _z_ge = bool((_i >> 0) & 1)
    _ew = "W" if _x_ge else "E"  # SBND: negative x is East
    _ns = "N" if _z_ge else "S"  # SBND: lower z is South
    _tb = "Top" if _y_ge else "Bottom"
    OCTANT_META.append(
        {
            "id": _i,
            "name": f"{_ew}-{_ns}-{_tb}",
            "slug": f"oct{_i}_{_ew}{_ns}{_tb}",
            "title": (
                f"x{'≥' if _x_ge else '<'}0, "
                f"y{'≥' if _y_ge else '<'}0, "
                f"z{'≥' if _z_ge else '<'}250"
            ),
            "x_ge0": _x_ge,
            "y_ge0": _y_ge,
            "z_ge250": _z_ge,
            "ew": _ew,
            "ns": _ns,
            "tb": _tb,
        }
    )
N_OCT = 8


def octant_ids_from_xyz(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    x0: float = X0,
    y0: float = Y0,
    z0: float = Z0,
) -> np.ndarray:
    """Return octant id 0..7 per row; -1 where any coordinate is non-finite."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    xid = (x >= x0).astype(np.int8)
    yid = (y >= y0).astype(np.int8)
    zid = (z >= z0).astype(np.int8)
    out = ((xid << 2) | (yid << 1) | zid).astype(np.int8)
    out = np.where(ok, out, np.int8(-1))
    return out


def _vertex_xyz(evt: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = multicol_get_series(evt, ("slc", "vertex", "x", "", "")).to_numpy(dtype=float)
    y = multicol_get_series(evt, ("slc", "vertex", "y", "", "")).to_numpy(dtype=float)
    z = multicol_get_series(evt, ("slc", "vertex", "z", "", "")).to_numpy(dtype=float)
    return x, y, z


def _empty_oct_hists(var_defs: Dict[str, dict]) -> List[Dict[str, np.ndarray]]:
    return [
        {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}
        for _ in range(N_OCT)
    ]


def _fill_hist_by_octant(
    hists_oct: List[Dict[str, np.ndarray]],
    var_name: str,
    values: np.ndarray,
    oct_ids: np.ndarray,
    bins: np.ndarray,
) -> None:
    values = np.asarray(values, dtype=float)
    oct_ids = np.asarray(oct_ids)
    if len(values) != len(oct_ids):
        raise ValueError(
            f"{var_name}: values ({len(values)}) vs oct_ids ({len(oct_ids)}) length mismatch"
        )
    for oid in range(N_OCT):
        m = oct_ids == oid
        if not np.any(m):
            continue
        hists_oct[oid][var_name] += histogram_var(values[m], bins)


def _store_keyed_by_octant(
    keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]],
    var_name: str,
    keys: Sequence[Optional[dc.PairKey]],
    values: np.ndarray,
    oct_ids: np.ndarray,
    *,
    trk_slot: Optional[int] = None,
) -> None:
    for key, val, oid in zip(keys, values, oct_ids):
        if key is None or oid < 0 or oid >= N_OCT:
            continue
        if not np.isfinite(val):
            continue
        pk: dc.PairKey = key if trk_slot is None else key + (int(trk_slot),)
        keyed_oct[int(oid)].setdefault(var_name, {})[pk] = float(val)


def process_sel_all_octants(
    df_file: str,
    *,
    hists_oct: List[Dict[str, np.ndarray]],
    keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]],
    var_defs: Dict[str, dict],
    stage_specs: Dict[str, List[Tuple[str, Any, str]]],
    collect_keyed: bool,
) -> None:
    n_split = get_n_split(df_file)
    for i in range(n_split):
        split: Dict[str, Optional[pd.DataFrame]] = {}
        for key in ("evt", "trk", "hdr"):
            try:
                split[key] = pd.read_hdf(df_file, key=f"{key}_{i}")
            except Exception:
                split[key] = None
        hdr = split.get("hdr")
        evt = split.get("evt")
        trk = split.get("trk")
        if evt is None or len(evt) == 0:
            continue

        entry_table = (
            dc._sel_all_entry_key_table(hdr, evt) if collect_keyed else None
        )
        attach_intrinsic_weights(evt, trk, "mc", use_mc_genweight=False)
        evt, _ = ensure_phi_and_kinematics_cols(evt, trk, None)
        state = {"evt": evt, "trk": trk, "hdr": hdr, "mcnu": None}

        for stage_key, cur in walk_pipeline(state, sample="mc"):
            cur_evt = cur.get("evt")
            if cur_evt is None or len(cur_evt) == 0:
                continue
            try:
                vx, vy, vz = _vertex_xyz(cur_evt)
                oct_ids = octant_ids_from_xyz(vx, vy, vz)
            except Exception:
                continue

            base_keys: Optional[List[Optional[dc.PairKey]]] = None
            if collect_keyed:
                base_keys = dc._evt_base_keys(cur_evt, entry_table)

            for var_name, vc, target in stage_specs.get(stage_key, []):
                if var_name not in var_defs:
                    continue
                bins = var_defs[var_name]["bins"]
                if target == "evt":
                    got = get_var_series(cur, vc, target)
                    if got is None:
                        continue
                    vals, _ = got
                    _fill_hist_by_octant(hists_oct, var_name, vals, oct_ids, bins)
                    if collect_keyed and base_keys is not None:
                        _store_keyed_by_octant(
                            keyed_oct, var_name, base_keys, vals, oct_ids
                        )
                elif target == "trk":
                    try:
                        v1 = multicol_get_series(
                            cur_evt.trk1, vc.var_evt_reco_col
                        ).to_numpy(dtype=float)
                        v2 = multicol_get_series(
                            cur_evt.trk2, vc.var_evt_reco_col
                        ).to_numpy(dtype=float)
                    except Exception:
                        continue
                    # Each track inherits its parent slice octant.
                    vals = np.concatenate([v1, v2])
                    oct_both = np.concatenate([oct_ids, oct_ids])
                    _fill_hist_by_octant(hists_oct, var_name, vals, oct_both, bins)
                    if collect_keyed and base_keys is not None:
                        _store_keyed_by_octant(
                            keyed_oct, var_name, base_keys, v1, oct_ids, trk_slot=0
                        )
                        _store_keyed_by_octant(
                            keyed_oct, var_name, base_keys, v2, oct_ids, trk_slot=1
                        )

        del split, state
        gc.collect()


def process_sel_mup_octants(
    df_file: str,
    *,
    hists_oct: List[Dict[str, np.ndarray]],
    keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]],
    var_defs: Dict[str, dict],
    collect_keyed: bool,
) -> None:
    n_split = get_n_split(df_file)
    for i in range(n_split):
        evt = None
        for evt_key in ("evt_cv", "evt"):
            try:
                evt = pd.read_hdf(df_file, key=f"{evt_key}_{i}")
                break
            except Exception:
                continue
        meta = None
        try:
            meta = pd.read_hdf(df_file, key=f"meta_{i}")
        except Exception:
            meta = None
        if evt is None or len(evt) == 0:
            continue

        attach_intrinsic_weights(evt, None, "mc", use_mc_genweight=False)
        evt, _ = ensure_phi_and_kinematics_cols(evt, None, None)
        try:
            vx, vy, vz = _vertex_xyz(evt)
            oct_ids = octant_ids_from_xyz(vx, vy, vz)
        except Exception:
            continue

        base_keys: Optional[List[Optional[dc.PairKey]]] = None
        if collect_keyed and meta is not None and len(meta) > 0:
            mf = dc._flatten_meta_event_cols(meta)
            mf = mf.copy()
            mf["E"] = mf["E"].astype(np.float32)
            left = pd.DataFrame(
                {
                    "__ntuple": evt.index.get_level_values("__ntuple").to_numpy(),
                    "entry": evt.index.get_level_values("entry").to_numpy(),
                    "E": multicol_get_series(
                        evt, ("mc", "E", "", "", "", "")
                    ).to_numpy(dtype=np.float32),
                    "row_i": np.arange(len(evt), dtype=np.int64),
                }
            )
            m = left.merge(mf, on=["__ntuple", "entry", "E"], how="left")
            m = m.sort_values("row_i")
            base_keys = [None] * len(evt)
            E = m["E"].to_numpy()
            run = m["run"].to_numpy()
            subrun = m["subrun"].to_numpy()
            ev = m["evt"].to_numpy()
            row_i = m["row_i"].to_numpy(dtype=np.int64)
            for j, ri in enumerate(row_i):
                if not np.isfinite(run[j]):
                    continue
                base_keys[int(ri)] = (
                    float(E[j]),
                    int(run[j]),
                    int(subrun[j]),
                    int(ev[j]),
                )

        for var_name, cfg in var_defs.items():
            if var_name == "integrated":
                continue
            try:
                vals = cfg["extract"](evt)
            except Exception:
                continue
            if len(vals) == 0:
                continue
            _fill_hist_by_octant(hists_oct, var_name, vals, oct_ids, cfg["bins"])
            if collect_keyed and base_keys is not None:
                _store_keyed_by_octant(keyed_oct, var_name, base_keys, vals, oct_ids)

        gc.collect()


def process_variation_octants(
    variation: str,
    matched_dir: str,
    filename_str: str,
    *,
    input_format: str,
    var_defs: Dict[str, dict],
    max_files: Optional[int],
    stage_specs: Optional[Dict[str, List[Tuple[str, Any, str]]]] = None,
    collect_keyed: bool = True,
) -> Tuple[List[Dict[str, np.ndarray]], List[Dict[str, Dict[dc.PairKey, float]]]]:
    files = dc.list_matched_files(matched_dir, filename_str)
    if not files:
        files = sorted(glob.glob(path.join(matched_dir, f"*{filename_str}*.df")))
        files = [f for f in files if "_matched" not in path.basename(f)]
    if max_files is not None:
        files = files[:max_files]
    print(f"[{variation}/{input_format}] {len(files)} files in {matched_dir}", flush=True)

    hists_oct = _empty_oct_hists(var_defs)
    keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]] = [{} for _ in range(N_OCT)]

    for fpath in tqdm(files, desc=f"oct {variation} {input_format}"):
        if input_format == "sel_all":
            process_sel_all_octants(
                fpath,
                hists_oct=hists_oct,
                keyed_oct=keyed_oct,
                var_defs=var_defs,
                stage_specs=stage_specs or {},
                collect_keyed=collect_keyed,
            )
        else:
            process_sel_mup_octants(
                fpath,
                hists_oct=hists_oct,
                keyed_oct=keyed_oct,
                var_defs=var_defs,
                collect_keyed=collect_keyed,
            )
    return hists_oct, keyed_oct


def _merge_keyed(
    keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]],
) -> Dict[str, Dict[dc.PairKey, float]]:
    """Union of per-octant keyed maps (later octants overwrite duplicates)."""
    out: Dict[str, Dict[dc.PairKey, float]] = {}
    for od in keyed_oct:
        for var_name, mp in od.items():
            out.setdefault(var_name, {}).update(mp)
    return out


def _frac_vals_for_octant(
    keyed_cv_oct: Dict[str, Dict[dc.PairKey, float]],
    keyed_dent_all: Dict[str, Dict[dc.PairKey, float]],
    var_name: str,
) -> np.ndarray:
    """Frac-diff for CV-octant keys; DENT value looked up globally (CV defines octant)."""
    return dc.compute_frac_diffs(
        keyed_cv_oct.get(var_name, {}),
        keyed_dent_all.get(var_name, {}),
    )


def plot_and_package_stage(
    stage: str,
    var_defs: Dict[str, dict],
    hists_cv: List[Dict[str, np.ndarray]],
    hists_dent: List[Dict[str, np.ndarray]],
    keyed_cv: List[Dict[str, Dict[dc.PairKey, float]]],
    keyed_dent: List[Dict[str, Dict[dc.PairKey, float]]],
    *,
    fig_dir: str,
    cache_path: str,
) -> dict:
    """Write per-octant plots and a recombination pickle."""
    keyed_dent_all = _merge_keyed(keyed_dent)
    payload: dict = {
        "stage": stage,
        "boundaries": {"x0": X0, "y0": Y0, "z0": Z0},
        "octant_meta": OCTANT_META,
        "var_defs": {
            k: {"label": v["label"], "bins": np.asarray(v["bins"])}
            for k, v in var_defs.items()
            if k != "integrated"
        },
        "frac_bins": dc.FRAC_DIFF_BINS,
        # hists[oct_id][variation][var] -> counts
        "hists": {},
        # frac_vals[oct_id][var] -> 1d array (concatenate to merge octants)
        # Octant assignment for frac-diff uses the *CV* reco vertex.
        "frac_vals": {},
        "frac_stats": {},
        "n_events_proxy": {},
    }

    for oid, meta in enumerate(OCTANT_META):
        slug = meta["slug"]
        odir = path.join(fig_dir, slug)
        odir_frac = path.join(odir, "fracdiff")
        makedirs(odir_frac, exist_ok=True)

        all_hists = {"cv": hists_cv[oid], "dent": hists_dent[oid]}
        payload["hists"][oid] = {
            "cv": {k: np.asarray(v) for k, v in hists_cv[oid].items()},
            "dent": {k: np.asarray(v) for k, v in hists_dent[oid].items()},
            "meta": meta,
        }
        payload["frac_vals"][oid] = {}
        payload["frac_stats"][oid] = {}
        payload["n_events_proxy"][oid] = {
            "cv": float(sum(np.sum(a) for a in hists_cv[oid].values())),
            "dent": float(sum(np.sum(a) for a in hists_dent[oid].values())),
        }

        for var_name, cfg in var_defs.items():
            if var_name == "integrated":
                continue
            if var_name not in hists_cv[oid] or var_name not in hists_dent[oid]:
                continue
            dc.plot_var_comparison(
                var_name,
                cfg,
                all_hists,
                fig_dir=odir,
                pot_scales=None,
            )

            fracs = _frac_vals_for_octant(keyed_cv[oid], keyed_dent_all, var_name)
            payload["frac_vals"][oid][var_name] = fracs
            payload["frac_stats"][oid][var_name] = {
                "n": int(len(fracs)),
                "mean": float(np.mean(fracs)) if len(fracs) else float("nan"),
                "std": float(np.std(fracs)) if len(fracs) else float("nan"),
                "median": float(np.median(fracs)) if len(fracs) else float("nan"),
            }
            cfg_oct = dict(cfg)
            cfg_oct["label"] = f"{cfg['label']}  [{meta['name']}]"
            print(
                f"[fracdiff/{stage}/{slug}] {var_name}: N={len(fracs)}",
                flush=True,
            )
            dc.plot_frac_diff(var_name, cfg_oct, fracs, fig_dir=odir_frac)

    payload["recombine_notes"] = (
        "Distribution overlays for a union of octants: sum hists[oid]['cv'|'dent'][var] "
        "over the desired oid list, then call dent_compare.plot_var_comparison. "
        "Frac-diff for a union: concatenate frac_vals[oid][var] arrays (CV-octant tagged), "
        "then dent_compare.plot_frac_diff / np.histogram with frac_bins. "
        "Quadrant examples (by ew/tb, ignoring ns): filter octant_meta ew/tb fields. "
        "N/S halves: filter ns; E/W: filter ew; Top/Bottom: filter tb."
    )

    dc.save_hists(cache_path, payload)
    return payload


def parse_args():
    p = argparse.ArgumentParser(description="DENT CV vs DENT per-octant comparison")
    p.add_argument("--out-base", default=_OUT_BASE)
    p.add_argument("--cv-all-dir", default=dc.DEFAULT_DIRS["cv_all"])
    p.add_argument("--dent-all-dir", default=dc.DEFAULT_DIRS["dent_all"])
    p.add_argument("--cv-mup-dir", default=dc.DEFAULT_DIRS["cv_mup"])
    p.add_argument("--dent-mup-dir", default=dc.DEFAULT_DIRS["dent_mup"])
    p.add_argument("--max-files", type=int, default=None)
    p.add_argument("--skip-match", action="store_true", default=True)
    p.add_argument("--do-match", action="store_true", help="Run matching first")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = path.join(args.out_base, "cache")
    fig_dir = path.join(args.out_base, "plots")
    fig_dir_mup = path.join(fig_dir, "sel_mup")
    makedirs(cache_dir, exist_ok=True)
    makedirs(fig_dir, exist_ok=True)

    if args.do_match:
        # reuse dent_compare matching
        class _A:
            pass
        a = _A()
        a.cv_all_dir = args.cv_all_dir
        a.dent_all_dir = args.dent_all_dir
        a.cv_mup_dir = args.cv_mup_dir
        a.dent_mup_dir = args.dent_mup_dir
        a.max_files = args.max_files
        a.out_base = args.out_base
        dc.run_matching(a)

    print("Octant boundaries: x=0, y=0, z=250", flush=True)
    for m in OCTANT_META:
        print(f"  {m['id']}: {m['name']:12s}  {m['title']}", flush=True)

    # ── sel_all ────────────────────────────────────────────────────────────
    sel_all_var_defs = dc.build_sel_all_var_defs()
    stage_specs = dc._stage_specs_by_key()

    hists_cv, keyed_cv = process_variation_octants(
        "cv", args.cv_all_dir, "sel_all",
        input_format="sel_all", var_defs=sel_all_var_defs,
        max_files=args.max_files, stage_specs=stage_specs, collect_keyed=True,
    )
    hists_dent, keyed_dent = process_variation_octants(
        "dent", args.dent_all_dir, "sel_all",
        input_format="sel_all", var_defs=sel_all_var_defs,
        max_files=args.max_files, stage_specs=stage_specs, collect_keyed=True,
    )
    plot_and_package_stage(
        "sel_all",
        sel_all_var_defs,
        hists_cv,
        hists_dent,
        keyed_cv,
        keyed_dent,
        fig_dir=fig_dir,
        cache_path=path.join(cache_dir, "dent_sel_all_octants.pkl"),
    )
    # Free keyed maps before sel_mup
    del keyed_cv, keyed_dent
    gc.collect()

    # ── sel_mup ────────────────────────────────────────────────────────────
    mup_var_defs = dc.build_mup_var_defs()
    hists_cv, keyed_cv = process_variation_octants(
        "cv", args.cv_mup_dir, "sel_mup",
        input_format="sel_mup", var_defs=mup_var_defs,
        max_files=args.max_files, collect_keyed=True,
    )
    hists_dent, keyed_dent = process_variation_octants(
        "dent", args.dent_mup_dir, "sel_mup",
        input_format="sel_mup", var_defs=mup_var_defs,
        max_files=args.max_files, collect_keyed=True,
    )
    plot_and_package_stage(
        "sel_mup",
        mup_var_defs,
        hists_cv,
        hists_dent,
        keyed_cv,
        keyed_dent,
        fig_dir=fig_dir_mup,
        cache_path=path.join(cache_dir, "dent_sel_mup_octants.pkl"),
    )

    print(f"\nDone. Outputs under {args.out_base}", flush=True)
    print(f"  plots:  {fig_dir}/oct*/ and {fig_dir_mup}/oct*/", flush=True)
    print(f"  cache:  {cache_dir}/dent_sel_*_octants.pkl", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
