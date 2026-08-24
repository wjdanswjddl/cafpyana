#!/usr/bin/env python3
"""
DENT vs CV by detector region: full sample, octants, and x=0 halves.

Regions (slice vertex):
  * all      — no spatial cut
  * octants  — x=0, y=0, z=250 (8 octants; SBND E/W–N/S–Top/Bottom)
  * halves   — x=0 only (E: x<0, W: x>=0)

For each region writes:
  * distribution overlays (diagnostic_*)
  * per-event (DENT-CV)/CV  (fracdiff/)
  * per-event DENT/CV       (ratio/)

Caches under ``cache/`` store histogram counts **and** paired CV/DENT arrays
so plots can be remade later without re-reading matched dfs
(see ``dent_replot_regions.py``).

Example:
    python dent_compare_regions.py --max-files 3 \\
        --out-base .../systematics-final/DENT-regions
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
    "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/DENT-regions"
)

X0, Y0, Z0 = 0.0, 0.0, 250.0

OCTANT_META: List[dict] = []
for _i in range(8):
    _x_ge = bool((_i >> 2) & 1)
    _y_ge = bool((_i >> 1) & 1)
    _z_ge = bool((_i >> 0) & 1)
    _ew = "W" if _x_ge else "E"
    _ns = "N" if _z_ge else "S"
    _tb = "Top" if _y_ge else "Bottom"
    OCTANT_META.append(
        {
            "id": _i,
            "kind": "octant",
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

HALF_META = [
    {
        "id": 0,
        "kind": "half",
        "name": "E",
        "slug": "half_E",
        "title": "x<0",
        "x_ge0": False,
        "ew": "E",
    },
    {
        "id": 1,
        "kind": "half",
        "name": "W",
        "slug": "half_W",
        "title": "x≥0",
        "x_ge0": True,
        "ew": "W",
    },
]
N_HALF = 2


def octant_ids_from_xyz(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    x0: float = X0,
    y0: float = Y0,
    z0: float = Z0,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    out = (
        ((x >= x0).astype(np.int8) << 2)
        | ((y >= y0).astype(np.int8) << 1)
        | (z >= z0).astype(np.int8)
    ).astype(np.int8)
    return np.where(ok, out, np.int8(-1))


def half_ids_from_x(x: np.ndarray, *, x0: float = X0) -> np.ndarray:
    """0 = E (x<0), 1 = W (x>=0); -1 if non-finite."""
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    out = np.where(x >= x0, np.int8(1), np.int8(0))
    return np.where(ok, out, np.int8(-1))


def _vertex_xyz(evt: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = multicol_get_series(evt, ("slc", "vertex", "x", "", "")).to_numpy(dtype=float)
    y = multicol_get_series(evt, ("slc", "vertex", "y", "", "")).to_numpy(dtype=float)
    z = multicol_get_series(evt, ("slc", "vertex", "z", "", "")).to_numpy(dtype=float)
    return x, y, z


def _empty_hists(var_defs: Dict[str, dict]) -> Dict[str, np.ndarray]:
    return {v: np.zeros(len(cfg["bins"]) - 1, dtype=float) for v, cfg in var_defs.items()}


def _fill_hist(
    hists: Dict[str, np.ndarray],
    var_name: str,
    values: np.ndarray,
    bins: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> None:
    values = np.asarray(values, dtype=float)
    if mask is not None:
        values = values[mask]
    if len(values) == 0:
        return
    hists[var_name] += histogram_var(values, bins)


def _store_keyed(
    keyed: Dict[str, Dict[dc.PairKey, float]],
    var_name: str,
    keys: Sequence[Optional[dc.PairKey]],
    values: np.ndarray,
    mask: Optional[np.ndarray] = None,
    *,
    trk_slot: Optional[int] = None,
) -> None:
    values = np.asarray(values, dtype=float)
    for i, (key, val) in enumerate(zip(keys, values)):
        if mask is not None and not mask[i]:
            continue
        if key is None or not np.isfinite(val):
            continue
        pk: dc.PairKey = key if trk_slot is None else key + (int(trk_slot),)
        keyed.setdefault(var_name, {})[pk] = float(val)


class RegionAccum:
    """Histograms + keyed maps for all / octants / halves."""

    def __init__(self, var_defs: Dict[str, dict]):
        self.var_defs = var_defs
        self.hists_all = _empty_hists(var_defs)
        self.hists_oct = [_empty_hists(var_defs) for _ in range(N_OCT)]
        self.hists_half = [_empty_hists(var_defs) for _ in range(N_HALF)]
        self.keyed_all: Dict[str, Dict[dc.PairKey, float]] = {}
        self.keyed_oct: List[Dict[str, Dict[dc.PairKey, float]]] = [{} for _ in range(N_OCT)]
        self.keyed_half: List[Dict[str, Dict[dc.PairKey, float]]] = [{} for _ in range(N_HALF)]

    def fill_evt_var(
        self,
        var_name: str,
        values: np.ndarray,
        oct_ids: np.ndarray,
        half_ids: np.ndarray,
        keys: Optional[Sequence[Optional[dc.PairKey]]],
        *,
        collect_keyed: bool,
        trk_slot: Optional[int] = None,
    ) -> None:
        bins = self.var_defs[var_name]["bins"]
        values = np.asarray(values, dtype=float)
        _fill_hist(self.hists_all, var_name, values, bins)
        if collect_keyed and keys is not None:
            _store_keyed(self.keyed_all, var_name, keys, values, trk_slot=trk_slot)

        for oid in range(N_OCT):
            m = oct_ids == oid
            if not np.any(m):
                continue
            _fill_hist(self.hists_oct[oid], var_name, values, bins, m)
            if collect_keyed and keys is not None:
                _store_keyed(
                    self.keyed_oct[oid], var_name, keys, values, m, trk_slot=trk_slot
                )

        for hid in range(N_HALF):
            m = half_ids == hid
            if not np.any(m):
                continue
            _fill_hist(self.hists_half[hid], var_name, values, bins, m)
            if collect_keyed and keys is not None:
                _store_keyed(
                    self.keyed_half[hid], var_name, keys, values, m, trk_slot=trk_slot
                )


def _process_sel_all_file(
    df_file: str,
    *,
    accum: RegionAccum,
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

        entry_table = dc._sel_all_entry_key_table(hdr, evt) if collect_keyed else None
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
                half_ids = half_ids_from_x(vx)
            except Exception:
                n = len(cur_evt)
                oct_ids = np.full(n, -1, dtype=np.int8)
                half_ids = np.full(n, -1, dtype=np.int8)

            base_keys: Optional[List[Optional[dc.PairKey]]] = None
            if collect_keyed:
                base_keys = dc._evt_base_keys(cur_evt, entry_table)

            for var_name, vc, target in stage_specs.get(stage_key, []):
                if var_name not in accum.var_defs:
                    continue
                if target == "evt":
                    got = get_var_series(cur, vc, target)
                    if got is None:
                        continue
                    vals, _ = got
                    accum.fill_evt_var(
                        var_name, vals, oct_ids, half_ids, base_keys,
                        collect_keyed=collect_keyed,
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
                    oct_both = np.concatenate([oct_ids, oct_ids])
                    half_both = np.concatenate([half_ids, half_ids])
                    # Fill each track with its slot in the pair key.
                    if collect_keyed and base_keys is not None:
                        accum.fill_evt_var(
                            var_name, v1, oct_ids, half_ids, base_keys,
                            collect_keyed=True, trk_slot=0,
                        )
                        accum.fill_evt_var(
                            var_name, v2, oct_ids, half_ids, base_keys,
                            collect_keyed=True, trk_slot=1,
                        )
                    else:
                        vals = np.concatenate([v1, v2])
                        accum.fill_evt_var(
                            var_name, vals, oct_both, half_both, None,
                            collect_keyed=False,
                        )

        del split, state
        gc.collect()


def _process_sel_mup_file(
    df_file: str,
    *,
    accum: RegionAccum,
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
            half_ids = half_ids_from_x(vx)
        except Exception:
            n = len(evt)
            oct_ids = np.full(n, -1, dtype=np.int8)
            half_ids = np.full(n, -1, dtype=np.int8)

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

        for var_name, cfg in accum.var_defs.items():
            if var_name == "integrated":
                continue
            try:
                vals = cfg["extract"](evt)
            except Exception:
                continue
            if len(vals) == 0:
                continue
            accum.fill_evt_var(
                var_name, vals, oct_ids, half_ids, base_keys,
                collect_keyed=collect_keyed,
            )
        gc.collect()


def process_variation_regions(
    variation: str,
    matched_dir: str,
    filename_str: str,
    *,
    input_format: str,
    var_defs: Dict[str, dict],
    max_files: Optional[int],
    stage_specs: Optional[Dict[str, List[Tuple[str, Any, str]]]] = None,
    collect_keyed: bool = True,
) -> RegionAccum:
    files = dc.list_matched_files(matched_dir, filename_str)
    if not files:
        files = sorted(glob.glob(path.join(matched_dir, f"*{filename_str}*.df")))
        files = [f for f in files if "_matched" not in path.basename(f)]
    if max_files is not None:
        files = files[:max_files]
    print(f"[{variation}/{input_format}] {len(files)} files", flush=True)

    accum = RegionAccum(var_defs)
    for fpath in tqdm(files, desc=f"regions {variation} {input_format}"):
        if input_format == "sel_all":
            _process_sel_all_file(
                fpath,
                accum=accum,
                stage_specs=stage_specs or {},
                collect_keyed=collect_keyed,
            )
        else:
            _process_sel_mup_file(
                fpath, accum=accum, collect_keyed=collect_keyed,
            )
    return accum


def _plot_one_region(
    region_key: str,
    meta: dict,
    var_defs: Dict[str, dict],
    hists_cv: Dict[str, np.ndarray],
    hists_dent: Dict[str, np.ndarray],
    keyed_cv: Dict[str, Dict[dc.PairKey, float]],
    keyed_dent_all: Dict[str, Dict[dc.PairKey, float]],
    *,
    fig_dir: str,
) -> dict:
    """Overlays + fracdiff + ratio; return cache blob for this region."""
    odir = path.join(fig_dir, meta["slug"]) if meta.get("slug") != "all" else fig_dir
    odir_frac = path.join(odir, "fracdiff")
    odir_ratio = path.join(odir, "ratio")
    makedirs(odir_frac, exist_ok=True)
    makedirs(odir_ratio, exist_ok=True)

    all_hists = {"cv": hists_cv, "dent": hists_dent}
    blob: dict = {
        "meta": meta,
        "hists": {
            "cv": {k: np.asarray(v) for k, v in hists_cv.items()},
            "dent": {k: np.asarray(v) for k, v in hists_dent.items()},
        },
        "paired": {},
        "vars": {},
    }

    label_tag = meta.get("name", region_key)
    for var_name, cfg in var_defs.items():
        if var_name == "integrated":
            continue
        if var_name not in hists_cv or var_name not in hists_dent:
            continue
        dc.plot_var_comparison(
            var_name, cfg, all_hists, fig_dir=odir, pot_scales=None,
        )

        cv_a, dent_a = dc.compute_paired_values(
            keyed_cv.get(var_name, {}),
            keyed_dent_all.get(var_name, {}),
        )
        fracs = (dent_a - cv_a) / cv_a if len(cv_a) else np.array([], dtype=float)
        ratios = dc.compute_ratios_from_paired(cv_a, dent_a)
        cfg_r = dict(cfg)
        cfg_r["label"] = f"{cfg['label']}  [{label_tag}]"
        print(
            f"[paired/{meta.get('slug', region_key)}] {var_name}: N={len(cv_a)}",
            flush=True,
        )
        hist_f = dc.plot_frac_diff(var_name, cfg_r, fracs, fig_dir=odir_frac)
        hist_r = dc.plot_ratio_dist(var_name, cfg_r, ratios, fig_dir=odir_ratio)
        blob["paired"][var_name] = {
            "cv": cv_a,
            "dent": dent_a,
            "frac": fracs,
            "ratio": ratios,
        }
        blob["vars"][var_name] = {
            "label": cfg.get("label", var_name),
            "n_paired": int(len(cv_a)),
            "frac_stats": dc._stats_dict(fracs),
            "ratio_stats": dc._stats_dict(ratios),
            "hist_frac": hist_f,
            "hist_ratio": hist_r,
        }
    return blob


def plot_and_package_stage(
    stage: str,
    var_defs: Dict[str, dict],
    accum_cv: RegionAccum,
    accum_dent: RegionAccum,
    *,
    fig_dir: str,
    cache_path: str,
) -> dict:
    keyed_dent_all = accum_dent.keyed_all
    payload: dict = {
        "stage": stage,
        "boundaries": {"x0": X0, "y0": Y0, "z0": Z0},
        "octant_meta": OCTANT_META,
        "half_meta": HALF_META,
        "frac_bins": dc.FRAC_DIFF_BINS,
        "ratio_bins": dc.RATIO_BINS,
        "cv_abs_eps": dc.CV_ABS_EPS,
        "var_defs": {
            k: {"label": v["label"], "bins": np.asarray(v["bins"])}
            for k, v in var_defs.items()
            if k != "integrated"
        },
        "regions": {},
        "recombine_notes": (
            "Each regions[key] has hists[cv|dent][var] and paired[var]={cv,dent,frac,ratio}. "
            "Union overlays: sum hists. Union frac/ratio: concatenate paired arrays "
            "(or recompute from paired cv/dent). "
            "E half ↔ octants with ew=='E'; W ↔ ew=='W'. "
            "Replot helper: dent_replot_regions.py"
        ),
    }

    # Full sample
    print(f"\n=== [{stage}] region=all ===", flush=True)
    payload["regions"]["all"] = _plot_one_region(
        "all",
        {"id": -1, "kind": "all", "name": "all", "slug": "all", "title": "full sample"},
        var_defs,
        accum_cv.hists_all,
        accum_dent.hists_all,
        accum_cv.keyed_all,
        keyed_dent_all,
        fig_dir=path.join(fig_dir, "all"),
    )

    # Octants
    for oid, meta in enumerate(OCTANT_META):
        print(f"\n=== [{stage}] region={meta['slug']} ===", flush=True)
        payload["regions"][meta["slug"]] = _plot_one_region(
            meta["slug"],
            meta,
            var_defs,
            accum_cv.hists_oct[oid],
            accum_dent.hists_oct[oid],
            accum_cv.keyed_oct[oid],
            keyed_dent_all,
            fig_dir=path.join(fig_dir, "octants"),
        )

    # Halves
    for hid, meta in enumerate(HALF_META):
        print(f"\n=== [{stage}] region={meta['slug']} ===", flush=True)
        payload["regions"][meta["slug"]] = _plot_one_region(
            meta["slug"],
            meta,
            var_defs,
            accum_cv.hists_half[hid],
            accum_dent.hists_half[hid],
            accum_cv.keyed_half[hid],
            keyed_dent_all,
            fig_dir=path.join(fig_dir, "halves"),
        )

    dc.save_hists(cache_path, payload)
    return payload


def parse_args():
    p = argparse.ArgumentParser(
        description="DENT CV vs DENT: full + octants + x=0 halves"
    )
    p.add_argument("--out-base", default=_OUT_BASE)
    p.add_argument("--cv-all-dir", default=dc.DEFAULT_DIRS["cv_all"])
    p.add_argument("--dent-all-dir", default=dc.DEFAULT_DIRS["dent_all"])
    p.add_argument("--cv-mup-dir", default=dc.DEFAULT_DIRS["cv_mup"])
    p.add_argument("--dent-mup-dir", default=dc.DEFAULT_DIRS["dent_mup"])
    p.add_argument("--max-files", type=int, default=None)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = path.join(args.out_base, "cache")
    fig_dir = path.join(args.out_base, "plots")
    fig_dir_mup = path.join(fig_dir, "sel_mup")
    makedirs(cache_dir, exist_ok=True)

    print("Boundaries: x=0, y=0, z=250 (octants); x=0 (halves E/W)", flush=True)
    for m in OCTANT_META:
        print(f"  oct {m['id']}: {m['name']:12s}  {m['title']}", flush=True)
    for m in HALF_META:
        print(f"  half: {m['name']:12s}  {m['title']}", flush=True)

    # ── sel_all ────────────────────────────────────────────────────────────
    sel_all_var_defs = dc.build_sel_all_var_defs()
    stage_specs = dc._stage_specs_by_key()
    accum_cv = process_variation_regions(
        "cv", args.cv_all_dir, "sel_all",
        input_format="sel_all", var_defs=sel_all_var_defs,
        max_files=args.max_files, stage_specs=stage_specs, collect_keyed=True,
    )
    accum_dent = process_variation_regions(
        "dent", args.dent_all_dir, "sel_all",
        input_format="sel_all", var_defs=sel_all_var_defs,
        max_files=args.max_files, stage_specs=stage_specs, collect_keyed=True,
    )
    plot_and_package_stage(
        "sel_all",
        sel_all_var_defs,
        accum_cv,
        accum_dent,
        fig_dir=fig_dir,
        cache_path=path.join(cache_dir, "dent_sel_all_regions.pkl"),
    )
    del accum_cv, accum_dent
    gc.collect()

    # ── sel_mup ────────────────────────────────────────────────────────────
    mup_var_defs = dc.build_mup_var_defs()
    accum_cv = process_variation_regions(
        "cv", args.cv_mup_dir, "sel_mup",
        input_format="sel_mup", var_defs=mup_var_defs,
        max_files=args.max_files, collect_keyed=True,
    )
    accum_dent = process_variation_regions(
        "dent", args.dent_mup_dir, "sel_mup",
        input_format="sel_mup", var_defs=mup_var_defs,
        max_files=args.max_files, collect_keyed=True,
    )
    plot_and_package_stage(
        "sel_mup",
        mup_var_defs,
        accum_cv,
        accum_dent,
        fig_dir=fig_dir_mup,
        cache_path=path.join(cache_dir, "dent_sel_mup_regions.pkl"),
    )

    print(f"\nDone. Outputs under {args.out_base}", flush=True)
    print(f"  plots/all/, plots/octants/, plots/halves/  (+ sel_mup/…)", flush=True)
    print(f"  cache/dent_sel_*_regions.pkl  (hists + paired cv/dent/frac/ratio)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
