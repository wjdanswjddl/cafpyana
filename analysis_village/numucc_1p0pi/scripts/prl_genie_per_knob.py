#!/usr/bin/env python3
"""PRL GENIE pipeline helper: chunk-map / merge / extract per-knob rate & xsec unc.

Used by ``run_prl_genie_per_knob.sh``. Modes are processed one at a time; within a
mode, ``.df`` files are mapped in parallel. Final products are **per-knob**
binned fractional uncertainties (rate and xsec), not per-mode totals.
"""
from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import os
import pickle
import re
import sys
import time
import traceback
from os import path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

_REPO = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    GENIE_GROUP_KNOBS,
    GENIE_GROUP_ORDER,
)
from analysis_village.numucc_1p0pi.final_selected_evt_vars import (  # noqa: E402
    CORE_SELECTED_EVT_VARIABLE_CONFIGS,
)
from analysis_village.numucc_1p0pi.scripts import get_systematics_genie as genie_mod  # noqa: E402
from analysis_village.numucc_1p0pi.scripts.syst_genie_aggregate import (  # noqa: E402
    GENIE_MERGE_COMBINED_KEY,
    run_genie_syst_aggregate,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402

# Variables requested for the PRL rate-vs-xsec comparison.
PRL_VAR_SAVE_NAMES: tuple[str, ...] = (
    "integrated",
    "muon-p",
    "proton-p",
    "muon-dir_z",
    "proton-dir_z",
    "tki-del_Tp",
    "tki-del_p",
    "tki-del_Tp_x",
    "tki-del_Tp_y",
    "tki-del_alpha",
    "tki-del_phi",
)

DF_GLOB_TEMPLATE = (
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_08_24_*__sel_mup-wgts_genie_{mode}"
)
MODE_RE = re.compile(r"__sel_mup-wgts_genie_([A-Za-z0-9]+)$")


def prl_var_configs() -> List[VariableConfig]:
    """``VariableConfig`` objects for the PRL variable list (order preserved)."""
    by_name = {vc.var_save_name: vc for vc in CORE_SELECTED_EVT_VARIABLE_CONFIGS}
    missing = [n for n in PRL_VAR_SAVE_NAMES if n not in by_name]
    if missing:
        raise SystemExit("missing VariableConfig for: %s" % missing)
    return [by_name[n] for n in PRL_VAR_SAVE_NAMES]


def discover_mode_dirs(
    df_root_glob: str = "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_08_24_*__sel_mup-wgts_genie_*",
    modes: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    """Return ``{mode: directory}`` for matching GENIE weight sample dirs."""
    allowed = set(modes) if modes else None
    found: Dict[str, str] = {}
    for d in sorted(glob.glob(df_root_glob)):
        if not path.isdir(d):
            continue
        m = MODE_RE.search(path.basename(d))
        if not m:
            continue
        mode = m.group(1)
        if allowed is not None and mode not in allowed:
            continue
        if mode not in GENIE_GROUP_KNOBS:
            print("[PRL] skip unknown mode dir %s (not in GENIE_GROUP_KNOBS)" % d)
            continue
        # Prefer first match; warn on duplicates.
        if mode in found:
            print("[PRL] WARNING: duplicate mode %s; keeping %s, ignoring %s" % (mode, found[mode], d))
            continue
        found[mode] = d
    # Stable order: GENIE_GROUP_ORDER first, then leftovers.
    ordered: Dict[str, str] = {}
    for mode in GENIE_GROUP_ORDER:
        if mode in found:
            ordered[mode] = found[mode]
    for mode in sorted(found):
        if mode not in ordered:
            ordered[mode] = found[mode]
    return ordered


def select_df_files(df_dir: str, fraction: float, max_files: int = 0) -> List[str]:
    files = sorted(glob.glob(path.join(df_dir, "*.df")))
    if not files:
        raise SystemExit("no .df files under %s" % df_dir)
    if fraction < 1.0:
        n_take = max(1, int(round(len(files) * fraction)))
        idx = np.linspace(0, len(files) - 1, n_take, dtype=int)
        files = [files[i] for i in idx]
    if max_files > 0:
        files = files[:max_files]
    return files


def _expected_pkl(chunks_dir: str, mode: str, df_file: str) -> str:
    stem = path.splitext(path.basename(df_file))[0]
    return path.join(chunks_dir, "genie__%s__%s.pkl" % (mode, stem))


def _worker_map(job: Dict[str, Any]) -> Dict[str, Any]:
    mode = job["mode"]
    df_file = job["df_file"]
    out_dir = job["out_dir"]
    out_pkl = _expected_pkl(out_dir, mode, df_file)
    if path.exists(out_pkl) and path.getsize(out_pkl) > 0:
        return {"ok": True, "status": "skipped", "df_file": df_file, "out_pkl": out_pkl}
    t0 = time.time()
    try:
        knobs = list(GENIE_GROUP_KNOBS[mode])
        syst_names = [("mc", k) for k in knobs]
        var_configs = prl_var_configs()
        out_path = genie_mod.run_chunk_map(
            df_file=df_file,
            out_dir=out_dir,
            genie_group=mode,
            var_configs=var_configs,
            syst_names=syst_names,
            input_stage="final",
        )
        return {
            "ok": True,
            "status": "ok",
            "df_file": df_file,
            "out_pkl": out_path,
            "elapsed": time.time() - t0,
        }
    except BaseException as e:
        return {
            "ok": False,
            "status": "failed",
            "df_file": df_file,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=8)),
            "elapsed": time.time() - t0,
        }


def run_map_phase(
    mode: str,
    df_files: Sequence[str],
    chunks_dir: str,
    workers: int,
) -> None:
    os.makedirs(chunks_dir, exist_ok=True)
    jobs = [{"mode": mode, "df_file": f, "out_dir": chunks_dir} for f in df_files]
    print("[PRL] %s chunk-map: %d files, workers=%d" % (mode, len(jobs), workers))
    if workers <= 1:
        results = [_worker_map(j) for j in jobs]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=workers) as pool:
            results = pool.map(_worker_map, jobs, chunksize=1)
    ok = sum(1 for r in results if r.get("ok"))
    fail = [r for r in results if not r.get("ok")]
    skipped = sum(1 for r in results if r.get("status") == "skipped")
    print("[PRL] %s chunk-map done: ok=%d skipped=%d failed=%d" % (mode, ok, skipped, len(fail)))
    if fail:
        fail_log = path.join(path.dirname(chunks_dir), "failed_chunk_map_%s.log" % mode)
        with open(fail_log, "w") as fh:
            for r in fail:
                fh.write("%s\n%s\n\n" % (r.get("df_file"), r.get("err", "")))
        raise SystemExit("[PRL] %s chunk-map failures → %s" % (mode, fail_log))


def frac_unc_pct(cov_frac: np.ndarray) -> np.ndarray:
    d = np.maximum(np.diag(np.asarray(cov_frac, dtype=float)), 0.0)
    return 100.0 * np.sqrt(d)


def extract_per_knob_uncertainties(
    cov_mat_dict: Dict[str, Dict[str, Any]],
    knob_to_mode: Dict[str, str],
    var_configs: Sequence[VariableConfig],
) -> Dict[str, Any]:
    """Build ``{knob: {mode, variables: {vsn: {bin_centers, rate_*, xsec_*}}}}``."""
    out: Dict[str, Any] = {
        "schema": "prl_genie_per_knob_unc_v1",
        "variables": list(PRL_VAR_SAVE_NAMES),
        "knobs": {},
    }
    vc_by = {vc.var_save_name: vc for vc in var_configs}

    # Collect knob names from cov dict (exclude aggregate keys).
    knob_names: set[str] = set()
    for vsn, row in cov_mat_dict.items():
        if not isinstance(row, dict):
            continue
        for k in row:
            if k in ("genie", "genie_rate") or k.endswith("_rate"):
                continue
            if k == GENIE_MERGE_COMBINED_KEY:
                continue
            knob_names.add(str(k))

    for kn in sorted(knob_names):
        kn_rate = "%s_rate" % kn
        entry: Dict[str, Any] = {
            "mode": knob_to_mode.get(kn, "unknown"),
            "variables": {},
        }
        for vsn in PRL_VAR_SAVE_NAMES:
            row = cov_mat_dict.get(vsn)
            if not isinstance(row, dict):
                continue
            xsec_cf = row.get(kn)
            rate_cf = row.get(kn_rate)
            if xsec_cf is None and rate_cf is None:
                continue
            vc = vc_by[vsn]
            centers = np.asarray(vc.bin_centers, dtype=float)
            var_entry: Dict[str, Any] = {
                "bin_centers": centers.tolist(),
                "bin_edges": np.asarray(vc.bins, dtype=float).tolist(),
            }
            if rate_cf is not None:
                rate_cf = np.asarray(rate_cf, dtype=float)
                var_entry["rate_cov_frac"] = rate_cf.tolist()
                var_entry["rate_frac_unc_pct"] = frac_unc_pct(rate_cf).tolist()
            if xsec_cf is not None:
                xsec_cf = np.asarray(xsec_cf, dtype=float)
                var_entry["xsec_cov_frac"] = xsec_cf.tolist()
                var_entry["xsec_frac_unc_pct"] = frac_unc_pct(xsec_cf).tolist()
            entry["variables"][vsn] = var_entry
        if entry["variables"]:
            out["knobs"][kn] = entry
    return out


def write_per_knob_outputs(payload: Dict[str, Any], out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    json_path = path.join(out_dir, "per_knob_rate_xsec_unc.json")
    with open(json_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print("[PRL] wrote", json_path)

    # Compact NPZ: one array per (knob, var, kind).
    npz_kw: Dict[str, np.ndarray] = {}
    for kn, kn_entry in payload["knobs"].items():
        safe_kn = kn.replace("/", "_")
        for vsn, vent in kn_entry["variables"].items():
            prefix = "%s__%s" % (safe_kn, vsn)
            npz_kw["%s__bin_centers" % prefix] = np.asarray(vent["bin_centers"], dtype=float)
            if "rate_frac_unc_pct" in vent:
                npz_kw["%s__rate_frac_unc_pct" % prefix] = np.asarray(
                    vent["rate_frac_unc_pct"], dtype=float
                )
            if "xsec_frac_unc_pct" in vent:
                npz_kw["%s__xsec_frac_unc_pct" % prefix] = np.asarray(
                    vent["xsec_frac_unc_pct"], dtype=float
                )
            if "rate_cov_frac" in vent:
                npz_kw["%s__rate_cov_frac" % prefix] = np.asarray(
                    vent["rate_cov_frac"], dtype=float
                )
            if "xsec_cov_frac" in vent:
                npz_kw["%s__xsec_cov_frac" % prefix] = np.asarray(
                    vent["xsec_cov_frac"], dtype=float
                )
    npz_path = path.join(out_dir, "per_knob_rate_xsec_unc.npz")
    np.savez_compressed(npz_path, **npz_kw)
    print("[PRL] wrote", npz_path)

    # Flat CSV-like summary for quick inspection.
    csv_path = path.join(out_dir, "per_knob_rate_xsec_unc_summary.csv")
    with open(csv_path, "w") as fh:
        fh.write("knob,mode,var_save_name,bin,bin_center,rate_frac_unc_pct,xsec_frac_unc_pct\n")
        for kn, kn_entry in payload["knobs"].items():
            mode = kn_entry["mode"]
            for vsn, vent in kn_entry["variables"].items():
                centers = vent["bin_centers"]
                rate = vent.get("rate_frac_unc_pct") or [float("nan")] * len(centers)
                xsec = vent.get("xsec_frac_unc_pct") or [float("nan")] * len(centers)
                for i, (c, r, x) in enumerate(zip(centers, rate, xsec)):
                    fh.write("%s,%s,%s,%d,%.6g,%.6g,%.6g\n" % (kn, mode, vsn, i, c, r, x))
    print("[PRL] wrote", csv_path)


def run_mode(
    mode: str,
    df_dir: str,
    work_dir: str,
    *,
    workers: int,
    fraction: float,
    max_files: int,
    xsec_unit: float,
    skip_map: bool,
) -> None:
    chunks_dir = path.join(work_dir, "chunks", mode)
    merge_dir = path.join(work_dir, "merged", mode)
    os.makedirs(chunks_dir, exist_ok=True)
    os.makedirs(merge_dir, exist_ok=True)

    df_files = select_df_files(df_dir, fraction=fraction, max_files=max_files)
    man = {
        "mode": mode,
        "df_dir": df_dir,
        "n_total": len(sorted(glob.glob(path.join(df_dir, "*.df")))),
        "n_used": len(df_files),
        "fraction": fraction,
        "max_files": max_files,
        "knobs": list(GENIE_GROUP_KNOBS[mode]),
        "df_files": df_files,
    }
    with open(path.join(work_dir, "manifest_%s.json" % mode), "w") as fh:
        json.dump(man, fh, indent=2)
    print("[PRL] %s: using %d / %d files" % (mode, man["n_used"], man["n_total"]))

    if not skip_map:
        run_map_phase(mode, df_files, chunks_dir, workers)

    var_configs = prl_var_configs()
    print("[PRL] %s chunk-merge …" % mode)
    genie_mod.run_chunk_merge(
        chunks_dir,
        merge_dir,
        mode,
        var_configs,
        xsec_unit=float(xsec_unit),
        bkgd_subtract=True,
        save_figs=False,
        npz_path=path.join(merge_dir, "genie_syst_%s.npz" % mode),
    )


def run_aggregate_and_extract(
    work_dir: str,
    modes: Sequence[str],
    *,
    xsec_unit: float,
) -> str:
    chunks_root = path.join(work_dir, "chunks")
    # Aggregate expects all genie__MODE__*.pkl under one chunks dir. Our layout is
    # chunks/<MODE>/genie__MODE__*.pkl — flatten via symlinks into a staging dir.
    staging = path.join(work_dir, "chunks_flat")
    os.makedirs(staging, exist_ok=True)
    for mode in modes:
        for pkl in glob.glob(path.join(chunks_root, mode, "genie__%s__*.pkl" % mode)):
            dest = path.join(staging, path.basename(pkl))
            if path.lexists(dest):
                os.remove(dest)
            os.symlink(path.abspath(pkl), dest)

    syst_disk = path.join(work_dir, "syst_disk")
    pkl_path = run_genie_syst_aggregate(
        staging,
        syst_disk,
        mc_df_stage="final",
        xsec_unit=float(xsec_unit),
        genie_groups=list(modes),
    )

    with open(pkl_path, "rb") as fh:
        cov_mat_dict = pickle.load(fh)

    knob_to_mode: Dict[str, str] = {}
    for mode in modes:
        for kn in GENIE_GROUP_KNOBS.get(mode, []):
            knob_to_mode[kn] = mode

    payload = extract_per_knob_uncertainties(cov_mat_dict, knob_to_mode, prl_var_configs())
    write_per_knob_outputs(payload, path.join(work_dir, "per_knob"))
    return pkl_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pd = sub.add_parser("discover", help="Print mode→directory mapping as JSON")
    pd.add_argument("--df-glob", default="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_08_24_*__sel_mup-wgts_genie_*")
    pd.add_argument("--modes", default="", help="Comma-separated mode filter")

    pm = sub.add_parser("run-mode", help="Chunk-map + merge one GENIE mode")
    pm.add_argument("--mode", required=True)
    pm.add_argument("--df-dir", required=True)
    pm.add_argument("--work-dir", required=True)
    pm.add_argument("--workers", type=int, default=8)
    pm.add_argument("--fraction", type=float, default=1.0)
    pm.add_argument("--max-files", type=int, default=0)
    pm.add_argument("--xsec-unit", type=float, default=1.0)
    pm.add_argument("--skip-map", action="store_true")

    pa = sub.add_parser("aggregate", help="Aggregate all modes + write per-knob unc")
    pa.add_argument("--work-dir", required=True)
    pa.add_argument("--modes", required=True, help="Comma-separated modes present under chunks/")
    pa.add_argument("--xsec-unit", type=float, default=1.0)

    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.cmd == "discover":
        modes = [s.strip() for s in args.modes.split(",") if s.strip()] or None
        found = discover_mode_dirs(args.df_glob, modes=modes)
        print(json.dumps(found, indent=2))
        return 0
    if args.cmd == "run-mode":
        run_mode(
            args.mode,
            args.df_dir,
            args.work_dir,
            workers=args.workers,
            fraction=args.fraction,
            max_files=args.max_files,
            xsec_unit=args.xsec_unit,
            skip_map=args.skip_map,
        )
        return 0
    if args.cmd == "aggregate":
        modes = [s.strip() for s in args.modes.split(",") if s.strip()]
        run_aggregate_and_extract(args.work_dir, modes, xsec_unit=args.xsec_unit)
        return 0
    raise SystemExit("unknown cmd")


if __name__ == "__main__":
    raise SystemExit(main())
