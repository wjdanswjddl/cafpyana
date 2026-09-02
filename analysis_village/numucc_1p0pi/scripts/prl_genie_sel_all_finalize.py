#!/usr/bin/env python3
"""Finalize sel_all GENIE syst: per-mode NPZ merges + combined cov_mat_dict + per-knob unc.

After the chunked map phase (``run_prl_genie_sel_all.sh``) finishes for each mode:

1. Wait until every group's chunk pickles are present (optional).
2. ``chunk-merge`` each group → ``merged/<GROUP>/genie_syst_<GROUP>.npz`` (per-knob cov).
3. Fold **all** groups into one ``syst_disk/GENIE/cov_mat_dict.pkl``.
4. Write ``per_knob/per_knob_all_vars.json`` (+ compact NPZ) with rate/xsec fractional
   uncertainties for every sel_all variable slug.

Used by ``run_prl_genie_sel_all_finalize.sh`` (typically inside tmux).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import subprocess
import sys
import time
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
    iter_genie_group_df_paths,
)
from analysis_village.numucc_1p0pi.scripts import get_systematics_genie as genie_mod  # noqa: E402
from analysis_village.numucc_1p0pi.scripts.get_systematics_genie import (  # noqa: E402
    GENIE_MERGE_COMBINED_KEY,
    genie_all_var_configs,
    run_chunk_merge,
)
from analysis_village.numucc_1p0pi.scripts.syst_genie_aggregate import (  # noqa: E402
    run_genie_syst_aggregate,
)
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig  # noqa: E402

DEFAULT_WORK = "/exp/sbnd/data/users/munjung/PRL_data/PRL_genie_sel_all"
DEFAULT_GROUPS = ("slim", "CCQE", "MEC", "RES", "nonRES", "Other")


def frac_unc_pct(cov_frac: np.ndarray) -> np.ndarray:
    c = np.asarray(cov_frac, dtype=np.float64)
    return 100.0 * np.sqrt(np.maximum(np.diag(c), 0.0))


def _count_chunks(chunks_dir: str, group: str) -> int:
    return len(glob.glob(path.join(chunks_dir, "genie__%s__*.pkl" % group)))


def _expected_df_count(group: str) -> int:
    return sum(1 for _ in iter_genie_group_df_paths(group, mc_df_stage="sel_all"))


def _pipeline_running(work_base: str) -> bool:
    needle = path.abspath(work_base)
    try:
        out = subprocess.check_output(["pgrep", "-af", "syst_genie_parallel.py"], text=True)
    except subprocess.CalledProcessError:
        return False
    for line in out.splitlines():
        if needle in line:
            return True
    return False


def wait_for_chunks(
    chunks_dir: str,
    groups: Sequence[str],
    *,
    poll_sec: int,
    work_base: str,
) -> None:
    """Block until each group has one chunk pickle per input .df and map workers exit."""
    print("[finalize] waiting for chunk-map to finish …", flush=True)
    while True:
        counts = {g: _count_chunks(chunks_dir, g) for g in groups}
        expected = {g: _expected_df_count(g) for g in groups}
        running = _pipeline_running(work_base)
        parts = ["%s=%d/%d" % (g, counts[g], expected[g]) for g in groups]
        print(
            "[finalize] %s  parallel_running=%s"
            % ("  ".join(parts), running),
            flush=True,
        )
        ready = all(counts[g] >= expected[g] for g in groups)
        if ready and not running:
            print("[finalize] all chunk pickles present", flush=True)
            return
        time.sleep(poll_sec)


def merge_group_npz(
    chunks_dir: str,
    merge_root: str,
    group: str,
    var_configs: List[VariableConfig],
    xsec_unit: float,
    *,
    force: bool,
) -> str:
    out_dir = path.join(merge_root, group)
    os.makedirs(out_dir, exist_ok=True)
    npz_path = path.join(out_dir, "genie_syst_%s.npz" % group)
    if path.isfile(npz_path) and not force:
        print("[finalize] skip merge group=%s (exists %s)" % (group, npz_path), flush=True)
        return npz_path
    n_chunks = _count_chunks(chunks_dir, group)
    if n_chunks == 0:
        raise SystemExit("[finalize] no chunks for group %s under %s" % (group, chunks_dir))
    print("[finalize] chunk-merge group=%s (%d pickles) …" % (group, n_chunks), flush=True)
    run_chunk_merge(
        chunks_dir,
        out_dir,
        group,
        var_configs,
        xsec_unit,
        bkgd_subtract=True,
        save_figs=False,
        npz_path=npz_path,
    )
    return npz_path


def knob_to_mode_map(groups: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for grp in groups:
        for kn in GENIE_GROUP_KNOBS.get(grp, []):
            out[kn] = grp
    return out


def extract_per_knob_uncertainties(
    cov_mat_dict: Dict[str, Dict[str, Any]],
    knob_to_mode: Dict[str, str],
    var_configs: Sequence[VariableConfig],
) -> Dict[str, Any]:
    """Per-knob rate/xsec fractional cov + diag unc for all sel_all variables."""
    var_slugs = [vc.var_save_name for vc in var_configs]
    out: Dict[str, Any] = {
        "schema": "prl_genie_sel_all_per_knob_unc_v1",
        "input_stage": "sel_all",
        "variables": var_slugs,
        "knobs": {},
    }

    knob_names: set[str] = set()
    for _vsn, row in cov_mat_dict.items():
        if not isinstance(row, dict):
            continue
        for k in row:
            if k in ("genie", "genie_rate") or k == GENIE_MERGE_COMBINED_KEY:
                continue
            if str(k).endswith("_rate"):
                continue
            knob_names.add(str(k))

    vc_by = {vc.var_save_name: vc for vc in var_configs}

    for kn in sorted(knob_names):
        kn_rate = "%s_rate" % kn
        entry: Dict[str, Any] = {
            "mode": knob_to_mode.get(kn, "unknown"),
            "variables": {},
        }
        for vsn in var_slugs:
            row = cov_mat_dict.get(vsn)
            if not isinstance(row, dict):
                continue
            xsec_cf = row.get(kn)
            rate_cf = row.get(kn_rate)
            if xsec_cf is None and rate_cf is None:
                continue
            vc = vc_by[vsn]
            var_entry: Dict[str, Any] = {
                "bin_centers": np.asarray(vc.bin_centers, dtype=float).tolist(),
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
    json_path = path.join(out_dir, "per_knob_all_vars.json")
    with open(json_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print("[finalize] wrote", json_path, flush=True)

    npz_kw: Dict[str, np.ndarray] = {}
    for kn, kn_entry in payload["knobs"].items():
        safe_kn = kn.replace("/", "_")
        for vsn, vent in kn_entry["variables"].items():
            prefix = "%s__%s" % (safe_kn, vsn)
            npz_kw["%s__bin_centers" % prefix] = np.asarray(vent["bin_centers"], dtype=float)
            for kind in ("rate_frac_unc_pct", "xsec_frac_unc_pct", "rate_cov_frac", "xsec_cov_frac"):
                if kind in vent:
                    npz_kw["%s__%s" % (prefix, kind)] = np.asarray(vent[kind], dtype=float)
    npz_path = path.join(out_dir, "per_knob_all_vars.npz")
    np.savez_compressed(npz_path, **npz_kw)
    print("[finalize] wrote", npz_path, flush=True)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-dir", default=DEFAULT_WORK)
    p.add_argument(
        "--groups",
        default=",".join(DEFAULT_GROUPS),
        help="Comma-separated GENIE groups to merge and fold (default: slim + interaction modes).",
    )
    p.add_argument("--xsec-unit", type=float, default=1.0)
    p.add_argument("--poll-sec", type=int, default=300, help="Wait poll interval when --wait.")
    p.add_argument("--wait", action="store_true", default=True, help="Wait for chunk map (default).")
    p.add_argument("--no-wait", action="store_false", dest="wait")
    p.add_argument("--force-merge", action="store_true", help="Re-run chunk-merge NPZ even if present.")
    p.add_argument("--skip-npz", action="store_true", help="Skip per-group NPZ; only aggregate + extract.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    work = path.abspath(path.expanduser(args.work_dir))
    chunks_dir = path.join(work, "chunks")
    merge_root = path.join(work, "merged")
    syst_disk = path.join(work, "syst_disk")
    per_knob_dir = path.join(work, "per_knob")

    raw_groups = [s.strip() for s in args.groups.split(",") if s.strip()]
    groups: List[str] = []
    for g in GENIE_GROUP_ORDER:
        if g in raw_groups:
            groups.append(g)
    for g in raw_groups:
        if g not in groups:
            groups.append(g)

    var_configs = genie_all_var_configs("sel_all")
    print("[finalize] work=%s groups=%s n_vars=%d" % (work, groups, len(var_configs)), flush=True)

    if args.wait:
        wait_for_chunks(chunks_dir, groups, poll_sec=args.poll_sec, work_base=work)

    if not args.skip_npz:
        for grp in groups:
            merge_group_npz(
                chunks_dir,
                merge_root,
                grp,
                var_configs,
                args.xsec_unit,
                force=args.force_merge,
            )

    pkl_path = run_genie_syst_aggregate(
        chunks_dir,
        syst_disk,
        mc_df_stage="sel_all",
        xsec_unit=float(args.xsec_unit),
        genie_groups=groups,
    )

    with open(pkl_path, "rb") as fh:
        cov_mat_dict = pickle.load(fh)

    payload = extract_per_knob_uncertainties(
        cov_mat_dict, knob_to_mode_map(groups), var_configs
    )
    write_per_knob_outputs(payload, per_knob_dir)

    summary = {
        "schema": "prl_genie_sel_all_finalize_v1",
        "work_dir": work,
        "groups": groups,
        "n_variables": len(var_configs),
        "n_knobs": len(payload["knobs"]),
        "cov_mat_dict_pkl": pkl_path,
        "per_knob_json": path.join(per_knob_dir, "per_knob_all_vars.json"),
        "per_knob_npz": path.join(per_knob_dir, "per_knob_all_vars.npz"),
        "merged_npz": {
            g: path.join(merge_root, g, "genie_syst_%s.npz" % g) for g in groups
        },
    }
    summ_path = path.join(work, "finalize_summary.json")
    with open(summ_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print("[finalize] wrote", summ_path, flush=True)
    print("[finalize] DONE combined cov →", pkl_path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
