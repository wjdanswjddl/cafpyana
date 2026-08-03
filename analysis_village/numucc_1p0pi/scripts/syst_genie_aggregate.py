#!/usr/bin/env python
"""Reduce phase for chunked GENIE reweight covariances.

Reads ``genie__<GROUP>__*.pkl`` from ``--chunks-dir`` (same tree as
``run_syst_genie_chunked.sh`` / ``syst_genie_parallel.py`` map phase), runs
:func:`get_systematics_genie.run_chunk_merge` once per knob group, then folds
per-knob **fractional** covariance matrices into the legacy ``cov_mat_dict``
pickle layout expected by :func:`analysis_village.numucc_1p0pi.utils.get_syst_unc`
and the breakdown notebook:

* ``<var_save_name>["genie"]`` / ``["genie_rate"]`` — independent sums of per-knob
  ``cov_frac`` (xsec / rate), matching the old notebook recipe.
* ``<var_save_name>[<knob>]`` / ``[<knob>_rate]`` — per-knob ``cov_frac`` blocks.

Writes ``<syst-disk-root>/GENIE/cov_mat_dict.pkl`` (see ``syst_disk_layout``).

Example::

    python syst_genie_aggregate.py \\
        --chunks-dir /path/to/genie_syst-chunked-20260101/chunks \\
        --syst-disk-root "$NUMUCC_SYST_DISK_ROOT"
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import pickle
import sys
import tempfile
from os import path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.dataset_locations import (  # noqa: E402
    GENIE_GROUP_ORDER,
    _genie_glob_map,
)
from analysis_village.numucc_1p0pi.scripts.get_systematics_genie import (  # noqa: E402
    GENIE_MERGE_COMBINED_KEY,
    genie_all_var_configs,
    run_chunk_merge,
)
from analysis_village.numucc_1p0pi.syst_disk_layout import (  # noqa: E402
    FILE_GENIE,
    SUB_GENIE,
    category_out_dir,
    normalized_root,
)


def _ordered_groups_with_chunks(
    chunks_dir: str,
    stage: str,
    allowed: Optional[Sequence[str]],
) -> List[str]:
    """``GENIE_GROUP_ORDER`` first, then remaining glob-map keys; only groups with pickles."""
    gmap = _genie_glob_map(stage)
    allowed_set = set(allowed) if allowed else None
    seen: set[str] = set()
    out: List[str] = []
    for g in GENIE_GROUP_ORDER:
        if g not in gmap:
            continue
        if allowed_set is not None and g not in allowed_set:
            continue
        seen.add(g)
        if glob.glob(path.join(chunks_dir, "genie__%s__*.pkl" % g)):
            out.append(g)
    for g in sorted(k for k in gmap if k not in seen):
        if allowed_set is not None and g not in allowed_set:
            continue
        if glob.glob(path.join(chunks_dir, "genie__%s__*.pkl" % g)):
            out.append(g)
    return out


def _detect_input_stage(chunks_dir: str) -> str:
    paths = sorted(glob.glob(path.join(chunks_dir, "genie__*__*.pkl")))
    if not paths:
        raise SystemExit("[genie-aggregate] no genie__*__*.pkl under %s" % chunks_dir)
    with open(paths[0], "rb") as f:
        d = pickle.load(f)
    st = d.get("input_stage", "final")
    if st not in ("final", "sel_all"):
        return "final"
    return st


def _zero_cov(vc) -> np.ndarray:
    n = len(vc.bin_centers)
    return np.zeros((n, n), dtype=np.float64)


def _init_cov_mat_dict(var_configs: list) -> Dict[str, Dict[str, np.ndarray]]:
    """Per-variable dict with zeroed ``genie`` / ``genie_rate`` totals."""
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for vc in var_configs:
        z = _zero_cov(vc)
        out[vc.var_save_name] = {
            "genie": z.copy(),
            "genie_rate": z.copy(),
        }
    return out


def _accumulate_group_syst(
    cov_mat_dict: Dict[str, Dict[str, np.ndarray]],
    group_syst: Dict[str, Any],
) -> None:
    """Add per-knob packs from one group's ``run_chunk_merge`` dict into ``cov_mat_dict``."""
    for knob, by_slug in group_syst.items():
        if knob == GENIE_MERGE_COMBINED_KEY:
            continue
        if not isinstance(by_slug, dict):
            continue
        for slug, packs in by_slug.items():
            if slug not in cov_mat_dict:
                continue
            row = cov_mat_dict[slug]
            if not isinstance(packs, dict):
                continue
            if "rate" in packs and isinstance(packs["rate"], dict):
                cf = np.asarray(packs["rate"].get("cov_frac"), dtype=np.float64)
                if cf.ndim == 2:
                    rk = "%s_rate" % knob
                    if rk not in row:
                        row[rk] = np.zeros_like(cf)
                    row[rk] += cf
                    row["genie_rate"] += cf
            if "xsec" in packs and isinstance(packs["xsec"], dict):
                cf = np.asarray(packs["xsec"].get("cov_frac"), dtype=np.float64)
                if cf.ndim == 2:
                    if knob not in row:
                        row[knob] = np.zeros_like(cf)
                    row[knob] += cf
                    row["genie"] += cf


def run_genie_syst_aggregate(
    chunks_dir: str,
    syst_disk_root: str,
    *,
    mc_df_stage: str = "final",
    xsec_unit: float = 1.0,
    genie_groups: Optional[Sequence[str]] = None,
) -> str:
    """Merge GENIE chunk pickles and write ``GENIE/cov_mat_dict.pkl``. Returns output path."""
    chunks_dir = path.abspath(path.expanduser(chunks_dir))
    root = normalized_root(syst_disk_root)
    genie_dir = category_out_dir(root, SUB_GENIE)
    os.makedirs(genie_dir, exist_ok=True)
    out_pkl = path.join(genie_dir, FILE_GENIE)

    if mc_df_stage not in ("final", "sel_all"):
        raise SystemExit("[genie-aggregate] mc_df_stage must be final or sel_all")

    detected = _detect_input_stage(chunks_dir)
    if detected != mc_df_stage:
        print(
            "[genie-aggregate] WARNING: chunks input_stage=%r differs from --mc-df-stage=%r "
            "(using chunk value for var registry and group map)"
            % (detected, mc_df_stage)
        )
    stage = detected
    var_configs = genie_all_var_configs(stage)

    allowed = [s.strip() for s in (genie_groups or []) if s.strip()] or None
    if allowed:
        gmap = _genie_glob_map(stage)
        bad = set(allowed) - set(gmap.keys())
        if bad:
            raise SystemExit(
                "[genie-aggregate] unknown GENIE group(s) %r; valid: %s"
                % (sorted(bad), tuple(gmap.keys()))
            )

    groups = _ordered_groups_with_chunks(chunks_dir, stage, allowed)
    if not groups:
        raise SystemExit(
            "[genie-aggregate] no genie__<GROUP>__*.pkl under %s (check MC_DF_STAGE / group filter)"
            % chunks_dir
        )

    print(
        "[genie-aggregate] chunks_dir=%s syst_disk_root=%s stage=%s groups=%s"
        % (chunks_dir, root, stage, ",".join(groups))
    )

    cov_mat_dict = _init_cov_mat_dict(var_configs)

    with tempfile.TemporaryDirectory(prefix="genie_agg_merge_") as tmp_root:
        for grp in groups:
            tmp_out = path.join(tmp_root, grp)
            os.makedirs(tmp_out, exist_ok=True)
            print("[genie-aggregate] chunk-merge group=%s …" % grp)
            group_syst = run_chunk_merge(
                chunks_dir,
                tmp_out,
                grp,
                var_configs,
                xsec_unit,
                bkgd_subtract=True,
                save_figs=False,
                npz_path=None,
            )
            _accumulate_group_syst(cov_mat_dict, group_syst)

    with open(out_pkl, "wb") as f:
        pickle.dump(cov_mat_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("[genie-aggregate] wrote", out_pkl)

    manifest = {
        "schema": "numucc_genie_cov_mat_dict_v1",
        "chunks_dir": chunks_dir,
        "mc_df_stage_cli": mc_df_stage,
        "input_stage": stage,
        "genie_groups": groups,
        "xsec_unit": float(xsec_unit),
        "output_pkl": out_pkl,
    }
    man_path = path.join(genie_dir, "genie_covariance_manifest.json")
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print("[genie-aggregate] wrote", man_path)

    return out_pkl


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--chunks-dir",
        "--chunks_dir",
        dest="chunks_dir",
        required=True,
        help="Directory with genie__<GROUP>__*.pkl from the GENIE map phase.",
    )
    p.add_argument(
        "--syst-disk-root",
        "--syst_disk_root",
        "--out-dir",
        dest="syst_disk_root",
        required=True,
        help="Syst disk root (writes GENIE/cov_mat_dict.pkl).",
    )
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--xsec-unit", type=float, default=1.0)
    p.add_argument(
        "--genie-groups",
        default="",
        help="Comma-separated subset of GENIE groups (default: all groups with chunk pickles).",
    )
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    groups = [s.strip() for s in args.genie_groups.split(",") if s.strip()] or None
    run_genie_syst_aggregate(
        args.chunks_dir,
        args.syst_disk_root,
        mc_df_stage=args.mc_df_stage,
        xsec_unit=float(args.xsec_unit),
        genie_groups=groups,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
