#!/usr/bin/env python3
"""Map phase: one MC ``.df`` → pickle with **joint-bin** GENIE **rate** histograms.

Default (``--mode stack``): inclusive selected-rate stack
(``get_univ_rates(..., cov_type=rate, bkgd_subtract=False)``) matching
:mod:`syst_cc_joint_multisim_chunk`. GENIE knobs come from
:func:`dataset_locations.joint_cc_genie_knobs_for_group` (FSI_compare →
``GENIE_slim_v3`` only). Pickle prefix ``nu__joint_cc_genie_stack__``.

Aggregate with :mod:`analysis_village.numucc_1p0pi.scripts.syst_cc_joint_genie_aggregate`.
"""

from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from os import path
from typing import Any, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x


os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.split_df_helpers import get_n_split  # noqa: E402

from analysis_village.numucc_1p0pi.dataset_locations import joint_cc_genie_knobs_for_group  # noqa: E402
from analysis_village.numucc_1p0pi.scripts.get_systematics_genie import (  # noqa: E402
    RATE_ACC_KEY,
    SystName,
    _annotate_topo_genie_phi,
    _prefix_mcnu_columns,
    add_mc_cc1p0pi_tki_mcnu,
    add_reco_cc1p0pi_tki_evtdf,
    add_truth_cc1p0pi_tki_evtdf,
    ensure_mc_level_phi_mcnu,
    normalize_and_infer_n_univ,
    validate_genie_dataframes,
    validate_split_pair,
)
from analysis_village.numucc_1p0pi.evt_derived_kinematics import ensure_derived_trk_kinematics_cols  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    joint_cc_genie_chunk_basename,
    parse_pair_slugs_csv,
    stack_jobs,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates  # noqa: E402

# Log stacked (X|Y) universe matrix shapes once per (knob, pair) per process unless
# ``NUMUCC_JOINT_GENIE_SHAPES=all`` (every HDF split) or ``NUMUCC_JOINT_GENIE_SHAPES=0`` (off).
_LOGGED_JOINT_GENIE_SHAPES: set[tuple[str, str]] = set()
_LOGGED_JOINT_GENIE_SKIPPED_KNOBS: set[str] = set()


def accumulate_joint_genie_rate_split(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    blob_root: dict[str, Any],
    jobs: tuple[tuple[str, tuple], ...],
    syst_names: Sequence[SystName],
    *,
    bkgd_subtract: bool = False,
) -> None:
    validate_split_pair(mc_evt_df, mc_nu_df, -1)
    rate_blk = blob_root.setdefault(RATE_ACC_KEY, {})

    for syst_name in syst_names:
        knob = syst_name[1]
        n_univ = normalize_and_infer_n_univ(
            mc_evt_df, mc_nu_df, syst_name, raise_if_missing=False
        )
        if n_univ <= 0:
            if knob not in _LOGGED_JOINT_GENIE_SKIPPED_KNOBS:
                _LOGGED_JOINT_GENIE_SKIPPED_KNOBS.add(knob)
                print(
                    "[cc-joint-genie-chunk] skip missing knob %r (not on this DF)" % knob,
                    flush=True,
                )
            continue

        for _pair_slug, vars_ in jobs:
            us = []
            cs = []
            names = []
            for var in vars_:
                u, c = get_univ_rates(
                    cov_type="rate",
                    syst_type="GENIE",
                    evtdf=mc_evt_df,
                    nudf=mc_nu_df,
                    var_config=var,
                    syst_name=syst_name,
                    n_univ=n_univ,
                    bkgd_subtract=bkgd_subtract,
                    plot=False,
                )
                u = np.asarray(u, dtype=np.float64)
                c = np.asarray(c, dtype=np.float64).reshape(-1)
                if u.ndim != 2:
                    raise ValueError("joint GENIE univ must be 2-D, got %s" % (u.shape,))
                if c.size != u.shape[1]:
                    raise ValueError(
                        "joint GENIE cv vs univ width for %s: %d vs %d"
                        % (var.var_save_name, c.size, u.shape[1])
                    )
                us.append(u)
                cs.append(c)
                names.append(var.var_save_name)
            n_univ_set = {u.shape[0] for u in us}
            if len(n_univ_set) != 1:
                raise ValueError("joint GENIE n_univ mismatch: %s" % [u.shape for u in us])
            u_j = np.hstack(us)
            c_j = np.concatenate(cs)
            _shape_log = os.environ.get("NUMUCC_JOINT_GENIE_SHAPES", "1").strip().lower()
            _sk = (knob, _pair_slug)
            if _shape_log != "0":
                if _shape_log == "all" or _sk not in _LOGGED_JOINT_GENIE_SHAPES:
                    if _shape_log != "all":
                        _LOGGED_JOINT_GENIE_SHAPES.add(_sk)
                    print(
                        "[cc-joint-genie-chunk] stacked shapes knob=%r slug=%r vars=%s: u_j=%s cv_j=%s"
                        % (_sk[0], _sk[1], names, u_j.shape, c_j.shape),
                        flush=True,
                    )
            slot = rate_blk.setdefault(knob, {}).setdefault(
                _pair_slug,
                {"univ": np.zeros_like(u_j, dtype=np.float64), "cv": np.zeros_like(c_j, dtype=np.float64)},
            )
            slot["univ"] += u_j
            slot["cv"] += c_j


def run_joint_genie_chunk_map(
    df_file: str,
    out_dir: str,
    genie_group: str,
    jobs: tuple[tuple[str, tuple], ...],
    *,
    max_splits: int = 0,
    input_stage: str = "final",
    bkgd_subtract: bool = False,
    mode: str = "stack",
) -> str:
    if input_stage != "final":
        raise SystemExit("[cc-joint-genie-chunk] only --input-stage final is implemented")
    knobs = joint_cc_genie_knobs_for_group(genie_group)
    if not knobs:
        raise SystemExit("[cc-joint-genie-chunk] no knobs for GENIE group %r" % genie_group)
    syst_names: list[SystName] = [("mc", k) for k in knobs]

    os.makedirs(out_dir, exist_ok=True)
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[cc-joint-genie-chunk] no HDF splits")

    blob_root: dict[str, Any] = {
        "kind": "joint_genie_cc_stack_chunk" if mode == "stack" else "joint_genie_cc_chunk",
        "input_stage": input_stage,
        "mode": mode,
        "bkgd_subtract": bool(bkgd_subtract),
        "meta": {
            "df_file": df_file,
            "splits_processed": n_use,
            "genie_group": genie_group,
            "genie_knobs": knobs,
            "input_stage": input_stage,
            "pairs": [p[0] for p in jobs],
            "mode": mode,
            "bkgd_subtract": bool(bkgd_subtract),
        },
        RATE_ACC_KEY: {},
    }

    for i in tqdm(range(n_use), desc="HDF splits (joint GENIE)"):
        mc_evt_df = pd.read_hdf(df_file, key=f"evt_{i}")
        mc_nu_df = pd.read_hdf(df_file, key=f"mcnu_{i}")
        validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context=f"split {i}:")
        mc_evt_df = mc_evt_df.copy()
        mc_nu_df = mc_nu_df.copy()
        _prefix_mcnu_columns(mc_nu_df)
        mc_evt_df = ensure_derived_trk_kinematics_cols(mc_evt_df)
        mc_evt_df = add_reco_cc1p0pi_tki_evtdf(mc_evt_df)
        mc_evt_df = add_truth_cc1p0pi_tki_evtdf(mc_evt_df)
        mc_nu_df = add_mc_cc1p0pi_tki_mcnu(mc_nu_df)
        _annotate_topo_genie_phi(mc_evt_df, mc_nu_df)
        accumulate_joint_genie_rate_split(
            mc_evt_df,
            mc_nu_df,
            blob_root,
            jobs,
            syst_names,
            bkgd_subtract=bkgd_subtract,
        )
        del mc_evt_df, mc_nu_df
        gc.collect()

    stem = path.splitext(path.basename(df_file))[0]
    out_path = path.join(out_dir, joint_cc_genie_chunk_basename(genie_group, stem, mode=mode))
    tmp_path = out_path + ".tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(blob_root, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, out_path)
    print("[cc-joint-genie-chunk] wrote", out_path)
    return out_path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--genie_group", required=True)
    p.add_argument("--input-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument(
        "--mode",
        choices=("stack", "pairs"),
        default="stack",
        help="``stack`` (default): inclusive 4-var vector. ``pairs``: legacy pairwise shards.",
    )
    p.add_argument(
        "--pairs",
        default=None,
        help="Comma-separated preset pair slugs (implies --mode pairs).",
    )
    p.add_argument(
        "--bkgd-subtract",
        action="store_true",
        help="Legacy signal-subtracted universes. Default is inclusive selected rate.",
    )
    return p.parse_args()


def _resolve_mode(args) -> str:
    if getattr(args, "pairs", None):
        return "pairs"
    return str(getattr(args, "mode", "stack") or "stack")


def compute_out_path(args) -> str:
    stem = path.splitext(path.basename(args.df_file))[0]
    return path.join(args.out_dir, joint_cc_genie_chunk_basename(args.genie_group, stem, mode=_resolve_mode(args)))


def run_with_args(args, *, skip_existing: bool = False) -> Tuple[str, str]:
    os.makedirs(args.out_dir, exist_ok=True)
    mode = _resolve_mode(args)
    pair_slugs = parse_pair_slugs_csv(args.pairs) if mode == "pairs" else None
    jobs = stack_jobs(mode=mode, pair_slugs=pair_slugs)
    out_path = compute_out_path(args)
    if skip_existing and path.exists(out_path):
        return out_path, "skipped"
    run_joint_genie_chunk_map(
        args.df_file,
        args.out_dir,
        args.genie_group,
        jobs,
        max_splits=int(args.max_splits),
        input_stage=args.input_stage,
        bkgd_subtract=bool(getattr(args, "bkgd_subtract", False)),
        mode=mode,
    )
    return out_path, "ok"


def main():
    args = parse_args()
    run_with_args(args)


if __name__ == "__main__":
    main()
