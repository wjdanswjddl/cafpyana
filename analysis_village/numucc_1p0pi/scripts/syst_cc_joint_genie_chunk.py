#!/usr/bin/env python3
"""Map phase: one MC ``.df`` → pickle with **joint-bin** GENIE **rate** histograms per kinematic pair.

For each GENIE knob in the active group, stacks per-universe **event-rate** histograms
(``get_univ_rates(..., cov_type=\"rate\")``) for variable **X** (proton) then **Y** (muon),
matching :mod:`syst_cc_joint_multisim_chunk` bin ordering. Output files are named
``nu__joint_cc_genie__<GROUP>__<stem>.pkl`` (see :mod:`syst_cc_joint_multisim_common`).

Reuses accumulator layout key ``rate_univ_cv`` from :mod:`get_systematics_genie` so
:func:`get_systematics_genie.merge_genie_chunk_pickles` can merge chunk files across ``.df``
stems within a group.

* ``input_stage=final`` only (same HDF keys as marginal GENIE chunk-map).
* ``sel_all`` is not implemented here — use the marginal GENIE pipeline for cut-stage tensors.
* Optional: ``NUMUCC_JOINT_GENIE_SHAPES`` — ``1`` (default): print stacked ``u_j`` / ``cv_j`` shapes
  once per (knob, kinematic pair) per worker; ``all``: every HDF split; ``0``: off.

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

from analysis_village.numucc_1p0pi.dataset_locations import GENIE_GROUP_KNOBS  # noqa: E402
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
    default_kinematic_joint_pairs,
    joint_cc_genie_chunk_basename,
    parse_pair_slugs_csv,
    select_joint_pairs,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates  # noqa: E402

# Log stacked (X|Y) universe matrix shapes once per (knob, pair) per process unless
# ``NUMUCC_JOINT_GENIE_SHAPES=all`` (every HDF split) or ``NUMUCC_JOINT_GENIE_SHAPES=0`` (off).
_LOGGED_JOINT_GENIE_SHAPES: set[tuple[str, str]] = set()


def accumulate_joint_genie_rate_split(
    mc_evt_df: pd.DataFrame,
    mc_nu_df: pd.DataFrame,
    blob_root: dict[str, Any],
    joint_pairs: tuple[tuple[str, object, object], ...],
    syst_names: Sequence[SystName],
    *,
    bkgd_subtract: bool = True,
) -> None:
    validate_split_pair(mc_evt_df, mc_nu_df, -1)
    rate_blk = blob_root.setdefault(RATE_ACC_KEY, {})

    for syst_name in syst_names:
        knob = syst_name[1]
        n_univ = normalize_and_infer_n_univ(mc_evt_df, mc_nu_df, syst_name)

        for _pair_slug, vx, vy in joint_pairs:
            u_x, c_x = get_univ_rates(
                cov_type="rate",
                syst_type="GENIE",
                evtdf=mc_evt_df,
                nudf=mc_nu_df,
                var_config=vx,
                syst_name=syst_name,
                n_univ=n_univ,
                bkgd_subtract=bkgd_subtract,
                plot=False,
            )
            u_y, c_y = get_univ_rates(
                cov_type="rate",
                syst_type="GENIE",
                evtdf=mc_evt_df,
                nudf=mc_nu_df,
                var_config=vy,
                syst_name=syst_name,
                n_univ=n_univ,
                bkgd_subtract=bkgd_subtract,
                plot=False,
            )
            u_x = np.asarray(u_x, dtype=np.float64)
            u_y = np.asarray(u_y, dtype=np.float64)
            c_x = np.asarray(c_x, dtype=np.float64).reshape(-1)
            c_y = np.asarray(c_y, dtype=np.float64).reshape(-1)
            if u_x.ndim != 2 or u_y.ndim != 2:
                raise ValueError("joint GENIE univ arrays must be 2-D, got %s and %s" % (u_x.shape, u_y.shape))
            if u_x.shape[0] != u_y.shape[0]:
                raise ValueError(
                    "joint GENIE n_univ mismatch X vs Y: %d vs %d (shapes %s vs %s)"
                    % (u_x.shape[0], u_y.shape[0], u_x.shape, u_y.shape)
                )
            if c_x.size != u_x.shape[1] or c_y.size != u_y.shape[1]:
                raise ValueError(
                    "joint GENIE cv length vs univ width: len(cv_X)=%d vs n_bins_X=%d; len(cv_Y)=%d vs n_bins_Y=%d"
                    % (c_x.size, u_x.shape[1], c_y.size, u_y.shape[1])
                )
            u_j = np.hstack([u_x, u_y])
            c_j = np.concatenate([c_x, c_y])
            _shape_log = os.environ.get("NUMUCC_JOINT_GENIE_SHAPES", "1").strip().lower()
            _sk = (knob, _pair_slug)
            if _shape_log != "0":
                if _shape_log == "all" or _sk not in _LOGGED_JOINT_GENIE_SHAPES:
                    if _shape_log != "all":
                        _LOGGED_JOINT_GENIE_SHAPES.add(_sk)
                    print(
                        "[cc-joint-genie-chunk] joint stacked shapes knob=%r pair=%r: "
                        "u_x=%s u_y=%s -> u_j=%s, cv_j=%s"
                        % (_sk[0], _sk[1], u_x.shape, u_y.shape, u_j.shape, c_j.shape),
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
    joint_pairs: tuple[tuple[str, object, object], ...],
    *,
    max_splits: int = 0,
    input_stage: str = "final",
    bkgd_subtract: bool = True,
) -> str:
    if input_stage != "final":
        raise SystemExit("[cc-joint-genie-chunk] only --input-stage final is implemented")
    knobs = list(GENIE_GROUP_KNOBS.get(genie_group) or [])
    if not knobs:
        raise SystemExit("[cc-joint-genie-chunk] no knobs for GENIE group %r" % genie_group)
    syst_names: list[SystName] = [("mc", k) for k in knobs]

    os.makedirs(out_dir, exist_ok=True)
    n_keys = int(get_n_split(df_file))
    n_use = n_keys if max_splits <= 0 else min(max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[cc-joint-genie-chunk] no HDF splits")

    blob_root: dict[str, Any] = {
        "kind": "joint_genie_cc_chunk",
        "input_stage": input_stage,
        "meta": {
            "df_file": df_file,
            "splits_processed": n_use,
            "genie_group": genie_group,
            "input_stage": input_stage,
            "pairs": [p[0] for p in joint_pairs],
        },
        RATE_ACC_KEY: {},
    }

    for i in tqdm(range(n_use), desc="HDF splits (joint GENIE)"):
        mc_evt_df = pd.read_hdf(df_file, key=f"evt_{i}")
        mc_nu_df = pd.read_hdf(df_file, key=f"mcnu_{i}")
        validate_genie_dataframes({"evt": mc_evt_df, "mcnu": mc_nu_df}, context=f"split {i}:")
        mc_evt_df = mc_evt_df.copy()
        mc_nu_df = mc_nu_df.copy()
        # _prefix_mcnu_columns(mc_nu_df)
        # mc_nu_df = ensure_mc_level_phi_mcnu(mc_nu_df)
        mc_evt_df = ensure_derived_trk_kinematics_cols(mc_evt_df)
        mc_evt_df = add_reco_cc1p0pi_tki_evtdf(mc_evt_df)
        mc_evt_df = add_truth_cc1p0pi_tki_evtdf(mc_evt_df)
        mc_nu_df = add_mc_cc1p0pi_tki_mcnu(mc_nu_df)
        _annotate_topo_genie_phi(mc_evt_df, mc_nu_df)
        accumulate_joint_genie_rate_split(
            mc_evt_df,
            mc_nu_df,
            blob_root,
            joint_pairs,
            syst_names,
            bkgd_subtract=bkgd_subtract,
        )
        del mc_evt_df, mc_nu_df
        gc.collect()

    stem = path.splitext(path.basename(df_file))[0]
    out_path = path.join(out_dir, joint_cc_genie_chunk_basename(genie_group, stem))
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
        "--pairs",
        default=None,
        help="Comma-separated preset pair slugs (default: all from default_kinematic_joint_pairs).",
    )
    return p.parse_args()


def compute_out_path(args) -> str:
    stem = path.splitext(path.basename(args.df_file))[0]
    return path.join(args.out_dir, joint_cc_genie_chunk_basename(args.genie_group, stem))


def run_with_args(args, *, skip_existing: bool = False) -> Tuple[str, str]:
    os.makedirs(args.out_dir, exist_ok=True)
    pair_slugs = parse_pair_slugs_csv(args.pairs)
    joint_pairs = select_joint_pairs(pair_slugs)
    out_path = compute_out_path(args)
    if skip_existing and path.exists(out_path):
        return out_path, "skipped"
    run_joint_genie_chunk_map(
        args.df_file,
        args.out_dir,
        args.genie_group,
        joint_pairs,
        max_splits=int(args.max_splits),
        input_stage=args.input_stage,
    )
    return out_path, "ok"


def main():
    args = parse_args()
    run_with_args(args)


if __name__ == "__main__":
    main()
