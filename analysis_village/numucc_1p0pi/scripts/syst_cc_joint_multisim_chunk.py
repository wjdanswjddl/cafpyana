#!/usr/bin/env python3
"""Map phase: one MC ``.df`` → pickle with **joint-bin** multisim histograms per kinematic pair.

For each requested systematic (MCstat / Flux / G4 — same column layout as
``syst_multisim_chunk.py``) and each kinematic pair ``(var_X, var_Y)``, this script
histograms **both** variables under the same universe weights into a single stacked
vector of length ``n_X + n_Y`` (**constrained / proton channel X first**, then **constraining /
muon channel Y**). Chunk pickles mirror the marginal multisim chunk schema except variable keys
are *preset* pair slugs (e.g. ``muon_p__proton_costheta``). Output basename prefix is
``nu__joint_cc__`` (see :mod:`syst_cc_joint_multisim_common`) so parallel ``skip_existing`` does
not reuse pre-expansion ``nu__joint__*`` shards.

* Optional: ``NUMUCC_JOINT_MULTISIM_SHAPES`` — ``1`` (default): print stacked ``u_j`` / ``cv_j``
  shapes once per ``(syst_type, syst_key, var_X, var_Y)`` per worker; ``all``: every call; ``0``: off.

Aggregate with :mod:`analysis_village.numucc_1p0pi.scripts.syst_cc_joint_multisim_aggregate`.
"""

from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from os import path
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x


os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))

from pyanalib.split_df_helpers import get_n_split  # noqa: E402

from analysis_village.numucc_1p0pi.categories import get_topo_category  # noqa: E402
from analysis_village.numucc_1p0pi.evt_derived_kinematics import ensure_derived_trk_kinematics_cols  # noqa: E402
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key  # noqa: E402
from analysis_village.numucc_1p0pi.syst_cc_joint_multisim_common import (  # noqa: E402
    default_kinematic_joint_pairs,
    joint_cc_multisim_chunk_basename,
    parse_pair_slugs_csv,
    select_joint_pairs,
)
from analysis_village.numucc_1p0pi.syst_multisim_common import (  # noqa: E402
    NEUTRINO_SYST_ORDER,
    DEFAULT_MULTISIM_CHUNK_SYST_NAMES,
    drop_bad_flux_knob_weights,
    drop_bad_g4_knob_weights,
    drop_bad_g4_weights,
    flux_mc_knob_names,
    g4_mc_knob_names,
    knob_nested_syst_block,
    syst_key_for_name,
)
from analysis_village.numucc_1p0pi.utils import get_univ_rates  # noqa: E402

# See module docstring: ``NUMUCC_JOINT_MULTISIM_SHAPES`` = 1 | all | 0
_LOGGED_JOINT_MULTISIM_SHAPES: set[tuple[str, str, str, str]] = set()


def _parse_syst_names(spec: str | None):
    if not spec or not str(spec).strip():
        return tuple(DEFAULT_MULTISIM_CHUNK_SYST_NAMES)
    raw = [x.strip() for x in str(spec).split(",") if x.strip()]
    if len(raw) == 1 and raw[0].lower() == "full":
        return tuple(NEUTRINO_SYST_ORDER)
    unk = set(raw) - set(NEUTRINO_SYST_ORDER)
    if unk:
        raise SystemExit("[cc-joint-multisim-chunk] unknown --syst-names entries: %s" % sorted(unk))
    ordered = tuple(sn for sn in NEUTRINO_SYST_ORDER if sn in raw)
    if len(ordered) != len(raw):
        raise SystemExit("[cc-joint-multisim-chunk] duplicate syst name in --syst-names")
    return ordered


def _count_univ_columns(evt_df: pd.DataFrame, syst_col_key: Any, cap: int = 512) -> int:
    n = 0
    for i in range(cap):
        if isinstance(syst_col_key, tuple):
            probe = syst_col_key + (f"univ_{i}",)
        else:
            probe = (syst_col_key, f"univ_{i}")
        if multicol_resolve_column_key(evt_df, probe) is None:
            break
        n += 1
    return n


def _univ_syst_key_for_df(evt_df: pd.DataFrame, sname: str) -> Any:
    sk = syst_key_for_name(sname)
    if sname != "MCstat":
        return sk
    if _count_univ_columns(evt_df, ("mc", "MCstat")) > 0:
        return ("mc", "MCstat")
    if _count_univ_columns(evt_df, "MCstat") > 0:
        return "MCstat"
    return sk


def _merge_joint_pack(dst: dict, pair_slug: str, univ: np.ndarray, cv: np.ndarray, err: str) -> None:
    u = np.asarray(univ, dtype=float)
    c = np.asarray(cv, dtype=float)
    if pair_slug not in dst:
        dst[pair_slug] = {"univ_events": u.copy(), "cv_events": c.copy()}
    else:
        cur = dst[pair_slug]
        if cur["univ_events"].shape != u.shape:
            raise ValueError(
                "joint shape mismatch %s: %s vs %s" % (err, cur["univ_events"].shape, u.shape)
            )
        cur["univ_events"] += u
        cur["cv_events"] += c


def _merge_flux_or_g4_joint_block(merged_block: dict, raw_block: dict, label: str) -> None:
    if not raw_block:
        return
    if knob_nested_syst_block(raw_block):
        if merged_block and not knob_nested_syst_block(merged_block):
            raise ValueError("[cc-joint] cannot mix flat %s with knob-nested chunks" % label)
        for knob, kb in raw_block.items():
            tgt = merged_block.setdefault(knob, {})
            for pair_slug, pack in kb.items():
                _merge_joint_pack(tgt, pair_slug, pack["univ_events"], pack["cv_events"], "%s/%s" % (knob, pair_slug))
        return
    if merged_block and knob_nested_syst_block(merged_block):
        raise ValueError("[cc-joint] cannot mix knob-nested %s with flat chunks" % label)
    for pair_slug, pack in raw_block.items():
        _merge_joint_pack(merged_block, pair_slug, pack["univ_events"], pack["cv_events"], pair_slug)


def _joint_univ_rates_one_syst(
    mc_evt_df: pd.DataFrame,
    var_X,
    var_Y,
    syst_name_key,
    n_univ: int,
    syst_type: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Stacked per-universe histograms: shape ``(n_univ, n_X + n_Y)``."""
    u_x, c_x = get_univ_rates(
        cov_type="rate",
        syst_type=syst_type,
        evtdf=mc_evt_df,
        nudf=None,
        var_config=var_X,
        syst_name=syst_name_key,
        n_univ=n_univ,
        bkgd_subtract=True,
    )
    u_y, c_y = get_univ_rates(
        cov_type="rate",
        syst_type=syst_type,
        evtdf=mc_evt_df,
        nudf=None,
        var_config=var_Y,
        syst_name=syst_name_key,
        n_univ=n_univ,
        bkgd_subtract=True,
    )
    if u_x.ndim != 2 or u_y.ndim != 2:
        raise ValueError("joint univ arrays must be 2-D, got %s and %s" % (u_x.shape, u_y.shape))
    if u_x.shape[0] != u_y.shape[0]:
        raise ValueError(
            "joint n_univ mismatch X vs Y: %d vs %d (shapes %s vs %s)"
            % (u_x.shape[0], u_y.shape[0], u_x.shape, u_y.shape)
        )
    c_x = np.asarray(c_x, dtype=np.float64).reshape(-1)
    c_y = np.asarray(c_y, dtype=np.float64).reshape(-1)
    if c_x.size != u_x.shape[1] or c_y.size != u_y.shape[1]:
        raise ValueError(
            "joint cv length vs univ width: len(cv_X)=%d vs n_bins_X=%d; len(cv_Y)=%d vs n_bins_Y=%d"
            % (c_x.size, u_x.shape[1], c_y.size, u_y.shape[1])
        )
    u_j = np.hstack([u_x, u_y])
    c_j = np.concatenate([c_x, c_y])
    _sk = (syst_type, repr(syst_name_key), var_X.var_save_name, var_Y.var_save_name)
    _shape_log = os.environ.get("NUMUCC_JOINT_MULTISIM_SHAPES", "1").strip().lower()
    if _shape_log != "0":
        if _shape_log == "all" or _sk not in _LOGGED_JOINT_MULTISIM_SHAPES:
            if _shape_log != "all":
                _LOGGED_JOINT_MULTISIM_SHAPES.add(_sk)
            print(
                "[cc-joint-multisim-chunk] joint stacked shapes syst=%r syst_key=%s var_X=%r var_Y=%r: "
                "u_x=%s u_y=%s -> u_j=%s cv_j=%s"
                % (
                    syst_type,
                    syst_name_key,
                    var_X.var_save_name,
                    var_Y.var_save_name,
                    u_x.shape,
                    u_y.shape,
                    u_j.shape,
                    c_j.shape,
                ),
                flush=True,
            )
    return u_j, c_j


def _accum_univ_joint_for_mc_knobs(
    mc_evt_df: pd.DataFrame,
    joint_pairs: tuple[tuple[str, object, object], ...],
    knobs: tuple[str, ...],
    acc_cat: dict,
    log_tag: str,
    n_univ_requested: int,
) -> None:
    for knob in knobs:
        sk = ("mc", knob)
        n_u = min(int(n_univ_requested), _count_univ_columns(mc_evt_df, sk))
        if n_u <= 0:
            continue
        knob_acc = acc_cat.setdefault(knob, {})
        for pair_slug, vx, vy in joint_pairs:
            try:
                uj, cj = _joint_univ_rates_one_syst(mc_evt_df, vx, vy, sk, n_u, syst_type=log_tag)
            except Exception as ex:
                print(f"[cc-joint-multisim-chunk] skip pair={pair_slug} {log_tag} knob={knob}: {ex}")
                continue
            _merge_joint_pack(knob_acc, pair_slug, uj, cj, "%s/%s/%s" % (log_tag, knob, pair_slug))


def _accumulate_joint_final(
    args,
    syst_names: tuple[str, ...],
    joint_pairs: tuple[tuple[str, object, object], ...],
) -> Dict[str, Any]:
    syst_active = frozenset(syst_names)
    acc_syst: Dict[str, Any] = {sn: {} for sn in NEUTRINO_SYST_ORDER}
    g4_knobs = g4_mc_knob_names() if getattr(args, "g4_mode", "knobs") == "knobs" else ()
    flux_knobs = (
        flux_mc_knob_names(args.flux_knob_groups)
        if getattr(args, "flux_mode", "knobs") == "knobs"
        else ()
    )

    def flush_evt(mc_evt_df: pd.DataFrame) -> None:
        if "topo_categ" not in mc_evt_df.columns:
            mc_evt_df = mc_evt_df.copy()
            mc_evt_df.loc[:, "topo_categ"] = get_topo_category(mc_evt_df)
        mc_evt_df = ensure_derived_trk_kinematics_cols(mc_evt_df)
        if "G4" in syst_active:
            if args.g4_mode == "knobs":
                mc_evt_df = drop_bad_g4_knob_weights(mc_evt_df, knobs=g4_knobs, n_univ=args.n_universe)
            else:
                mc_evt_df = drop_bad_g4_weights(mc_evt_df, n_univ=args.n_universe)
        if "Flux" in syst_active and args.flux_mode == "knobs" and flux_knobs:
            mc_evt_df = drop_bad_flux_knob_weights(mc_evt_df, knobs=flux_knobs, n_univ=args.n_universe)
        if len(mc_evt_df) == 0:
            return
        for sname in syst_names:
            if sname == "Flux" and args.flux_mode == "knobs":
                if flux_knobs:
                    _accum_univ_joint_for_mc_knobs(
                        mc_evt_df, joint_pairs, flux_knobs, acc_syst["Flux"], "Flux", args.n_universe
                    )
                continue
            if sname == "Flux" and args.flux_mode == "bundled":
                sk = _univ_syst_key_for_df(mc_evt_df, "Flux")
                n_u = min(int(args.n_universe), _count_univ_columns(mc_evt_df, sk))
                if n_u > 0:
                    for pair_slug, vx, vy in joint_pairs:
                        try:
                            uj, cj = _joint_univ_rates_one_syst(mc_evt_df, vx, vy, sk, n_u, syst_type="Flux")
                        except Exception as ex:
                            print(f"[cc-joint-multisim-chunk] skip pair={pair_slug} Flux bundled: {ex}")
                            continue
                        _merge_joint_pack(acc_syst["Flux"], pair_slug, uj, cj, "Flux/%s" % pair_slug)
                continue
            if sname == "G4" and args.g4_mode == "knobs":
                if g4_knobs:
                    _accum_univ_joint_for_mc_knobs(
                        mc_evt_df, joint_pairs, g4_knobs, acc_syst["G4"], "G4", args.n_universe
                    )
                continue
            if sname == "G4" and args.g4_mode == "bundled":
                sk = _univ_syst_key_for_df(mc_evt_df, "G4")
                n_u = min(int(args.n_universe), _count_univ_columns(mc_evt_df, sk))
                if n_u > 0:
                    for pair_slug, vx, vy in joint_pairs:
                        try:
                            uj, cj = _joint_univ_rates_one_syst(mc_evt_df, vx, vy, sk, n_u, syst_type="G4")
                        except Exception as ex:
                            print(f"[cc-joint-multisim-chunk] skip pair={pair_slug} G4 bundled: {ex}")
                            continue
                        _merge_joint_pack(acc_syst["G4"], pair_slug, uj, cj, "G4/%s" % pair_slug)
                continue
            sk = _univ_syst_key_for_df(mc_evt_df, sname)
            n_u = min(int(args.n_universe), _count_univ_columns(mc_evt_df, sk))
            if n_u <= 0:
                continue
            stype = "MCstat" if sname == "MCstat" else sname
            for pair_slug, vx, vy in joint_pairs:
                try:
                    uj, cj = _joint_univ_rates_one_syst(mc_evt_df, vx, vy, sk, n_u, syst_type=stype)
                except Exception as ex:
                    print(f"[cc-joint-multisim-chunk] skip pair={pair_slug} syst={sname}: {ex}")
                    continue
                _merge_joint_pack(acc_syst[sname], pair_slug, uj, cj, "%s/%s" % (sname, pair_slug))

    n_keys = int(get_n_split(args.df_file))
    n_use = n_keys if args.max_splits <= 0 else min(args.max_splits, n_keys)
    if n_use <= 0:
        raise SystemExit("[cc-joint-multisim-chunk] no splits")
    for i in tqdm(range(n_use), desc="HDF splits (joint)"):
        mc_evt_df = pd.read_hdf(args.df_file, key=f"evt_{i}")
        flush_evt(mc_evt_df)
        del mc_evt_df
        gc.collect()
    return acc_syst


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument(
        "--input-stage",
        choices=("final", "sel_all"),
        default="final",
        help="``final``: legacy evt-only path. ``sel_all``: not implemented here (use marginal chunk).",
    )
    p.add_argument(
        "--syst-names",
        default=None,
        help="Comma-separated subset of MCstat,Flux,G4 (default: Flux,G4). Use ``full`` for all three.",
    )
    p.add_argument("--n-universe", type=int, default=100)
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument("--g4-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument("--flux-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument(
        "--flux-knob-groups",
        default="all",
        help="Passed to ``flux_mc_knob_names`` when --flux-mode knobs.",
    )
    p.add_argument(
        "--pairs",
        default=None,
        help="Comma-separated pair slugs (default: all kinematic pairs). Example: muon_p__proton_costheta",
    )
    return p.parse_args()


def compute_out_path(args, syst_names) -> str:
    stem = path.splitext(path.basename(args.df_file))[0]
    tag = "_".join(syst_names)
    return path.join(args.out_dir, joint_cc_multisim_chunk_basename(tag, stem))


def run_with_args(args, *, skip_existing: bool = False) -> Tuple[str, str]:
    os.makedirs(args.out_dir, exist_ok=True)
    syst_names = _parse_syst_names(args.syst_names)
    if args.input_stage != "final":
        raise SystemExit("[cc-joint-multisim-chunk] only --input-stage final is implemented")
    pair_slugs = parse_pair_slugs_csv(args.pairs)
    joint_pairs = select_joint_pairs(pair_slugs)

    out_path = compute_out_path(args, syst_names)
    if skip_existing and path.exists(out_path):
        return out_path, "skipped"

    acc_syst = _accumulate_joint_final(args, syst_names, joint_pairs)
    blob = {
        "kind": "joint_multisim_cc_chunk",
        "df_file": args.df_file,
        "syst_names_computed": list(syst_names),
        "pairs": [p[0] for p in joint_pairs],
        "splits_processed": int(get_n_split(args.df_file)) if args.max_splits <= 0 else int(args.max_splits),
        "n_univ_requested": args.n_universe,
        "syst": acc_syst,
    }
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(blob, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, out_path)
    print("[cc-joint-multisim-chunk] wrote", out_path)
    return out_path, "ok"


def main():
    args = parse_args()
    run_with_args(args)


if __name__ == "__main__":
    main()
