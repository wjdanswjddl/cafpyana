#!/usr/bin/env python
"""Compare in-memory vs chunked event selection on one MC shard."""
from __future__ import annotations

import argparse
import sys
import tempfile
from os import path

import numpy as np
import pandas as pd

REPO = path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
sys.path.insert(0, REPO)

from analysis_village.numucc_1p0pi.categories import DETECTOR
from analysis_village.numucc_1p0pi.dataset_locations import iter_event_selection_df_paths
from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_runner
from analysis_village.numucc_1p0pi.selection_framework import (
    aggregate_chunk_files,
    apply_global_exposure_scales,
    merge_samples,
)
from pyanalib.split_df_helpers import get_n_split, load_dfs

SCRIPTS = path.join(REPO, "analysis_village", "numucc_1p0pi", "scripts")
sys.path.insert(0, SCRIPTS)
import event_selection_chunk as chunk_mod


def _load_one_file(df_file: str, max_splits: int) -> dict:
    n = int(get_n_split(df_file))
    cap = min(n, max_splits) if max_splits > 0 else n
    keys = ["evt", "trk", "hdr"]
    if chunk_mod._hdf_has_mcnu(df_file):
        keys.append("mcnu")
    return load_dfs(df_file, keys2load=keys, n_max_concat=cap)


def _purity_from_eff(eff_stage: dict) -> dict:
    var0 = next(iter(eff_stage)) if eff_stage else None
    out = {"final_stage": "2prong-mup", "var": var0}
    if not var0:
        return out
    ea = eff_stage[var0]
    n_raw = float(getattr(ea, "n_at_stage_int_raw", 0.0))
    out.update(
        {
            "n_at_stage_int": ea.n_at_stage_int,
            "n_at_stage_int_raw": n_raw,
            "n_signal_int_raw": ea.n_total_signal_int_raw,
            "purity_weighted": ea.n_total_signal_int / ea.n_at_stage_int * 100
            if ea.n_at_stage_int > 0
            else 0.0,
            "purity_raw": ea.n_total_signal_int_raw / n_raw * 100 if n_raw > 0 else 0.0,
        }
    )
    return out


def _purity_from_bar(bar_stage: dict) -> dict:
    """Notebook-style purity: signal topology bin / total MC at stage."""
    if "topology" not in bar_stage:
        return {}
    bb = bar_stage["topology"]
    total = float(np.sum(bb.mc_counts))
    if total <= 0:
        return {}
    # topology_list order: signal (1p0pi) is index 0
    signal = float(bb.mc_counts[0])
    return {
        "purity_from_bar": 100.0 * signal / total,
        "n_signal_bar": signal,
        "n_total_bar": total,
    }


def run_inmemory_parity(state: dict) -> dict:
    runner = build_runner("mc")
    chunk_mod.attach_intrinsic_weights(state["evt"], state["trk"], "mc")
    evt_df, mcnu_df = chunk_mod._ensure_phi_and_kinematics_cols(
        state["evt"], state["trk"], state.get("mcnu")
    )
    runner.run({"evt": evt_df, "trk": state["trk"], "hdr": state["hdr"], "mcnu": mcnu_df})
    out = {"engine": "inmemory_runner"}
    out.update(_purity_from_eff(runner.eff.get("2prong-mup", {})))
    out.update(_purity_from_bar(runner.bar.get("2prong-mup", {})))
    out["bar_mc_totals"] = {
        sk: float(bar["topology"].total_count("mc"))
        for sk, bar in runner.bar.items()
        if "topology" in bar
    }
    return out


def run_chunk_parity(df_file: str, out_dir: str, max_splits: int) -> dict:
    out_pkl = path.join(out_dir, f"mc__{path.splitext(path.basename(df_file))[0]}.pkl")
    if path.isfile(out_pkl):
        out_pkl.unlink()
    cmd = [
        sys.executable,
        path.join(SCRIPTS, "event_selection_chunk.py"),
        "--df_file",
        df_file,
        "--sample",
        "mc",
        "--out_dir",
        out_dir,
        "--load_mode",
        "splits",
    ]
    if max_splits > 0:
        cmd.extend(["--max_splits", str(max_splits)])
    import subprocess

    subprocess.run(cmd, check=True)
    agg = aggregate_chunk_files([out_pkl])
    merged = merge_samples({"mc": agg})
    totals = type(
        "T",
        (),
        {
            "data_pot": 1.0,
            "mc_pot": 1.0,
            "dirt_pot": 0.0,
            "data_gates_bnb": 1.0,
            "intime_gates": 0.0,
            "offbeam_gates": 0.0,
        },
    )()
    apply_global_exposure_scales(merged, totals, f_offbeam_coincident=0.08)
    out = {"engine": "chunk_pickle", "pkl": out_pkl}
    out.update(_purity_from_eff(merged["eff"].get("2prong-mup", {})))
    out.update(_purity_from_bar(merged["bar"].get("2prong-mup", {})))
    out["bar_mc_totals"] = {
        sk: float(by_bt["topology"].total_count("mc"))
        for sk, by_bt in merged.get("bar", {}).items()
        if "topology" in by_bt
    }
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--df_file", default=None)
    p.add_argument("--max_splits", type=int, default=2)
    args = p.parse_args()

    df_file = args.df_file or next(iter_event_selection_df_paths("mc"))
    print(f"[parity] DETECTOR={DETECTOR!r}")
    print(f"[parity] input={df_file}")
    print(f"[parity] max_splits={args.max_splits}")

    state = _load_one_file(df_file, args.max_splits)
    print(
        f"[parity] loaded evt={len(state['evt'])} trk={len(state['trk'])} "
        f"mcnu={'mcnu' in state}"
    )

    mem = run_inmemory_parity(state)
    with tempfile.TemporaryDirectory(prefix="es_parity_") as tmp:
        chk = run_chunk_parity(df_file, tmp, args.max_splits)

    print("\n=== Purity at final stage (2prong-mup) ===")
    for label, res in (("in-memory runner", mem), ("chunk pickle", chk)):
        print(f"  {label}:")
        print(f"    weighted purity: {res.get('purity_weighted', float('nan')):.4f}%")
        print(f"    raw purity:      {res.get('purity_raw', float('nan')):.4f}%")
        print(f"    n_signal raw:  {res.get('n_signal_int_raw', res.get('n_signal_bar', float('nan'))):.0f}")
        print(f"    n_total raw:   {res.get('n_at_stage_int_raw', res.get('n_total_bar', float('nan'))):.0f}")
        if res.get("purity_from_bar") is not None:
            print(f"    purity (bar):  {res.get('purity_from_bar', float('nan')):.4f}%")

    print("\n=== Per-stage MC totals (bar breakdown, unscaled weights) ===")
    stages = sorted(set(mem.get("bar_mc_totals", {})) | set(chk.get("bar_mc_totals", {})))
    print(f"{'stage':24s} {'in-memory':>12s} {'chunk':>12s} {'diff':>10s}")
    for sk in stages:
        a = mem.get("bar_mc_totals", {}).get(sk, float("nan"))
        b = chk.get("bar_mc_totals", {}).get(sk, float("nan"))
        diff = b - a if np.isfinite(a) and np.isfinite(b) else float("nan")
        print(f"{sk:24s} {a:12.1f} {b:12.1f} {diff:10.1f}")

    pw_diff = abs(mem.get("purity_weighted", 0) - chk.get("purity_weighted", 0))
    pr_diff = abs(mem.get("purity_raw", 0) - chk.get("purity_raw", 0))
    print(f"\n[parity] |Δ purity_weighted| = {pw_diff:.6f}%  |Δ purity_raw| = {pr_diff:.6f}%")
    if pw_diff > 0.01 or pr_diff > 0.01:
        print("[parity] WARN: purity mismatch between engines")
        sys.exit(1)
    print("[parity] OK: in-memory and chunk paths agree on this shard")


if __name__ == "__main__":
    main()
