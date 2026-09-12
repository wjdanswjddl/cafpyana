#!/usr/bin/env python
"""Stream-sum per-file histcounts → one full-stat hist (never load all files).

Two backends
------------
* ``dense-rate`` (default) — unpack each file to ``{cv, univ}`` arrays, add in
  place, pack once at the end. Best for Flux / G4 rate-only campaigns; keeps
  resident set near one file + dense acc (~tens of GB peak during load).
* ``long`` — pairwise ``sum_histcounts_dfs([acc, next])`` for any kind
  (GENIE xsec rows included). Peak ~2× one long table.

Also writes optional rate covariances (per-knob + combined) for dense-rate.

Examples
--------
Flux campaign (resume-friendly)::

    cd /exp/sbnd/app/users/munjung/xsec/freeze/cafpyana
    source envs/venv_py310_cafpyana/bin/activate
    export PYTHONPATH=\"$PWD:${PYTHONPATH:-}\"

    nohup python analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum.py \\
      --input-glob '/pnfs/sbnd/scratch/users/munjung/cafpyana_tmp/syst_histcounts_from_df_2026_09_07_140239/dfs/hist_mc_flux/*.df' \\
      --out-df /pnfs/sbnd/scratch/users/munjung/cafpyana_tmp/syst_histcounts_from_df_2026_09_07_140239/summed/hist_mc_flux__fullstat.df \\
      --family Flux --mode dense-rate --build-covs \\
      --checkpoint-every 25 \\
      > /tmp/syst_histcounts_stream_sum_flux.log 2>&1 &

Smoke (2 files)::

    python analysis_village/numucc_1p0pi/scripts/syst_histcounts_stream_sum.py \\
      --input-glob '.../hist_mc_flux/*.df' --max-files 2 \\
      --out-df /tmp/flux_sum_smoke.df --family Flux
"""
from __future__ import annotations

import argparse
import gc
import glob
import json
import os
import pickle
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime
from os import path
from typing import Any, Dict, List, Optional, Sequence, Set

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.syst_histcounts import (
    add_rate_dense_inplace,
    combine_indep_knob_frac_covs,
    histcounts_var_configs_df,
    load_syst_hists_from_df_file,
    load_var_configs_from_df_file,
    nbins_by_var_from_configs,
    rate_cov_from_univ_cv,
    rate_dense_to_hist_df,
    sum_histcounts_dfs,
    unpack_rate_from_df,
)


# Slim / product knobs excluded from combined Flux/G4 totals (same spirit as notebook).
DEFAULT_SLIM_SKIP = frozenset(
    {
        "flux_total",
        "g4_total",
        "Flux",
        "G4",
    }
)


def parse_cli(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--input-glob",
        required=True,
        help="Glob of per-file histcounts .df (quote the glob)",
    )
    p.add_argument(
        "--out-df",
        required=True,
        help="Output path for summed histcounts HDF (.df)",
    )
    p.add_argument(
        "--family",
        default="Flux",
        help="Family label for packing / filtering (Flux, G4, GENIE, …)",
    )
    p.add_argument(
        "--mode",
        choices=("dense-rate", "long"),
        default="dense-rate",
        help="Accumulation backend (default dense-rate)",
    )
    p.add_argument("--max-files", type=int, default=0, help="If >0, only first N files")
    p.add_argument(
        "--checkpoint-every",
        type=int,
        default=50,
        help="Write resume state every N files (0=disable)",
    )
    p.add_argument(
        "--state-file",
        default="",
        help="JSON resume state (default: <out-df>.stream_state.json)",
    )
    p.add_argument(
        "--no-resume",
        action="store_true",
        help="Ignore existing state file and start fresh",
    )
    p.add_argument(
        "--build-covs",
        action="store_true",
        help="After dense-rate sum, write per-var rate cov pickle next to out-df",
    )
    p.add_argument(
        "--skip-hist-df",
        action="store_true",
        help="Do not pack/write the long syst_hists HDF (dense pickle is enough for covs)",
    )
    p.add_argument(
        "--cov-pkl",
        default="",
        help="Cov pickle path (default: <out-df>.rate_covs.pkl)",
    )
    p.add_argument(
        "--slim-skip",
        default="",
        help="Comma-separated knobs to exclude from combined family cov",
    )
    p.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Log progress every N files",
    )
    return p.parse_args(argv)


def _state_path(args: argparse.Namespace) -> str:
    if args.state_file:
        return args.state_file
    return args.out_df + ".stream_state.json"


def _load_state(state_path: str) -> Dict[str, Any]:
    if not path.isfile(state_path):
        return {"done": [], "mode": None, "family": None}
    with open(state_path, "r") as fh:
        return json.load(fh)


def _save_state(state_path: str, state: Dict[str, Any]) -> None:
    os.makedirs(path.dirname(path.abspath(state_path)) or ".", exist_ok=True)
    tmp = state_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
    os.replace(tmp, state_path)


def _checkpoint_dense_path(state_path: str) -> str:
    return state_path + ".dense.pkl"


def _write_out_df(
    out_df: str,
    hist: pd.DataFrame,
    *,
    var_cfg_df: Optional[pd.DataFrame],
) -> None:
    os.makedirs(path.dirname(path.abspath(out_df)) or ".", exist_ok=True)
    local_dir = None
    try:
        # Prefer local temp then copy when targeting pnfs (same lesson as from-df).
        if out_df.startswith("/pnfs/"):
            import shutil
            import tempfile

            local_dir = tempfile.mkdtemp(prefix="histcounts_stream_sum_")
            tmp = path.join(local_dir, path.basename(out_df))
        else:
            tmp = out_df + ".tmp"

        with pd.HDFStore(tmp, mode="w", complevel=5, complib="zlib") as store:
            store.put("syst_hists_0", hist, format="table")
            if var_cfg_df is not None and len(var_cfg_df):
                store.put("var_configs_0", var_cfg_df, format="table")
            store.put("split", pd.DataFrame({"n_split": [1]}), format="table")

        if tmp != out_df:
            import shutil

            shutil.copy2(tmp, out_df)
            try:
                os.remove(tmp)
            except OSError:
                pass
        else:
            os.replace(tmp, out_df)
    finally:
        if local_dir is not None:
            import shutil

            shutil.rmtree(local_dir, ignore_errors=True)


def _build_rate_covs(
    rate: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    *,
    family: str,
    slim_skip: Set[str],
) -> Dict[str, Dict[str, Any]]:
    by_var: Dict[str, Dict[str, Any]] = defaultdict(dict)
    packs_by_var: Dict[str, List[Any]] = defaultdict(list)
    cv_by_var: Dict[str, np.ndarray] = {}
    n_slug = 0
    for knob, vars_d in rate.items():
        for slug, pack in vars_d.items():
            univ = np.asarray(pack["univ"], dtype=float)
            if univ.size == 0:
                continue
            rp = rate_cov_from_univ_cv(univ, pack["cv"])
            by_var[slug][knob] = rp
            if knob not in slim_skip:
                packs_by_var[slug].append(rp)
            cv_by_var[slug] = np.asarray(pack["cv"], dtype=float)
            n_slug += 1
            if n_slug % 50 == 0:
                print(
                    "[stream-sum] cov progress  knob-var pairs=%d" % n_slug,
                    flush=True,
                )
    key = family.lower()
    for slug, packs in packs_by_var.items():
        by_var[slug][key] = combine_indep_knob_frac_covs(packs, cv_by_var[slug])
    return dict(by_var)


def _run_dense(args: argparse.Namespace, files: List[str]) -> int:
    state_path = _state_path(args)
    dense_ckpt = _checkpoint_dense_path(state_path)
    slim_skip = set(DEFAULT_SLIM_SKIP)
    if args.slim_skip.strip():
        slim_skip |= {s.strip() for s in args.slim_skip.split(",") if s.strip()}

    done: Set[str] = set()
    acc: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    if not args.no_resume and path.isfile(state_path) and path.isfile(dense_ckpt):
        state = _load_state(state_path)
        if state.get("mode") in (None, "dense-rate") and state.get("family") in (
            None,
            args.family,
        ):
            done = set(state.get("done", []))
            with open(dense_ckpt, "rb") as fh:
                acc = pickle.load(fh)
            print(
                "[stream-sum] resume dense-rate  done=%d  acc_knobs=%d"
                % (len(done), len(acc)),
                flush=True,
            )

    # nbins from first available file's var_configs
    nbins_by_var = None
    var_cfg_df = None
    probe = next((p for p in files if True), None)
    if probe is not None:
        try:
            cfg = load_var_configs_from_df_file(probe)
            if cfg:
                nbins_by_var = nbins_by_var_from_configs(cfg)
                var_cfg_df = histcounts_var_configs_df(list(cfg.values()))
        except Exception as e:
            print("[stream-sum] warn: var_configs from %s: %s" % (probe, e), flush=True)

    pending = [p for p in files if path.abspath(p) not in done]
    print(
        "[stream-sum] dense-rate  total=%d  pending=%d  family=%s"
        % (len(files), len(pending), args.family),
        flush=True,
    )

    t0 = time.time()
    for i, p in enumerate(pending, start=1):
        t1 = time.time()
        df = load_syst_hists_from_df_file(p)
        piece = unpack_rate_from_df(
            df, family=args.family, nbins_by_var=nbins_by_var
        )
        del df
        add_rate_dense_inplace(acc, piece)
        del piece
        gc.collect()
        done.add(path.abspath(p))
        dt = time.time() - t1
        if args.progress_every > 0 and (
            i % args.progress_every == 0 or i == len(pending)
        ):
            rate = i / max(time.time() - t0, 1e-6)
            eta_h = (len(pending) - i) / rate / 3600.0 if rate > 0 else float("inf")
            print(
                "[stream-sum] %d/%d  last=%.1fs  rate=%.2f/file/h  ETA~%.1fh  %s"
                % (
                    i,
                    len(pending),
                    dt,
                    rate * 3600.0,
                    eta_h,
                    path.basename(p),
                ),
                flush=True,
            )
        if args.checkpoint_every > 0 and (
            i % args.checkpoint_every == 0 or i == len(pending)
        ):
            with open(dense_ckpt + ".tmp", "wb") as fh:
                pickle.dump(acc, fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(dense_ckpt + ".tmp", dense_ckpt)
            _save_state(
                state_path,
                {
                    "mode": "dense-rate",
                    "family": args.family,
                    "done": sorted(done),
                    "updated": datetime.now().isoformat(timespec="seconds"),
                    "n_done": len(done),
                },
            )
            print(
                "[stream-sum] checkpoint  done=%d  %s" % (len(done), state_path),
                flush=True,
            )

    # Always persist dense accumulator (small vs long hist); resume/covs use this.
    final_dense = args.out_df + ".dense.pkl"
    with open(final_dense + ".tmp", "wb") as fh:
        pickle.dump(acc, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(final_dense + ".tmp", final_dense)
    print("[stream-sum] wrote dense %s  knobs=%d" % (final_dense, len(acc)), flush=True)

    if not args.skip_hist_df:
        print(
            "[stream-sum] packing long hist DF (%d knobs) — can take a while…"
            % len(acc),
            flush=True,
        )
        hist = rate_dense_to_hist_df(acc, family=args.family)
        if var_cfg_df is None:
            var_cfg_df = histcounts_var_configs_df()
        _write_out_df(args.out_df, hist, var_cfg_df=var_cfg_df)
        print(
            "[stream-sum] wrote %s  rows=%d  knobs=%d"
            % (args.out_df, len(hist), len(acc)),
            flush=True,
        )
        del hist
        gc.collect()
    else:
        print("[stream-sum] skip long hist DF (--skip-hist-df)", flush=True)
        # Touch a tiny marker so --out-df path is documented even without HDF.
        os.makedirs(path.dirname(path.abspath(args.out_df)) or ".", exist_ok=True)

    if args.build_covs:
        cov_path = args.cov_pkl or (args.out_df + ".rate_covs.pkl")
        print("[stream-sum] building rate covs…", flush=True)
        covs = _build_rate_covs(acc, family=args.family, slim_skip=slim_skip)
        os.makedirs(path.dirname(path.abspath(cov_path)) or ".", exist_ok=True)
        with open(cov_path, "wb") as fh:
            pickle.dump(
                {"family": args.family, "by_var": covs, "n_files": len(done)},
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        print(
            "[stream-sum] wrote covs %s  vars=%d" % (cov_path, len(covs)),
            flush=True,
        )
    return 0


def _run_long(args: argparse.Namespace, files: List[str]) -> int:
    state_path = _state_path(args)
    acc_path = state_path + ".long_acc.df"
    done: Set[str] = set()
    acc = None
    if not args.no_resume and path.isfile(state_path) and path.isfile(acc_path):
        state = _load_state(state_path)
        if state.get("mode") in (None, "long"):
            done = set(state.get("done", []))
            acc = load_syst_hists_from_df_file(acc_path)
            print(
                "[stream-sum] resume long  done=%d  rows=%d"
                % (len(done), len(acc)),
                flush=True,
            )

    pending = [p for p in files if path.abspath(p) not in done and p not in done]
    print(
        "[stream-sum] long  total=%d  pending=%d" % (len(files), len(pending)),
        flush=True,
    )
    var_cfg_df = None
    if files:
        try:
            cfg = load_var_configs_from_df_file(files[0])
            if cfg:
                var_cfg_df = histcounts_var_configs_df(list(cfg.values()))
        except Exception:
            pass

    t0 = time.time()
    for i, p in enumerate(pending, start=1):
        t1 = time.time()
        df = load_syst_hists_from_df_file(p)
        if acc is None:
            acc = df
        else:
            acc = sum_histcounts_dfs([acc, df])
            del df
            gc.collect()
        done.add(path.abspath(p))
        dt = time.time() - t1
        if args.progress_every > 0 and (
            i % args.progress_every == 0 or i == len(pending)
        ):
            rate = i / max(time.time() - t0, 1e-6)
            eta_h = (len(pending) - i) / rate / 3600.0 if rate > 0 else float("inf")
            print(
                "[stream-sum] %d/%d  rows=%d  last=%.1fs  ETA~%.1fh  %s"
                % (i, len(pending), len(acc), dt, eta_h, path.basename(p)),
                flush=True,
            )
        if args.checkpoint_every > 0 and (
            i % args.checkpoint_every == 0 or i == len(pending)
        ):
            _write_out_df(acc_path, acc, var_cfg_df=var_cfg_df)
            _save_state(
                state_path,
                {
                    "mode": "long",
                    "family": args.family,
                    "done": sorted(done),
                    "updated": datetime.now().isoformat(timespec="seconds"),
                    "n_done": len(done),
                },
            )
            print("[stream-sum] checkpoint  done=%d" % len(done), flush=True)

    if acc is None:
        raise SystemExit("[stream-sum] no input files / empty sum")
    if var_cfg_df is None:
        var_cfg_df = histcounts_var_configs_df()
    _write_out_df(args.out_df, acc, var_cfg_df=var_cfg_df)
    print("[stream-sum] wrote %s  rows=%d" % (args.out_df, len(acc)), flush=True)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_cli(argv)
    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit("[stream-sum] no files match %r" % args.input_glob)
    if args.max_files > 0:
        files = files[: int(args.max_files)]

    print(
        "[stream-sum] %s  files=%d  mode=%s  out=%s"
        % (
            datetime.now().isoformat(timespec="seconds"),
            len(files),
            args.mode,
            args.out_df,
        ),
        flush=True,
    )
    try:
        if args.mode == "dense-rate":
            return _run_dense(args, files)
        return _run_long(args, files)
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
