#!/usr/bin/env python
"""Parallel dense-rate histcounts stream-sum (shard → merge → covs).

The serial ``syst_histcounts_stream_sum.py`` is I/O bound at ~1 file/min.
This driver shards the file list across workers, each writing a partial
dense pickle, then merges with :func:`add_rate_dense_inplace` (associative).

Can seed from an in-progress serial checkpoint::

    --seed-dense PATH.dense.pkl --seed-done-json PATH.stream_state.json

so a slow serial run can be cut over without redoing finished files.

Example (Flux remaining + seed)::

    python syst_histcounts_stream_sum_parallel.py \\
      --input-glob '.../dfs/hist_mc_flux/*.df' \\
      --out-df '.../summed/hist_mc_flux__fullstat.df' \\
      --family Flux --workers 16 \\
      --seed-dense '.../summed/hist_mc_flux__fullstat.df.stream_state.json.dense.pkl' \\
      --seed-done-json '.../summed/hist_mc_flux__fullstat.df.stream_state.json' \\
      --build-covs --skip-hist-df
"""
from __future__ import annotations

import argparse
import gc
import glob
import json
import multiprocessing as mp
import os
import pickle
import shutil
import sys
import tempfile
import time
import traceback
from datetime import datetime
from os import path
from typing import Any, Dict, List, Optional, Sequence, Set

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.syst_histcounts import (  # noqa: E402
    add_rate_dense_inplace,
    combine_indep_knob_frac_covs,
    histcounts_var_configs_df,
    load_syst_hists_from_df_file,
    load_var_configs_from_df_file,
    nbins_by_var_from_configs,
    rate_cov_from_univ_cv,
    unpack_rate_from_df,
)

DEFAULT_SLIM_SKIP = frozenset(
    {
        "flux_total",
        "g4_total",
        "Flux",
        "G4",
        "slim",
        "slim_multisim",
        "Flux_slim",
        "Flux_slim_multisim",
        "G4_slim",
        "G4_slim_multisim",
    }
)


def parse_cli(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input-glob", required=True)
    p.add_argument("--out-df", required=True)
    p.add_argument("--family", default="Flux")
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--max-files", type=int, default=0)
    p.add_argument(
        "--seed-dense",
        default="",
        help="Existing dense.pkl to start from (e.g. serial checkpoint)",
    )
    p.add_argument(
        "--seed-done-json",
        default="",
        help="stream_state.json with 'done' list of already-summed abspaths",
    )
    p.add_argument("--build-covs", action="store_true")
    p.add_argument("--skip-hist-df", action="store_true", default=True)
    p.add_argument("--slim-skip", default="")
    p.add_argument("--progress-every", type=int, default=5)
    p.add_argument(
        "--partials-dir",
        default="",
        help="Directory for shard pickles (default: <out-df>.partials/)",
    )
    p.add_argument(
        "--gc-every",
        type=int,
        default=20,
        help="gc.collect every N files per worker (0=never; default 20)",
    )
    return p.parse_args(argv)


def _build_rate_covs(
    rate: Dict[str, Dict[str, Dict[str, Any]]],
    *,
    family: str,
    slim_skip: Set[str],
) -> Dict[str, Dict[str, Any]]:
    from collections import defaultdict

    by_var: Dict[str, Dict[str, Any]] = defaultdict(dict)
    packs_by_var: Dict[str, List[Any]] = defaultdict(list)
    cv_by_var: Dict[str, Any] = {}
    n_slug = 0
    for knob, vars_d in rate.items():
        for slug, pack in vars_d.items():
            univ = pack["univ"]
            if univ is None or getattr(univ, "size", 0) == 0:
                continue
            rp = rate_cov_from_univ_cv(univ, pack["cv"])
            by_var[slug][knob] = rp
            if knob not in slim_skip:
                packs_by_var[slug].append(rp)
            cv_by_var[slug] = pack["cv"]
            n_slug += 1
            if n_slug % 50 == 0:
                print("[stream-sum-par] cov progress pairs=%d" % n_slug, flush=True)
    key = family.lower()
    for slug, packs in packs_by_var.items():
        by_var[slug][key] = combine_indep_knob_frac_covs(packs, cv_by_var[slug])
    return dict(by_var)


def _worker(job: dict) -> dict:
    """Sum one shard of files → partial dense pickle."""
    shard_id = int(job["shard_id"])
    files: List[str] = job["files"]
    out_pkl = job["out_pkl"]
    family = job["family"]
    nbins_by_var = job["nbins_by_var"]
    progress_every = int(job["progress_every"])
    gc_every = int(job["gc_every"])
    t0 = time.time()
    acc: Dict[str, Dict[str, Dict[str, Any]]] = {}
    done: List[str] = []
    try:
        for i, p in enumerate(files, start=1):
            df = load_syst_hists_from_df_file(p)
            piece = unpack_rate_from_df(df, family=family, nbins_by_var=nbins_by_var)
            del df
            add_rate_dense_inplace(acc, piece)
            del piece
            done.append(path.abspath(p))
            if gc_every > 0 and (i % gc_every == 0):
                gc.collect()
            if progress_every > 0 and (i % progress_every == 0 or i == len(files)):
                rate = i / max(time.time() - t0, 1e-6)
                print(
                    "[stream-sum-par] shard %d  %d/%d  rate=%.1f/file/h  %s"
                    % (
                        shard_id,
                        i,
                        len(files),
                        rate * 3600.0,
                        path.basename(p),
                    ),
                    flush=True,
                )
        # Write locally first (pnfs pickle dumps are slow / flaky).
        local_dir = tempfile.mkdtemp(prefix="histcounts_partial_")
        local_pkl = path.join(local_dir, path.basename(out_pkl))
        try:
            with open(local_pkl, "wb") as fh:
                pickle.dump(
                    {"acc": acc, "done": done, "shard_id": shard_id},
                    fh,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            os.makedirs(path.dirname(out_pkl) or ".", exist_ok=True)
            shutil.copy2(local_pkl, out_pkl)
        finally:
            shutil.rmtree(local_dir, ignore_errors=True)
        return {
            "ok": True,
            "shard_id": shard_id,
            "out_pkl": out_pkl,
            "n_files": len(done),
            "elapsed": time.time() - t0,
        }
    except BaseException as e:
        return {
            "ok": False,
            "shard_id": shard_id,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=12)),
            "elapsed": time.time() - t0,
        }


def _load_seed(seed_dense: str, seed_done_json: str) -> tuple[dict, Set[str]]:
    acc: dict = {}
    done: Set[str] = set()
    if seed_dense and path.isfile(seed_dense):
        with open(seed_dense, "rb") as fh:
            acc = pickle.load(fh)
        print(
            "[stream-sum-par] loaded seed dense %s  knobs=%d"
            % (seed_dense, len(acc)),
            flush=True,
        )
    if seed_done_json and path.isfile(seed_done_json):
        with open(seed_done_json, "r") as fh:
            state = json.load(fh)
        done = {path.abspath(p) for p in state.get("done", [])}
        print(
            "[stream-sum-par] seed done list: %d files from %s"
            % (len(done), seed_done_json),
            flush=True,
        )
    return acc, done


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_cli(argv)
    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit("[stream-sum-par] no files match %r" % args.input_glob)
    if args.max_files > 0:
        files = files[: int(args.max_files)]

    # Load only the done-list before forking (NOT the dense seed) so workers
    # do not inherit a multi-GB CoW copy of the accumulator.
    done: Set[str] = set()
    if args.seed_done_json and path.isfile(args.seed_done_json):
        with open(args.seed_done_json, "r") as fh:
            state = json.load(fh)
        done = {path.abspath(p) for p in state.get("done", [])}
        print(
            "[stream-sum-par] seed done list: %d files from %s"
            % (len(done), args.seed_done_json),
            flush=True,
        )
    pending = [p for p in files if path.abspath(p) not in done]
    print(
        "[stream-sum-par] %s  total=%d  seed_done=%d  pending=%d  workers=%d  family=%s"
        % (
            datetime.now().isoformat(timespec="seconds"),
            len(files),
            len(done),
            len(pending),
            args.workers,
            args.family,
        ),
        flush=True,
    )

    nbins_by_var = None
    probe = next(iter(files), None)
    if probe is not None:
        try:
            cfg = load_var_configs_from_df_file(probe)
            if cfg:
                nbins_by_var = nbins_by_var_from_configs(cfg)
        except Exception as e:
            print("[stream-sum-par] warn var_configs: %s" % e, flush=True)

    partials_dir = args.partials_dir or (args.out_df + ".partials")
    os.makedirs(partials_dir, exist_ok=True)
    os.makedirs(path.dirname(path.abspath(args.out_df)) or ".", exist_ok=True)

    n_workers = max(1, int(args.workers))
    shards: List[List[str]] = [[] for _ in range(n_workers)]
    for i, p in enumerate(pending):
        shards[i % n_workers].append(p)
    # Drop empty shards
    shard_jobs = []
    for sid, flist in enumerate(shards):
        if not flist:
            continue
        shard_jobs.append(
            {
                "shard_id": sid,
                "files": flist,
                "out_pkl": path.join(partials_dir, "shard_%03d.dense.pkl" % sid),
                "family": args.family,
                "nbins_by_var": nbins_by_var,
                "progress_every": int(args.progress_every),
                "gc_every": int(args.gc_every),
            }
        )

    partial_paths: List[str] = []
    if shard_jobs:
        print(
            "[stream-sum-par] launching %d shards (%d pending files)"
            % (len(shard_jobs), len(pending)),
            flush=True,
        )
        # fork without seed dense in parent → workers stay lean
        ctx = mp.get_context("fork")
        fails = 0
        with ctx.Pool(processes=len(shard_jobs), maxtasksperchild=1) as pool:
            for res in pool.imap_unordered(_worker, shard_jobs):
                if not res["ok"]:
                    fails += 1
                    print(
                        "[stream-sum-par] FAIL shard %s: %s"
                        % (res.get("shard_id"), str(res.get("err", ""))[:400]),
                        flush=True,
                    )
                else:
                    partial_paths.append(res["out_pkl"])
                    print(
                        "[stream-sum-par] shard %s ok  n=%d  elapsed=%.1fs"
                        % (res["shard_id"], res["n_files"], res["elapsed"]),
                        flush=True,
                    )
        if fails:
            raise SystemExit("[stream-sum-par] %d shard(s) failed" % fails)
    else:
        print("[stream-sum-par] nothing pending — seed only", flush=True)

    # Load seed dense only now (after workers exit) and merge partials one-by-one.
    acc, _ = _load_seed(args.seed_dense, "")
    print(
        "[stream-sum-par] merging %d partials into seed…" % len(partial_paths),
        flush=True,
    )
    t_merge = time.time()
    for i, pp in enumerate(sorted(partial_paths), start=1):
        with open(pp, "rb") as fh:
            blob = pickle.load(fh)
        add_rate_dense_inplace(acc, blob["acc"])
        done.update(path.abspath(x) for x in blob.get("done", []))
        del blob
        gc.collect()
        print(
            "[stream-sum-par] merged %d/%d  %s" % (i, len(partial_paths), path.basename(pp)),
            flush=True,
        )
    print(
        "[stream-sum-par] merge done in %.1fs  knobs=%d  done_files=%d"
        % (time.time() - t_merge, len(acc), len(done)),
        flush=True,
    )

    final_dense = args.out_df + ".dense.pkl"
    local_dir = tempfile.mkdtemp(prefix="histcounts_final_")
    try:
        local_pkl = path.join(local_dir, path.basename(final_dense))
        with open(local_pkl, "wb") as fh:
            pickle.dump(acc, fh, protocol=pickle.HIGHEST_PROTOCOL)
        shutil.copy2(local_pkl, final_dense)
    finally:
        shutil.rmtree(local_dir, ignore_errors=True)
    print("[stream-sum-par] wrote dense %s" % final_dense, flush=True)

    state_path = args.out_df + ".stream_state.json"
    state_payload = {
        "mode": "dense-rate",
        "family": args.family,
        "done": sorted(done),
        "updated": datetime.now().isoformat(timespec="seconds"),
        "n_done": len(done),
        "parallel": True,
    }
    # dCache/pnfs often rejects open("w") on an existing file; write tmp then replace.
    state_tmp = state_path + ".tmp"
    with open(state_tmp, "w") as fh:
        json.dump(state_payload, fh, indent=2)
    try:
        os.replace(state_tmp, state_path)
    except OSError:
        try:
            os.remove(state_path)
        except OSError:
            pass
        shutil.copy2(state_tmp, state_path)
        try:
            os.remove(state_tmp)
        except OSError:
            pass
    print("[stream-sum-par] wrote state %s  n_done=%d" % (state_path, len(done)), flush=True)

    if args.build_covs:
        slim_skip = set(DEFAULT_SLIM_SKIP)
        if args.slim_skip.strip():
            slim_skip |= {s.strip() for s in args.slim_skip.split(",") if s.strip()}
        cov_path = args.out_df + ".rate_covs.pkl"
        print("[stream-sum-par] building rate covs…", flush=True)
        covs = _build_rate_covs(acc, family=args.family, slim_skip=slim_skip)
        with open(cov_path, "wb") as fh:
            pickle.dump(
                {"family": args.family, "by_var": covs, "n_files": len(done)},
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        print(
            "[stream-sum-par] wrote covs %s  vars=%d" % (cov_path, len(covs)),
            flush=True,
        )

    # Touch out-df path as a marker (skip long hist by default).
    if args.skip_hist_df:
        open(args.out_df + ".SKIP_HIST_DF", "w").close()
    print("[stream-sum-par] DONE family=%s" % args.family, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
