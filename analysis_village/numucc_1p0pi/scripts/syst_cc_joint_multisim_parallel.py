#!/usr/bin/env python
"""Parallel map dispatcher for :mod:`syst_cc_joint_multisim_chunk` (joint-bin multisim CC).

Mirrors :mod:`syst_multisim_parallel` job construction (same ``iter_multisim_chunk_tasks`` rows
and per-syst output routing) but each worker calls :func:`syst_cc_joint_multisim_chunk.run_with_args`
so one pickle carries stacked ``(n_X+n_Y)`` rates per kinematic pair.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import multiprocessing as mp
import sys
import time
import traceback
from datetime import datetime
from os import path
from typing import List, Optional, Tuple

sys.path.append(
    path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
)

from analysis_village.numucc_1p0pi.dataset_locations import iter_multisim_chunk_tasks
from analysis_village.numucc_1p0pi.syst_multisim_common import NEUTRINO_SYST_ORDER
from analysis_village.numucc_1p0pi.scripts import syst_cc_joint_multisim_chunk as joint_chunk_mod


def _parse_only_systs(spec: str) -> Tuple[str, ...]:
    spec = (spec or "all").strip().lower()
    if spec == "all" or spec == "":
        return ("Flux", "G4")
    if spec == "full":
        return tuple(NEUTRINO_SYST_ORDER)
    raw = [x.strip() for x in spec.split(",") if x.strip()]
    norm: List[str] = []
    for token in raw:
        u = token.upper().replace("-", "_")
        if u in ("MCSTAT", "MC_STAT", "MC"):
            norm.append("MCstat")
        elif u == "FLUX":
            norm.append("Flux")
        elif u == "G4":
            norm.append("G4")
        else:
            raise SystemExit(
                "[cc-joint-multisim-parallel] ERROR: unknown syst type %r in --syst-types %r"
                % (token, spec)
            )
    seen = set()
    uniq = []
    for s in norm:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return tuple(sn for sn in NEUTRINO_SYST_ORDER if sn in uniq)


def _chunk_out_for(
    syst: str,
    *,
    chunks_multisim: str,
    chunks_mcstat: str,
    chunks_flux: str,
    chunks_g4: str,
) -> str:
    if syst == "COMBINED":
        return path.join(chunks_multisim, "Combined")
    if syst == "MCstat":
        return chunks_mcstat
    if syst == "Flux":
        return chunks_flux
    if syst == "G4":
        return chunks_g4
    raise ValueError("unknown syst route %r" % syst)


def _syst_names_arg_for(syst: str, only_systs: Tuple[str, ...], all_three: bool) -> str:
    if syst == "COMBINED":
        if all_three:
            return "full"
        return ",".join(only_systs)
    return syst


def _build_jobs(
    *,
    mc_df_stage: str,
    only_systs: Tuple[str, ...],
    all_three: bool,
    max_files: int,
    chunks_multisim: str,
    chunks_mcstat: str,
    chunks_flux: str,
    chunks_g4: str,
) -> List[Tuple[str, str, str, str]]:
    only_set = set(only_systs)
    raw: List[Tuple[str, str, str, str]] = []
    for syst, df_path in iter_multisim_chunk_tasks(mc_df_stage):
        if syst != "COMBINED" and not all_three and syst not in only_set:
            continue
        out_dir = _chunk_out_for(
            syst,
            chunks_multisim=chunks_multisim,
            chunks_mcstat=chunks_mcstat,
            chunks_flux=chunks_flux,
            chunks_g4=chunks_g4,
        )
        sn_arg = _syst_names_arg_for(syst, only_systs, all_three)
        raw.append((syst, df_path, out_dir, sn_arg))
    if max_files > 0 and max_files < len(raw):
        raw = raw[:max_files]
    return raw


def _worker(job: dict) -> dict:
    syst = job["syst"]
    df_file = job["df_file"]
    started = time.time()
    try:
        ns = argparse.Namespace(
            df_file=df_file,
            out_dir=job["out_dir"],
            input_stage=job["input_stage"],
            syst_names=job["syst_names"] or None,
            n_universe=int(job["n_universe"]),
            max_splits=int(job["max_splits"]),
            g4_mode=job["g4_mode"],
            flux_mode=job["flux_mode"],
            flux_knob_groups=job["flux_knob_groups"],
            pairs=job.get("pairs"),
        )
        out_path, status = joint_chunk_mod.run_with_args(ns, skip_existing=True)
        return {
            "ok": True,
            "status": status,
            "syst": syst,
            "df_file": df_file,
            "out_path": out_path,
            "elapsed": time.time() - started,
        }
    except SystemExit as e:
        return {
            "ok": False,
            "status": "failed",
            "syst": syst,
            "df_file": df_file,
            "err": "SystemExit: %s" % e,
            "elapsed": time.time() - started,
        }
    except BaseException as e:
        return {
            "ok": False,
            "status": "failed",
            "syst": syst,
            "df_file": df_file,
            "err": "%s: %s\n%s" % (type(e).__name__, e, traceback.format_exc(limit=8)),
            "elapsed": time.time() - started,
        }


def parse_cli(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--syst-types", default="all", help="all | full | comma list (MCstat,Flux,G4).")
    p.add_argument("--max-files", type=int, default=0, help="Cap map jobs after syst filter (0 = all).")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--chunks-multisim", required=True)
    p.add_argument("--chunks-mcstat", required=True)
    p.add_argument("--chunks-flux", required=True)
    p.add_argument("--chunks-g4", required=True)
    p.add_argument("--failed-log", required=True)
    p.add_argument("--g4-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument("--flux-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument("--flux-knob-groups", default="all")
    p.add_argument("--n-universe", type=int, default=100)
    p.add_argument("--max-splits", type=int, default=0)
    p.add_argument(
        "--pairs",
        default=None,
        help="Optional CSV of pair slugs for syst_cc_joint_multisim_chunk (--pairs).",
    )
    return p.parse_args(argv)


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _failed_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main(argv: Optional[List[str]] = None) -> int:
    cli = parse_cli(argv)
    if cli.mc_df_stage != "final":
        print(
            "[cc-joint-multisim-parallel] ERROR: only --mc-df-stage final is supported "
            "(joint chunk does not implement sel_all).",
            file=sys.stderr,
        )
        return 2
    only_systs = _parse_only_systs(cli.syst_types)
    all_three = set(only_systs) == set(NEUTRINO_SYST_ORDER)

    for d in (
        path.join(cli.chunks_multisim, "Combined"),
        cli.chunks_mcstat,
        cli.chunks_flux,
        cli.chunks_g4,
    ):
        os.makedirs(d, exist_ok=True)
    os.makedirs(path.dirname(cli.failed_log) or ".", exist_ok=True)

    jobs = _build_jobs(
        mc_df_stage=cli.mc_df_stage,
        only_systs=only_systs,
        all_three=all_three,
        max_files=int(cli.max_files or 0),
        chunks_multisim=cli.chunks_multisim,
        chunks_mcstat=cli.chunks_mcstat,
        chunks_flux=cli.chunks_flux,
        chunks_g4=cli.chunks_g4,
    )
    total = len(jobs)
    print(
        "[cc-joint-multisim-parallel] mc_df_stage=%s only_systs=%s all_three=%s "
        "g4_mode=%s flux_mode=%s flux_knob_groups=%s"
        % (cli.mc_df_stage, only_systs, all_three, cli.g4_mode, cli.flux_mode, cli.flux_knob_groups)
    )
    print("[cc-joint-multisim-parallel] chunk-map queue: %d job(s); workers=%d" % (total, cli.workers))
    if total == 0:
        return 0

    workers = max(1, min(int(cli.workers), total))
    payloads = [
        {
            "syst": s,
            "df_file": p,
            "out_dir": d,
            "syst_names": sn,
            "input_stage": cli.mc_df_stage,
            "n_universe": int(cli.n_universe),
            "max_splits": int(cli.max_splits),
            "g4_mode": cli.g4_mode,
            "flux_mode": cli.flux_mode,
            "flux_knob_groups": cli.flux_knob_groups,
            "pairs": cli.pairs,
        }
        for (s, p, d, sn) in jobs
    ]

    done = 0
    failures = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers, maxtasksperchild=32) as pool, open(cli.failed_log, "a") as flog:
        for res in pool.imap_unordered(_worker, payloads, chunksize=1):
            done += 1
            syst = res["syst"]
            df = res["df_file"]
            if res["ok"]:
                tag = "(skip existing)" if res["status"] == "skipped" else "END"
                print(
                    "[cc-joint-multisim-parallel] progress map overall %d/%d  syst=%s  %s %s "
                    "elapsed=%.1fs df_file=%s -> %s"
                    % (done, total, syst, tag, _now_iso(), res["elapsed"], df, res.get("out_path"))
                )
            else:
                failures += 1
                print(
                    "[cc-joint-multisim-parallel] progress map overall %d/%d  syst=%s  FAILED %s "
                    "elapsed=%.1fs df_file=%s\n%s"
                    % (done, total, syst, _now_iso(), res["elapsed"], df, res["err"]),
                    file=sys.stderr,
                )
                flog.write("%s\t%s\t%s\t%s\n" % (_failed_ts(), cli.mc_df_stage, syst, df))
                flog.flush()
    print(
        "[cc-joint-multisim-parallel] map DONE total=%d failures=%d failed_log=%s"
        % (total, failures, cli.failed_log)
    )
    if failures and failures == total:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
