#!/usr/bin/env python
"""Parallel dispatcher for the multisim (MCstat / Flux / G4) chunk-map phase.

This replaces the serial ``for f in files; do python syst_multisim_chunk.py ...; done``
loop in :mod:`run_syst_multisim_chunked.sh`. The dispatcher:

* Imports ``syst_multisim_chunk`` **once** in the master process; workers ``fork``
  from it (cheap copy-on-write on Linux), so per-file Python-startup +
  import cost is paid only ``N_workers`` times, not ``N_files`` times.
* Calls :func:`syst_multisim_chunk.run_with_args` directly inside each worker
  (no subprocess re-launch per file).
* Writes the **same per-file pickles** as the serial path, so:
    - resume/skip-existing semantics are preserved bit-for-bit,
    - the aggregator (:mod:`syst_multisim_aggregate`) is unchanged,
    - one bad file cannot corrupt another (process-level isolation +
      atomic ``.tmp`` → rename inside the chunk writer).
* Logs the same ``[multisim-run] progress map overall ...`` lines as the bash
  loop, so existing log scrapers / dashboards keep working.

CLI flags mirror the relevant ``run_syst_multisim_chunked.sh`` env vars so the
shell driver can dispatch a single Python invocation for the entire map phase.
"""
from __future__ import annotations

# Cap BLAS / OpenMP thread fan-out BEFORE numpy / pandas / etc. import.
# Each worker process otherwise spawns ``N_CPU`` threads which thrashes the
# scheduler when we run many workers in parallel.
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
from analysis_village.numucc_1p0pi.scripts import syst_multisim_chunk as chunk_mod


# ---------------------------------------------------------------------------
# Job construction
# ---------------------------------------------------------------------------
def _parse_only_systs(spec: str) -> Tuple[str, ...]:
    """Mirror of the bash ``MULTISIM_SYST_TYPES`` normalization.

    ``all`` → Flux+G4 (legacy default — MCstat opt-in).
    ``full`` → MCstat+Flux+G4.
    Comma-separated → that subset (validated, re-ordered to NEUTRINO_SYST_ORDER).
    """
    spec = (spec or "all").strip().lower()
    if spec == "all" or spec == "":
        return ("Flux", "G4")
    if spec == "full":
        return tuple(NEUTRINO_SYST_ORDER)
    raw = [x.strip() for x in spec.split(",") if x.strip()]
    norm: List[str] = []
    for s in raw:
        u = s.upper().replace("-", "_")
        if u in ("MCSTAT", "MC_STAT", "MC"):
            norm.append("MCstat")
        elif u == "FLUX":
            norm.append("Flux")
        elif u == "G4":
            norm.append("G4")
        else:
            raise SystemExit(
                "[multisim-parallel] ERROR: unknown syst type %r in --syst-types %r" % (s, spec)
            )
    seen = set()
    uniq = []
    for s in norm:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return tuple(sn for sn in NEUTRINO_SYST_ORDER if sn in uniq)


def _chunk_out_for(syst: str, *, chunks_multisim: str, chunks_mcstat: str,
                   chunks_flux: str, chunks_g4: str) -> str:
    if syst == "COMBINED":
        return path.join(chunks_multisim, "Combined")
    if syst == "MCstat":
        return chunks_mcstat
    if syst == "Flux":
        return chunks_flux
    if syst == "G4":
        return chunks_g4
    raise ValueError("unknown syst route %r" % syst)


def _syst_names_arg_for(
    syst: str,
    only_systs: Tuple[str, ...],
    all_three: bool,
) -> str:
    """Value to pass for ``--syst-names`` / Namespace.syst_names.

    - ``COMBINED`` + all three  → ``"full"`` (chunk script expands to MCstat,Flux,G4).
    - ``COMBINED`` + subset     → the explicit CSV (e.g. ``"Flux,G4"``).
    - Per-syst job              → just that one name.
    """
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
    """Yield list of (syst, df_file, chunk_out_dir, syst_names_arg) jobs.

    Filtering matches the legacy shell logic exactly: ``COMBINED`` rows are
    always emitted; per-syst rows only when that systematic is in
    ``only_systs`` (or when running all three).
    """
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


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------
def _worker(job: dict) -> dict:
    """Run a single (syst, df_file) chunk-map job in a worker process.

    ``job`` keys: ``syst``, ``df_file``, ``out_dir``, ``syst_names``, plus the
    chunk-script knobs ``input_stage``, ``var_set``, ``n_universe``,
    ``max_splits``, ``g4_mode``, ``flux_mode``, ``flux_knob_groups``.
    Returns a dict the master uses for progress + failure logging.
    """
    syst = job["syst"]
    df_file = job["df_file"]
    started = time.time()
    try:
        ns = argparse.Namespace(
            df_file=df_file,
            out_dir=job["out_dir"],
            input_stage=job["input_stage"],
            var_set=job["var_set"],
            syst_names=job["syst_names"] or None,
            n_universe=job["n_universe"],
            max_splits=job["max_splits"],
            g4_mode=job["g4_mode"],
            flux_mode=job["flux_mode"],
            flux_knob_groups=job["flux_knob_groups"],
        )
        out_path, status = chunk_mod.run_with_args(ns, skip_existing=True)
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


# ---------------------------------------------------------------------------
# CLI / driver
# ---------------------------------------------------------------------------
def parse_cli(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--var-set", choices=("final", "intermediate", "both", "sel_all"), default="final")
    p.add_argument("--syst-types", default="all",
                   help="all | full | comma list (MCstat,Flux,G4). Same semantics as the bash driver.")
    p.add_argument("--max-files", type=int, default=0,
                   help="Cap total map jobs after syst filtering (0 = no cap).")
    p.add_argument("--workers", type=int, default=8,
                   help="Number of worker processes. Capped to len(jobs).")
    p.add_argument("--chunks-multisim", required=True,
                   help="Chunk root for COMBINED jobs (``<root>/Combined/`` will receive pickles).")
    p.add_argument("--chunks-mcstat", required=True)
    p.add_argument("--chunks-flux", required=True)
    p.add_argument("--chunks-g4", required=True)
    p.add_argument("--failed-log", required=True,
                   help="Append-only log; one row per failed (syst, df) pair.")
    p.add_argument("--g4-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument("--flux-mode", choices=("knobs", "bundled"), default="knobs")
    p.add_argument("--flux-knob-groups", default="all")
    p.add_argument("--n-universe", type=int, default=100)
    p.add_argument("--max-splits", type=int, default=0)
    return p.parse_args(argv)


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _failed_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main(argv: Optional[List[str]] = None) -> int:
    cli = parse_cli(argv)
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
        "[multisim-parallel] mc_df_stage=%s only_systs=%s all_three=%s "
        "var_set=%s g4_mode=%s flux_mode=%s flux_knob_groups=%s"
        % (cli.mc_df_stage, only_systs, all_three, cli.var_set,
           cli.g4_mode, cli.flux_mode, cli.flux_knob_groups)
    )
    print(
        "[multisim-parallel] chunk-map queue: %d (syst, .df) job(s); workers=%d"
        % (total, cli.workers)
    )
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
            "var_set": cli.var_set,
            "n_universe": int(cli.n_universe),
            "max_splits": int(cli.max_splits),
            "g4_mode": cli.g4_mode,
            "flux_mode": cli.flux_mode,
            "flux_knob_groups": cli.flux_knob_groups,
        }
        for (s, p, d, sn) in jobs
    ]

    done = 0
    failures = 0
    # ``maxtasksperchild`` prevents long-lived workers from accumulating pandas
    # memory growth over hundreds of files (RSS would otherwise climb steadily).
    # 32 is a balance: amortises import cost vs. periodic worker recycling.
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers, maxtasksperchild=32) as pool, open(cli.failed_log, "a") as flog:
        for res in pool.imap_unordered(_worker, payloads, chunksize=1):
            done += 1
            syst = res["syst"]
            df = res["df_file"]
            if res["ok"]:
                tag = "(skip existing)" if res["status"] == "skipped" else "END"
                print(
                    "[multisim-parallel] progress map overall %d/%d  syst=%s  %s %s "
                    "elapsed=%.1fs df_file=%s -> %s"
                    % (done, total, syst, tag, _now_iso(), res["elapsed"], df, res.get("out_path"))
                )
            else:
                failures += 1
                print(
                    "[multisim-parallel] progress map overall %d/%d  syst=%s  FAILED %s "
                    "elapsed=%.1fs df_file=%s\n%s"
                    % (done, total, syst, _now_iso(), res["elapsed"], df, res["err"]),
                    file=sys.stderr,
                )
                flog.write("%s\t%s\t%s\t%s\n" % (_failed_ts(), cli.mc_df_stage, syst, df))
                flog.flush()
    print(
        "[multisim-parallel] map DONE total=%d failures=%d failed_log=%s"
        % (total, failures, cli.failed_log)
    )
    # Same convention as the bash driver: non-zero exit only when **all** jobs failed,
    # so partial completions still let the aggregate phase run on what we have.
    if failures and failures == total:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
