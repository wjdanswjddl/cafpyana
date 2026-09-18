#!/usr/bin/env python
"""Parallel sel_all cosmics map + aggregate for the Sep-1 updated offbeam/intime DFs.

Writes ``<SYST_DISK>/Cosmics/cosmics_syst_dict.npz`` (Product A cut-stage + Product B
final vars in one NPZ). Skips existing chunk pickles.
"""
from __future__ import annotations

import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana/envs/venv_py310_cafpyana/bin/python")
CHUNK_PY = REPO / "analysis_village/numucc_1p0pi/scripts/syst_cosmics_chunk.py"
AGG_PY = REPO / "analysis_village/numucc_1p0pi/scripts/syst_cosmics_aggregate.py"

SYST_DISK = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final")
WORK = Path(
    os.environ.get(
        "COSMICS_WORK",
        "/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/Cosmics/work-20260918",
    )
)
CHUNKS = WORK / "chunks" / "sel_all"
SELECTED_MC = Path(
    "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/2026_09_01_195328__sel_mup-mc-BNB_cosmics-CV"
)
N_WORKERS = int(os.environ.get("COSMICS_WORKERS", "8"))
LOG = WORK / "launch.log"

ENV = os.environ.copy()
ENV["PYTHONPATH"] = str(REPO) + (os.pathsep + ENV["PYTHONPATH"] if ENV.get("PYTHONPATH") else "")
ENV.setdefault("MPLBACKEND", "Agg")
ENV.setdefault("OMP_NUM_THREADS", "1")
ENV.setdefault("MKL_NUM_THREADS", "1")
ENV.setdefault("OPENBLAS_NUM_THREADS", "1")
ENV.setdefault("NUMEXPR_NUM_THREADS", "1")


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def jobs():
    sys.path.insert(0, str(REPO))
    from analysis_village.numucc_1p0pi.dataset_locations import iter_cosmics_chunk_df_paths

    out = []
    for sample in ("offbeam", "intime"):
        for p in iter_cosmics_chunk_df_paths(sample, input_stage="sel_all"):
            out.append((sample, p))
    return out


def run_one(sample: str, df_file: str) -> tuple[str, str, int]:
    stem = os.path.splitext(os.path.basename(df_file))[0]
    out = CHUNKS / f"cosmics__{sample}__{stem}.pkl"
    if out.is_file():
        return "skip", str(out), 0
    r = subprocess.run(
        [
            str(PY),
            str(CHUNK_PY),
            "--input-stage",
            "sel_all",
            "--sample",
            sample,
            "--df_file",
            df_file,
            "--out_dir",
            str(CHUNKS),
        ],
        env=ENV,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        err = (r.stderr or r.stdout or "")[-2000:]
        return "fail", f"{df_file}\n{err}", r.returncode
    return "ok", str(out), 0


def main() -> int:
    CHUNKS.mkdir(parents=True, exist_ok=True)
    SYST_DISK.mkdir(parents=True, exist_ok=True)
    task_list = jobs()
    log(f"start n_jobs={len(task_list)} workers={N_WORKERS} chunks={CHUNKS}")
    n_ok = n_skip = n_fail = 0
    fails = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = [ex.submit(run_one, s, p) for s, p in task_list]
        done = 0
        for fut in as_completed(futs):
            status, detail, _rc = fut.result()
            done += 1
            if status == "ok":
                n_ok += 1
            elif status == "skip":
                n_skip += 1
            else:
                n_fail += 1
                fails.append(detail)
            if done % 25 == 0 or done == len(futs):
                log(f"map {done}/{len(futs)} ok={n_ok} skip={n_skip} fail={n_fail}")
    if fails:
        fail_log = WORK / "map_failures.log"
        fail_log.write_text("\n\n".join(fails))
        log(f"map failures written to {fail_log}")
        return 1
    log("aggregate begin")
    r = subprocess.run(
        [
            str(PY),
            str(AGG_PY),
            "--chunks_dir",
            str(CHUNKS),
            "--syst-disk-root",
            str(SYST_DISK),
            "--no-plots",
            "--selected-mc-df",
            str(SELECTED_MC),
        ],
        env=ENV,
    )
    log(f"aggregate exit={r.returncode}")
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
