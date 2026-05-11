#!/usr/bin/env python
"""End-to-end workflow test driver for numucc_1p0pi.

Runs, in order:

1. **Multisim map** — ``syst_multisim_chunk.py`` (capped).
2. **Multisim aggregate** — ``syst_multisim_aggregate.py``.
3. **DetVar map** — ``syst_detvar_chunk.py`` per WireMod tag / ``.df`` (see ``dataset_locations.DETVAR_WIREMOD_GLOBS``), capped per tag.
4. **DetVar aggregate** — ``syst_detvar_aggregate.py`` → ``detector_syst_dict.npz``.
5. **Event-selection map** — ``event_selection_chunk.py`` per sample (mc, data, …), capped per type.
6. **Event-selection reduce** — ``event_selection_aggregate.py`` with ``--syst-disk-root`` pointing at
   the unified syst tree (``syst_disk_layout``): ``MCstat/``, ``Flux/``, ``G4/``, ``GENIE/``,
   ``Cosmics/``, ``Detector/``. This driver passes ``--syst-disk-root`` only when multisim included
   cosmics **and** DetVar ran (complete disk covariance tree). Set ``NUMUCC_GENIE_COV_MAT_PKL`` so
   ``GENIE/cov_mat_dict.pkl`` exists before stage 6.

``utils.get_syst_unc`` requires **all** category files under that root (no legacy fallbacks).

Limits:

* ``--max-files N`` (``0`` = no limit): **per** event-selection sample and **per** WireMod tag for
  DetVar; for multisim, at most **N** distinct CAF paths in the multisim task list.

Outputs under ``--output-dir``::

    <output-dir>/
      logs/workflow_<timestamp>.log
      logs/manifest_<timestamp>.json
      multisim/chunks/           # nu__*.pkl shards
      syst_disk/                 # MCstat/, Flux/, G4/, Cosmics/, Detector/, GENIE/ (see syst_disk_layout)
      detvar/chunks/
      event_selection/chunks/  event_selection/plots/

Examples::

    python run_workflow_test.py --output-dir /tmp/nu_test --max-files 1
    python run_workflow_test.py -o ./wf_out --max-files 3 --skip-multisim
    python run_workflow_test.py -o ./wf_out --dry-run
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


# Repo root: analysis_village/numucc_1p0pi/scripts/<this_file> → parents × 4
_SCRIPTS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPTS_DIR.parent.parent.parent.parent


def _ensure_pythonpath() -> None:
    r = str(_REPO_ROOT)
    if r not in sys.path:
        sys.path.insert(0, r)


def _ts_tag() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _seed_genie_cov_pkl(syst_disk: Path, logger: WorkflowLog, dry_run: bool) -> None:
    """Copy ``NUMUCC_GENIE_COV_MAT_PKL`` into ``<syst_disk>/GENIE/cov_mat_dict.pkl`` if set."""
    src_env = os.environ.get("NUMUCC_GENIE_COV_MAT_PKL")
    if not src_env:
        return
    sp = Path(os.path.expanduser(src_env))
    if not sp.is_file():
        logger.line(
            "WARNING: NUMUCC_GENIE_COV_MAT_PKL=%s is not a regular file; skip GENIE seed."
            % src_env
        )
        return
    dst = syst_disk / "GENIE" / "cov_mat_dict.pkl"
    if dry_run:
        logger.line("DRY-RUN: would copy %s -> %s" % (sp, dst))
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(sp, dst)
    logger.line("Seeded %s from NUMUCC_GENIE_COV_MAT_PKL" % dst)


@dataclass
class Manifest:
    started_at: str
    output_dir: str
    max_files: int
    mc_df_stage: str
    var_set: str
    dry_run: bool
    stages: List[Dict[str, Any]] = field(default_factory=list)

    def add_stage(self, name: str, entries: List[Dict[str, Any]]) -> None:
        self.stages.append({"stage": name, "runs": entries})


class WorkflowLog:
    def __init__(self, log_path: Path, manifest_path: Path, manifest: Manifest):
        self.log_path = log_path
        self.manifest_path = manifest_path
        self.manifest = manifest
        self._fp = open(log_path, "w", encoding="utf-8")

    def close(self) -> None:
        self._fp.close()
        self.manifest_path.write_text(json.dumps(asdict(self.manifest), indent=2), encoding="utf-8")

    def line(self, msg: str) -> None:
        self._fp.write(msg.rstrip() + "\n")
        self._fp.flush()

    def banner(self, title: str) -> None:
        bar = "=" * 78
        self.line("")
        self.line(bar)
        self.line(title)
        self.line(bar)

    def subsection(self, title: str) -> None:
        self.line("")
        self.line("--- " + title + " ---")


def _limit_multisim_tasks(
    tasks: Sequence[Tuple[str, str]],
    max_files: int,
) -> List[Tuple[str, str]]:
    """Keep tasks whose ``.df`` path is among the first ``max_files`` unique paths (order preserved)."""
    if max_files <= 0:
        return list(tasks)
    seen_order: List[str] = []
    seen_set: set[str] = set()
    for _, p in tasks:
        if p not in seen_set:
            seen_set.add(p)
            seen_order.append(p)
        if len(seen_order) >= max_files:
            break
    allowed = set(seen_order[:max_files])
    return [(sn, p) for sn, p in tasks if p in allowed]


def _limit_paths(paths: List[str], max_files: int) -> List[str]:
    if max_files <= 0:
        return paths
    return paths[:max_files]


def _limit_detvar_jobs(jobs: List[Tuple[str, str]], max_files: int) -> List[Tuple[str, str]]:
    """Keep at most ``max_files`` CAFs per WireMod tag (see ``DETVAR_WIREMOD_GLOBS`` order)."""
    if max_files <= 0:
        return jobs
    per_tag: Dict[str, int] = {}
    out: List[Tuple[str, str]] = []
    for tag, p in jobs:
        n = per_tag.get(tag, 0)
        if n >= max_files:
            continue
        per_tag[tag] = n + 1
        out.append((tag, p))
    return out


def _run_python(
    script: Path,
    args: List[str],
    *,
    env: Optional[dict] = None,
    dry_run: bool,
    logger: WorkflowLog,
    stage: str,
    meta: Dict[str, Any],
) -> int:
    cmd = [sys.executable, str(script)] + args
    meta["command"] = cmd
    logger.line("")
    logger.line("SCRIPT:     %s" % script.name)
    logger.line("STAGE:      %s" % stage)
    logger.line("COMMAND:    %s" % " ".join(cmd))
    for k, v in meta.items():
        if k in ("command",):
            continue
        logger.line("META[%s]: %s" % (k, v))
    if dry_run:
        logger.line("ACTION:     DRY-RUN (not executed)")
        meta["exit_code"] = None
        meta["dry_run"] = True
        return 0
    r = subprocess.run(cmd, cwd=str(_REPO_ROOT), env=env or os.environ)
    meta["exit_code"] = r.returncode
    logger.line("EXIT_CODE:  %s" % r.returncode)
    return int(r.returncode)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Root directory for multisim + event_selection subtrees and logs/",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=2,
        help="Max .df files per event-selection sample (mc,data,...) and max distinct "
        "CAF paths for multisim map. Use 0 for no limit (full workflow). Default: 2",
    )
    p.add_argument("--mc-df-stage", choices=("final", "sel_all"), default="final")
    p.add_argument("--var-set", choices=("final", "intermediate", "both"), default="final")
    p.add_argument("--skip-multisim", action="store_true", help="Skip multisim (stages 1–2).")
    p.add_argument(
        "--skip-event-selection",
        action="store_true",
        help="Skip event-selection map/aggregate (stages 5–6). Implies no DetVar (nothing consumes it).",
    )
    p.add_argument(
        "--skip-detvar",
        action="store_true",
        help="Skip DetVar map/aggregate (stages 3–4). Omit disk syst covariances unless "
        "NUMUCC_SYST_DISK_ROOT already holds a complete syst_disk_layout tree.",
    )
    p.add_argument(
        "--multisim-include-cosmics",
        action="store_true",
        help="Pass cosmics loading to syst_multisim_aggregate (heavy). Default: skip cosmics.",
    )
    p.add_argument(
        "--multisim-no-plots",
        action="store_true",
        help="syst_multisim_aggregate.py --no-plots (faster).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands only; do not execute.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    _ensure_pythonpath()

    from analysis_village.numucc_1p0pi.dataset_locations import (
        DETVAR_WIREMOD_GLOBS,
        EVENT_SELECTION_GLOBS,
        iter_detvar_chunk_jobs,
        iter_event_selection_df_paths,
        iter_multisim_chunk_tasks,
    )

    out = Path(args.output_dir).resolve()
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    tag = _ts_tag()
    log_path = logs / ("workflow_%s.log" % tag)
    manifest_path = logs / ("manifest_%s.json" % tag)

    manifest = Manifest(
        started_at=datetime.now().isoformat(),
        output_dir=str(out),
        max_files=args.max_files,
        mc_df_stage=args.mc_df_stage,
        var_set=args.var_set,
        dry_run=args.dry_run,
    )

    logger = WorkflowLog(log_path, manifest_path, manifest)
    logger.banner("numucc_1p0pi workflow test")
    logger.line("Output directory: %s" % out)
    logger.line("Log file:        %s" % log_path)
    logger.line("Manifest:        %s" % manifest_path)
    logger.line("Repo root:       %s" % _REPO_ROOT)
    logger.line("max_files:       %s (0 = unlimited)" % args.max_files)
    logger.line("PYTHON:          %s" % sys.executable)

    multisim_chunks = out / "multisim" / "chunks"
    syst_disk = out / "syst_disk"
    detvar_chunks = out / "detvar" / "chunks"
    es_chunks = out / "event_selection" / "chunks"
    es_plots = out / "event_selection" / "plots"

    run_detvar = (not args.skip_detvar) and (not args.skip_event_selection)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(_REPO_ROOT) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")

    exit_code = 0

    # ----- Stage 1–2: multisim -----
    if not args.skip_multisim:
        logger.banner("STAGE 1 — Multisim map (syst_multisim_chunk.py)")
        raw_tasks = list(iter_multisim_chunk_tasks(args.mc_df_stage))
        tasks = _limit_multisim_tasks(raw_tasks, args.max_files)
        logger.line(
            "Multisim tasks: %d total after limit (from %d raw tasks; max_files=%s)"
            % (len(tasks), len(raw_tasks), args.max_files)
        )

        multisim_chunk_py = _SCRIPTS_DIR / "syst_multisim_chunk.py"
        multisim_chunks.mkdir(parents=True, exist_ok=True)
        runs_ms: List[Dict[str, Any]] = []

        if not tasks:
            logger.line("No multisim tasks — skipping map and aggregate (check globs / NUMUCC_SPRING_GEN1_ROOT).")
            manifest.add_stage("multisim_map", [])
            manifest.add_stage("multisim_aggregate", [{"skipped": True, "reason": "no_tasks"}])
        else:
            for idx, (syst_name, df_path) in enumerate(tasks, start=1):
                stem = Path(df_path).stem
                if syst_name == "COMBINED":
                    out_leaf = "nu__{}.pkl".format(stem)
                    chunk_args = [
                        "--df_file",
                        df_path,
                        "--out_dir",
                        str(multisim_chunks),
                        "--var-set",
                        args.var_set,
                    ]
                else:
                    out_leaf = "nu__{}__{}.pkl".format(syst_name, stem)
                    chunk_args = [
                        "--df_file",
                        df_path,
                        "--out_dir",
                        str(multisim_chunks),
                        "--var-set",
                        args.var_set,
                        "--syst-names",
                        syst_name,
                    ]
                meta: Dict[str, Any] = {
                    "invocation_index": idx,
                    "syst_task_label": syst_name,
                    "input_df": df_path,
                    "expected_pickle_leaf": out_leaf,
                }
                rc = _run_python(
                    multisim_chunk_py,
                    chunk_args,
                    env=env,
                    dry_run=args.dry_run,
                    logger=logger,
                    stage="multisim_map",
                    meta=meta,
                )
                runs_ms.append(meta)
                if rc != 0:
                    exit_code = rc
                    logger.line("WARNING: multisim map failed; continuing to allow inspection.")
            manifest.add_stage("multisim_map", runs_ms)

            logger.banner("STAGE 2 — Multisim aggregate (syst_multisim_aggregate.py)")
            syst_disk.mkdir(parents=True, exist_ok=True)
            agg_args = [
                "--chunks_dir",
                str(multisim_chunks),
                "--syst-disk-root",
                str(syst_disk),
                "--mc-df-stage",
                args.mc_df_stage,
                "--var-set",
                args.var_set,
            ]
            if not args.multisim_include_cosmics:
                agg_args.append("--skip-cosmics")
            if args.multisim_no_plots:
                agg_args.append("--no-plots")

            meta_agg: Dict[str, Any] = {
                "chunks_dir": str(multisim_chunks),
                "syst_disk_root": str(syst_disk),
                "skip_cosmics": not args.multisim_include_cosmics,
                "no_plots": args.multisim_no_plots,
            }
            rc = _run_python(
                _SCRIPTS_DIR / "syst_multisim_aggregate.py",
                agg_args,
                env=env,
                dry_run=args.dry_run,
                logger=logger,
                stage="multisim_aggregate",
                meta=meta_agg,
            )
            manifest.add_stage("multisim_aggregate", [meta_agg])
            if rc != 0:
                exit_code = rc
    else:
        logger.line("Skipping multisim (--skip-multisim).")

    # ----- Stage 3–4: DetVar -----
    if run_detvar:
        logger.banner("STAGE 3 — DetVar map (syst_detvar_chunk.py)")
        raw_dv = list(iter_detvar_chunk_jobs())
        dv_jobs = _limit_detvar_jobs(raw_dv, args.max_files)
        logger.line(
            "DetVar jobs: %d after limit (from %d raw; max_files=%s per WireMod tag); tags=%s"
            % (
                len(dv_jobs),
                len(raw_dv),
                args.max_files,
                list(dict.fromkeys(t for t, _ in DETVAR_WIREMOD_GLOBS)),
            )
        )
        detvar_chunks.mkdir(parents=True, exist_ok=True)
        runs_dv: List[Dict[str, Any]] = []

        if not dv_jobs:
            logger.line(
                "No DetVar inputs matched — skipping map/aggregate (extend DETVAR_WIREMOD_GLOBS / paths)."
            )
            manifest.add_stage("detvar_map", [])
            manifest.add_stage("detvar_aggregate", [{"skipped": True, "reason": "no_jobs"}])
        else:
            detvar_chunk_py = _SCRIPTS_DIR / "syst_detvar_chunk.py"
            for idx, (wm_tag, df_path) in enumerate(dv_jobs, start=1):
                stem = Path(df_path).stem
                meta_dv = {
                    "invocation_index": idx,
                    "wiremod_tag": wm_tag,
                    "input_df": df_path,
                    "expected_pickle_leaf": "%s__%s.pkl" % (wm_tag, stem),
                }
                rc = _run_python(
                    detvar_chunk_py,
                    [
                        "--df_file",
                        df_path,
                        "--wiremod_tag",
                        wm_tag,
                        "--out_dir",
                        str(detvar_chunks),
                    ],
                    env=env,
                    dry_run=args.dry_run,
                    logger=logger,
                    stage="detvar_map",
                    meta=meta_dv,
                )
                runs_dv.append(meta_dv)
                if rc != 0:
                    exit_code = rc
                    logger.line("WARNING: syst_detvar_chunk failed; continuing.")
            manifest.add_stage("detvar_map", runs_dv)

            logger.banner("STAGE 4 — DetVar aggregate (syst_detvar_aggregate.py)")
            syst_disk.mkdir(parents=True, exist_ok=True)
            meta_dv_agg: Dict[str, Any] = {
                "in_dir": str(detvar_chunks),
                "syst_disk_root": str(syst_disk),
            }
            rc = _run_python(
                _SCRIPTS_DIR / "syst_detvar_aggregate.py",
                ["--in_dir", str(detvar_chunks), "--syst-disk-root", str(syst_disk)],
                env=env,
                dry_run=args.dry_run,
                logger=logger,
                stage="detvar_aggregate",
                meta=meta_dv_agg,
            )
            manifest.add_stage("detvar_aggregate", [meta_dv_agg])
            if rc != 0:
                exit_code = rc
    elif args.skip_event_selection:
        logger.line("Skipping DetVar (no event-selection aggregate).")
    else:
        logger.line("Skipping DetVar (--skip-detvar).")

    # ----- Stage 5–6: event selection -----
    if not args.skip_event_selection:
        logger.banner("STAGE 5 — Event selection map (event_selection_chunk.py)")
        es_chunks.mkdir(parents=True, exist_ok=True)
        event_chunk_py = _SCRIPTS_DIR / "event_selection_chunk.py"
        runs_es: List[Dict[str, Any]] = []

        sample_order = ("mc", "data", "intime", "offbeam", "dirt")
        for sample in sample_order:
            paths = list(iter_event_selection_df_paths(sample))
            paths = _limit_paths(paths, args.max_files)
            logger.subsection("sample=%s glob=%s → %d file(s)" % (sample, EVENT_SELECTION_GLOBS[sample], len(paths)))
            for idx, df_path in enumerate(paths, start=1):
                meta = {
                    "sample": sample,
                    "invocation_index": idx,
                    "input_df": df_path,
                    "glob_pattern": EVENT_SELECTION_GLOBS[sample],
                }
                rc = _run_python(
                    event_chunk_py,
                    [
                        "--df_file",
                        df_path,
                        "--sample",
                        sample,
                        "--out_dir",
                        str(es_chunks),
                    ],
                    env=env,
                    dry_run=args.dry_run,
                    logger=logger,
                    stage="event_selection_map",
                    meta=meta,
                )
                runs_es.append(meta)
                if rc != 0:
                    exit_code = rc
                    logger.line("WARNING: event_selection_chunk failed for sample=%s path=%s" % (sample, df_path))

        manifest.add_stage("event_selection_map", runs_es)

        det_npz_src = syst_disk / "Detector" / "detector_syst_dict.npz"

        logger.banner("STAGE 6 — Event selection aggregate (event_selection_aggregate.py)")
        es_plots.mkdir(parents=True, exist_ok=True)
        agg_es_args = [
            "--in_dir",
            str(es_chunks),
            "--out_dir",
            str(es_plots),
        ]
        # Disk syst requires a complete syst_disk_layout tree (utils.get_syst_unc has no fallbacks).
        use_disk_syst = (
            not args.skip_multisim
            and args.multisim_include_cosmics
            and run_detvar
        )
        if use_disk_syst:
            _seed_genie_cov_pkl(syst_disk, logger, args.dry_run)
            agg_es_args.extend(["--syst-disk-root", str(syst_disk)])

        meta_es: Dict[str, Any] = {
            "in_dir": str(es_chunks),
            "out_dir": str(es_plots),
            "syst_disk_root": str(syst_disk) if use_disk_syst else None,
            "detector_npz_expected": str(det_npz_src) if run_detvar else None,
        }
        rc = _run_python(
            _SCRIPTS_DIR / "event_selection_aggregate.py",
            agg_es_args,
            env=env,
            dry_run=args.dry_run,
            logger=logger,
            stage="event_selection_aggregate",
            meta=meta_es,
        )
        manifest.add_stage("event_selection_aggregate", [meta_es])
        if rc != 0:
            exit_code = rc
    else:
        logger.line("Skipping event selection (--skip-event-selection).")

    logger.banner("DONE")
    logger.line("Exit code (worst subprocess): %s" % exit_code)
    logger.line("Plots (event selection): %s" % es_plots)
    logger.line("Syst disk tree:           %s" % syst_disk)
    logger.close()

    print("[workflow-test] Wrote %s" % log_path)
    print("[workflow-test] Wrote %s" % manifest_path)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
