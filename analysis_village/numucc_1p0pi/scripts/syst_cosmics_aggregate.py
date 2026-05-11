#!/usr/bin/env python
"""Reduce phase: merge ``cosmics__offbeam__*.pkl`` + ``cosmics__intime__*.pkl`` → ``cosmics_syst_dict.npz``.

Gate scaling ``sum(offbeam gates) / sum(intime gates)`` matches
``files_config.get_ana_dfs(option='cosmics_systs')`` applied to intime histograms.

Both input-stage flavours produced by ``syst_cosmics_chunk.py`` are accepted and
auto-detected from the first pickle:

* ``input_stage="final"``: per-var histograms live under ``d["hists"]``. Output
  NPZ keys are the final-variable ``var_save_name``s (legacy layout).
* ``input_stage="sel_all"``: per-stage histograms live under
  ``d["stage_hists"][stage_key][var_save_name]``. Because cut-stage and final
  ``var_save_name``s are disjoint, the merged dict can be flattened to the same
  flat layout, so the output NPZ still uses ``{var_save_name: {Cosmics: …}}``
  -- downstream loaders (e.g. ``selected_events_intermediate_cuts``) do not need
  to know that the dict now also contains cut-variable covariances.

Usage::

    python syst_cosmics_aggregate.py --chunks_dir CHUNKS --syst-disk-root ROOT

Requires ``NUMUCC_SYST_DISK_ROOT`` or ``--syst-disk-root`` (writes under ``Cosmics/``).
"""
from __future__ import annotations

import argparse
import glob
import logging
import os
import pickle
import sys
import traceback
from datetime import datetime
from os import makedirs, path
from typing import Any, Dict, List, MutableMapping, Sequence, Tuple

import numpy as np

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(x=None, **kwargs):
        return x


_REPO_ROOT = path.abspath(path.join(path.dirname(__file__), "..", "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_village.numucc_1p0pi.syst_cosmics_common import build_variable_configs
from analysis_village.numucc_1p0pi.syst_disk_layout import SUB_COSMICS, SYST_DISK_ENV
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (
    CUT_STAGE_VAR_SPECS,
    FINAL_STAGE_KEY,
    final_stage_var_configs,
)
from analysis_village.numucc_1p0pi.scripts.get_systematics_cosmics import (
    process_variable_cosmics_from_histograms,
    save_cosmics_npz,
)


# ---------------------------------------------------------------------------
# Detect & merge chunk pickles
# ---------------------------------------------------------------------------
def _detect_input_stage(chunks_dir: str) -> str:
    """Peek at the first pickle to decide the merge layout."""
    pattern = path.join(chunks_dir, "cosmics__*__*.pkl")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise SystemExit("[cosmics-aggregate] no pickles matching %s" % pattern)
    with open(paths[0], "rb") as f:
        d = pickle.load(f)
    if "stage_hists" in d:
        return "sel_all"
    if "hists" in d:
        return "final"
    stage = d.get("input_stage")
    if stage in ("final", "sel_all"):
        return stage
    raise SystemExit(
        "[cosmics-aggregate] cannot determine input_stage from %s: keys=%s"
        % (paths[0], sorted(d.keys()))
    )


def _merge_chunks_final(sample: str, chunks_dir: str) -> Tuple[float, Dict[str, np.ndarray]]:
    pattern = path.join(chunks_dir, "cosmics__%s__*.pkl" % sample)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise SystemExit("[cosmics-aggregate] no pickles matching %s" % pattern)
    gates_tot = 0.0
    hsum: Dict[str, np.ndarray] | None = None
    for fp in tqdm(paths, desc="merge %s" % sample):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        gates_tot += float(d["gates"])
        # Normally the final-stage chunk payload has ``hists``. If a user accidentally
        # mixes sel_all + final chunks under the same chunks_dir, tolerate it by
        # extracting the final-stage histograms from ``stage_hists``.
        hists = d.get("hists")
        if hists is None and "stage_hists" in d:
            hists = d["stage_hists"].get(FINAL_STAGE_KEY, {})
        if hists is None:
            raise KeyError(
                "Chunk pickle missing both 'hists' and 'stage_hists': %s (keys=%s)"
                % (fp, sorted(d.keys()))
            )
        if hsum is None:
            hsum = {k: np.array(v, dtype=np.float64, copy=True) for k, v in hists.items()}
        else:
            for k, v in hists.items():
                hsum[k] += np.asarray(v, dtype=np.float64)
    assert hsum is not None
    return gates_tot, hsum


def _merge_chunks_sel_all(sample: str, chunks_dir: str) -> Tuple[
    float,
    Dict[str, np.ndarray],
    Dict[str, str],
]:
    """Merge sel_all chunks. Returns (gates_total, flat_hists, vsn_to_stage_key)."""
    pattern = path.join(chunks_dir, "cosmics__%s__*.pkl" % sample)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise SystemExit("[cosmics-aggregate] no pickles matching %s" % pattern)
    gates_tot = 0.0
    flat: Dict[str, np.ndarray] = {}
    vsn_to_stage: Dict[str, str] = {}
    for fp in tqdm(paths, desc="merge %s" % sample):
        with open(fp, "rb") as f:
            d = pickle.load(f)
        gates_tot += float(d["gates"])
        stage_hists = d["stage_hists"]
        for stage_key, var_map in stage_hists.items():
            for vsn, h in var_map.items():
                arr = np.asarray(h, dtype=np.float64)
                if vsn not in flat:
                    flat[vsn] = arr.copy()
                    vsn_to_stage[vsn] = stage_key
                else:
                    if flat[vsn].shape != arr.shape:
                        raise SystemExit(
                            "[cosmics-aggregate] shape mismatch for %s: %s vs %s"
                            % (vsn, flat[vsn].shape, arr.shape)
                        )
                    flat[vsn] += arr
    return gates_tot, flat, vsn_to_stage


# ---------------------------------------------------------------------------
# Variable config catalogue (cut + final) for the sel_all path
# ---------------------------------------------------------------------------
def _sel_all_var_configs() -> List[Any]:
    """All VariableConfig objects produced by the sel_all chunk path."""
    out: List[Any] = []
    seen: set[str] = set()
    for spec in CUT_STAGE_VAR_SPECS:
        vc = spec.var_config
        if vc.var_save_name in seen:
            continue
        out.append(vc)
        seen.add(vc.var_save_name)
    for vc in final_stage_var_configs():
        if vc.var_save_name in seen:
            continue
        out.append(vc)
        seen.add(vc.var_save_name)
    return out


# ---------------------------------------------------------------------------
# Aggregate driver
# ---------------------------------------------------------------------------
def run_aggregate(args: argparse.Namespace) -> None:
    """Shared entry for CLI and ``get_systematics_cosmics.py aggregate``."""
    tag = getattr(args, "out_tag", None) or datetime.now().strftime("%Y%m%d")
    chunks_dir = args.chunks_dir
    save_fig_dir = path.join(args.syst_disk_root, SUB_COSMICS)
    makedirs(save_fig_dir, exist_ok=True)

    log_path = getattr(args, "error_log", None) or path.join(save_fig_dir, "cosmics_failures.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    logger = logging.getLogger("cosmics_syst_agg")
    logger.info("Run tag=%s chunks_dir=%s -> %s", tag, chunks_dir, save_fig_dir)

    input_stage = getattr(args, "input_stage", None) or _detect_input_stage(chunks_dir)
    if getattr(args, "input_stage", None):
        logger.info("Using input_stage=%s (CLI override)", input_stage)
    else:
        logger.info("Detected input_stage=%s", input_stage)

    if input_stage == "final":
        g_off, h_off = _merge_chunks_final("offbeam", chunks_dir)
        g_in, h_in_raw = _merge_chunks_final("intime", chunks_dir)
        var_configs = build_variable_configs(getattr(args, "vars", None))
        vsn_to_stage: Dict[str, str] = {}
    else:
        g_off, h_off, vsn_to_stage_off = _merge_chunks_sel_all("offbeam", chunks_dir)
        g_in, h_in_raw, vsn_to_stage_in = _merge_chunks_sel_all("intime", chunks_dir)
        vsn_to_stage = {**vsn_to_stage_in, **vsn_to_stage_off}
        sel_all_vcs = _sel_all_var_configs()
        requested = getattr(args, "vars", None)
        if requested:
            req = set(requested)
            var_configs = [vc for vc in sel_all_vcs if vc.var_save_name in req]
            missing = req - {vc.var_save_name for vc in var_configs}
            if missing:
                logger.warning("Requested vars not produced by sel_all chunk: %s", sorted(missing))
        else:
            var_configs = sel_all_vcs

    scale = float(g_off / g_in) if g_in else 1.0
    logger.info(
        "Gates offbeam=%.6e intime=%.6e scale(data/mc)=%.6f",
        g_off,
        g_in,
        scale,
    )

    h_in = {k: scale * v for k, v in h_in_raw.items()}

    save_plots = not getattr(args, "no_plots", False)
    cv_mode = getattr(args, "cv_mode", "offbeam")

    syst_dict: Dict[str, MutableMapping[str, Any]] = {}

    for var_config in tqdm(var_configs, desc="cosmics"):
        slug = var_config.var_save_name
        if slug not in h_off or slug not in h_in:
            logger.warning(
                "Skipping %s: missing histogram (offbeam=%s intime=%s)",
                slug,
                slug in h_off,
                slug in h_in,
            )
            continue
        try:
            pay = process_variable_cosmics_from_histograms(
                h_off[slug],
                h_in[slug],
                var_config,
                cv_mode,
                save_fig_dir,
                save_plots,
            )
            if vsn_to_stage:
                pay["stage_key"] = vsn_to_stage.get(slug, "")
            syst_dict.setdefault(slug, {})
            syst_dict[slug]["Cosmics"] = pay
        except Exception:
            logger.error(
                "FAILED variable=%s\n%s",
                slug,
                traceback.format_exc(),
            )

    if not getattr(args, "no_save_npz", False) and syst_dict:
        save_cosmics_npz(syst_dict, path.join(save_fig_dir, "cosmics_syst_dict.npz"))

    logger.info(
        "Done -> %s (log %s) input_stage=%s n_vars=%d",
        save_fig_dir,
        log_path,
        input_stage,
        len(syst_dict),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chunks_dir", required=True)
    p.add_argument(
        "--syst-disk-root",
        default=None,
        help="Writes to <root>/%s/ (default: env %s)." % (SUB_COSMICS, SYST_DISK_ENV),
    )
    p.add_argument(
        "--input-stage",
        dest="input_stage",
        choices=("final", "sel_all"),
        default=None,
        help="Override how chunk pickles are interpreted. If unset, auto-detects from pickles.",
    )
    p.add_argument("--out-tag", default=None)
    p.add_argument("--error-log", default=None)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--no-save-npz", action="store_true")
    p.add_argument(
        "--cv-mode",
        choices=("mean", "intime", "offbeam"),
        default="offbeam",
    )
    p.add_argument("--vars", nargs="*", default=None)
    args = p.parse_args(argv)
    root = args.syst_disk_root or os.environ.get(SYST_DISK_ENV)
    if not root:
        p.error("Pass --syst-disk-root or set %s." % SYST_DISK_ENV)
    args.syst_disk_root = path.abspath(path.expanduser(root.rstrip("/")))
    return args


def main(argv: list[str] | None = None) -> None:
    run_aggregate(parse_args(argv))


if __name__ == "__main__":
    main()
