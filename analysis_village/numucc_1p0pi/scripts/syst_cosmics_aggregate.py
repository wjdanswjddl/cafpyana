#!/usr/bin/env python
"""Reduce phase: merge ``cosmics__offbeam__*.pkl`` + ``cosmics__intime__*.pkl`` → ``cosmics_syst_dict.npz``.

Gate scaling ``sum(offbeam gates) / sum(intime gates)`` matches
``files_config.get_ana_dfs(option='cosmics_systs')`` applied to intime histograms.

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
from typing import Any, Dict, MutableMapping

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
from analysis_village.numucc_1p0pi.scripts.get_systematics_cosmics import (
    process_variable_cosmics_from_histograms,
    save_cosmics_npz,
)


def _merge_chunks(sample: str, chunks_dir: str) -> tuple[float, Dict[str, np.ndarray]]:
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
        hists = d["hists"]
        if hsum is None:
            hsum = {k: np.array(v, dtype=np.float64, copy=True) for k, v in hists.items()}
        else:
            for k, v in hists.items():
                hsum[k] += np.asarray(v, dtype=np.float64)
    assert hsum is not None
    return gates_tot, hsum


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

    g_off, h_off = _merge_chunks("offbeam", chunks_dir)
    g_in, h_in_raw = _merge_chunks("intime", chunks_dir)
    scale = float(g_off / g_in) if g_in else 1.0
    logger.info(
        "Gates offbeam=%.6e intime=%.6e scale(data/mc)=%.6f",
        g_off,
        g_in,
        scale,
    )

    h_in = {k: scale * v for k, v in h_in_raw.items()}

    var_configs = build_variable_configs(getattr(args, "vars", None))
    save_plots = not getattr(args, "no_plots", False)
    cv_mode = getattr(args, "cv_mode", "offbeam")

    syst_dict: Dict[str, MutableMapping[str, Any]] = {}

    for var_config in tqdm(var_configs, desc="cosmics"):
        slug = var_config.var_save_name
        try:
            pay = process_variable_cosmics_from_histograms(
                h_off[slug],
                h_in[slug],
                var_config,
                cv_mode,
                save_fig_dir,
                save_plots,
            )
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

    logger.info("Done -> %s (log %s)", save_fig_dir, log_path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chunks_dir", required=True)
    p.add_argument(
        "--syst-disk-root",
        default=None,
        help="Writes to <root>/%s/ (default: env %s)." % (SUB_COSMICS, SYST_DISK_ENV),
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
