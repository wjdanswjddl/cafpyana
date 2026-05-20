#!/usr/bin/env python3
"""
Find events common to **0xSCE**, **2xSCE**, and **CV** (calovar / nominal) MC
samples and write per-file ``_matched`` HDF5 outputs.

Same two-phase workflow as ``wiremod_match_common_events.py`` (meta scan, then
per-file write). Defaults to three variations and copies only ``meta`` and
``evt_cv`` (no calo-shift universes).

Usage
-----
    python sce_match_common_events.py --use-default-variations \\
        --common-keys-pkl /pnfs/.../sce_common_keys.pkl \\
        --summary-csv /pnfs/.../sce_matched_summary.csv

    python sce_match_common_events.py \\
        --variation 0xSCE /pnfs/.../0xSCE/merged_perTPC \\
        --variation 2xSCE /pnfs/.../2xSCE/merged_perTPC \\
        --variation cv   /pnfs/.../calovar/merged_perTPC

    python sce_match_common_events.py --phase meta ... --common-keys-pkl common.pkl
    python sce_match_common_events.py --phase write --common-keys-pkl common.pkl ...
"""

from __future__ import annotations

import sys
from os import path
from typing import Sequence

_SCRIPT_DIR = path.dirname(path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import wiremod_match_common_events as _core

# SCE comparison: cv universe only (see notebooks/sce.ipynb)
DEFAULT_KEYS2LOAD = ["meta", "evt_cv"]

DEFAULT_VARIATIONS = {
    "0xSCE": "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_05_12_113644__sel_2prong-mc-BNB_cosmics-0xSCE/merged_perTPC",
    "2xSCE": "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_05_12_114000__sel_2prong-mc-BNB_cosmics-2xSCE/merged_perTPC",
    "cv": "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs/"
    "2026_05_16_190744__sel_2prong-mc-BNB_cosmics-calovar/merged_perTPC",
}

_core.DEFAULT_KEYS2LOAD = DEFAULT_KEYS2LOAD
_core.DEFAULT_VARIATIONS = DEFAULT_VARIATIONS

_orig_build_parser = _core.build_parser


def build_parser():
    p = _orig_build_parser()
    p.description = (
        "Select events common to 0xSCE, 2xSCE, and CV MC samples and write _matched .df files."
    )
    for action in p._actions:
        if action.dest == "use_default_variations":
            action.help = (
                "Use bundled 0xSCE, 2xSCE, and cv (calovar) merged_perTPC directories "
                "(override with --variation)."
            )
        elif action.dest == "keys":
            action.default = ",".join(DEFAULT_KEYS2LOAD)
            action.help = "Comma-separated HDF keys to copy (default: meta,evt_cv)."
    return p


def main(argv: Sequence[str] | None = None) -> int:
    _core.build_parser = build_parser
    return _core.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
