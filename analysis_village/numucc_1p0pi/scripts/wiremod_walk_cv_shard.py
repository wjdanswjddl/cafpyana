#!/usr/bin/env python3
"""Matched CV walk. Same as ``wiremod_walk_shard.py --cv``."""
from __future__ import annotations

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from analysis_village.numucc_1p0pi.scripts.wiremod_walk_shard import main


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--cv" not in argv:
        argv = ["--cv", *argv]
    raise SystemExit(main(argv))
