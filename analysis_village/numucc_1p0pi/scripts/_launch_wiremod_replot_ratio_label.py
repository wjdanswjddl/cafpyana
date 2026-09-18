#!/usr/bin/env python3
"""Remake WireMod plots after ratio-panel ylabel update."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
SCRIPT = REPO / "analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py"


def main() -> int:
    logpath = CACHE / "wiremod_make_plots_actual_env.log"
    logpath.write_bytes(b"")
    lf = open(logpath, "ab")
    p = subprocess.Popen(
        [PY, "-u", str(SCRIPT), "--out-base", str(OUT), "--skip-inspect"],
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
        env={
            **os.environ,
            "BEARER_TOKEN_FILE": os.environ.get(
                "BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}"
            ),
        },
    )
    (CACHE / "wiremod_make_plots_actual_env.pid").write_text(f"{p.pid}\n")
    print(f"launched pid={p.pid} log={logpath}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
