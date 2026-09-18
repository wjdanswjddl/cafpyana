#!/usr/bin/env python3
"""Launch wiremod_make_plots.py (symmetric envelope) with log under WireMod/cache."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
SCRIPTS = REPO / "analysis_village" / "numucc_1p0pi" / "scripts"


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    logpath = CACHE / "wiremod_make_plots_symenv.log"
    pidpath = CACHE / "wiremod_make_plots_symenv.pid"
    logpath.write_bytes(b"")
    lf = open(logpath, "ab")
    p = subprocess.Popen(
        [PY, "-u", str(SCRIPTS / "wiremod_make_plots.py"), "--out-base", str(OUT)],
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
        env={**os.environ, "BEARER_TOKEN_FILE": os.environ.get("BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}")},
    )
    pidpath.write_text(f"{p.pid}\n")
    print(f"launched pid={p.pid} log={logpath}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
