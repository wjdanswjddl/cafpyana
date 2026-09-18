#!/usr/bin/env python3
"""Launch verify_wiremod_sel_all_eventcounts.py in background."""
from __future__ import annotations

import os
import signal
import subprocess
import time
from pathlib import Path

CACHE = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod/cache")
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPT = REPO / "analysis_village/numucc_1p0pi/scripts/verify_wiremod_sel_all_eventcounts.py"


def main() -> int:
    try:
        out = subprocess.check_output(["pgrep", "-f", "verify_wiremod_sel_all_eventcounts.py"], text=True)
        for tok in out.split():
            try:
                os.kill(int(tok), signal.SIGTERM)
            except OSError:
                pass
        time.sleep(1)
    except subprocess.CalledProcessError:
        pass

    log = CACHE / "wiremod_sel_all_eventcount_verify.log"
    log.write_bytes(b"")
    lf = open(log, "ab")
    p = subprocess.Popen(
        [PY, "-u", str(SCRIPT)],
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
    )
    (CACHE / "wiremod_sel_all_eventcount_verify.pid").write_text(f"{p.pid}\n")
    print(f"verify pid={p.pid}", flush=True)
    time.sleep(8)
    print(f"alive={Path(f'/proc/{p.pid}').exists()}", flush=True)
    print(log.read_text()[:1200], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
