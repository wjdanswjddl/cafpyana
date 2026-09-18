#!/usr/bin/env python3
"""Launch (or relaunch) XTXW matched artkey dedupe only."""
from __future__ import annotations

import os
import signal
import subprocess
import time
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPT = REPO / "analysis_village/numucc_1p0pi/scripts/dedupe_matched_artkeys.py"


def main() -> int:
    # kill existing dedupe only
    out = subprocess.check_output(["pgrep", "-f", "dedupe_matched_artkeys.py"], text=True, stderr=subprocess.DEVNULL) if False else ""
    try:
        out = subprocess.check_output(["pgrep", "-f", "dedupe_matched_artkeys.py"], text=True)
    except subprocess.CalledProcessError:
        out = ""
    for tok in out.split():
        try:
            os.kill(int(tok), signal.SIGTERM)
        except OSError:
            pass
    time.sleep(1)

    log = CACHE / "wiremod_dedupe_xtxw.log"
    log.write_bytes(b"")
    lf = open(log, "ab")
    p = subprocess.Popen(
        [PY, "-u", str(SCRIPT), "--matched-dir", str(OUT / "matched" / "xtxw")],
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
    )
    (CACHE / "wiremod_dedupe_xtxw.pid").write_text(f"{p.pid}\n")
    print(f"dedupe pid={p.pid}", flush=True)
    time.sleep(10)
    print(f"alive={Path(f'/proc/{p.pid}').exists()}", flush=True)
    print(log.read_text()[:1000], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
