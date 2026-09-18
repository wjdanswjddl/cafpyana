#!/usr/bin/env python3
"""Launch fast XTXW drop-map scan + replace supervisor (YZ walks already running)."""
from __future__ import annotations

import os
import shutil
import signal
import subprocess
import time
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"


def _kill_matching(pattern: str) -> None:
    try:
        out = subprocess.check_output(["pgrep", "-f", pattern], text=True)
    except subprocess.CalledProcessError:
        return
    for tok in out.split():
        try:
            os.kill(int(tok), signal.SIGTERM)
        except OSError:
            pass


def main() -> int:
    os.chdir(REPO)
    _kill_matching("dedupe_matched_artkeys.py")
    _kill_matching("_walk_supervisor_refix.sh")
    _kill_matching("_wiremod_walk_supervisor_refix.sh")
    time.sleep(2)

    drop_out = CACHE / "wiremod_xtxw_drop_map.pkl"
    if drop_out.exists():
        drop_out.unlink()

    log = CACHE / "wiremod_dedupe_xtxw.log"
    log.write_bytes(b"")
    lf = open(log, "ab")
    p = subprocess.Popen(
        [
            PY,
            "-u",
            str(SCRIPTS / "dedupe_matched_artkeys.py"),
            "--matched-dir",
            str(OUT / "matched" / "xtxw"),
            "--drop-map-out",
            str(drop_out),
        ],
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
    )
    (CACHE / "wiremod_dedupe_xtxw.pid").write_text(f"{p.pid}\n")
    print(f"drop-map scan pid={p.pid}", flush=True)

    shutil.copy2(SCRIPTS / "_wiremod_walk_supervisor_refix.sh", CACHE / "_walk_supervisor_refix.sh")
    (CACHE / "_walk_supervisor_refix.sh").chmod(0o755)
    # preserve YZ-waiting state: append note and restart supervisor fresh
    # (YZ pids file still valid)
    (CACHE / "wiremod_walk_supervisor_nohup.log").write_bytes(b"")
    sp = subprocess.Popen(
        ["bash", str(CACHE / "_walk_supervisor_refix.sh")],
        stdout=open(CACHE / "wiremod_walk_supervisor_nohup.log", "ab"),
        stderr=subprocess.STDOUT,
        start_new_session=True,
        cwd=str(REPO),
    )
    (CACHE / "wiremod_walk_supervisor.pid").write_text(f"{sp.pid}\n")
    print(f"supervisor pid={sp.pid}", flush=True)

    time.sleep(8)
    print(f"scan alive={Path(f'/proc/{p.pid}').exists()}", flush=True)
    print(f"sup alive={Path(f'/proc/{sp.pid}').exists()}", flush=True)
    print("--- dedupe log ---", flush=True)
    print(log.read_text()[:800], flush=True)
    print("--- supervisor ---", flush=True)
    print((CACHE / "wiremod_walk_supervisor.log").read_text()[:400], flush=True)
    yz = (CACHE / "wiremod_walk_yz_pids.txt").read_text().split()
    print("yz alive", sum(Path(f"/proc/{x}").exists() for x in yz), "/", len(yz), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
