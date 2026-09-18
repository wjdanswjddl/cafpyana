#!/usr/bin/env python3
"""Merge WireMod walk shards and remake plots (post XTXW drop-map rewalk)."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

CACHE = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod/cache")
OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"


def main() -> int:
    os.chdir(REPO)
    yz = [str(CACHE / f"wiremod_walk_yz_s{i}.pkl") for i in range(6)]
    xtxw = [str(CACHE / f"wiremod_walk_xtxw_s{i}.pkl") for i in range(6)]
    cv = [str(CACHE / f"wiremod_walk_cv_s{i}.pkl") for i in range(6)]
    for p in yz + xtxw + cv:
        assert Path(p).is_file(), p

    mlog = CACHE / "wiremod_merge_walk.log"
    mlog.write_bytes(b"")
    cmd = [
        PY,
        "-u",
        "analysis_village/numucc_1p0pi/scripts/wiremod_merge_walk_shards.py",
        "--out-base",
        str(OUT),
        "--yz-shard-ckpts",
        *yz,
        "--xtxw-shard-ckpts",
        *xtxw,
        "--cv-shard-ckpts",
        *cv,
        "--cv-campaign",
        "2026_09_04_172912__sel_all-mc-CV",
    ]
    print("MERGE ...", flush=True)
    r = subprocess.run(cmd, cwd=str(REPO), stdout=open(mlog, "ab"), stderr=subprocess.STDOUT)
    print(f"merge_rc={r.returncode}", flush=True)
    print(mlog.read_text()[-1200:], flush=True)
    if r.returncode != 0:
        return r.returncode

    plog = CACHE / "wiremod_make_plots.log"
    plog.write_bytes(b"")
    cmd2 = [
        PY,
        "-u",
        "analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py",
        "--out-base",
        str(OUT),
        "--skip-inspect",
    ]
    print("PLOTS ...", flush=True)
    r2 = subprocess.run(cmd2, cwd=str(REPO), stdout=open(plog, "ab"), stderr=subprocess.STDOUT)
    print(f"plots_rc={r2.returncode}", flush=True)
    print(plog.read_text()[-800:], flush=True)
    return r2.returncode


if __name__ == "__main__":
    raise SystemExit(main())
