#!/usr/bin/env python3
"""Remake WireMod NPZs+plots with WireMod in-file cv envelope baseline."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PY = "/exp/sbnd/app/users/munjung/env/bin/python"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
SCRIPT = REPO / "analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py"
OUTS = [
    Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod"),
    Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod-chi2mu25"),
]


def main() -> int:
    for out in OUTS:
        cache = out / "cache"
        if not (cache / "wiremod_sel_all_products.pkl").is_file():
            print(f"skip (no products): {out}", flush=True)
            continue
        logpath = Path(f"/tmp/wiremod_make_plots_infile_cv_{out.name}.log")
        pidpath = Path(f"/tmp/wiremod_make_plots_infile_cv_{out.name}.pid")
        logpath.write_bytes(b"")
        lf = open(logpath, "ab")
        p = subprocess.Popen(
            [PY, "-u", str(SCRIPT), "--out-base", str(out), "--skip-inspect"],
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
        pidpath.write_text(f"{p.pid}\n")
        # best-effort mirror under cache
        try:
            (cache / "wiremod_make_plots_infile_cv.pid").write_text(f"{p.pid}\n")
        except OSError:
            pass
        print(f"launched pid={p.pid} out={out.name} log={logpath}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
