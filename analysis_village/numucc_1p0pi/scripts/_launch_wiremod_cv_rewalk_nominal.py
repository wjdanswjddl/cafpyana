#!/usr/bin/env python3
"""Rewalk WireMod CV (cut=30) after cut-campaign pollution, then merge+plots.

YZ/XTXW shards are reused (Sep-18, nominal NU_SCORE_TH=0.45).
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"
OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
N_SHARDS = 6
BATCH = 25
CV_CAMPAIGN = "2026_09_19_035612__sel_all-mc-CV-updatecalo-cvonly"


def _log(msg: str) -> None:
    print(msg, flush=True)


def _pid_alive(pid: int) -> bool:
    try:
        with open(f"/proc/{pid}/stat") as fh:
            rest = fh.read().rsplit(")", 1)[-1].strip()
            state = rest.split(None, 1)[0] if rest else ""
        return state not in ("Z", "")
    except OSError:
        return False


def main() -> int:
    os.chdir(REPO)
    env = {
        **os.environ,
        "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}",
        "BEARER_TOKEN_FILE": os.environ.get(
            "BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}"
        ),
    }

    # Assert nominal thresholds before launch.
    rc = subprocess.call(
        [
            PY, "-c",
            "from analysis_village.numucc_1p0pi.makedf.selections import "
            "NU_SCORE_TH, MU_CHI2MU_TH; "
            "assert float(NU_SCORE_TH)==0.45, NU_SCORE_TH; "
            "assert int(MU_CHI2MU_TH)==30, MU_CHI2MU_TH; "
            "print('thresholds OK', NU_SCORE_TH, MU_CHI2MU_TH)",
        ],
        cwd=str(REPO),
        env=env,
    )
    if rc != 0:
        return rc

    # Kill leftover CV walks for this tree.
    try:
        out = subprocess.check_output(["pgrep", "-af", "wiremod_walk_cv_shard"], text=True)
    except subprocess.CalledProcessError:
        out = ""
    for line in out.splitlines():
        if str(CACHE) not in line:
            continue
        try:
            os.kill(int(line.split(None, 1)[0]), signal.SIGTERM)
        except (OSError, ValueError):
            pass
    time.sleep(2)

    pids = []
    for sid in range(N_SHARDS):
        ck = CACHE / f"wiremod_walk_cv_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        logp = CACHE / f"wiremod_walk_cv_s{sid}.log"
        logp.write_bytes(b"")
        lf = open(logp, "ab")
        p = subprocess.Popen(
            [
                PY, "-u", str(SCRIPTS / "wiremod_walk_cv_shard.py"),
                "--files-pkl", str(CACHE / "wiremod_cv_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid),
                "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH),
                "--rss-limit-gb", "20",
                "--mu-chi2mu-th", "30",
            ],
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=str(REPO),
            start_new_session=True,
            env=env,
        )
        _log(f"launched cv walk {p.pid} s{sid}")
        pids.append(p.pid)

    while any(_pid_alive(p) for p in pids):
        alive = sum(_pid_alive(p) for p in pids)
        _log(f"waiting cv walks {alive}/{len(pids)}")
        time.sleep(60)
    for p in pids:
        try:
            os.waitpid(p, os.WNOHANG)
        except ChildProcessError:
            pass

    for sid in range(N_SHARDS):
        txt = (CACHE / f"wiremod_walk_cv_s{sid}.log").read_text()
        if "done pot=" not in txt:
            raise SystemExit(f"incomplete cv s{sid}")
        if "NU_SCORE_TH=0.45" not in txt:
            raise SystemExit(f"cv s{sid} missing nominal NU_SCORE_TH log line")

    yz = [str(CACHE / f"wiremod_walk_yz_s{i}.pkl") for i in range(N_SHARDS)]
    xtxw = [str(CACHE / f"wiremod_walk_xtxw_s{i}.pkl") for i in range(N_SHARDS)]
    cv = [str(CACHE / f"wiremod_walk_cv_s{i}.pkl") for i in range(N_SHARDS)]
    merge_log = CACHE / "wiremod_merge_walk_chi2new_fix.log"
    plots_log = CACHE / "wiremod_make_plots_chi2new_fix.log"
    _log("merge")
    with open(merge_log, "w") as lf:
        rc = subprocess.call(
            [
                PY, "-u", str(SCRIPTS / "wiremod_merge_walk_shards.py"),
                "--out-base", str(OUT),
                "--yz-shard-ckpts", *yz,
                "--xtxw-shard-ckpts", *xtxw,
                "--cv-shard-ckpts", *cv,
                "--cv-campaign", CV_CAMPAIGN,
                "--mu-chi2mu-th", "30",
            ],
            cwd=str(REPO), stdout=lf, stderr=subprocess.STDOUT, env=env,
        )
    if rc != 0:
        raise SystemExit(f"merge failed rc={rc}")
    _log("plots")
    with open(plots_log, "w") as lf:
        rc = subprocess.call(
            [
                PY, "-u", str(SCRIPTS / "wiremod_make_plots.py"),
                "--out-base", str(OUT), "--skip-inspect",
            ],
            cwd=str(REPO), stdout=lf, stderr=subprocess.STDOUT, env=env,
        )
    if rc != 0:
        raise SystemExit(f"plots failed rc={rc}")
    _log(f"done OUT={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
