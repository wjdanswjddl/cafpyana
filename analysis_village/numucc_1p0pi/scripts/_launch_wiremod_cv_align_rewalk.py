#!/usr/bin/env python3
"""Drop CV-missing events from WireMod walks and recalculate envelopes.

1. Build CV-align drop maps (YZ + XTXW).
2. Rewalk YZ/XTXW (cut=30) with those maps; reuse existing CV walk.
3. Merge + remake WireMod plots.
4. Same for WireMod-chi2mu25 (full YZ/XTXW/CV rewalk at mu_chi2mu_th=25).
"""
from __future__ import annotations

import os
import pickle
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"
OUT30 = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE30 = OUT30 / "cache"
OUT25 = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod-chi2mu25")
CACHE25 = OUT25 / "cache"
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


def _kill(pattern: str, must_contain: str | None = None) -> None:
    try:
        out = subprocess.check_output(["pgrep", "-af", pattern], text=True)
    except subprocess.CalledProcessError:
        return
    for line in out.splitlines():
        if must_contain and must_contain not in line:
            continue
        tok = line.split(None, 1)[0]
        try:
            os.kill(int(tok), signal.SIGTERM)
        except (OSError, ValueError):
            pass


def _launch(cmd: list[str], logpath: Path) -> int:
    logpath.parent.mkdir(parents=True, exist_ok=True)
    logpath.write_bytes(b"")
    lf = open(logpath, "ab")
    p = subprocess.Popen(
        cmd,
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
        env={
            **os.environ,
            "BEARER_TOKEN_FILE": os.environ.get(
                "BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}"
            ),
            "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}",
        },
    )
    _log(f"launched {p.pid} {' '.join(cmd[2:6])}...")
    return p.pid


def _wait_pids(pids: list[int], label: str, poll_s: int = 60) -> None:
    while True:
        alive = [p for p in pids if _pid_alive(p)]
        if not alive:
            for p in pids:
                try:
                    os.waitpid(p, os.WNOHANG)
                except ChildProcessError:
                    pass
            _log(f"{label} all done")
            return
        _log(f"{label} waiting {len(alive)}/{len(pids)} alive")
        time.sleep(poll_s)


def _assert_done(cache: Path, prefix: str, n: int) -> None:
    for sid in range(n):
        logp = cache / f"{prefix}_s{sid}.log"
        txt = logp.read_text() if logp.is_file() else ""
        if "done pot=" not in txt:
            raise SystemExit(f"incomplete: {logp}")


def build_drop_maps() -> None:
    logp = CACHE30 / "wiremod_cv_align_drop_maps.log"
    rc = subprocess.run(
        [
            PY, "-u", str(SCRIPTS / "build_wiremod_cv_align_drop_maps.py"),
            "--out-base", str(OUT30),
            "--workers", "12",
        ],
        cwd=str(REPO),
        stdout=open(logp, "wb"),
        stderr=subprocess.STDOUT,
        env={
            **os.environ,
            "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}",
        },
    )
    _log(logp.read_text()[-2000:])
    if rc.returncode != 0:
        raise SystemExit(f"drop-map build failed rc={rc.returncode}")


def rewalk_wm(cache: Path, *, mu_chi2mu_th: float | None, label: str) -> None:
    """Rewalk YZ + XTXW (and CV if mu cut overridden)."""
    yz_drop = CACHE30 / "wiremod_yz_cv_align_drop_map.pkl"
    xtxw_drop = CACHE30 / "wiremod_xtxw_cv_align_drop_map.pkl"
    if not yz_drop.is_file() or not xtxw_drop.is_file():
        raise SystemExit("missing cv-align drop maps")

    # Ensure file lists exist in this cache.
    for lab in ("yz", "xtxw", "cv"):
        src = CACHE30 / f"wiremod_{lab}_all_files.pkl"
        dst = cache / f"wiremod_{lab}_all_files.pkl"
        if src.is_file() and (not dst.is_file() or cache != CACHE30):
            shutil.copy2(src, dst)

    if cache != CACHE30:
        shutil.copy2(yz_drop, cache / yz_drop.name)
        shutil.copy2(xtxw_drop, cache / "wiremod_xtxw_drop_map.pkl")
        shutil.copy2(xtxw_drop, cache / xtxw_drop.name)

    _kill("wiremod_walk_shard.py", must_contain=str(cache))
    if mu_chi2mu_th is not None:
        _kill("wiremod_walk_cv_shard.py", must_contain=str(cache))

    pids = []
    for lab, drop in (("yz", yz_drop if cache == CACHE30 else cache / yz_drop.name),
                      ("xtxw", xtxw_drop if cache == CACHE30 else cache / xtxw_drop.name)):
        for sid in range(N_SHARDS):
            ck = cache / f"wiremod_walk_{lab}_s{sid}.pkl"
            if ck.exists():
                ck.unlink()
            cmd = [
                PY, "-u", str(SCRIPTS / "wiremod_walk_shard.py"),
                "--files-pkl", str(cache / f"wiremod_{lab}_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid),
                "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH),
                "--rss-limit-gb", "20",
                "--label", f"{lab}_{label}",
                "--drop-map-pkl", str(drop),
            ]
            if mu_chi2mu_th is not None:
                cmd += ["--mu-chi2mu-th", str(mu_chi2mu_th)]
            pids.append(_launch(cmd, cache / f"wiremod_walk_{lab}_s{sid}.log"))

    if mu_chi2mu_th is not None:
        for sid in range(N_SHARDS):
            ck = cache / f"wiremod_walk_cv_s{sid}.pkl"
            if ck.exists():
                ck.unlink()
            # wiremod_walk_cv_shard has no --label flag
            cmd = [
                PY, "-u", str(SCRIPTS / "wiremod_walk_cv_shard.py"),
                "--files-pkl", str(cache / "wiremod_cv_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid),
                "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH),
                "--rss-limit-gb", "20",
                "--mu-chi2mu-th", str(mu_chi2mu_th),
            ]
            pids.append(_launch(cmd, cache / f"wiremod_walk_cv_s{sid}.log"))

    _wait_pids(pids, f"walk_{label}")
    _assert_done(cache, "wiremod_walk_yz", N_SHARDS)
    _assert_done(cache, "wiremod_walk_xtxw", N_SHARDS)
    if mu_chi2mu_th is not None:
        _assert_done(cache, "wiremod_walk_cv", N_SHARDS)


def merge_and_plot(out: Path, cache: Path, *, mu_chi2mu_th: float | None) -> None:
    yz = [str(cache / f"wiremod_walk_yz_s{i}.pkl") for i in range(N_SHARDS)]
    xtxw = [str(cache / f"wiremod_walk_xtxw_s{i}.pkl") for i in range(N_SHARDS)]
    cv = [str(cache / f"wiremod_walk_cv_s{i}.pkl") for i in range(N_SHARDS)]
    mlog = cache / "wiremod_merge_walk_cvalign.log"
    cmd = [
        PY, "-u", str(SCRIPTS / "wiremod_merge_walk_shards.py"),
        "--out-base", str(out),
        "--yz-shard-ckpts", *yz,
        "--xtxw-shard-ckpts", *xtxw,
        "--cv-shard-ckpts", *cv,
        "--cv-campaign", CV_CAMPAIGN,
    ]
    if mu_chi2mu_th is not None:
        cmd += ["--mu-chi2mu-th", str(mu_chi2mu_th)]
    r = subprocess.run(
        cmd, cwd=str(REPO), stdout=open(mlog, "wb"), stderr=subprocess.STDOUT
    )
    _log(mlog.read_text()[-1500:])
    if r.returncode != 0:
        raise SystemExit(f"merge failed rc={r.returncode}")

    plog = cache / "wiremod_make_plots_cvalign.log"
    r2 = subprocess.run(
        [
            PY, "-u", str(SCRIPTS / "wiremod_make_plots.py"),
            "--out-base", str(out),
            "--skip-inspect",
        ],
        cwd=str(REPO),
        stdout=open(plog, "wb"),
        stderr=subprocess.STDOUT,
    )
    _log(plog.read_text()[-1200:])
    if r2.returncode != 0:
        raise SystemExit(f"plots failed rc={r2.returncode}")


def main() -> int:
    os.chdir(REPO)
    CACHE30.mkdir(parents=True, exist_ok=True)
    note = CACHE30 / "README_cv_align.txt"
    note.write_text(
        "CV updatecalo-cvonly is missing ~3.3% of the original common keys.\n"
        "WireMod YZ/XTXW walks are filtered via cv-align drop maps so the\n"
        "event set matches matched/cv before envelopes are recomputed.\n"
    )

    _log("=== build drop maps ===")
    build_drop_maps()

    _log("=== rewalk cut30 (YZ+XTXW only; CV reused) ===")
    rewalk_wm(CACHE30, mu_chi2mu_th=None, label="cvalign30")
    _log("=== merge+plots cut30 ===")
    merge_and_plot(OUT30, CACHE30, mu_chi2mu_th=None)

    _log("=== setup chi2mu25 ===")
    CACHE25.mkdir(parents=True, exist_ok=True)
    matched = OUT25 / "matched"
    if not matched.exists() and not matched.is_symlink():
        matched.symlink_to(OUT30 / "matched")
    for lab in ("yz", "xtxw", "cv"):
        shutil.copy2(CACHE30 / f"wiremod_{lab}_all_files.pkl", CACHE25 / f"wiremod_{lab}_all_files.pkl")

    _log("=== rewalk cut25 (YZ+XTXW+CV) ===")
    rewalk_wm(CACHE25, mu_chi2mu_th=25.0, label="cvalign25")
    _log("=== merge+plots cut25 ===")
    merge_and_plot(OUT25, CACHE25, mu_chi2mu_th=25.0)

    _log("ALL DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
