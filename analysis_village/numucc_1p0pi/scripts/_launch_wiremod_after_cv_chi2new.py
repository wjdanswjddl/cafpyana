#!/usr/bin/env python3
"""After CV updatecalo-cvonly finishes: rematch CV, rewalk cut=30, then cut=25.

Reuse existing YZ/XTXW walk shards for mu_chi2mu_th=30 (already chi2_*_new).
Only rematch + rewalk CV for cut 30, then full chi2mu25 launch.
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
DFS = Path("/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs")
CV_RAW = DFS / "2026_09_19_035612__sel_all-mc-CV-updatecalo-cvonly"
CV_CAMPAIGN = CV_RAW.name
OUT30 = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE30 = OUT30 / "cache"
OUT25 = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod-chi2mu25")
N_SHARDS = 6
BATCH = 25
EXPECTED_DFS = 2000


def _log(msg: str) -> None:
    print(msg, flush=True)


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


def _pid_alive(pid: int) -> bool:
    """True if pid is a live (non-zombie) process."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            # comm may contain spaces/parens; state is after last ')'
            rest = fh.read().rsplit(")", 1)[-1].strip()
            state = rest.split(None, 1)[0] if rest else ""
        return state not in ("Z", "")
    except OSError:
        return False


def _wait_pids(pids: list[int], label: str, poll_s: int = 60) -> None:
    while True:
        alive = [p for p in pids if _pid_alive(p)]
        if not alive:
            # Reap zombies so they don't linger under this orchestrator.
            for p in pids:
                try:
                    os.waitpid(p, os.WNOHANG)
                except ChildProcessError:
                    pass
            _log(f"{label} all done")
            return
        _log(f"{label} waiting {len(alive)}/{len(pids)} alive")
        time.sleep(poll_s)


def _assert_done_logs(cache: Path, prefix: str, n: int) -> None:
    for sid in range(n):
        logp = cache / f"{prefix}_s{sid}.log"
        txt = logp.read_text() if logp.is_file() else ""
        if "done pot=" not in txt:
            raise SystemExit(f"incomplete: {logp}")


def cv_jobs_ready() -> tuple[bool, str]:
    n_df = len(list(CV_RAW.glob("*.df"))) if CV_RAW.is_dir() else 0
    n_cluster = -1  # unknown if jobsub_q unavailable
    try:
        q = subprocess.check_output(
            ["jobsub_q", "-G", "sbnd", "--user", "munjung"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
        n_cluster = sum(1 for ln in q.splitlines() if "29995450" in ln)
    except Exception as exc:
        # Token / jobsub hiccups: still proceed if enough dfs landed.
        msg = f"dfs={n_df}/{EXPECTED_DFS} cluster_jobs=? ({exc})"
        return n_df >= EXPECTED_DFS * 0.98, msg
    msg = f"dfs={n_df}/{EXPECTED_DFS} cluster_jobs={n_cluster}"
    ready = n_cluster == 0 and n_df >= EXPECTED_DFS * 0.98
    return ready, msg


def rematch_cv(*, force: bool = False) -> None:
    keys = CACHE30 / "wiremod_common_keys_sel_all.pkl"
    if not keys.is_file():
        raise SystemExit(f"missing common keys: {keys}")
    matched_cv = OUT30 / "matched" / "cv"
    out_pkl = CACHE30 / "wiremod_cv_all_files.pkl"

    if not force and matched_cv.is_dir():
        files = sorted(
            str(p)
            for p in matched_cv.glob("*_matched.df")
            if "sel_all" in p.name and "updatecalo-cvonly" in p.name
        )
        if len(files) >= EXPECTED_DFS * 0.98:
            with open(out_pkl, "wb") as fh:
                pickle.dump(files, fh, protocol=pickle.HIGHEST_PROTOCOL)
            _log(f"reuse matched cv n={len(files)} -> {out_pkl}")
            return

    archive = OUT30 / "matched" / f"cv_sep4_pre_chi2new_{time.strftime('%Y%m%d_%H%M%S')}"
    if matched_cv.exists() and not matched_cv.is_symlink():
        if archive.exists():
            shutil.rmtree(archive)
        matched_cv.rename(archive)
        _log(f"archived old matched cv -> {archive}")
    matched_cv.mkdir(parents=True, exist_ok=True)

    _kill("wiremod_match_shard_write.py", must_contain="updatecalo-cvonly")
    pids = []
    for sid in range(N_SHARDS):
        logp = CACHE30 / f"wiremod_match_cv_chi2new_s{sid}.log"
        pid = _launch(
            [
                PY, "-u", str(SCRIPTS / "wiremod_match_shard_write.py"),
                "--raw-dir", str(CV_RAW),
                "--matched-out-dir", str(matched_cv),
                "--common-keys-pkl", str(keys),
                "--variation-name", "cv",
                "--shard-id", str(sid),
                "--n-shards", str(N_SHARDS),
                "--drop-partials",
            ],
            logp,
        )
        pids.append(pid)
    _wait_pids(pids, "match_cv")
    files = sorted(
        str(p)
        for p in matched_cv.glob("*_matched.df")
        if "sel_all" in p.name
    )
    if len(files) < EXPECTED_DFS * 0.98:
        raise SystemExit(f"matched cv too few: {len(files)}")
    with open(out_pkl, "wb") as fh:
        pickle.dump(files, fh, protocol=pickle.HIGHEST_PROTOCOL)
    _log(f"wrote {out_pkl} n={len(files)}")


def rewalk_cv_cut30() -> None:
    _kill("wiremod_walk_cv_shard.py", must_contain=str(CACHE30))
    pids = []
    for sid in range(N_SHARDS):
        ck = CACHE30 / f"wiremod_walk_cv_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        pid = _launch(
            [
                PY, "-u", str(SCRIPTS / "wiremod_walk_cv_shard.py"),
                "--files-pkl", str(CACHE30 / "wiremod_cv_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid),
                "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH),
                "--rss-limit-gb", "20",
            ],
            CACHE30 / f"wiremod_walk_cv_s{sid}.log",
        )
        pids.append(pid)
    (CACHE30 / "wiremod_walk_cv_pids.txt").write_text(
        " " + " ".join(map(str, pids)) + "\n"
    )
    _wait_pids(pids, "walk_cv30")
    _assert_done_logs(CACHE30, "wiremod_walk_cv", N_SHARDS)


def merge_plots_cut30() -> None:
    for lab in ("yz", "xtxw"):
        _assert_done_logs(CACHE30, f"wiremod_walk_{lab}", N_SHARDS)
    yz = [str(CACHE30 / f"wiremod_walk_yz_s{i}.pkl") for i in range(N_SHARDS)]
    xtxw = [str(CACHE30 / f"wiremod_walk_xtxw_s{i}.pkl") for i in range(N_SHARDS)]
    cv = [str(CACHE30 / f"wiremod_walk_cv_s{i}.pkl") for i in range(N_SHARDS)]
    merge_log = CACHE30 / "wiremod_merge_walk_chi2new.log"
    plots_log = CACHE30 / "wiremod_make_plots_chi2new.log"
    _log("merge cut=30")
    with open(merge_log, "w") as lf:
        rc = subprocess.call(
            [
                PY, "-u", str(SCRIPTS / "wiremod_merge_walk_shards.py"),
                "--out-base", str(OUT30),
                "--yz-shard-ckpts", *yz,
                "--xtxw-shard-ckpts", *xtxw,
                "--cv-shard-ckpts", *cv,
                "--cv-campaign", CV_CAMPAIGN,
            ],
            cwd=str(REPO),
            stdout=lf,
            stderr=subprocess.STDOUT,
            env={**os.environ, "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}"},
        )
    if rc != 0:
        raise SystemExit(f"merge cut30 failed rc={rc}; see {merge_log}")
    _log("plots cut=30")
    with open(plots_log, "w") as lf:
        rc = subprocess.call(
            [
                PY, "-u", str(SCRIPTS / "wiremod_make_plots.py"),
                "--out-base", str(OUT30),
                "--skip-inspect",
            ],
            cwd=str(REPO),
            stdout=lf,
            stderr=subprocess.STDOUT,
            env={**os.environ, "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}"},
        )
    if rc != 0:
        raise SystemExit(f"plots cut30 failed rc={rc}; see {plots_log}")
    _log(f"cut30 done OUT={OUT30}")


def launch_cut25() -> None:
    _log(f"launching chi2mu25 into {OUT25}")
    # _launch_wiremod_chi2mu25 already points at new CV campaign in merge.
    rc = subprocess.call(
        [PY, "-u", str(SCRIPTS / "_launch_wiremod_chi2mu25.py")],
        cwd=str(REPO),
        env={
            **os.environ,
            "PYTHONPATH": f"{REPO}:{os.environ.get('PYTHONPATH', '')}",
            "BEARER_TOKEN_FILE": os.environ.get(
                "BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}"
            ),
        },
    )
    if rc != 0:
        raise SystemExit(f"chi2mu25 launcher failed rc={rc}")
    _log("chi2mu25 walks+supervisor launched")


def main() -> int:
    os.chdir(REPO)
    if not CV_RAW.is_dir():
        raise SystemExit(f"missing CV raw dir: {CV_RAW}")
    ready, msg = cv_jobs_ready()
    if not ready:
        raise SystemExit(f"CV jobs not ready: {msg}")
    _log(f"CV ready: {msg}")
    rematch_cv()
    rewalk_cv_cut30()
    merge_plots_cut30()
    launch_cut25()
    _log("orchestrator finished (cut25 still walking under supervisor)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
