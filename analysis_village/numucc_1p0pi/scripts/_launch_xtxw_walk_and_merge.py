#!/usr/bin/env python3
"""Launch fresh XTXW walk shards with drop-map (after YZ + drop-map already done)."""
from __future__ import annotations

import os
import pickle
import signal
import subprocess
import time
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"
DROP_MAP = CACHE / "wiremod_xtxw_drop_map.pkl"
N_SHARDS = 6
BATCH = 25


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
    assert DROP_MAP.is_file(), f"missing {DROP_MAP}"

    _kill_matching("_walk_supervisor_refix.sh")
    _kill_matching("_wiremod_walk_supervisor_refix.sh")
    _kill_matching("wiremod_walk_shard.py --label XTXW")
    _kill_matching("wiremod_walk_shard.py.*XTXW")
    time.sleep(2)

    # refresh file list without importing syst_detvar_common (matplotlib)
    matched = sorted(
        str(p)
        for p in (OUT / "matched" / "xtxw").glob("*_matched.df")
        if "sel_all" in p.name
    )
    with open(CACHE / "wiremod_xtxw_all_files.pkl", "wb") as fh:
        pickle.dump(matched, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"xtxw files={len(matched)}", flush=True)

    pids = []
    for sid in range(N_SHARDS):
        log = CACHE / f"wiremod_walk_xtxw_s{sid}.log"
        ck = CACHE / f"wiremod_walk_xtxw_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        log.write_bytes(b"")
        lf = open(log, "ab")
        p = subprocess.Popen(
            [
                PY,
                "-u",
                str(SCRIPTS / "wiremod_walk_shard.py"),
                "--files-pkl",
                str(CACHE / "wiremod_xtxw_all_files.pkl"),
                "--checkpoint",
                str(ck),
                "--shard-id",
                str(sid),
                "--n-shards",
                str(N_SHARDS),
                "--batch-size",
                str(BATCH),
                "--rss-limit-gb",
                "20",
                "--label",
                "XTXW",
                "--drop-map-pkl",
                str(DROP_MAP),
            ],
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
        pids.append(p.pid)
        print(f"XTXW s{sid} pid={p.pid}", flush=True)

    (CACHE / "wiremod_walk_xtxw_pids.txt").write_text(" " + " ".join(map(str, pids)) + "\n")

    # merge+plots supervisor that only waits on XTXW
    sup = CACHE / "_xtxw_merge_plots_supervisor.sh"
    sup.write_text(
        f"""#!/bin/bash
set -e
CACHE="{CACHE}"
OUT="{OUT}"
PY="{PY}"
N_SHARDS={N_SHARDS}
SUP=$CACHE/wiremod_walk_supervisor.log
cd {REPO}
export BEARER_TOKEN_FILE=${{BEARER_TOKEN_FILE:-/tmp/bt_u$(id -u)}}
echo "xtxw_merge_supervisor start $(date -Is)" >> "$SUP"

wait_pids() {{
  local pf=$1
  while true; do
    alive=0
    for pid in $(cat "$pf" 2>/dev/null); do
      kill -0 "$pid" 2>/dev/null && alive=1
    done
    [ "$alive" -eq 0 ] && return 0
    sleep 60
  done
}}

echo "waiting XTXW $(date -Is)" >> "$SUP"
wait_pids "$CACHE/wiremod_walk_xtxw_pids.txt"
echo "XTXW done $(date -Is)" >> "$SUP"
for sid in $(seq 0 $((N_SHARDS-1))); do
  grep -q 'done pot=' "$CACHE/wiremod_walk_xtxw_s${{sid}}.log" || {{ echo "XTXW $sid incomplete" >> "$SUP"; exit 1; }}
done

YZ_CKPTS=""; XTXW_CKPTS=""; CV_CKPTS=""
for sid in $(seq 0 $((N_SHARDS-1))); do
  YZ_CKPTS="$YZ_CKPTS $CACHE/wiremod_walk_yz_s${{sid}}.pkl"
  XTXW_CKPTS="$XTXW_CKPTS $CACHE/wiremod_walk_xtxw_s${{sid}}.pkl"
  CV_CKPTS="$CV_CKPTS $CACHE/wiremod_walk_cv_s${{sid}}.pkl"
done
echo "merge $(date -Is)" >> "$SUP"
: > "$CACHE/wiremod_merge_walk.log"
# shellcheck disable=SC2086
$PY -u analysis_village/numucc_1p0pi/scripts/wiremod_merge_walk_shards.py \\
  --out-base "$OUT" \\
  --yz-shard-ckpts $YZ_CKPTS \\
  --xtxw-shard-ckpts $XTXW_CKPTS \\
  --cv-shard-ckpts $CV_CKPTS \\
  --cv-campaign 2026_09_04_172912__sel_all-mc-CV >> "$CACHE/wiremod_merge_walk.log" 2>&1
echo "merge_rc=$? " >> "$SUP"
tail -10 "$CACHE/wiremod_merge_walk.log" >> "$SUP"

echo "plots $(date -Is)" >> "$SUP"
: > "$CACHE/wiremod_make_plots.log"
$PY -u analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py \\
  --out-base "$OUT" --skip-inspect >> "$CACHE/wiremod_make_plots.log" 2>&1
echo "plots_rc=$? $(date -Is)" >> "$SUP"
tail -8 "$CACHE/wiremod_make_plots.log" >> "$SUP"
"""
    )
    sup.chmod(0o755)
    sp = subprocess.Popen(
        ["bash", str(sup)],
        stdout=open(CACHE / "wiremod_walk_supervisor_nohup.log", "ab"),
        stderr=subprocess.STDOUT,
        start_new_session=True,
        cwd=str(REPO),
    )
    (CACHE / "wiremod_walk_supervisor.pid").write_text(f"{sp.pid}\n")
    print(f"merge/plots supervisor pid={sp.pid}", flush=True)

    time.sleep(8)
    for sid, pid in enumerate(pids):
        print(f"  s{sid} alive={Path(f'/proc/{pid}').exists()}", flush=True)
    print((CACHE / f"wiremod_walk_xtxw_s0.log").read_text()[:600], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
