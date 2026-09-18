#!/usr/bin/env python3
"""Launch XTXW dedupe + full YZ walk + supervisor (merge/plots after XTXW)."""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
CACHE = OUT / "cache"
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
SCRIPTS = REPO / "analysis_village" / "numucc_1p0pi" / "scripts"


def launch(cmd, logpath: Path, pidpath: Path | None = None) -> int:
    logpath.write_bytes(b"")
    lf = open(logpath, "ab")
    p = subprocess.Popen(
        cmd,
        stdout=lf,
        stderr=subprocess.STDOUT,
        cwd=str(REPO),
        start_new_session=True,
        env={**os.environ, "BEARER_TOKEN_FILE": os.environ.get("BEARER_TOKEN_FILE", f"/tmp/bt_u{os.getuid()}")},
    )
    if pidpath is not None:
        pidpath.write_text(f"{p.pid}\n")
    print(f"launched pid={p.pid} {' '.join(cmd[:4])} ...", flush=True)
    return p.pid


def main() -> int:
    os.chdir(REPO)
    subprocess.run(["pkill", "-f", "dedupe_matched_artkeys.py"], check=False)
    subprocess.run(["pkill", "-f", "wiremod_walk_shard.py"], check=False)
    subprocess.run(["pkill", "-f", "_wiremod_walk_supervisor_refix.sh"], check=False)
    time.sleep(2)

    mon = CACHE / "_mem_monitor_walk_shards.sh"
    mon.write_text(
        "#!/bin/bash\n"
        "PF=$1; OUTLOG=$2\n"
        "while true; do\n"
        "  alive=0\n"
        '  line="$(date -Is)"\n'
        '  for pid in $(cat "$PF" 2>/dev/null); do\n'
        "    if kill -0 $pid 2>/dev/null; then\n"
        "      alive=1\n"
        '      rss=$(ps -o rss= -p $pid | tr -d " ")\n'
        '      gb=$(awk -v k=${rss:-0} \'BEGIN{printf "%.2f", k/1024/1024}\')\n'
        '      line="$line pid=$pid rss=${gb}GiB"\n'
        "    fi\n"
        "  done\n"
        '  echo "$line" >> "$OUTLOG"\n'
        "  [ $alive -eq 0 ] && break\n"
        "  sleep 120\n"
        "done\n"
    )
    mon.chmod(0o755)

    dpid = launch(
        [
            PY,
            "-u",
            str(SCRIPTS / "dedupe_matched_artkeys.py"),
            "--matched-dir",
            str(OUT / "matched" / "xtxw"),
        ],
        CACHE / "wiremod_dedupe_xtxw.log",
        CACHE / "wiremod_dedupe_xtxw.pid",
    )

    yz_pids = []
    for sid in range(6):
        ck = CACHE / f"wiremod_walk_yz_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        pid = launch(
            [
                PY,
                "-u",
                str(SCRIPTS / "wiremod_walk_shard.py"),
                "--files-pkl",
                str(CACHE / "wiremod_yz_all_files.pkl"),
                "--checkpoint",
                str(ck),
                "--shard-id",
                str(sid),
                "--n-shards",
                "6",
                "--batch-size",
                "25",
                "--rss-limit-gb",
                "20",
                "--label",
                "YZ",
            ],
            CACHE / f"wiremod_walk_yz_s{sid}.log",
        )
        yz_pids.append(pid)
    (CACHE / "wiremod_walk_yz_pids.txt").write_text(" " + " ".join(map(str, yz_pids)) + "\n")

    (CACHE / "wiremod_mem_monitor_walk.log").write_bytes(b"")
    subprocess.Popen(
        ["bash", str(mon), str(CACHE / "wiremod_walk_yz_pids.txt"), str(CACHE / "wiremod_mem_monitor_walk.log")],
        start_new_session=True,
    )

    # copy supervisor into cache for local paths, keep repo script as source
    import shutil

    shutil.copy2(SCRIPTS / "_wiremod_walk_supervisor_refix.sh", CACHE / "_walk_supervisor_refix.sh")
    (CACHE / "_walk_supervisor_refix.sh").chmod(0o755)
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

    time.sleep(6)
    for label, pid in [("dedupe", dpid), ("sup", sp.pid)] + [(f"yz{i}", p) for i, p in enumerate(yz_pids)]:
        print(f"  {label} {pid} alive={Path(f'/proc/{pid}').exists()}", flush=True)
    print("--- dedupe ---", flush=True)
    print((CACHE / "wiremod_dedupe_xtxw.log").read_text()[:800], flush=True)
    print("--- yz0 ---", flush=True)
    print("\n".join((CACHE / "wiremod_walk_yz_s0.log").read_text().splitlines()[-8:]), flush=True)
    print("--- supervisor ---", flush=True)
    print((CACHE / "wiremod_walk_supervisor.log").read_text()[:400], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
