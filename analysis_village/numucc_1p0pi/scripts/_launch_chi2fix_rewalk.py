#!/usr/bin/env python3
"""Launch YZ+XTXW full rewalks (chi2_new remap) then merge+plots."""
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


def _kill(pattern: str) -> None:
    try:
        out = subprocess.check_output(["pgrep", "-f", pattern], text=True)
    except subprocess.CalledProcessError:
        return
    for tok in out.split():
        try:
            os.kill(int(tok), signal.SIGTERM)
        except OSError:
            pass


def _launch(cmd, logpath: Path) -> int:
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
        },
    )
    print(f"launched {p.pid} {' '.join(cmd[2:4])}...", flush=True)
    return p.pid


def main() -> int:
    os.chdir(REPO)
    assert DROP_MAP.is_file(), DROP_MAP
    _kill("wiremod_walk_shard.py")
    _kill("_xtxw_merge_plots_supervisor")
    _kill("_chi2fix_supervisor")
    time.sleep(2)

    # refresh file lists without matplotlib
    for lab in ("yz", "xtxw"):
        files = sorted(
            str(p)
            for p in (OUT / "matched" / lab).glob("*_matched.df")
            if "sel_all" in p.name
        )
        with open(CACHE / f"wiremod_{lab}_all_files.pkl", "wb") as fh:
            pickle.dump(files, fh, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"{lab} files={len(files)}", flush=True)

    yz_pids = []
    for sid in range(N_SHARDS):
        ck = CACHE / f"wiremod_walk_yz_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        pid = _launch(
            [
                PY, "-u", str(SCRIPTS / "wiremod_walk_shard.py"),
                "--files-pkl", str(CACHE / "wiremod_yz_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid), "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH), "--rss-limit-gb", "20",
                "--label", "YZ",
            ],
            CACHE / f"wiremod_walk_yz_s{sid}.log",
        )
        yz_pids.append(pid)
    (CACHE / "wiremod_walk_yz_pids.txt").write_text(" " + " ".join(map(str, yz_pids)) + "\n")

    xtxw_pids = []
    for sid in range(N_SHARDS):
        ck = CACHE / f"wiremod_walk_xtxw_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        pid = _launch(
            [
                PY, "-u", str(SCRIPTS / "wiremod_walk_shard.py"),
                "--files-pkl", str(CACHE / "wiremod_xtxw_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid), "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH), "--rss-limit-gb", "20",
                "--label", "XTXW",
                "--drop-map-pkl", str(DROP_MAP),
            ],
            CACHE / f"wiremod_walk_xtxw_s{sid}.log",
        )
        xtxw_pids.append(pid)
    (CACHE / "wiremod_walk_xtxw_pids.txt").write_text(
        " " + " ".join(map(str, xtxw_pids)) + "\n"
    )

    # supervisor waits for both, then merge+plots
    sup = CACHE / "_chi2fix_supervisor.sh"
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
echo "chi2fix_supervisor start $(date -Is)" > "$SUP"

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

echo "waiting YZ $(date -Is)" >> "$SUP"
wait_pids "$CACHE/wiremod_walk_yz_pids.txt"
echo "YZ done $(date -Is)" >> "$SUP"
for sid in $(seq 0 $((N_SHARDS-1))); do
  grep -q 'done pot=' "$CACHE/wiremod_walk_yz_s${{sid}}.log" || {{ echo "YZ $sid incomplete" >> "$SUP"; exit 1; }}
done

echo "waiting XTXW $(date -Is)" >> "$SUP"
wait_pids "$CACHE/wiremod_walk_xtxw_pids.txt"
echo "XTXW done $(date -Is)" >> "$SUP"
for sid in $(seq 0 $((N_SHARDS-1))); do
  grep -q 'done pot=' "$CACHE/wiremod_walk_xtxw_s${{sid}}.log" || {{ echo "XTXW $sid incomplete" >> "$SUP"; exit 1; }}
done

# quick univ-diff sanity from products after merge would be better; check shard0
$PY - <<'PY'
import pickle, numpy as np
from pathlib import Path
ck=Path("{CACHE}")/"wiremod_walk_yz_s0.pkl"
with open(ck,"rb") as fh: st=pickle.load(fh)
acc=st["acc"]
bu=acc["by_universe"]
var="muon-p"
cv=np.asarray(bu["cv"]["hists_final"][var])
nd=sum(1 for u,p in bu.items() if u!="cv" and not np.array_equal(np.asarray(p["hists_final"][var]),cv))
print(f"YZ s0 muon-p univs differing from cv: {{nd}}/{{len(bu)-1}}")
if nd<1:
    raise SystemExit("YZ universes still identical after chi2 remap")
PY

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
tail -12 "$CACHE/wiremod_merge_walk.log" >> "$SUP"

echo "plots $(date -Is)" >> "$SUP"
: > "$CACHE/wiremod_make_plots.log"
$PY -u analysis_village/numucc_1p0pi/scripts/wiremod_make_plots.py \\
  --out-base "$OUT" --skip-inspect >> "$CACHE/wiremod_make_plots.log" 2>&1
echo "plots_rc=$? $(date -Is)" >> "$SUP"
tail -8 "$CACHE/wiremod_make_plots.log" >> "$SUP"
"""
    )
    sup.chmod(0o755)
    (CACHE / "wiremod_walk_supervisor_nohup.log").write_bytes(b"")
    sp = subprocess.Popen(
        ["bash", str(sup)],
        stdout=open(CACHE / "wiremod_walk_supervisor_nohup.log", "ab"),
        stderr=subprocess.STDOUT,
        start_new_session=True,
        cwd=str(REPO),
    )
    (CACHE / "wiremod_walk_supervisor.pid").write_text(f"{sp.pid}\n")
    print(f"supervisor {sp.pid}", flush=True)
    time.sleep(10)
    for label, pids in (("yz", yz_pids), ("xtxw", xtxw_pids)):
        alive = sum(Path(f"/proc/{p}").exists() for p in pids)
        print(f"{label} alive {alive}/{len(pids)}", flush=True)
    print((CACHE / "wiremod_walk_yz_s0.log").read_text()[:500], flush=True)
    print((CACHE / "wiremod_walk_supervisor.log").read_text()[:400], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
