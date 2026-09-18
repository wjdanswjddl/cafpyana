#!/usr/bin/env python3
"""Launch WireMod selection walk with mu chi2 cut = 25 into a parallel out dir.

Reuses matched sel_all files from the chi2mu=30 WireMod tree (matching is
pre-selection). Writes products/NPZs/plots under::

  .../systematics-final/WireMod-chi2mu25

without touching the original ``WireMod`` outputs.
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

SRC = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod")
OUT = Path("/exp/sbnd/data/users/munjung/plots/numucc1p0pi/systematics-final/WireMod-chi2mu25")
CACHE = OUT / "cache"
REPO = Path("/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")
PY = "/exp/sbnd/app/users/munjung/env/bin/python"
SCRIPTS = REPO / "analysis_village/numucc_1p0pi/scripts"
DROP_MAP_SRC = SRC / "cache" / "wiremod_xtxw_drop_map.pkl"
N_SHARDS = 6
BATCH = 25
MU_CHI2MU_TH = 25.0


def _kill(pattern: str) -> None:
    try:
        out = subprocess.check_output(["pgrep", "-af", pattern], text=True)
    except subprocess.CalledProcessError:
        return
    for line in out.splitlines():
        # only kill jobs targeting this OUT tree
        if str(OUT) not in line and "chi2mu25" not in line:
            continue
        tok = line.split(None, 1)[0]
        try:
            os.kill(int(tok), signal.SIGTERM)
        except (OSError, ValueError):
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
    print(f"launched {p.pid} {' '.join(cmd[2:5])}...", flush=True)
    return p.pid


def _setup_tree() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    matched = OUT / "matched"
    if matched.is_symlink() or matched.exists():
        pass
    else:
        matched.symlink_to(SRC / "matched")
        print(f"symlinked {matched} -> {SRC / 'matched'}", flush=True)

    note = OUT / "README_chi2mu25.txt"
    note.write_text(
        "WireMod envelope products with muon chi2 cut mu_chi2mu_th=25\n"
        f"(canonical default MU_CHI2MU_TH=30 lives under {SRC}).\n"
        "Matched sel_all files are reused via symlink; only the selection walk\n"
        "and envelope/NPZ/plots are recomputed here.\n"
    )

    if not DROP_MAP_SRC.is_file():
        raise SystemExit(f"missing drop map: {DROP_MAP_SRC}")
    drop_dst = CACHE / "wiremod_xtxw_drop_map.pkl"
    if not drop_dst.is_file() or drop_dst.stat().st_mtime < DROP_MAP_SRC.stat().st_mtime:
        shutil.copy2(DROP_MAP_SRC, drop_dst)
        print(f"copied drop map -> {drop_dst}", flush=True)

    for lab in ("yz", "xtxw", "cv"):
        src_pkl = SRC / "cache" / f"wiremod_{lab}_all_files.pkl"
        dst_pkl = CACHE / f"wiremod_{lab}_all_files.pkl"
        if src_pkl.is_file():
            shutil.copy2(src_pkl, dst_pkl)
            with open(dst_pkl, "rb") as fh:
                n = len(pickle.load(fh))
            print(f"copied {lab} file list n={n}", flush=True)
        else:
            files = sorted(
                str(p)
                for p in (OUT / "matched" / lab).glob("*_matched.df")
                if "sel_all" in p.name
            )
            with open(dst_pkl, "wb") as fh:
                pickle.dump(files, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"built {lab} file list n={len(files)}", flush=True)


def main() -> int:
    os.chdir(REPO)
    # Avoid importing selections.py (heavy deps); confirm canonical default from source text.
    sel_py = REPO / "analysis_village/numucc_1p0pi/makedf/selections.py"
    sel_txt = sel_py.read_text()
    if "MU_CHI2MU_TH  = 30" not in sel_txt and "MU_CHI2MU_TH = 30" not in sel_txt:
        raise SystemExit(f"expected MU_CHI2MU_TH=30 in {sel_py}")
    print(f"canonical MU_CHI2MU_TH=30 (selections.py); this run uses {MU_CHI2MU_TH}", flush=True)

    _kill("wiremod_walk_shard.py")
    _kill("wiremod_walk_cv_shard.py")
    _kill("_chi2mu25_supervisor")
    time.sleep(2)
    _setup_tree()

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
                "--mu-chi2mu-th", str(MU_CHI2MU_TH),
            ],
            CACHE / f"wiremod_walk_yz_s{sid}.log",
        )
        yz_pids.append(pid)
    (CACHE / "wiremod_walk_yz_pids.txt").write_text(
        " " + " ".join(map(str, yz_pids)) + "\n"
    )

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
                "--drop-map-pkl", str(CACHE / "wiremod_xtxw_drop_map.pkl"),
                "--mu-chi2mu-th", str(MU_CHI2MU_TH),
            ],
            CACHE / f"wiremod_walk_xtxw_s{sid}.log",
        )
        xtxw_pids.append(pid)
    (CACHE / "wiremod_walk_xtxw_pids.txt").write_text(
        " " + " ".join(map(str, xtxw_pids)) + "\n"
    )

    cv_pids = []
    for sid in range(N_SHARDS):
        ck = CACHE / f"wiremod_walk_cv_s{sid}.pkl"
        if ck.exists():
            ck.unlink()
        pid = _launch(
            [
                PY, "-u", str(SCRIPTS / "wiremod_walk_cv_shard.py"),
                "--files-pkl", str(CACHE / "wiremod_cv_all_files.pkl"),
                "--checkpoint", str(ck),
                "--shard-id", str(sid), "--n-shards", str(N_SHARDS),
                "--batch-size", str(BATCH), "--rss-limit-gb", "20",
                "--mu-chi2mu-th", str(MU_CHI2MU_TH),
            ],
            CACHE / f"wiremod_walk_cv_s{sid}.log",
        )
        cv_pids.append(pid)
    (CACHE / "wiremod_walk_cv_pids.txt").write_text(
        " " + " ".join(map(str, cv_pids)) + "\n"
    )

    sup = CACHE / "_chi2mu25_supervisor.sh"
    sup.write_text(
        f"""#!/bin/bash
set -e
CACHE="{CACHE}"
OUT="{OUT}"
PY="{PY}"
N_SHARDS={N_SHARDS}
MU_TH={MU_CHI2MU_TH}
SUP=$CACHE/wiremod_walk_supervisor.log
cd {REPO}
export BEARER_TOKEN_FILE=${{BEARER_TOKEN_FILE:-/tmp/bt_u$(id -u)}}
echo "chi2mu25_supervisor start $(date -Is) mu_chi2mu_th=$MU_TH" > "$SUP"

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

for lab in yz xtxw cv; do
  echo "waiting $lab $(date -Is)" >> "$SUP"
  wait_pids "$CACHE/wiremod_walk_${{lab}}_pids.txt"
  echo "$lab done $(date -Is)" >> "$SUP"
  for sid in $(seq 0 $((N_SHARDS-1))); do
    grep -q 'done pot=' "$CACHE/wiremod_walk_${{lab}}_s${{sid}}.log" \\
      || {{ echo "$lab $sid incomplete" >> "$SUP"; exit 1; }}
  done
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
  --cv-campaign 2026_09_04_172912__sel_all-mc-CV \\
  --mu-chi2mu-th "$MU_TH" >> "$CACHE/wiremod_merge_walk.log" 2>&1
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
    time.sleep(8)
    for label, pids in (("yz", yz_pids), ("xtxw", xtxw_pids), ("cv", cv_pids)):
        alive = sum(Path(f"/proc/{p}").exists() for p in pids)
        print(f"{label} alive {alive}/{len(pids)}", flush=True)
    print((CACHE / "wiremod_walk_yz_s0.log").read_text()[:600], flush=True)
    print((CACHE / "wiremod_walk_supervisor.log").read_text()[:400], flush=True)
    print(f"OUT={OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
