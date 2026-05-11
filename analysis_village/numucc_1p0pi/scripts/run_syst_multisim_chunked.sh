#!/usr/bin/env bash
#
# Chunked Flux / G4 / MCstat covariance driver (map → optional reduce).
# -----------------------------------------------------------------------------
# Phase 1: for each (syst, ``.df``) from ``dataset_locations.iter_multisim_chunk_tasks``,
#          run syst_multisim_chunk.py with ``--syst-names syst`` → ``nu__<syst>__<stem>.pkl``
#          (per-syst input dirs are configured in MULTISIM_SYST_GLOBS_FINAL / _SEL_ALL).
# Phase 2: syst_multisim_aggregate.py merges chunks and writes plots + legacy NPZs
#          (unless SKIP_AGGREGATE=1).
#
# Paths are owned by ``analysis_village.numucc_1p0pi.dataset_locations``:
#   MC_DF_STAGE=final|sel_all  → which glob map (MULTISIM_SYST_GLOBS_*).
#
# Environment:
#   WORK_BASE       Output root (default from dataset_locations.default_multisim_syst_work_root)
#   MC_DF_STAGE     final (default) or sel_all
#   VAR_SET         final | intermediate | both (passed to chunk + aggregate)
#   SKIP_AGGREGATE  1 to run phase 1 only
#   SKIP_COSMICS    1 → aggregate with --skip-cosmics (nu uncertainties only)
#   NO_PLOTS        1 → aggregate with --no-plots
#   NO_LEGACY_NPZ   1 → aggregate with --no-legacy-npz
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_multisim_syst_work_root
print(default_multisim_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR="$WORK_BASE/chunks"
# Neutrino multisim + cosmics NPZs use syst_disk_layout directly under WORK_BASE (MCstat/, Flux/, …).
SYST_DISK_ROOT="$WORK_BASE"
FAILED_LOG="$WORK_BASE/failed_multisim_df_files.log"

MC_DF_STAGE=${MC_DF_STAGE:-final}
VAR_SET=${VAR_SET:-final}

mkdir -p "$CHUNKS_DIR" "$SYST_DISK_ROOT"

chunk_py="$THIS_DIR/syst_multisim_chunk.py"
agg_py="$THIS_DIR/syst_multisim_aggregate.py"

echo "[multisim-run] WORK_BASE=$WORK_BASE  MC_DF_STAGE=$MC_DF_STAGE  VAR_SET=$VAR_SET"
echo "[multisim-run] logging failures to $FAILED_LOG"

while IFS=$'\t' read -r syst f; do
    [[ -z "${syst:-}" ]] && continue
    out_stem="$(basename "$f" .df)"
    if [[ "$syst" == "COMBINED" ]]; then
        out_pkl="$CHUNKS_DIR/nu__${out_stem}.pkl"
    else
        out_pkl="$CHUNKS_DIR/nu__${syst}__${out_stem}.pkl"
    fi
    if [[ -f "$out_pkl" ]]; then
        echo "[multisim-run] skip existing $out_pkl"
        continue
    fi
    ok=0
    if [[ "$syst" == "COMBINED" ]]; then
        if python "$chunk_py" \
            --df_file "$f" \
            --out_dir "$CHUNKS_DIR" \
            --var-set "$VAR_SET"; then
            ok=1
        fi
    else
        if python "$chunk_py" \
            --df_file "$f" \
            --out_dir "$CHUNKS_DIR" \
            --var-set "$VAR_SET" \
            --syst-names "$syst"; then
            ok=1
        fi
    fi
    if [[ "$ok" -ne 1 ]]; then
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[multisim-run] FAILED: syst=$syst file=$f" >&2
        printf '%s\t%s\t%s\t%s\n' "$ts" "$MC_DF_STAGE" "$syst" "$f" >> "$FAILED_LOG"
        continue
    fi
done < <(python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_multisim_chunk_tasks
stage = os.environ.get('MC_DF_STAGE', 'final')
for sn, p in iter_multisim_chunk_tasks(stage):
    print('%s\t%s' % (sn, p))
")

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[multisim-run] SKIP_AGGREGATE=1 — chunks only → $CHUNKS_DIR"
    exit 0
fi

agg_cmd=(
    python "$agg_py"
    --chunks_dir "$CHUNKS_DIR"
    --syst-disk-root "$SYST_DISK_ROOT"
    --mc-df-stage "$MC_DF_STAGE"
    --var-set "$VAR_SET"
)
[[ "${SKIP_COSMICS:-0}" == "1" ]] && agg_cmd+=(--skip-cosmics)
[[ "${NO_PLOTS:-0}" == "1" ]] && agg_cmd+=(--no-plots)
[[ "${NO_LEGACY_NPZ:-0}" == "1" ]] && agg_cmd+=(--no-legacy-npz)

echo "[multisim-run] aggregating: ${agg_cmd[*]}"
"${agg_cmd[@]}"

echo "[multisim-run] DONE syst_disk_layout tree → $SYST_DISK_ROOT (MCstat/, Flux/, G4/, Cosmics/ unless SKIP_COSMICS=1)"
echo "[multisim-run] Point consumers at NUMUCC_SYST_DISK_ROOT=$SYST_DISK_ROOT (needs GENIE/, Detector/ filled for utils.get_syst_unc)."
