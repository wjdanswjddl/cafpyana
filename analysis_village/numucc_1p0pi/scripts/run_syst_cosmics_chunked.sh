#!/usr/bin/env bash
#
# Cosmics systematic driver (offbeam + intime map → aggregate).
# -----------------------------------------------------------------------------
# Phase 1: ``syst_cosmics_chunk.py`` for every offbeam / intime ``.df``.
# Phase 2: ``syst_cosmics_aggregate.py`` → ``<SYST_DISK_ROOT>/Cosmics/cosmics_syst_dict.npz``.
#
# INPUT_STAGE controls which input-files glob the chunk job sees AND what the
# chunk script does:
#   sel_all (default) → ``dataset_locations.EVENT_SELECTION_GLOBS``: raw
#                       evt/trk/hdr dfs. Chunk re-runs the full event selection
#                       and saves per-stage cut-variable + final-variable hists.
#   final             → ``dataset_locations.SELECTED_EVENTS_GLOBS``: already
#                       final-selected dfs. Chunk just histograms the final
#                       variables (legacy behaviour).
#
# Usage:
#   run_syst_cosmics_chunked.sh [--input-stage final|sel_all] [SYST_DISK_ROOT]
#   run_syst_cosmics_chunked.sh [--input_stage final|sel_all] [SYST_DISK_ROOT]
#
# ``--input-stage`` / ``--input_stage`` overrides the ``INPUT_STAGE`` environment
# variable for this invocation (default when unset: ``sel_all``).
#
# Environment:
#   WORK_BASE          Chunk dir root (``default_cosmics_syst_work_root``)
#   CHUNKS_DIR         (default ``$WORK_BASE/chunks``)
#   NUMUCC_SYST_DISK_ROOT  Override for aggregate output (else first positional, else
#                          ``dataset_locations.default_syst_disk_root()``)
#   SKIP_AGGREGATE     1 → map phase only (no ``syst_cosmics_aggregate.py``)
#   INPUT_STAGE        sel_all (default) | final  (CLI flag wins if given)
#   VERBOSE            If non-empty: ``set -x``, ``python -u``, chunk ``--verbose``,
#                      EXIT trap with status, extra timestamps around each chunk.
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

# Parse optional flags; remaining args become ``$@`` (first positional = syst disk root).
INPUT_STAGE_CLI=""
_positional=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-stage|--input_stage)
            if [[ $# -lt 2 ]]; then
                echo "[cosmics-run] ERROR: $1 requires a value: final or sel_all" >&2
                exit 2
            fi
            INPUT_STAGE_CLI="$2"
            shift 2
            ;;
        -h|--help)
            cat <<'EOF'
Usage: run_syst_cosmics_chunked.sh [OPTIONS] [SYST_DISK_ROOT]

Options:
  --input-stage VALUE   final | sel_all (same as --input_stage)
  -h, --help            Show this message

Environment: INPUT_STAGE, NUMUCC_SYST_DISK_ROOT, SKIP_AGGREGATE, WORK_BASE, CHUNKS_DIR, VERBOSE.
Optional positional SYST_DISK_ROOT overrides env and ``default_syst_disk_root()`` from
``dataset_locations``.
EOF
            exit 0
            ;;
        -*)
            echo "[cosmics-run] ERROR: unknown option: $1 (use --help)" >&2
            exit 2
            ;;
        *)
            _positional+=("$1")
            shift
            ;;
    esac
done
if ((${#_positional[@]} > 0)); then
    set -- "${_positional[@]}"
else
    set --
fi

if [[ -n "${VERBOSE:-}" ]]; then
    set -x
    export PS4='+ [${BASH_SOURCE##*/}:${LINENO}] '
    trap 'echo "[cosmics-run] EXIT trap: status=$? at $(date -Is)" >&2' EXIT
fi

INPUT_STAGE="${INPUT_STAGE_CLI:-${INPUT_STAGE:-sel_all}}"
export INPUT_STAGE
case "$INPUT_STAGE" in
    sel_all|final) ;;
    *)
        echo "[cosmics-run] ERROR: INPUT_STAGE must be 'sel_all' or 'final', got '$INPUT_STAGE'" >&2
        exit 2
        ;;
esac

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_cosmics_syst_work_root
print(default_cosmics_syst_work_root('${TODAY}'))
")}
CHUNKS_DIR=${CHUNKS_DIR:-"$WORK_BASE/chunks"}
SYST_DISK_ROOT="${NUMUCC_SYST_DISK_ROOT:-${1:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_root
print(default_syst_disk_root())
")}}"

# Avoid mixing ``final`` and ``sel_all`` pickles under the same chunks dir.
# (Mixing breaks aggregation because payload schemas differ.)
CHUNKS_DIR="${CHUNKS_DIR%/}/${INPUT_STAGE}"

mkdir -p "$CHUNKS_DIR" "$SYST_DISK_ROOT"
echo "[cosmics-run] SYST_DISK_ROOT=$SYST_DISK_ROOT  (Cosmics/ aggregate target)"

chunk_py="$THIS_DIR/syst_cosmics_chunk.py"
agg_py="$THIS_DIR/syst_cosmics_aggregate.py"

echo "[cosmics-run] CHUNKS_DIR=$CHUNKS_DIR"
echo "[cosmics-run] REPO_ROOT=$REPO_ROOT"
echo "[cosmics-run] WORK_BASE=$WORK_BASE"
echo "[cosmics-run] INPUT_STAGE=$INPUT_STAGE"
echo "[cosmics-run] date=$(date -Is) host=$(hostname) pid=$$"
echo "[cosmics-run] ulimit -a (resource limits):"
ulimit -a 2>&1 | sed 's/^/[cosmics-run] /' || true
echo "[cosmics-run] python=$(command -v python3 || command -v python) $(python3 --version 2>&1 || true)"
if [[ -n "${VERBOSE:-}" ]]; then
    echo "[cosmics-run] VERBOSE is set: bash xtrace on, python -u, chunk --verbose"
    export PYTHONUNBUFFERED=1
fi

_cosmics_total_all="$(
    INPUT_STAGE="$INPUT_STAGE" python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_cosmics_chunk_df_paths
stage = os.environ.get('INPUT_STAGE', 'sel_all')
n = 0
for sample in ('offbeam', 'intime'):
    n += sum(1 for _ in iter_cosmics_chunk_df_paths(sample, input_stage=stage))
print(n)
"
)"
echo "[cosmics-run] input_stage=$INPUT_STAGE: ${_cosmics_total_all} input .df file(s) total (offbeam + intime)"

_cosmics_done_all=0
for sample in offbeam intime; do
    mapfile -t _cosmics_dfs < <(INPUT_STAGE="$INPUT_STAGE" python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_cosmics_chunk_df_paths
stage = os.environ.get('INPUT_STAGE', 'sel_all')
for p in iter_cosmics_chunk_df_paths('${sample}', input_stage=stage):
    print(p)
")
    _cosmics_n="${#_cosmics_dfs[@]}"
    echo "[cosmics-run] sample=$sample input_stage=$INPUT_STAGE: ${_cosmics_n} input .df file(s) queued"
    # Use a word-style ``for`` loop (not C-style ``for ((…))``) so older Bash never
    # mispairs ``done`` with the outer ``for sample`` loop.
    _cosmics_i=0
    for f in "${_cosmics_dfs[@]}"; do
        ((_cosmics_i++)) || true
        _cur=$_cosmics_i
        [[ -z "$f" ]] && continue
        out_stem="$(basename "$f" .df)"
        out_pkl="$CHUNKS_DIR/cosmics__${sample}__${out_stem}.pkl"
        if [[ -f "$out_pkl" ]]; then
            ((_cosmics_done_all++)) || true
            echo "[cosmics-run] progress map overall ${_cosmics_done_all}/${_cosmics_total_all}  sample=$sample ${_cur}/${_cosmics_n} (skip existing) $out_pkl"
            continue
        fi
        echo "[cosmics-run] progress map overall $((_cosmics_done_all + 1))/${_cosmics_total_all}  sample=$sample ${_cur}/${_cosmics_n} BEGIN $(date -Is) df_file=$f"
        if [[ -n "${VERBOSE:-}" ]]; then
            python -u "$chunk_py" --verbose --input-stage "$INPUT_STAGE" --sample "$sample" --df_file "$f" --out_dir "$CHUNKS_DIR"
        else
            python "$chunk_py" --input-stage "$INPUT_STAGE" --sample "$sample" --df_file "$f" --out_dir "$CHUNKS_DIR"
        fi
        ((_cosmics_done_all++)) || true
        echo "[cosmics-run] progress map overall ${_cosmics_done_all}/${_cosmics_total_all}  sample=$sample ${_cur}/${_cosmics_n} END $(date -Is) df_file=$f"
    done
done

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[cosmics-run] SKIP_AGGREGATE=1 — map only → $CHUNKS_DIR"
    exit 0
fi

echo "[cosmics-run] progress aggregate 1/1  BEGIN $(date -Is)"
python "$agg_py" --chunks_dir "$CHUNKS_DIR" --syst-disk-root "$SYST_DISK_ROOT"
echo "[cosmics-run] progress aggregate 1/1  END $(date -Is)"
echo "[cosmics-run] $(date -Is) DONE -> ${SYST_DISK_ROOT}/Cosmics/"
