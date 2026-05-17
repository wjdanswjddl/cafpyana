#!/usr/bin/env bash
#
# Chunked Flux / G4 / MCstat covariance driver (map → optional reduce).
# -----------------------------------------------------------------------------
# Phase 1: for each (syst, ``.df``) from ``dataset_locations.iter_multisim_chunk_tasks``,
#          run syst_multisim_chunk.py with ``--out_dir``:
#          ``$WORK_BASE/chunks/Combined/`` (multisim_syst-chunked-*),
#          ``$MCSTAT_WORK_BASE/chunks/`` (mcstat_syst-chunked-*),
#          ``$G4_WORK_BASE/chunks/`` (g4_syst-chunked-*), ``$FLUX_WORK_BASE/chunks/`` (flux_syst-chunked-*).
#          → ``nu__<syst>__<stem>.pkl`` (per-syst input dirs: MULTISIM_SYST_GLOBS_*).
# Phase 2: syst_multisim_aggregate.py merges chunks and writes plots + legacy NPZs
#          (unless SKIP_AGGREGATE=1).
#
# Paths are owned by ``analysis_village.numucc_1p0pi.dataset_locations``:
#   MC_DF_STAGE=final|sel_all  → which glob map (MULTISIM_SYST_GLOBS_*).
#
# MC_DF_STAGE also drives ``--input-stage`` on the chunk script:
#   final   → chunk reads only ``evt_{i}`` and runs the legacy
#             ``get_univ_rates``-based covariance on final-selected dfs.
#   sel_all → chunk reads ``evt+trk+hdr`` and re-runs the full numuCC 1p0pi
#             event-selection pipeline, recording per-universe rates at every
#             cut stage AND at the final stage. The aggregator auto-detects
#             this from the chunk pickles and switches its var catalogue.
#
# Environment:
#   WORK_BASE           Combined-only chunk root (default: default_multisim_syst_work_root)
#   MCSTAT_WORK_BASE    MCstat chunk root (default: default_mcstat_syst_work_root); NUMUCC_MCSTAT_SYST_WORK_BASE
#   G4_WORK_BASE        G4 chunk root (default: default_g4_syst_work_root); NUMUCC_G4_SYST_WORK_BASE
#   FLUX_WORK_BASE      Flux chunk root (default: default_flux_syst_work_root); NUMUCC_FLUX_SYST_WORK_BASE
#   MC_DF_STAGE     final (default) or sel_all
#   VAR_SET         final | intermediate | both | sel_all
#                   (ignored when MC_DF_STAGE=sel_all; aggregator forces sel_all)
#   MULTISIM_SYST_TYPES  Subset of {MCstat,Flux,G4}. Default ``all`` → Flux,G4 only (MCstat opt-in).
#                       Use ``full`` for MCstat+Flux+G4. Examples: Flux | MCstat,Flux | full
#   SKIP_AGGREGATE  1 to run phase 1 only (default: run aggregate)
#   NUMUCC_SYST_DISK_ROOT  Syst disk root (``default_syst_disk_root()`` reads this env var;
#                          aggregate always resolves to an absolute path).
#   NO_PLOTS        1 → aggregate with --no-plots
#   NO_LEGACY_NPZ   1 → aggregate with --no-legacy-npz
#   G4_MODE         knobs (default) | bundled — passed to syst_multisim_chunk.py --g4-mode.
#                   knobs uses ``makedf.g4syst.g4_systematics`` → (mc, knob, univ_i) columns.
#   FLUX_MODE       knobs (default) | bundled — passed as --flux-mode (makedf.bnbsyst flux knobs).
#   FLUX_KNOB_GROUPS  Default all → bnbsyst.regen_systematics; else e.g. beam,hadron,xsec
#                     (comma-separated BNB_FLUX_GROUPS keys), passed as --flux-knob-groups.
#   MAX_FILES       Cap map jobs after syst filter (0 = all). Overridden by
#                   ``--max-files N`` / ``-n N`` on the command line.
#   WORKERS         Map-phase process-pool size for ``syst_multisim_parallel.py``
#                   (default: ``min(nproc, 8)``). Overridden by ``--workers N`` /
#                   ``-j N``. Each worker calls ``syst_multisim_chunk.run_with_args``
#                   directly inside a forked process — no Python re-launch per file —
#                   so imports / module setup are paid ``WORKERS`` times instead of
#                   once per ``.df``. Per-file pickle outputs are identical to the
#                   serial path; ``skip-existing`` + the aggregator continue to work.
#
# Logs ``progress map overall`` / ``progress aggregate`` (k/total, BEGIN/END timestamps).
# -----------------------------------------------------------------------------
set -euo pipefail
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

# -----------------------------------------------------------------------------
# Optional systematic subset configuration (Flux/G4/MCstat only)
# -----------------------------------------------------------------------------
_usage() {
    cat <<'EOF'
run_syst_multisim_chunked.sh [OPTIONS]

Options:
  --syst-types TYPES   Subset of {MCstat,Flux,G4}. Default ``all`` → Flux,G4; ``full`` → all three.
                        Examples: all | full | Flux | MCstat,Flux
  --max-files|-n N     Run at most N map jobs (after syst filter). 0 = no limit.
  --workers|-j N       Number of parallel worker processes for the map phase
                       (default: min(nproc, 8)). Imports are shared via fork so this
                       is the main speed-up over the legacy per-file Python subprocess.
  -h, --help           Show this message

Environment:
  MULTISIM_SYST_TYPES  Same as --syst-types (default: all → Flux,G4; use full or MCstat,... for MCstat).
  MAX_FILES            Same as --max-files (CLI wins).
  WORKERS              Same as --workers (CLI wins).
  WORK_BASE            multisim chunk root (Combined/ only); see script header.
  MCSTAT_WORK_BASE     MCstat chunk root (default mcstat_syst-chunked-*); or NUMUCC_MCSTAT_SYST_WORK_BASE.
  G4_WORK_BASE         G4 chunk root (default g4_syst-chunked-*); or NUMUCC_G4_SYST_WORK_BASE.
  FLUX_WORK_BASE       Flux chunk root (default flux_syst-chunked-*); or NUMUCC_FLUX_SYST_WORK_BASE.
EOF
}

_CLI_SYST_TYPES=""
_CLI_MAX_FILES=""
_CLI_WORKERS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --syst-types)
            if [[ -z "${2:-}" ]]; then
                echo "[multisim-run] ERROR: --syst-types requires a value" >&2
                exit 2
            fi
            _CLI_SYST_TYPES="$2"
            shift 2
            ;;
        --max-files|-n)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "[multisim-run] ERROR: $1 requires a non-negative integer" >&2
                exit 2
            fi
            _CLI_MAX_FILES="$2"
            shift 2
            ;;
        --workers|-j)
            if [[ -z "${2:-}" ]] || ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]]; then
                echo "[multisim-run] ERROR: $1 requires a positive integer" >&2
                exit 2
            fi
            _CLI_WORKERS="$2"
            shift 2
            ;;
        -h|--help)
            _usage
            exit 0
            ;;
        *)
            echo "[multisim-run] ERROR: unknown option: $1 (use --help)" >&2
            exit 2
            ;;
    esac
done

MAX_FILES="${_CLI_MAX_FILES:-${MAX_FILES:-0}}"
_DEFAULT_WORKERS=$(python3 -c "import os; print(min(os.cpu_count() or 8, 8))")
WORKERS="${_CLI_WORKERS:-${WORKERS:-$_DEFAULT_WORKERS}}"

MULTISIM_SYST_TYPES="${_CLI_SYST_TYPES:-${MULTISIM_SYST_TYPES:-all}}"
MULTISIM_SYST_TYPES="$(echo "${MULTISIM_SYST_TYPES}" | tr '[:upper:]' '[:lower:]')"
_NEUTRINO_ORDER=('MCstat' 'Flux' 'G4')

if [[ "$MULTISIM_SYST_TYPES" == "all" || -z "$MULTISIM_SYST_TYPES" ]]; then
    ONLY_SYSTS=('Flux' 'G4')
elif [[ "${MULTISIM_SYST_TYPES}" == "full" ]]; then
    ONLY_SYSTS=("${_NEUTRINO_ORDER[@]}")
else
    _raw_norm=()
    IFS=',' read -ra _raw_systs <<< "${MULTISIM_SYST_TYPES}"
    for _s in "${_raw_systs[@]}"; do
        _sn="$(echo "${_s}" | tr '[:lower:]' '[:upper:]' )"
        case "${_sn}" in
            MCSTAT|MC_STAT|MC) _raw_norm+=('MCstat') ;;
            FLUX) _raw_norm+=('Flux') ;;
            G4) _raw_norm+=('G4') ;;
            *)
                echo "[multisim-run] ERROR: unknown syst type '${_s}' (expected MCstat,Flux,G4, all, or full)" >&2
                exit 2
                ;;
        esac
    done

    # De-dup while preserving set membership
    _uniq=()
    for _sn in "${_raw_norm[@]}"; do
        if [[ " ${_uniq[*]} " != *" ${_sn} "* ]]; then
            _uniq+=("${_sn}")
        fi
    done

    # Re-order to match syst_multisim_chunk.py (NEUTRINO_SYST_ORDER)
    ONLY_SYSTS=()
    for _sn in "${_NEUTRINO_ORDER[@]}"; do
        if [[ " ${_uniq[*]} " == *" ${_sn} "* ]]; then
            ONLY_SYSTS+=("${_sn}")
        fi
    done
fi

ONLY_SYSTS_CSV="$(IFS=','; echo "${ONLY_SYSTS[*]}")"
ONLY_SYSTS_TAG="$(IFS=_; echo "${ONLY_SYSTS[*]}")"

_ALL_SYSTS=0
if [[ "${ONLY_SYSTS_CSV}" == "MCstat,Flux,G4" ]]; then
    _ALL_SYSTS=1
fi

export FILTER_ONLY_SYSTS_CSV="${ONLY_SYSTS_CSV}"

TODAY=$(date +%Y%m%d)
WORK_BASE=${WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_multisim_syst_work_root
print(default_multisim_syst_work_root('${TODAY}'))
")}
G4_WORK_BASE=${G4_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_g4_syst_work_root
print(default_g4_syst_work_root('${TODAY}'))
")}
FLUX_WORK_BASE=${FLUX_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_flux_syst_work_root
print(default_flux_syst_work_root('${TODAY}'))
")}
MCSTAT_WORK_BASE=${MCSTAT_WORK_BASE:-$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_mcstat_syst_work_root
print(default_mcstat_syst_work_root('${TODAY}'))
")}
CHUNKS_MULTISIM="$WORK_BASE/chunks"
CHUNKS_MCSTAT="$MCSTAT_WORK_BASE/chunks"
CHUNKS_G4="$G4_WORK_BASE/chunks"
CHUNKS_FLUX="$FLUX_WORK_BASE/chunks"
SYST_DISK_ROOT="$(python3 -c "
import sys
from pathlib import Path
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_root
print(Path(default_syst_disk_root()).resolve())
")"
FAILED_LOG="$WORK_BASE/failed_multisim_df_files.log"

MC_DF_STAGE=${MC_DF_STAGE:-final}
VAR_SET=${VAR_SET:-final}
case "$MC_DF_STAGE" in
    final|sel_all) ;;
    *)
        echo "[multisim-run] ERROR: MC_DF_STAGE must be 'final' or 'sel_all', got '$MC_DF_STAGE'" >&2
        exit 2
        ;;
esac
export MC_DF_STAGE
INPUT_STAGE="$MC_DF_STAGE"

G4_MODE="${G4_MODE:-knobs}"
case "$G4_MODE" in
    knobs|bundled) ;;
    *)
        echo "[multisim-run] ERROR: G4_MODE must be 'knobs' or 'bundled', got '$G4_MODE'" >&2
        exit 2
        ;;
esac

FLUX_MODE="${FLUX_MODE:-knobs}"
case "$FLUX_MODE" in
    knobs|bundled) ;;
    *)
        echo "[multisim-run] ERROR: FLUX_MODE must be 'knobs' or 'bundled', got '$FLUX_MODE'" >&2
        exit 2
        ;;
esac

FLUX_KNOB_GROUPS="${FLUX_KNOB_GROUPS:-all}"

mkdir -p "${CHUNKS_MULTISIM}/Combined" "$CHUNKS_MCSTAT" "$CHUNKS_G4" "$CHUNKS_FLUX" "$SYST_DISK_ROOT"

parallel_py="$THIS_DIR/syst_multisim_parallel.py"
agg_py="$THIS_DIR/syst_multisim_aggregate.py"

echo "[multisim-run] WORK_BASE=$WORK_BASE  CHUNKS_MULTISIM=$CHUNKS_MULTISIM (Combined/)"
echo "[multisim-run] MCSTAT_WORK_BASE=$MCSTAT_WORK_BASE  CHUNKS_MCSTAT=$CHUNKS_MCSTAT"
echo "[multisim-run] G4_WORK_BASE=$G4_WORK_BASE  CHUNKS_G4=$CHUNKS_G4"
echo "[multisim-run] FLUX_WORK_BASE=$FLUX_WORK_BASE  CHUNKS_FLUX=$CHUNKS_FLUX"
echo "[multisim-run] SYST_DISK_ROOT=$SYST_DISK_ROOT  (NPZs → MCstat/, Flux/, G4/ under this tree)"
echo "[multisim-run] MC_DF_STAGE=$MC_DF_STAGE  VAR_SET=$VAR_SET  INPUT_STAGE=$INPUT_STAGE  G4_MODE=$G4_MODE  FLUX_MODE=$FLUX_MODE  FLUX_KNOB_GROUPS=$FLUX_KNOB_GROUPS"
echo "[multisim-run] logging failures to $FAILED_LOG"
echo "[multisim-run] MAX_FILES=${MAX_FILES} (0 = no cap on map jobs)"
echo "[multisim-run] WORKERS=${WORKERS} parallel chunk-map worker processes"

# ``syst_multisim_parallel.py`` imports the chunk code once in master and forks
# workers, dispatching one (syst, .df) job per worker at a time via
# multiprocessing.Pool. Pickle outputs and skip-existing semantics match the
# legacy per-file invocation exactly, so the aggregator step below is unchanged.
parallel_args=(
    --mc-df-stage "$MC_DF_STAGE"
    --var-set "$VAR_SET"
    --syst-types "$MULTISIM_SYST_TYPES"
    --max-files "$MAX_FILES"
    --workers "$WORKERS"
    --chunks-multisim "$CHUNKS_MULTISIM"
    --chunks-mcstat "$CHUNKS_MCSTAT"
    --chunks-flux "$CHUNKS_FLUX"
    --chunks-g4 "$CHUNKS_G4"
    --failed-log "$FAILED_LOG"
    --g4-mode "$G4_MODE"
    --flux-mode "$FLUX_MODE"
    --flux-knob-groups "$FLUX_KNOB_GROUPS"
)

echo "[multisim-run] progress map BEGIN $(date -Is) (workers=${WORKERS})"
set +e
python3 "$parallel_py" "${parallel_args[@]}"
_ms_map_rc=$?
set -e
echo "[multisim-run] progress map END $(date -Is) rc=${_ms_map_rc}"
if [[ "$_ms_map_rc" -ne 0 ]]; then
    # Non-zero only when ALL jobs failed — see syst_multisim_parallel.main().
    echo "[multisim-run] map phase failed entirely (rc=${_ms_map_rc}); aborting before aggregate" >&2
    exit "$_ms_map_rc"
fi

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[multisim-run] SKIP_AGGREGATE=1 — chunks: $CHUNKS_MULTISIM (Combined) | $CHUNKS_MCSTAT (MCstat) | $CHUNKS_G4 (G4) | $CHUNKS_FLUX (Flux)"
    exit 0
fi

# Aggregate unions all roots (dedupes by path). Optional legacy typed dirs under multisim.
_agg_chunks=(
    "$CHUNKS_MULTISIM"
    "$CHUNKS_MCSTAT"
    "$CHUNKS_G4"
    "$CHUNKS_FLUX"
)
[[ -d "${CHUNKS_MULTISIM}/G4" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/G4")
[[ -d "${CHUNKS_MULTISIM}/Flux" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/Flux")
[[ -d "${CHUNKS_MULTISIM}/MCstat" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/MCstat")

agg_cmd=(
    python "$agg_py"
    --chunks_dir "${_agg_chunks[@]}"
    --syst-disk-root "$SYST_DISK_ROOT"
    --mc-df-stage "$MC_DF_STAGE"
    --var-set "$VAR_SET"
    --syst-types "$ONLY_SYSTS_CSV"
)
[[ "${NO_PLOTS:-0}" == "1" ]] && agg_cmd+=(--no-plots)
[[ "${NO_LEGACY_NPZ:-0}" == "1" ]] && agg_cmd+=(--no-legacy-npz)

echo "[multisim-run] progress aggregate 1/1  BEGIN $(date -Is) ${agg_cmd[*]}"
"${agg_cmd[@]}"
echo "[multisim-run] progress aggregate 1/1  END $(date -Is)"

echo "[multisim-run] DONE syst_disk_layout tree → $SYST_DISK_ROOT (MCstat/, Flux/, G4/)"
if [[ "${_ALL_SYSTS}" -eq 1 ]]; then
    echo "[multisim-run] computed sources: MCstat, Flux, G4"
else
    echo "[multisim-run] computed sources subset: ${ONLY_SYSTS[*]}"
fi
echo "[multisim-run] Point consumers at NUMUCC_SYST_DISK_ROOT=$SYST_DISK_ROOT (export to match this run, or rely on the same default from dataset_locations.default_syst_disk_root)."
