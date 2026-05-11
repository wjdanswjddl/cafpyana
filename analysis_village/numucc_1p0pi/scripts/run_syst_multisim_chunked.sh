#!/usr/bin/env bash
#
# Chunked Flux / G4 / MCstat covariance driver (map → optional reduce).
# -----------------------------------------------------------------------------
# Phase 1: for each (syst, ``.df``) from ``dataset_locations.iter_multisim_chunk_tasks``,
#          run syst_multisim_chunk.py with ``--out_dir``:
#          ``$WORK_BASE/chunks/{Combined,MCstat}/`` (multisim_syst-chunked-*),
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
#   WORK_BASE         MCstat+Combined chunk root (default: default_multisim_syst_work_root)
#   G4_WORK_BASE      G4 chunk root (default: default_g4_syst_work_root); env NUMUCC_G4_SYST_WORK_BASE
#   FLUX_WORK_BASE    Flux chunk root (default: default_flux_syst_work_root); NUMUCC_FLUX_SYST_WORK_BASE
#   MC_DF_STAGE     final (default) or sel_all
#   VAR_SET         final | intermediate | both | sel_all
#                   (ignored when MC_DF_STAGE=sel_all; aggregator forces sel_all)
#   MULTISIM_SYST_TYPES  Comma-separated subset of {MCstat,Flux,G4} to run (default: all).
#                       Examples: all | Flux | Flux,G4 | MCstat,Flux
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
  --syst-types TYPES   Comma-separated subset of {MCstat,Flux,G4} (default: all).
                        Examples: all | Flux | Flux,G4 | MCstat,Flux
  --max-files|-n N     Run at most N map jobs (after syst filter). 0 = no limit.
  -h, --help           Show this message

Environment:
  MULTISIM_SYST_TYPES  Same as --syst-types (if provided, overrides default).
  MAX_FILES            Same as --max-files (CLI wins).
  WORK_BASE            multisim chunk root (Combined+MCstat); see script header for defaults.
  G4_WORK_BASE         G4 chunk root (default g4_syst-chunked-*); or set NUMUCC_G4_SYST_WORK_BASE for the Python default.
  FLUX_WORK_BASE       Flux chunk root (default flux_syst-chunked-*); or NUMUCC_FLUX_SYST_WORK_BASE.
EOF
}

_CLI_SYST_TYPES=""
_CLI_MAX_FILES=""
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

MULTISIM_SYST_TYPES="${_CLI_SYST_TYPES:-${MULTISIM_SYST_TYPES:-all}}"
_NEUTRINO_ORDER=('MCstat' 'Flux' 'G4')

if [[ "$MULTISIM_SYST_TYPES" == "all" || -z "$MULTISIM_SYST_TYPES" ]]; then
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
                echo "[multisim-run] ERROR: unknown syst type '${_s}' (expected MCstat,Flux,G4 or 'all')" >&2
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
CHUNKS_MULTISIM="$WORK_BASE/chunks"
CHUNKS_G4="$G4_WORK_BASE/chunks"
CHUNKS_FLUX="$FLUX_WORK_BASE/chunks"
_multisim_chunk_out_dir() {
    case "$1" in
        COMBINED) echo "${CHUNKS_MULTISIM}/Combined" ;;
        MCstat) echo "${CHUNKS_MULTISIM}/MCstat" ;;
        Flux) echo "${CHUNKS_FLUX}" ;;
        G4) echo "${CHUNKS_G4}" ;;
        *) echo "${CHUNKS_MULTISIM}" ;;
    esac
}
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

mkdir -p "${CHUNKS_MULTISIM}/Combined" "${CHUNKS_MULTISIM}/MCstat" "$CHUNKS_G4" "$CHUNKS_FLUX" "$SYST_DISK_ROOT"

chunk_py="$THIS_DIR/syst_multisim_chunk.py"
agg_py="$THIS_DIR/syst_multisim_aggregate.py"

echo "[multisim-run] WORK_BASE=$WORK_BASE  CHUNKS_MULTISIM=$CHUNKS_MULTISIM (Combined/, MCstat/)"
echo "[multisim-run] G4_WORK_BASE=$G4_WORK_BASE  CHUNKS_G4=$CHUNKS_G4"
echo "[multisim-run] FLUX_WORK_BASE=$FLUX_WORK_BASE  CHUNKS_FLUX=$CHUNKS_FLUX"
echo "[multisim-run] SYST_DISK_ROOT=$SYST_DISK_ROOT  (NPZs → MCstat/, Flux/, G4/ under this tree)"
echo "[multisim-run] MC_DF_STAGE=$MC_DF_STAGE  VAR_SET=$VAR_SET  INPUT_STAGE=$INPUT_STAGE  G4_MODE=$G4_MODE  FLUX_MODE=$FLUX_MODE  FLUX_KNOB_GROUPS=$FLUX_KNOB_GROUPS"
echo "[multisim-run] logging failures to $FAILED_LOG"
echo "[multisim-run] MAX_FILES=${MAX_FILES} (0 = no cap on map jobs)"

mapfile -t _ms_map_jobs < <(python3 -c "
import os, sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import iter_multisim_chunk_tasks
stage = os.environ.get('MC_DF_STAGE', 'final')
only_csv = os.environ.get('FILTER_ONLY_SYSTS_CSV', 'MCstat,Flux,G4')
only_systs = set(x.strip() for x in only_csv.split(',') if x.strip())
all_three = len(only_systs) == 3 and only_systs == {'MCstat','Flux','G4'}
for sn, p in iter_multisim_chunk_tasks(stage):
    if sn == 'COMBINED' or all_three or sn in only_systs:
        print('%s\t%s' % (sn, p))
")
_ms_map_n="${#_ms_map_jobs[@]}"
_ms_full_n="${_ms_map_n}"
if [[ "$MAX_FILES" =~ ^[0-9]+$ ]] && [[ "$MAX_FILES" -gt 0 ]] && [[ "$MAX_FILES" -lt "$_ms_map_n" ]]; then
    echo "[multisim-run] MAX_FILES=${MAX_FILES}: trimming queue from ${_ms_map_n} job(s)"
    _ms_map_jobs=( "${_ms_map_jobs[@]:0:$MAX_FILES}" )
    _ms_map_n="${#_ms_map_jobs[@]}"
fi
_ms_map_total="${_ms_map_n}"
echo "[multisim-run] chunk-map queue: ${_ms_map_total} (.df, syst) job(s) for MC_DF_STAGE=$MC_DF_STAGE"
if [[ "${_ms_full_n}" -ne "${_ms_map_total}" ]]; then
    echo "[multisim-run] (full filtered queue without cap would be ${_ms_full_n} job(s))"
fi
_ms_map_done=0

for ((_msi = 0; _msi < _ms_map_n; _msi++)); do
    line="${_ms_map_jobs[_msi]}"
    [[ -z "$line" ]] && continue
    IFS=$'\t' read -r syst f <<<"$line"
    if [[ -z "${syst:-}" ]] || [[ -z "${f:-}" ]]; then
        continue
    fi
    _ms_cur=$((_msi + 1))
    out_stem="$(basename "$f" .df)"
    CHUNK_OUT="$(_multisim_chunk_out_dir "$syst")"
    if [[ "$syst" == "COMBINED" ]]; then
        if [[ "${_ALL_SYSTS}" -eq 1 ]]; then
            out_pkl="$CHUNK_OUT/nu__${out_stem}.pkl"
        else
            # When --syst-names is a proper subset, chunk output becomes:
            #   nu__<MCstat_Flux_...>__<stem>.pkl
            out_pkl="$CHUNK_OUT/nu__${ONLY_SYSTS_TAG}__${out_stem}.pkl"
        fi
    else
        out_pkl="$CHUNK_OUT/nu__${syst}__${out_stem}.pkl"
    fi
    if [[ -f "$out_pkl" ]]; then
        ((_ms_map_done++)) || true
        echo "[multisim-run] progress map overall ${_ms_map_done}/${_ms_map_total}  job ${_ms_cur}/${_ms_map_n}  syst=$syst  (skip existing) $out_pkl"
        continue
    fi
    echo "[multisim-run] progress map overall $((_ms_map_done + 1))/${_ms_map_total}  job ${_ms_cur}/${_ms_map_n}  syst=$syst  BEGIN $(date -Is) df_file=$f"
    ok=0
    if [[ "$syst" == "COMBINED" ]]; then
        if [[ "${_ALL_SYSTS}" -eq 1 ]]; then
            if python "$chunk_py" \
                --df_file "$f" \
                --out_dir "$CHUNK_OUT" \
                --input-stage "$INPUT_STAGE" \
                --var-set "$VAR_SET" \
                --g4-mode "$G4_MODE" \
                --flux-mode "$FLUX_MODE" \
                --flux-knob-groups "$FLUX_KNOB_GROUPS"; then
                ok=1
            fi
        else
            if python "$chunk_py" \
                --df_file "$f" \
                --out_dir "$CHUNK_OUT" \
                --input-stage "$INPUT_STAGE" \
                --var-set "$VAR_SET" \
                --g4-mode "$G4_MODE" \
                --flux-mode "$FLUX_MODE" \
                --flux-knob-groups "$FLUX_KNOB_GROUPS" \
                --syst-names "$ONLY_SYSTS_CSV"; then
                ok=1
            fi
        fi
    else
        if python "$chunk_py" \
            --df_file "$f" \
            --out_dir "$CHUNK_OUT" \
            --input-stage "$INPUT_STAGE" \
            --var-set "$VAR_SET" \
            --g4-mode "$G4_MODE" \
            --flux-mode "$FLUX_MODE" \
            --flux-knob-groups "$FLUX_KNOB_GROUPS" \
            --syst-names "$syst"; then
            ok=1
        fi
    fi
    if [[ "$ok" -ne 1 ]]; then
        ((_ms_map_done++)) || true
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[multisim-run] progress map overall ${_ms_map_done}/${_ms_map_total}  job ${_ms_cur}/${_ms_map_n}  syst=$syst  FAILED $(date -Is) df_file=$f" >&2
        printf '%s\t%s\t%s\t%s\n' "$ts" "$MC_DF_STAGE" "$syst" "$f" >> "$FAILED_LOG"
        continue
    fi
    ((_ms_map_done++)) || true
    echo "[multisim-run] progress map overall ${_ms_map_done}/${_ms_map_total}  job ${_ms_cur}/${_ms_map_n}  syst=$syst  END $(date -Is) df_file=$f"
done

if [[ "${SKIP_AGGREGATE:-0}" == "1" ]]; then
    echo "[multisim-run] SKIP_AGGREGATE=1 — chunks: $CHUNKS_MULTISIM (Combined, MCstat) | $CHUNKS_G4 (G4) | $CHUNKS_FLUX (Flux)"
    exit 0
fi

# Aggregate unions all roots (dedupes by path). Optional legacy typed dirs under multisim.
_agg_chunks=(
    "$CHUNKS_MULTISIM"
    "$CHUNKS_G4"
    "$CHUNKS_FLUX"
)
[[ -d "${CHUNKS_MULTISIM}/G4" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/G4")
[[ -d "${CHUNKS_MULTISIM}/Flux" ]] && _agg_chunks+=("${CHUNKS_MULTISIM}/Flux")

agg_cmd=(
    python "$agg_py"
    --chunks_dir "${_agg_chunks[@]}"
    --syst-disk-root "$SYST_DISK_ROOT"
    --mc-df-stage "$MC_DF_STAGE"
    --var-set "$VAR_SET"
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
