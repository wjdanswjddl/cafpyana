# Shared helpers for chunked event-selection map/reduce drivers.
# Source from run_event_selection_chunked.sh or run_event_selection_paths.sh — do not execute directly.

: "${THIS_DIR:?THIS_DIR must be set before sourcing event_selection_chunked_lib.sh}"
REPO_ROOT="$(cd "$THIS_DIR/../../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:$PYTHONPATH}"

# ---- defaults (override via env before sourcing) ---------------------------------
TODAY="${TODAY:-$(date +%Y%m%d)}"
WORK_BASE="${WORK_BASE:-/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-chunked-${TODAY}}"
CHUNKS_DIR="${CHUNKS_DIR:-$WORK_BASE/chunks}"
PLOTS_DIR="${PLOTS_DIR:-$WORK_BASE/plots}"
FAILED_LOG="${FAILED_LOG:-$WORK_BASE/failed_df_files.log}"

# MC universe weights accumulated in chunk pickles → hatched syst bands in aggregate.
MC_UNIV_SYST="${MC_UNIV_SYST:-Flux,G4,GENIE}"

# NPZ fallback when a variable is missing from chunked universes (utils.get_syst_unc).
if [[ -z "${NUMUCC_SYST_DISK_ROOT:-}" ]]; then
    NUMUCC_SYST_DISK_ROOT="$(python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import default_syst_disk_root
print(default_syst_disk_root().expanduser().resolve())
")"
fi
export NUMUCC_SYST_DISK_ROOT

# Extra flags forwarded to event_selection_aggregate.py (e.g. --cosmic_estimate offbeam).
AGG_EXTRA_ARGS="${AGG_EXTRA_ARGS:-}"

# Expand EVENT_SELECTION_GLOBS from dataset_locations.py (all *.df per sample directory).
# Override directories in dataset_locations.py or NUMUCC_SPRING_GEN1_ROOT before running.
# Prints one "sample|/path/to/file.df" per line (stable sort via sorted_glob).
event_selection_discover_jobs_from_dataset_locations() {
    python3 -c "
import sys
sys.path.insert(0, '${REPO_ROOT}')
from analysis_village.numucc_1p0pi.dataset_locations import EVENT_SELECTION_GLOBS, sorted_glob
for sample in ('mc', 'data', 'intime', 'offbeam', 'dirt'):
    pattern = EVENT_SELECTION_GLOBS[sample]
    paths = sorted_glob(pattern)
    if not paths:
        print('[discover] sample=%s no files matched %s' % (sample, pattern), file=sys.stderr)
        continue
    print('[discover] sample=%s  %d files  (%s)' % (sample, len(paths), pattern), file=sys.stderr)
    for p in paths:
        print('%s|%s' % (sample, p))
"
}

event_selection_chunk_one() {
    local sample="$1"
    local f="$2"
    local out_pkl_base out_pkl chunk_args=()

    out_pkl_base="$(basename "$f" .df)"
    out_pkl="$CHUNKS_DIR/${sample}__${out_pkl_base}.pkl"

    if [[ -f "$out_pkl" ]]; then
        echo "[run] $out_pkl exists; skipping"
        return 0
    fi

    if [[ -n "$MC_UNIV_SYST" && "$sample" == "mc" ]]; then
        chunk_args+=(--mc-univ-syst "$MC_UNIV_SYST")
    fi

    if ! python "$THIS_DIR/event_selection_chunk.py" \
        --df_file "$f" \
        --sample "$sample" \
        --out_dir "$CHUNKS_DIR" \
        "${chunk_args[@]}"; then
        local ts
        ts="$(date '+%Y-%m-%d %H:%M:%S')"
        echo "[run] FAILED chunk (see $FAILED_LOG): $f" >&2
        printf '%s\t%s\t%s\n' "$ts" "$sample" "$f" >>"$FAILED_LOG"
        return 1
    fi
    return 0
}

event_selection_run_map() {
    mkdir -p "$CHUNKS_DIR"
    echo "[run] WORK_BASE=$WORK_BASE"
    echo "[run] CHUNKS_DIR=$CHUNKS_DIR"
    echo "[run] MC_UNIV_SYST=${MC_UNIV_SYST:-<disabled>}"
    echo "[run] NUMUCC_SYST_DISK_ROOT=$NUMUCC_SYST_DISK_ROOT"
    echo "[run] logging chunk failures to $FAILED_LOG"

    local sample f
    for entry in "$@"; do
        sample="${entry%%|*}"
        f="${entry#*|}"
        if [[ ! -f "$f" ]]; then
            echo "[run] sample=$sample missing file: $f -- skipping" >&2
            continue
        fi
        echo "[run] sample=$sample file=$f"
        event_selection_chunk_one "$sample" "$f" || true
    done
}

event_selection_run_aggregate() {
    mkdir -p "$PLOTS_DIR"
    echo "[run] aggregating + plotting -> $PLOTS_DIR"
    local agg_syst_args=()
    if [[ -d "$NUMUCC_SYST_DISK_ROOT" ]]; then
        agg_syst_args=(--syst-disk-root "$NUMUCC_SYST_DISK_ROOT")
    else
        echo "[run] WARN: syst disk not found ($NUMUCC_SYST_DISK_ROOT); overlay bands from MC universes in chunks only"
    fi
    # shellcheck disable=SC2086
    python "$THIS_DIR/event_selection_aggregate.py" \
        --in_dir "$CHUNKS_DIR" \
        --out_dir "$PLOTS_DIR" \
        "${agg_syst_args[@]}" \
        $AGG_EXTRA_ARGS
    echo "[run] DONE plots -> $PLOTS_DIR"
    echo "[run] merged_histdata.pkl -> $PLOTS_DIR/merged_histdata.pkl"
}
