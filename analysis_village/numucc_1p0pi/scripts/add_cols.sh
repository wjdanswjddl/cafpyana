#!/usr/bin/env bash
# Enrich split .df files for one sample subdirectory under cafpyana_out/dfs.
#
# Usage:
#   ./add_cols.sh <this_subdir>
#   ./add_cols.sh 2026_09_03_032705__sel_all-mc-BNB_cosmics-detvar_DENT
#
# Extra args after the subdir are forwarded to update_dfs_add_columns.py
# (e.g. --limit 1 --dry-run --process-only).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

base_dir="/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <this_subdir> [extra args for update_dfs_add_columns.py...]" >&2
  echo "Example: $0 2026_09_03_032705__sel_all-mc-BNB_cosmics-detvar_DENT" >&2
  exit 1
fi

this_subdir="$1"
shift

dir="${base_dir}/${this_subdir}"
if [[ ! -d "$dir" ]]; then
  echo "ERROR: directory not found: $dir" >&2
  exit 1
fi

# Prefer local/data log path if pnfs log writes fail; default next to the input dir.
log_file="${dir}/update.log"

echo "dirs: ${dir}"
echo "log:  ${log_file}"

python "${SCRIPT_DIR}/update_dfs_add_columns.py" \
  --dirs "${dir}" \
  --log-file "${log_file}" \
  -v \
  "$@"
