#!/usr/bin/env bash
# Batch gsimple flux + raytrace + z-slab config sweep.
# Uses a merged unique-z accumulation (one project_to_z per plane per file).
set -euo pipefail

fluxdir="/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/fluxFiles/bnb/BooNEtoGSimple/configK-v1/july2023/neutrinoMode/"
outdir="/exp/sbnd/data/users/munjung/flux/SBND_gsimple"

# Optional: limit for a quick test
# extra_args=(--n-files 50)
extra_args=()

python gsimple_batch.py \
  --gsimple-dir "$fluxdir" \
  --batch-size 500 \
  --out-dir "$outdir" \
  --save-npz "$outdir/flux_histograms.npz" \
  --no-plots \
  "${extra_args[@]}"
# Default --slab-z-ns: 5 11 21 26 51 101 201 + batch z0_500_n26
# Writes slab_zconfig_flux_histograms.npz and raytrace_flux_histograms.npz
# Avoid very large n (e.g. 441): ~800+ z planes/file is very slow even with union-z.
