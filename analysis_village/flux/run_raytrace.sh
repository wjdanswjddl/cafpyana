fluxdir="/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/fluxFiles/bnb/BooNEtoGSimple/configK-v1/july2023/neutrinoMode/"
python gsimple_raytrace_batch.py \
  --gsimple-dir $fluxdir \
  --batch-size 500 \
  --out-dir /exp/sbnd/data/users/munjung/flux/SBND_gsimple_raytrace \
  --save-npz /exp/sbnd/data/users/munjung/flux/SBND_gsimple_raytrace/raytrace_flux_histograms.npz
