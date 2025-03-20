#!/bin/bash 

export machine=${HOSTNAME}
if [[ $machine == *sbnd* || $machine == *jupyter* ]]; then
  echo "working on a sbnd machine"
  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
  export CAFPYANA_GRID_OUT_DIR="/pnfs/sbnd/scratch/users/$USER/cafpyana_out"
  mkdir -p $CAFPYANA_GRID_OUT_DIR
fi
if [[ $machine == *icarus* ]]; then
  echo "working on a icarus machine"
  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
  export CAFPYANA_GRID_OUT_DIR="/pnfs/icarus/scratch/users/$USER/cafpyana_out"
  mkdir -p $CAFPYANA_GRID_OUT_DIR
fi
spack load hdf5@1.14.3
spack load xrootd@5.6.1

VENV_NAME=venv_py39_cafpyana
source ./envs/$VENV_NAME/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/xrootd-5.6.1-py3.9-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH
export CAFPYANA_WD=`pwd`
