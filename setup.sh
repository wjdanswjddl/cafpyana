#!/bin/bash 

export machine=${HOSTNAME}
if [[ $machine == *sbnd* ]]; then
  echo "working on a sbnd machine"
  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
fi
if [[ $machine == *icarus* ]]; then
  echo "working on a icarus machine"
  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
fi
spack load hdf5@1.14.3

VENV_NAME=venv_py39_cafpyana
source $VENV_NAME/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd/lib64
