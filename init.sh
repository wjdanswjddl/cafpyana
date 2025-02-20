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

######################################################
#### setup virtual python env if it is not already set
######################################################
# Define the Python version and virtual environment name
#!/bin/bash

# Define the Python version and virtual environment name
PYTHON_VERSION=3.9.21
VENV_NAME=venv_py39_cafpyana

# Check if virtual environment already exists
if [ -d "$VENV_NAME" ]; then
    echo "Virtual environment '$VENV_NAME' already exists. Activating it."
else
    # Create the virtual environment
    python3 -m venv $VENV_NAME
    echo "Virtual environment '$VENV_NAME' created."
fi

# Activate the virtual environment
source $VENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip
pip install wheel setuptools
pip install -r ./data/requirements.txt

# Deactivate virtual environment
echo "Virtual environment '$VENV_NAME' set up successfully with Python $PYTHON_VERSION and required packages installed."
echo "If you want to exit from this virtual env, do $ deactivate"
echo "If you want to activate this virtual evn again, do $ source '$VENV_NAME'/bin/activate"

######################################################
export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd/lib64
