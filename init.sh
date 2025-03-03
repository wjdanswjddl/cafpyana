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
spack load xrootd@5.6.1

######################################################
#### setup virtual python env if it is not already set
######################################################
# Define the Python version and virtual environment name
#!/bin/bash
cd envs
# Define the Python version and virtual environment name
PYTHON_VERSION=3.9.21
VENV_NAME=venv_py39_cafpyana

# Check if virtual environment already exists
if [ -d "$VENV_NAME" ]; then
    echo "Virtual environment '$VENV_NAME' already exists. Activating it."
else
    # Create the virtual environment
    python -m venv $VENV_NAME
    echo "Virtual environment '$VENV_NAME' created."
fi

# Activate the virtual environment
source $VENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip
pip install wheel setuptools
pip install -r ./envs/pip_requirements.txt

# Deactivate virtual environment
echo "Virtual environment '$VENV_NAME' set up successfully with Python $PYTHON_VERSION and required packages installed."
echo "If you want to exit from this virtual env, do $ deactivate"
echo "If you want to activate this virtual evn again, do $ source '$VENV_NAME'/bin/activate"
cd ..
######################################################

######################################################
# Needed to install xrootd -- which, by the way, is super annoying
###################################################### 
OLDPATH=$PATH
PATH=$PATH:$PWD
cd envs
ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
which cmake3
wget https://files.pythonhosted.org/packages/96/e9/32107ac154c33c6bafd53a5f8444290938c3557210276e5deabb82f74b8f/xrootd-5.6.1.tar.gz
tar -zxvf xrootd-5.6.1.tar.gz
rm xrootd-5.6.1.tar.gz
cd xrootd-5.6.1
python setup.py install
cd ../..
PATH=$OLDPATH

###################################################### 

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/xrootd-5.6.1-py3.9-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH
