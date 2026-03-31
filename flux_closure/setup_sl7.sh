source ~/.profile
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup dk2nudata v01_11_00 -q e26:prof
setup dk2nugenie v01_11_00sbn2 -q e26:prof

source /exp/sbnd/app/users/kplows/NUISANCE/nuisance/build/Linux/setup.sh
echo Using `which nuiscomp`

source /exp/sbnd/app/users/kplows/NUISANCE/python_analysis/bin/activate
