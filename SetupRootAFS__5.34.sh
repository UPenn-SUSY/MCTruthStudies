source /afs/cern.ch/sw/lcg/external/gcc/4.7.2/x86_64-slc5/setup.sh

export PYTHONDIR=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc5-gcc47-opt/root/

export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# export PATH=$PATH:$ROOTSYS/bin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib

export PYTHONPATH=$ROOTSYS/pyroot:$ROOTSYS/lib:$PYTHONPATH
