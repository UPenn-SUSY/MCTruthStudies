cd LHAPDF
./configure --prefix=${PWD}/..
make
make install
cd ..

LD_LIBRARY_PATH=${PWD}/lib/:$LD_LIBRARY_PATH

