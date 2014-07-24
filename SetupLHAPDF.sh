tar xzvf lhapdf-5.9.1.tar.gz LHAPDF
cd LHAPDF
./configure --prefix=${PWD}/..
make
make install
cd ..

LD_LIBRARY_PATH=${PWD}/lib/:$LD_LIBRARY_PATH

