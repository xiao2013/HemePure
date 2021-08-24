#!/bin/bash
## Compilation/build script for HEMELB
## Run from found location

## MODULE loads
##GCC compilers
MODULES(){

#Module environment on ARCHER2
#module restore PrgEnv-gnu 

#Module environment on SuperMUC-NG, default compilers are fine
module load cmake

module list
#Export compiler shortcuts as named on given machine
export CC=mpiicc && export CXX=mpiicpc

}

## HEMELB build
# 1) Dependencies
DEPbuild(){
cd dep
rm -rf build
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} ..
make -j  && echo "Done HemeLB Dependencies"

cd ../..
}

SRCbuild_ARCHER2(){
cd src
rm -rf build
mkdir build
cd build

cmake -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCTEMPLATE_LIBRARY=<PathToRepoOnARCHER2>/dep/install/lib/libctemplate.a -DHEMELB_USE_MPI_WIN=OFF ..


make -j && echo "Done HemeLB Source"

cd ../..
}


SRCbuild_SNG(){
cd src
rm -rf build_DEFAULT
mkdir build_DEFAULT
cd build_DEFAULT

cmake -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DHEMELB_USE_GMYPLUS=OFF -DHEMELB_USE_MPI_WIN=ON -DHEMELB_USE_SSE3=ON -DHEMELB_USE_AVX2=OFF ..

make -j && echo "Done HemeLB Source"

cd ../..
}

MODULES
DEPbuild
SRCbuild_SNG
