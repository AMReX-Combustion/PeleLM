#!/bin/bash -l

module load llvm
module load cmake
module load python
module load py-matplotlib
module load py-six
module load py-numpy
module load py-cycler
module load py-bottleneck
module load py-python-dateutil
module load py-cython
module load py-nose
module load py-numexpr
module load py-packaging
module load py-pandas
module load py-pillow
module load py-pytz
module load py-setuptools
module load py-kiwisolver
module load py-pyparsing
module load cppcheck

ln -s ${CPPCHECK_ROOT_DIR}/cfg/std.cfg

#export CXXFLAGS="-fsanitize=address -fno-omit-frame-pointer"

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=clang++ \
      -DCMAKE_C_COMPILER:STRING=clang \
      -DCMAKE_Fortran_COMPILER:STRING=gfortran \
      \ #-DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DPELELM_DIM:STRING=3 \
      -DPELELM_ENABLE_AMREX_EB:BOOL=ON \
      -DPELELM_ENABLE_MPI:BOOL=OFF \
      -DPELELM_ENABLE_TESTS:BOOL=ON \
      -DPELELM_ENABLE_FCOMPARE:BOOL=OFF \
      -DPELELM_ENABLE_FCOMPARE_FOR_TESTS:BOOL=OFF \
      -DPELELM_ENABLE_SUNDIALS:BOOL=ON \
      -DPELELM_ENABLE_MASA:BOOL=OFF \
      -DMASA_DIR:STRING=$(spack location -i masa) \
      -DPELELM_ENABLE_ALL_WARNINGS:BOOL=OFF \
      -DPELELM_ENABLE_CPPCHECK:BOOL=OFF \
      -DPELELM_ENABLE_CLANG_TIDY:BOOL=OFF \
      -DPELELM_ENABLE_CUDA:BOOL=OFF \
      -DAMReX_CUDA_ARCH=Volta \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DPELELM_PRECISION:STRING=DOUBLE \
      -DPELELM_ENABLE_FPE_TRAP_FOR_TESTS:BOOL=OFF \
      .. 
cmake --build . --parallel $(nproc)
ctest -j $(nproc)
