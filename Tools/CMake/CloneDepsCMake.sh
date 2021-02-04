#!/usr/bin/env bash

echo "Getting PeleLM external repos for CMake..."
DEPS_HOME=$1

set -x
set -e

git clone https://github.com/AMReX-Codes/amrex.git ${DEPS_HOME}/amrex
git clone https://github.com/AMReX-Codes/IAMR.git ${DEPS_HOME}/IAMR
git clone https://github.com/AMReX-Combustion/PelePhysics.git ${DEPS_HOME}/PelePhysics
git clone https://github.com/LLNL/sundials.git ${DEPS_HOME}/sundials
git clone https://github.com/google/googletest.git ${DEPS_HOME}/googletest

#IAMR needs a specific commit
#(cd ${DEPS_HOME}/IAMR && git checkout 11a13e5159c5085d3a221437b9d79eb6e3467fc3)
