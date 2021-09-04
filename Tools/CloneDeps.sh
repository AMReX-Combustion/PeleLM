#!/usr/bin/env bash

echo "Getting PeleLM dependencies - tests ... "
export PELELM_HOME=${PWD}/..
mkdir build
git clone https://github.com/AMReX-Codes/amrex.git build/amrex
export AMREX_HOME=${PWD}/build/amrex
git clone -b development https://github.com/AMReX-Codes/IAMR.git build/IAMR
export IAMR_HOME=${PWD}/build/IAMR
git clone -b main https://github.com/AMReX-Codes/AMReX-Hydro.git build/AMReX-Hydro
export AMREX_HYDRO_HOME=${PWD}/build/AMReX-Hydro
git clone -b development https://github.com/AMReX-Combustion/PelePhysics.git build/PelePhysics
export PELE_PHYSICS_HOME=${PWD}/build/PelePhysics
