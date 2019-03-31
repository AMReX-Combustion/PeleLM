#!/bin/bash

#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=32
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -A m2860
#SBATCH -J FIAB_00016_proc_Intel_HSW
#SBATCH -t 00:10:00
#SBATCH -L SCRATCH
#SBATCH -o FIAB_00016_proc_Intel_HSW.%j.out
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90,TIME_LIMIT

EXE_DIR="${PWD}"
EXE="PeleLM3d.intel.haswell.MPI.ex"

INPUTS_DIR=${PWD}
# A relatively fair comparison of HSW vs KNL is "constant number of nodes", or
# "constant fraction of a node" if a run uses less than the full node. Although
# the KNL nodes have 68 cores, we assume they have only 64 (and we only use 64
# per node). The inputs file names correspond to the number of MPI processes on
# a KNL node. So, e.g., "inputs.3d-regt_8proc" corresponds to "1/4 of a single
# KNL node".  So for the corresponding Haswell run, we use the same inputs
# file, but use 1/4 of a Haswell node, which in the case of
# "inputs.3d-regt_8proc" corresponds to 4 MPI procs for the Haswell node, not 8
# procs. So the Haswell run scripts will have a factor of 2 fewer MPI procs
# than the inputs file names would suggest.
INPUTS="inputs.3d-regt_32proc"

PROBIN_DIR=${PWD}
PROBIN="probin.3d.test"

EXTRA="amr.plot_files_output=0 max_step=1 amrex.signal_handling=0 init_shrink=0.1"

WORKDIR=${SCRATCH}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}
mkdir -p ${WORKDIR}
cp ${EXE_DIR}/${EXE} ${INPUTS_DIR}/${INPUTS} ${PROBIN_DIR}/${PROBIN} ${WORKDIR}
cd ${WORKDIR}

sbcast --compress=lz4 ./${EXE} /tmp/${EXE}
srun -n 16 -c 2 --cpu-bind=cores /tmp/${EXE} ./${INPUTS} ${EXTRA}
