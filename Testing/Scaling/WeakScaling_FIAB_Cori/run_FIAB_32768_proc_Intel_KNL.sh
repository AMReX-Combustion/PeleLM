#!/bin/bash

#SBATCH --ntasks=32768
#SBATCH --ntasks-per-node=64
#SBATCH -C knl
#SBATCH -q debug
#SBATCH -A m2860
#SBATCH -J FIAB_32768_proc_Intel_KNL
#SBATCH -t 00:15:00
#SBATCH -L SCRATCH
#SBATCH -o FIAB_32768_proc_Intel_KNL.%j.out
#SBATCH -S 4
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90,TIME_LIMIT

EXE_DIR="${HOME}/PeleLM/Exec/FlameInABox"
EXE="PeleLM3d.intel.mic-knl.MPI.ex"

INPUTS_DIR=${PWD}
INPUTS="inputs.3d-regt_32768proc"

PROBIN_DIR=${PWD}
PROBIN="probin.3d.test"

EXTRA="amr.plot_files_output=0 max_step=1 amrex.signal_handling=0 init_shrink=0.1"

WORKDIR=${SCRATCH}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}
mkdir -p ${WORKDIR}
cp ${EXE_DIR}/${EXE} ${INPUTS_DIR}/${INPUTS} ${PROBIN_DIR}/${PROBIN} ${WORKDIR}
cd ${WORKDIR}

sbcast --compress=lz4 ./${EXE} /tmp/${EXE}

# If your job fails and you see errors like this in the output:
#
# Wed Aug 29 08:09:39 2018: [PE_14400]:inet_recv:inet_recv: recv error (fd=25) Connection reset by peer
# Wed Aug 29 08:09:39 2018: [PE_14400]:_pmi_network_barrier:_pmi_inet_recv from target 15 failed pmi errno -1
# Wed Aug 29 08:09:39 2018: [PE_14400]:_pmi_init:network_barrier failed
#
# then your srun may be timing out. It shouldn't do that, and you should report
# it to NERSC consulting. But if it does and you want to work around it ASAP,
# uncomment the following line to increase the timeout that Slurm allows before
# terminating the srun.
#export PMI_MMAP_SYNC_WAIT_TIME=300

srun -n 32768 -c 4 --cpu-bind=cores /tmp/${EXE} ${INPUTS} ${EXTRA}
