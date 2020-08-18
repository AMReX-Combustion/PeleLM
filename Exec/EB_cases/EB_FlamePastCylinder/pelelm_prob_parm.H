#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real MeanFlow;
    extern AMREX_GPU_DEVICE_MANAGED  int FlowDir;
}

#endif
