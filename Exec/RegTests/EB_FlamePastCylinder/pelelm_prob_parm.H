#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
    amrex::Real P_mean    = 101325.0;
    amrex::Real T_mean    = 298.0;
    amrex::Real MeanFlow  = 3.0;
    int FlowDir           = 0;
};
#endif
