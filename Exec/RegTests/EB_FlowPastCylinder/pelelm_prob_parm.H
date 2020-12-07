#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
   // Shared params
    int probType = 0;  
    amrex::Real meanFlowMag = 3.0;
    amrex::Real T_mean      = 298.0;
    amrex::Real P_mean      = 101325.0;

    // Wave params
    int wave_type = 0;
    amrex::Real delta_wave  = 2.0e-1;
    amrex::Real xwave       = 0.0;
    amrex::Real RC          = 0.02e-1;
};
#endif
