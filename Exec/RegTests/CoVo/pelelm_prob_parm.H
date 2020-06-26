#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
    : amrex::Gpu::Managed
{
    amrex::Real meanFlowMag = 50.0_rt;
    amrex::Real T_mean = 298.0_rt;
    amrex::Real P_mean = 101325.0_rt;
    amrex::Real rvort  = 0.07_rt;
    amrex::Real xvort  = 0.5_rt;
    amrex::Real yvort  = 0.5_rt;
    amrex::Real forcevort  = 6.0_rt;

    int  meanFlowDir = 1;

};

#endif
