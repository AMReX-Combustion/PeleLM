#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real T_mean = 298.0_rt;
   amrex::Real P_mean = 101325.0_rt;
   amrex::Real meanFlowMag = 0.0;
   int  meanFlowDir = 1;
};
#endif
