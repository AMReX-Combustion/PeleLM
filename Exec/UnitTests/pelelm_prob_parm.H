#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
    amrex::Real P_mean   = 101325.0_rt;
    amrex::Real standoff = 0.0_rt;
    amrex::Real pertmag  = 0.0004_rt;
};

#endif
