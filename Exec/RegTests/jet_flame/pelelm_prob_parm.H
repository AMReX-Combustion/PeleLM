#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <turbinflow.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real P_mean   = 101325.0_rt;

    TurbParm tp;
    amrex::GpuArray<bool,2*AMREX_SPACEDIM> turb_ok = {{false}};

    bool do_turb = false;

};
#endif
