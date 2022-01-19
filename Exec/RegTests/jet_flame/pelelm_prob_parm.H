#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <turbinflow.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real P_mean   = 101325.0_rt;
    amrex::Real standoff = 0.0_rt;
    amrex::Real pertmag  = 0.0004_rt;
    amrex::Real meanFlowMag = 0.0;
    int  meanFlowDir = 1;

    TurbParm tp;
    amrex::GpuArray<bool,2*AMREX_SPACEDIM> turb_ok = {{false}};

    bool do_turb = false;
 
};
#endif
