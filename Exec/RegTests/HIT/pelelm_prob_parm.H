#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean = 101325.0;
   amrex::Real T_mean = 300.0;

   int input_resolution = 0;
   amrex::Real urms0 = 1.0;
   amrex::Real uin_norm = 1.0;

   amrex::Real Linput = 0.0;

   amrex::Real* d_uinput = nullptr;
   amrex::Real* d_vinput = nullptr;
   amrex::Real* d_winput = nullptr;
   amrex::Real* d_xarray = nullptr;
   amrex::Real* d_xdiff = nullptr;
};
#endif
