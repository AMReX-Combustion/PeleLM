#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   // Shared params
   int probType = 0;  
   amrex::Real T_mean = 298.0;
   amrex::Real P_mean = 101325.0;
   amrex::Real meanFlowMag = 0.0;
   int  meanFlowDir = 1;

   // CoVo params
   amrex::Real rvort = 0.07;
   amrex::Real xvort = 0.5;
   amrex::Real yvort = 0.5;
   amrex::Real forcevort = 6.0;

   // CoGau & DifGau params
   amrex::Real rgauss = 0.1;
   amrex::Real xgauss = 0.5;
   amrex::Real ygauss = 0.5;
   amrex::Real ampgauss = 0.1;
   int gauss_type = 0;
};
#endif
