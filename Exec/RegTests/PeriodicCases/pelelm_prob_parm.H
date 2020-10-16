#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
   // Shared params
   extern AMREX_GPU_DEVICE_MANAGED  int probType;  
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag;
   extern AMREX_GPU_DEVICE_MANAGED  int  meanFlowDir;

   // CoVo params
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real rvort;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real xvort;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real yvort;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real forcevort;

   // CoGau & DifGau params
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real rgauss;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real xgauss;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real ygauss;
   extern AMREX_GPU_DEVICE_MANAGED  amrex::Real ampgauss;
   extern AMREX_GPU_DEVICE_MANAGED  int gauss_type;
}

#endif
