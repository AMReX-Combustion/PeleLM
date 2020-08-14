#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Utility.H>

namespace ProbParm
{
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real rgauss;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real xgauss;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real ygauss;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real ampgauss;

    extern AMREX_GPU_DEVICE_MANAGED  std::string gauss_type;

    extern AMREX_GPU_DEVICE_MANAGED  int  meanFlowDir;
}

#endif
