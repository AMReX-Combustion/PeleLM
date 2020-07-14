#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real standoff;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real pertmag;

extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N;
extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M;
extern AMREX_GPU_DEVICE_MANAGED bool pmf_do_average;

extern amrex::Gpu::ManagedVector<amrex::Real>* pmf_X;
extern amrex::Gpu::ManagedVector<amrex::Real>* pmf_Y;

extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_X;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_Y;

extern std::string pmf_datafile;
extern amrex::Vector<std::string> pmf_names;


}

#endif
