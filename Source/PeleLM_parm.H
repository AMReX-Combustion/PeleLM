#ifndef PELELM_AC_PARM_H_
#define PELELM_AC_PARM_H_
namespace ACParm {
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real  ctrl_dV;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real  ctrl_V_in;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real  ctrl_tBase;
    extern AMREX_GPU_DEVICE_MANAGED  unsigned int ctrl_active;
}
#endif
