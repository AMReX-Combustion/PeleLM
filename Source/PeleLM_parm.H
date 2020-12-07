#ifndef PELELM_AC_PARM_H_
#define PELELM_AC_PARM_H_
struct ACParm
   : amrex::Gpu::Managed
{
    amrex::Real  ctrl_dV;
    amrex::Real  ctrl_V_in;
    amrex::Real  ctrl_tBase;
    unsigned int ctrl_active;
};
#endif
