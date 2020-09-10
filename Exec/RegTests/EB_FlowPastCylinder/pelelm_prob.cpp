
#include <AMReX_PROB_AMR_F.H>
#include "pelelm_prob_parm.H"
#include <AMReX_ParmParse.H>

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag = 3.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean = 298.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real xblob  = 0.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real delta_blob  = 2e-1;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real RC  = 0.02e-1;
}

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("meanFlowMag", ProbParm::meanFlowMag);
        pp.query("T_mean", ProbParm::T_mean);
        pp.query("P_mean", ProbParm::P_mean);
        pp.query("xblob", ProbParm::xblob);
        pp.query("RC", ProbParm::RC);
        pp.query("delta_blob", ProbParm::delta_blob);

    }
}
