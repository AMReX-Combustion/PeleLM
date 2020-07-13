
#include <AMReX_PROB_AMR_F.H>
#include "pelelm_prob_parm.H"
#include <AMReX_ParmParse.H>

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag = 50.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean = 298.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real rvort  = 0.07;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real xvort  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real yvort  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real forcevort  = 6.0;

    AMREX_GPU_DEVICE_MANAGED  int  meanFlowDir = 1;
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
        pp.query("rvort", ProbParm::rvort);
        pp.query("xvort", ProbParm::xvort);
        pp.query("yvort", ProbParm::yvort);
        pp.query("forcevort", ProbParm::forcevort);
        pp.query("meanFlowDir", ProbParm::meanFlowDir);

    }
}
