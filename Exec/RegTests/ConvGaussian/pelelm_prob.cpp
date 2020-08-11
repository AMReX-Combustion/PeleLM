
#include <AMReX_PROB_AMR_F.H>
#include "pelelm_prob_parm.H"
#include <AMReX_ParmParse.H>

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag = 1.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean = 298.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real rgauss  = 0.1;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real xgauss  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real ygauss  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real ampgauss  = 0.1;

    AMREX_GPU_DEVICE_MANAGED std::string gauss_type = "Spec";

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
        pp.query("gaussian_rad", ProbParm::rgauss);
        pp.query("gaussian_x0", ProbParm::xgauss);
        pp.query("gaussian_y0", ProbParm::ygauss);
        pp.query("gaussian_ampl", ProbParm::ampgauss);
        pp.query("gaussian_type", ProbParm::gauss_type);
        pp.query("meanFlowDir", ProbParm::meanFlowDir);

    }
}
