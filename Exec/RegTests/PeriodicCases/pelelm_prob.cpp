#include <pelelm_prob.H>

namespace ProbParm
{
    // Shared parameters
    AMREX_GPU_DEVICE_MANAGED  int probType = 0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean = 298.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real meanFlowMag = 0.0;
    AMREX_GPU_DEVICE_MANAGED  int  meanFlowDir = 1;

    // CoVo params
    AMREX_GPU_DEVICE_MANAGED  amrex::Real rvort  = 0.07;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real xvort  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real yvort  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real forcevort  = 6.0;

    // CoGau & DifGau params
    AMREX_GPU_DEVICE_MANAGED  amrex::Real rgauss  = 0.1;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real xgauss  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real ygauss  = 0.5;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real ampgauss  = 0.1;
    AMREX_GPU_DEVICE_MANAGED  int gauss_type = 0;
}

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        std::string type;
        pp.query("type", type);
        pp.query("T_mean", ProbParm::T_mean);
        pp.query("P_mean", ProbParm::P_mean);

        if ( type == "ConvectedVortex" ) {
           ProbParm::probType = 0;
           pp.query("rvort", ProbParm::rvort);
           pp.query("xvort", ProbParm::xvort);
           pp.query("yvort", ProbParm::yvort);
           pp.query("forcevort", ProbParm::forcevort);
           pp.query("meanFlowDir", ProbParm::meanFlowDir);
           pp.query("meanFlowMag", ProbParm::meanFlowMag);
        } else if ( type == "ConvectedGaussian" ) {
           ProbParm::probType = 1;
           pp.query("gaussian_rad", ProbParm::rgauss);
           pp.query("gaussian_x0", ProbParm::xgauss);
           pp.query("gaussian_y0", ProbParm::ygauss);
           pp.query("gaussian_ampl", ProbParm::ampgauss);
           std::string gtype;
           pp.query("gaussian_type", gtype);
           if ( gtype == "Spec" ) {
              ProbParm::gauss_type = 0;
           } else if ( gtype == "Temp" ) { 
              ProbParm::gauss_type = 1;
           } else {
              amrex::Print() << " Unknown prob.gaussian_type ! Should be Spec or Temp \n";
              amrex::Abort();
           }
           pp.query("meanFlowDir", ProbParm::meanFlowDir);
           pp.query("meanFlowMag", ProbParm::meanFlowMag);
        } else if ( type == "DiffusedGaussian" ) {
           ProbParm::probType = 2;
           pp.query("gaussian_rad", ProbParm::rgauss);
           pp.query("gaussian_x0", ProbParm::xgauss);
           pp.query("gaussian_y0", ProbParm::ygauss);
           pp.query("gaussian_ampl", ProbParm::ampgauss);
           std::string gtype;
           pp.query("gaussian_type", gtype);
           if ( gtype == "Spec" ) {
              ProbParm::gauss_type = 0;
           } else if ( gtype == "Temp" ) { 
              ProbParm::gauss_type = 1;
           } else {
              amrex::Print() << " Unknown prob.gaussian_type ! Should be Spec or Temp \n";
              amrex::Abort();
           }
        } else { 
            amrex::Print() << " Unknown prob.type ! Should be ConvectedVortex, ConvectedGaussian or DiffusedGaussian \n";
            amrex::Abort();
        }
    }
}
