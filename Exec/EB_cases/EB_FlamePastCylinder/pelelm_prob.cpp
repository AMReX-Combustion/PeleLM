#include <pelelm_prob.H>
#include <pelelm_prob_parm.H>
#include <pmf.H>

namespace ProbParm
{
   AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real T_mean = 298.0;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real MeanFlow = 3.0;
   AMREX_GPU_DEVICE_MANAGED  int         FlowDir = AMREX_SPACEDIM - 1;
} // namespace ProbParm

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                        const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");
        pp.query("P_mean", ProbParm::P_mean);
        pp.query("T_mean", ProbParm::T_mean);
        pp.query("MeanFlow", ProbParm::MeanFlow);
        pp.query("FlowDir",  ProbParm::FlowDir);

        amrex::Print() << " Flow direction : "<< ProbParm::FlowDir << "\n";    

        std::string pmf_datafile;
        pp.query("pmf_datafile", pmf_datafile);
        int pmf_do_average = 1; 
        PMF::read_pmf(pmf_datafile, pmf_do_average); 
    }
}
