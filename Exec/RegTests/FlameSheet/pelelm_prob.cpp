#include <PeleLM.H>
#include <pelelm_prob.H>
#include <pmf.H>


namespace ProbParm
{
   AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real standoff = - 0.022;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real pertmag = 0.0004;
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
        pp.query("standoff", ProbParm::standoff);
        pp.query("pertmag", ProbParm::pertmag);


        std::string pmf_datafile;
        pp.query("pmf_datafile", pmf_datafile);
        int pmf_do_average = 1;
        PMF::read_pmf(pmf_datafile, pmf_do_average);
    }
}
