#include <PeleLM.H>
#include <pelelm_prob.H>

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        std::string type;
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("P_mean", PeleLM::prob_parm->P_mean);
    }
}
