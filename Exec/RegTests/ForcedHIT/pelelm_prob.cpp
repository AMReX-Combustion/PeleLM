#include <PeleLM.H>
#include <pelelm_prob.H>
#include <pmf.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // amrex::ParmParse pp("prob");

        // pp.query("T_in",   PeleLM::prob_parm->T_in);
        // pp.query("P_mean", PeleLM::prob_parm->P_mean);
        // pp.query("phi_in", PeleLM::prob_parm->phi_in);
        // pp.query("vn_in", PeleLM::prob_parm->vn_in);

    }
}
