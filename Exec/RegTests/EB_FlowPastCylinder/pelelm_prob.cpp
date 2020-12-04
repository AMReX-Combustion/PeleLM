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

        pp.query("meanFlowMag", PeleLM::prob_parm->meanFlowMag);
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("P_mean", PeleLM::prob_parm->P_mean);
        pp.query("xblob", PeleLM::prob_parm->xblob);
        pp.query("RC", PeleLM::prob_parm->RC);
        pp.query("delta_blob", PeleLM::prob_parm->delta_blob);

    }
}
