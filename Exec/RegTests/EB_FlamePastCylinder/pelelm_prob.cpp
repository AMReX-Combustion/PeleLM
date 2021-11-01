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

        pp.query("P_mean", PeleLM::prob_parm->P_mean);
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("MeanFlow", PeleLM::prob_parm->MeanFlow);
        pp.query("FlowDir",  PeleLM::prob_parm->FlowDir);

        amrex::Print() << " Flow direction : "<< PeleLM::prob_parm->FlowDir << "\n";
        PeleLM::pmf_data.initialize();
    }
}
