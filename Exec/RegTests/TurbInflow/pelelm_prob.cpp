#include <PeleLM.H>
#include <pelelm_prob.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");
       
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("P_mean", PeleLM::prob_parm->P_mean);
        pp.query("flowDir", PeleLM::prob_parm->meanFlowDir);
        pp.query("flowMag", PeleLM::prob_parm->meanFlowMag);
    
       /*
       if (!m_incompressible) {
          auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
          amrex::ParmParse pptr("transport");
          pp.query("const_viscosity", trans_parm.const_viscosity);
          pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
          pp.query("const_conductivity", trans_parm.const_conductivity);
          pp.query("const_diffusivity", trans_parm.const_diffusivity);
          PeleLM::trans_parms.sync_to_device();
       }
       */
    }
}
