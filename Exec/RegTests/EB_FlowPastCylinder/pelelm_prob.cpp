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
        pp.query("type", type);
        pp.query("meanFlowMag", PeleLM::prob_parm->meanFlowMag);
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("P_mean", PeleLM::prob_parm->P_mean);

        if ( type == "VortexShedding" ) {
           PeleLM::prob_parm->probType = 0;
        } else if ( type == "Wave" ) {
           PeleLM::prob_parm->probType = 1;
           pp.query("xwave", PeleLM::prob_parm->xwave);
           pp.query("RC",    PeleLM::prob_parm->RC);
           pp.query("delta_wave", PeleLM::prob_parm->delta_wave);
           std::string gtype;
           pp.query("wave_type", gtype);
           if ( gtype == "Spec" ) {
              PeleLM::prob_parm->wave_type = 0;
           } else if ( gtype == "Temp" ) { 
              PeleLM::prob_parm->wave_type = 1;
           } else {
              amrex::Print() << " Unknown prob.wave_type ! Should be Spec or Temp \n";
              amrex::Abort();
           }
        } else { 
            amrex::Print() << " Unknown prob.type ! Should be VortexShedding or Wave \n";
            amrex::Abort();
        }

        auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
        amrex::ParmParse pptr("transport");
        pp.query("const_viscosity", trans_parm.const_viscosity);
        pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
        pp.query("const_conductivity", trans_parm.const_conductivity);
        pp.query("const_diffusivity", trans_parm.const_diffusivity);
        PeleLM::trans_parms.sync_to_device();
    }
}
