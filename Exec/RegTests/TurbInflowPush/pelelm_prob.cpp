#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include "turbinflow.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
         amrex::ParmParse pp("prob");
         
         pp.query("T_mean", PeleLM::prob_parm->T_mean);
         pp.query("P_mean", PeleLM::prob_parm->P_mean);

         PeleLM::prob_parm->do_turb = false;
#ifdef PELE_USE_TURBINFLOW
         if (pp.countval("turb_file") > 0) {
#if AMREX_SPACEDIM==2
            amrex::Abort("Turbulence inflow unsupported in 2D.");
#endif
            std::string turb_file = "";
            pp.query("turb_file", turb_file);
            amrex::Real turb_scale_loc = 1.0;
            pp.query("turb_scale_loc", turb_scale_loc);
            amrex::Real turb_scale_vel = 1.0;
            pp.query("turb_scale_vel", turb_scale_vel);

            PeleLM::prob_parm->do_turb = true;

            // Hold nose here - required because of dynamically allocated data in tp
            AMREX_ASSERT_WITH_MESSAGE(PeleLM::prob_parm->tp.tph == nullptr,"Can only be one TurbParmHost");
            PeleLM::prob_parm->tp.tph = new TurbParmHost;

            amrex::Vector<amrex::Real> turb_center = {
              {0.5 * (probhi[0] + problo[0]), 0.5 * (probhi[1] + problo[1])}};
            pp.queryarr("turb_center", turb_center);
            AMREX_ASSERT_WITH_MESSAGE(turb_center.size() == 2, "turb_center must have two elements");
            for (int n = 0; n < turb_center.size(); ++n) {
              turb_center[n] *= turb_scale_loc;
            }

            int turb_nplane = 32;
            pp.query("turb_nplane", turb_nplane);
            AMREX_ASSERT(turb_nplane > 0);
            amrex::Real turb_conv_vel = 1;
            pp.query("turb_conv_vel", turb_conv_vel);
            AMREX_ASSERT(turb_conv_vel > 0);

            init_turbinflow(turb_file, turb_scale_loc, turb_scale_vel, turb_center, turb_conv_vel,
                            turb_nplane, PeleLM::prob_parm->tp);
        }
#endif
        auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
        amrex::ParmParse pptr("transport");
        pp.query("const_viscosity", trans_parm.const_viscosity);
        pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
        pp.query("const_conductivity", trans_parm.const_conductivity);
        pp.query("const_diffusivity", trans_parm.const_diffusivity);
        PeleLM::trans_parms.sync_to_device();

    } 
}
