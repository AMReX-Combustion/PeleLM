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
        pp.query("T_mean", PeleLM::prob_parm->T_mean);
        pp.query("P_mean", PeleLM::prob_parm->P_mean);

        if ( type == "ConvectedVortex" ) {
           PeleLM::prob_parm->probType = 0;
           pp.query("rvort", PeleLM::prob_parm->rvort);
           pp.query("xvort", PeleLM::prob_parm->xvort);
           pp.query("yvort", PeleLM::prob_parm->yvort);
           pp.query("forcevort", PeleLM::prob_parm->forcevort);
           pp.query("meanFlowDir", PeleLM::prob_parm->meanFlowDir);
           pp.query("meanFlowMag", PeleLM::prob_parm->meanFlowMag);
        } else if ( type == "ConvectedGaussian" ) {
           PeleLM::prob_parm->probType = 1;
           pp.query("gaussian_rad", PeleLM::prob_parm->rgauss);
           pp.query("gaussian_x0", PeleLM::prob_parm->xgauss);
           pp.query("gaussian_y0", PeleLM::prob_parm->ygauss);
           pp.query("gaussian_ampl", PeleLM::prob_parm->ampgauss);
           std::string gtype;
           pp.query("gaussian_type", gtype);
           if ( gtype == "Spec" ) {
              PeleLM::prob_parm->gauss_type = 0;
           } else if ( gtype == "Temp" ) { 
              PeleLM::prob_parm->gauss_type = 1;
           } else {
              amrex::Print() << " Unknown prob.gaussian_type ! Should be Spec or Temp \n";
              amrex::Abort();
           }
           pp.query("meanFlowDir", PeleLM::prob_parm->meanFlowDir);
           pp.query("meanFlowMag", PeleLM::prob_parm->meanFlowMag);
        } else if ( type == "DiffusedGaussian" ) {
           PeleLM::prob_parm->probType = 2;
           pp.query("gaussian_rad", PeleLM::prob_parm->rgauss);
           pp.query("gaussian_x0", PeleLM::prob_parm->xgauss);
           pp.query("gaussian_y0", PeleLM::prob_parm->ygauss);
           pp.query("gaussian_ampl", PeleLM::prob_parm->ampgauss);
           std::string gtype;
           pp.query("gaussian_type", gtype);
           if ( gtype == "Spec" ) {
              PeleLM::prob_parm->gauss_type = 0;
           } else if ( gtype == "Temp" ) { 
              PeleLM::prob_parm->gauss_type = 1;
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
