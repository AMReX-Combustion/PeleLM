#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include "pelelm_prob_parm.H"
#include <pmf_data.H>
#include <IndexDefines.H>
#include <EOS.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 PmfData const * /*pmf_data*/)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    state(i,j,k,DEF_Temp) = ProbParm::T_mean;

    amrex::Real Yt[NUM_SPECIES];
    Yt[O2_ID] = 0.233;
    Yt[N2_ID] = 0.767;

    if ( ProbParm::probType == 0 ) { // CoVo
       const amrex::Real deltax = x - ProbParm::xvort;
       const amrex::Real deltay = y - ProbParm::yvort;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = ProbParm::rvort * ProbParm::rvort;
       const amrex::Real u_vort = -ProbParm::forcevort*deltay/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real v_vort = ProbParm::forcevort*deltax/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real w_vort = 0.;

       switch(ProbParm::meanFlowDir) {
         case 1  :
           state(i,j,k,Xvel) = ProbParm::meanFlowMag + u_vort;
           state(i,j,k,Yvel) = v_vort;
           break;
         case -1  :
           state(i,j,k,Xvel) = -ProbParm::meanFlowMag + u_vort;
           state(i,j,k,Yvel) = v_vort;
           break;
         case 2  :
           state(i,j,k,Xvel) = u_vort;
           state(i,j,k,Yvel) = ProbParm::meanFlowMag + v_vort;
           break;
         case -2  :
           state(i,j,k,Xvel) = u_vort;
           state(i,j,k,Yvel) = -ProbParm::meanFlowMag + v_vort;
           break;
         case 3  :
           state(i,j,k,Xvel) = ProbParm::meanFlowMag + u_vort;
           state(i,j,k,Yvel) = ProbParm::meanFlowMag + v_vort;
           break; 
         case -3  :
           state(i,j,k,Xvel) = -ProbParm::meanFlowMag + u_vort;
           state(i,j,k,Yvel) = -ProbParm::meanFlowMag + v_vort;
           break;
       }
#if (AMREX_SPACEDIM == 3)
       state(i,j,k,Zvel) = w_vort;
#endif
    } else if ( ProbParm::probType == 1 ) { // CoGau
       const amrex::Real deltax = x - ProbParm::xgauss;
       const amrex::Real deltay = y - ProbParm::ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = ProbParm::rgauss * ProbParm::rgauss;
       if ( ProbParm::gauss_type == 0 ) { // Spec
          Yt[O2_ID] += Yt[O2_ID] * ProbParm::ampgauss * std::exp(-d_sq/r_sq);
          Yt[N2_ID] = 1.0 - Yt[O2_ID];
       } else if ( ProbParm::gauss_type == 1 ) { // Temp
          state(i,j,k,DEF_Temp) = ProbParm::T_mean * ( 1.0 + ProbParm::ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
       switch(ProbParm::meanFlowDir) {
         case 1  :
           state(i,j,k,Xvel) = ProbParm::meanFlowMag;
           state(i,j,k,Yvel) = 0.0;
           break;
         case -1  :
           state(i,j,k,Xvel) = -ProbParm::meanFlowMag;
           state(i,j,k,Yvel) = 0.0;
           break;
         case 2  :
           state(i,j,k,Xvel) = 0.0;
           state(i,j,k,Yvel) = ProbParm::meanFlowMag;
           break;
         case -2  :
           state(i,j,k,Xvel) = 0.0;
           state(i,j,k,Yvel) = -ProbParm::meanFlowMag;
           break;
         case 3  :
           state(i,j,k,Xvel) = ProbParm::meanFlowMag;
           state(i,j,k,Yvel) = ProbParm::meanFlowMag;
           break; 
         case -3  :
           state(i,j,k,Xvel) = -ProbParm::meanFlowMag;
           state(i,j,k,Yvel) = -ProbParm::meanFlowMag;
           break;
       }
#if (AMREX_SPACEDIM == 3)
       state(i,j,k,Zvel) = 0.0;
#endif
    } else if ( ProbParm::probType == 2 ) { // DifGau
       const amrex::Real deltax = x - ProbParm::xgauss;
       const amrex::Real deltay = y - ProbParm::ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = ProbParm::rgauss * ProbParm::rgauss;
       if ( ProbParm::gauss_type == 0 ) { // Spec
          Yt[O2_ID] += Yt[O2_ID] * ProbParm::ampgauss * std::exp(-d_sq/r_sq);
          Yt[N2_ID] = 1.0 - Yt[O2_ID];
       } else if ( ProbParm::gauss_type == 1 ) { //Temp
          state(i,j,k,DEF_Temp) = ProbParm::T_mean * ( 1.0 + ProbParm::ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
       state(i,j,k,Xvel) = 0.0;
       state(i,j,k,Yvel) = 0.0;
#if (AMREX_SPACEDIM == 3)
       state(i,j,k,Zvel) = 0.0;
#endif
    }

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, Yt, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), Yt, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = Yt[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  PmfData const * /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

}
#endif
