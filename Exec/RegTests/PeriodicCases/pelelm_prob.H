#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include "pelelm_prob_parm.H"
#include <PMFData.H>
#include <PeleLM_parm.H>
#include <IndexDefines.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    state(i,j,k,DEF_Temp) = prob_parm.T_mean;

    amrex::Real Yt[NUM_SPECIES];
    Yt[O2_ID] = 0.233;
    Yt[N2_ID] = 0.767;

    if ( prob_parm.probType == 0 ) { // CoVo
       const amrex::Real deltax = x - prob_parm.xvort;
       const amrex::Real deltay = y - prob_parm.yvort;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = prob_parm.rvort * prob_parm.rvort;
       const amrex::Real u_vort = -prob_parm.forcevort*deltay/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real v_vort = prob_parm.forcevort*deltax/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real w_vort = 0.;

       switch(prob_parm.meanFlowDir) {
         case 1  :
           state(i,j,k,Xvel) = prob_parm.meanFlowMag + u_vort;
           state(i,j,k,Yvel) = v_vort;
           break;
         case -1  :
           state(i,j,k,Xvel) = -prob_parm.meanFlowMag + u_vort;
           state(i,j,k,Yvel) = v_vort;
           break;
         case 2  :
           state(i,j,k,Xvel) = u_vort;
           state(i,j,k,Yvel) = prob_parm.meanFlowMag + v_vort;
           break;
         case -2  :
           state(i,j,k,Xvel) = u_vort;
           state(i,j,k,Yvel) = -prob_parm.meanFlowMag + v_vort;
           break;
         case 3  :
           state(i,j,k,Xvel) = prob_parm.meanFlowMag + u_vort;
           state(i,j,k,Yvel) = prob_parm.meanFlowMag + v_vort;
           break; 
         case -3  :
           state(i,j,k,Xvel) = -prob_parm.meanFlowMag + u_vort;
           state(i,j,k,Yvel) = -prob_parm.meanFlowMag + v_vort;
           break;
       }
#if (AMREX_SPACEDIM == 3)
       state(i,j,k,Zvel) = w_vort;
#endif
    } else if ( prob_parm.probType == 1 ) { // CoGau
       const amrex::Real deltax = x - prob_parm.xgauss;
       const amrex::Real deltay = y - prob_parm.ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
       if ( prob_parm.gauss_type == 0 ) { // Spec
          Yt[O2_ID] += Yt[O2_ID] * prob_parm.ampgauss * std::exp(-d_sq/r_sq);
          Yt[N2_ID] = 1.0 - Yt[O2_ID];
       } else if ( prob_parm.gauss_type == 1 ) { // Temp
          state(i,j,k,DEF_Temp) = prob_parm.T_mean * ( 1.0 + prob_parm.ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
       switch(prob_parm.meanFlowDir) {
         case 1  :
           state(i,j,k,Xvel) = prob_parm.meanFlowMag;
           state(i,j,k,Yvel) = 0.0;
           break;
         case -1  :
           state(i,j,k,Xvel) = -prob_parm.meanFlowMag;
           state(i,j,k,Yvel) = 0.0;
           break;
         case 2  :
           state(i,j,k,Xvel) = 0.0;
           state(i,j,k,Yvel) = prob_parm.meanFlowMag;
           break;
         case -2  :
           state(i,j,k,Xvel) = 0.0;
           state(i,j,k,Yvel) = -prob_parm.meanFlowMag;
           break;
         case 3  :
           state(i,j,k,Xvel) = prob_parm.meanFlowMag;
           state(i,j,k,Yvel) = prob_parm.meanFlowMag;
           break; 
         case -3  :
           state(i,j,k,Xvel) = -prob_parm.meanFlowMag;
           state(i,j,k,Yvel) = -prob_parm.meanFlowMag;
           break;
       }
#if (AMREX_SPACEDIM == 3)
       state(i,j,k,Zvel) = 0.0;
#endif
    } else if ( prob_parm.probType == 2 ) { // DifGau
       const amrex::Real deltax = x - prob_parm.xgauss;
       const amrex::Real deltay = y - prob_parm.ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
       if ( prob_parm.gauss_type == 0 ) { // Spec
          Yt[O2_ID] += Yt[O2_ID] * prob_parm.ampgauss * std::exp(-d_sq/r_sq);
          Yt[N2_ID] = 1.0 - Yt[O2_ID];
       } else if ( prob_parm.gauss_type == 1 ) { //Temp
          state(i,j,k,DEF_Temp) = prob_parm.T_mean * ( 1.0 + prob_parm.ampgauss * std::exp(-d_sq/r_sq) );
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
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, Yt, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,DEF_Temp), Yt, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = Yt[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& /*prob_parm*/,
  ACParm const& /*ac_parm*/,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
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
