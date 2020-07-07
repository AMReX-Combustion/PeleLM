#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "pelelm_prob_parm.H"

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include "EOS.H"

AMREX_GPU_DEVICE
inline
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& pressure,
                 amrex::GeometryData const& geomdata)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    state(i,j,k,DEF_Temp) = ProbParm::T_mean;

    amrex::Real Yt[NUM_SPECIES];
    Yt[0] = 0.233;
    Yt[1] = 0.767;


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

    pressure(i,j,k,0) = ProbParm::P_mean;

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

/*

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real vel_ext[AMREX_SPACEDIM],
  amrex::Real rho_ext,
  amrex::Real Y_ext[NUM_SPECIES],
  amrex::Real Temp_ext,
  amrex::Real enth_ext,
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();




  if (sgn == -1) {

vel_ext[0] = ProbParm::meanFlowMag;
vel_ext[1] = 0.0;




  } else {

    amrex::Abort("Inflow not implemented for hi_x");

  }




}


*/



#endif
