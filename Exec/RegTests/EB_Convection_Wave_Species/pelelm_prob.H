#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "pelelm_prob_parm.H"

//#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>

AMREX_GPU_DEVICE
inline
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    amrex::Real dist = sqrt(pow(x-ProbParm::xblob,2.) );

    state(i,j,k,DEF_Temp) = ProbParm::T_mean;

    amrex::Real Yt[NUM_SPECIES];

    Yt[0] = 0.233*(1.0 + ProbParm::delta_blob*exp(-dist/(2.0 * ProbParm::RC)));
    Yt[1] = 1.0 - Yt[0];


    state(i,j,k,Xvel) = ProbParm::meanFlowMag;
    state(i,j,k,Yvel) = 0.0;
#if (AMREX_SPACEDIM == 3)
    state(i,j,k,Zvel) = 0.0;
#endif

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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();

  if (sgn == 1) {

s_ext[Xvel] = ProbParm::meanFlowMag;
s_ext[Yvel] = 0.0;
#if (AMREX_SPACEDIM == 3)
s_ext[Zvel] = 0.0;
#endif

s_ext[DEF_Temp] = ProbParm::T_mean;

    amrex::Real Yt[NUM_SPECIES];
    Yt[0] = 0.233;
    Yt[1] = 0.767;


    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, Yt, ProbParm::T_mean, rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(ProbParm::T_mean, Yt, RhoH_temp);
    s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = Yt[n] * s_ext[Density];
    }



  } else {

    amrex::Abort("Inflow not implemented for hi_x");

  }


}






#endif
