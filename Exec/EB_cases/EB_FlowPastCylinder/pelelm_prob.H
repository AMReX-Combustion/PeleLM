#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
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

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real L_y = prob_hi[1] - prob_lo[1];

    amrex::Real massfrac[2] = {0.233,0.767};

    state(i,j,k,DEF_Temp) = ProbParm::T_mean;

    state(i,j,k,Xvel) = 0.0;
    state(i,j,k,Yvel) = 0.0;
#if ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Zvel) = 0.0;
#endif

    switch(ProbParm::FlowDir) {
      case 0:
         state(i,j,k,Xvel) = ProbParm::MeanFlow;
         break;
      case 1:
         state(i,j,k,Yvel) = ProbParm::MeanFlow;
         break;
#if ( AMREX_SPACEDIM == 3 ) 
      case 2:
         state(i,j,k,Zvel) = ProbParm::MeanFlow;
         break;
#endif
    }

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
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
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();

  if (sgn == 1) {
    s_ext[Xvel] = 0.0;
    s_ext[Yvel] = 0.0;
#if ( AMREX_SPACEDIM == 3 ) 
    s_ext[Zvel] = 0.0;
#endif
    switch(ProbParm::FlowDir) {
      case 0:
         s_ext[Xvel] = ProbParm::MeanFlow;
         break;
      case 1:
         s_ext[Yvel] = ProbParm::MeanFlow;
         break;
#if ( AMREX_SPACEDIM == 3 ) 
      case 2:
         s_ext[Zvel] = ProbParm::MeanFlow;
         break;
#endif
    }

    s_ext[DEF_Temp] = ProbParm::T_mean;

    amrex::Real massfrac[2] = {0.233,0.767};

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    EOS::TY2H(s_ext[DEF_Temp], massfrac, s_ext[DEF_RhoH]);
    s_ext[DEF_RhoH] *= 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
    }
  }

}

AMREX_GPU_DEVICE
inline
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
