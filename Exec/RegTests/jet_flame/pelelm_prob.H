#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    state(i,j,k,DEF_Temp) = prob_parm.T_ox;

    for (int n = 0; n < NUM_SPECIES; n++)
    {
        massfrac[n] = prob_parm.Y_ox[n];
    }

    state(i,j,k,Xvel) = 0;
    state(i,j,k,Yvel) = 0;
    state(i,j,k,Zvel) = 0.;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_func(amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start)/(duration*tau))*
    std::tanh((-time + start + duration)/(duration*tau))+1.0);
  eta_time = std::min(std::max(eta_time,0.0),1.0);
  return eta_time;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real *x, /*x[AMREX_SPACEDIM]*/
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data,
  amrex::Real& eta_out)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real radius = {0.0};
  amrex::Real enthalpy;
  int ierr;

  if (sgn == 1 && idir == 2) {

    radius = std::sqrt(std::pow(x[0]-prob_parm.center_xy[0],2) + std::pow(x[1]-prob_parm.center_xy[1],2));

    s_ext[Xvel] = 0.0;
    s_ext[Yvel] = 0.0;

    amrex::Real eta = {0};
    eta = 0.5*(1.0 - std::tanh(2.0*(radius-prob_parm.D*0.5)/(0.1*prob_parm.D*0.5)));
    eta *= time_func(time,prob_parm.inj_start,prob_parm.inj_dur,prob_parm.tau);
    eta_out = eta;

    
    for (int n = 0; n < NUM_SPECIES; n++)
    {
      massfrac[n] = eta * prob_parm.Y_fuel[n] + (1.0-eta) * prob_parm.Y_ox[n];
    }
    enthalpy = eta * prob_parm.H_fuel + (1.0-eta)*prob_parm.H_ox;

    GET_T_GIVEN_HY(&enthalpy,massfrac,&s_ext[DEF_Temp],&ierr);

    s_ext[Zvel] = eta*prob_parm.v_in;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
    s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif
