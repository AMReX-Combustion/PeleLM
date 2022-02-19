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

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
    auto eos = pele::physics::PhysicsType::eos();

    const amrex::Real L_x = prob_hi[0] - prob_lo[0];

    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert;
    amrex::Real y1, y2;

    pert = 0.0;
    if (prob_parm.pertmag > 0.0) 
    {

#if ( AMREX_SPACEDIM == 2 ) 
       pert = prob_parm.pertmag * 
             (1.0 * std::sin(2 * Pi * 4 * x / L_x) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / L_x) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / L_x) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / L_x) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / L_x));


    }

    y1 = (y - prob_parm.standoff - 0.5*dx[1] + pert)*100;
    y2 = (y - prob_parm.standoff + 0.5*dx[1] + pert)*100;

#elif ( AMREX_SPACEDIM == 3 ) 
       const amrex::Real L_y = prob_hi[1] - prob_lo[1];
       pert = prob_parm.pertmag *
             (1.0 * std::sin(2 * Pi * 4 * x / L_x) *
                std::sin(2 * Pi * 5 * y / L_y) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / L_x) *
                std::sin(2 * Pi * 4 * (y - 0.0053765) / L_y) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / L_x) *
                std::sin(2 * Pi * 3 * (y - 0.02137) / L_y) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / L_x) *
                std::sin(2 * Pi * 6 * (y - 0.018) / L_y) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / L_x));
    }

    y1 = (z - prob_parm.standoff - 0.5*dx[2] + pert)*100;
    y2 = (z - prob_parm.standoff + 0.5*dx[2] + pert)*100;
#endif

    pele::physics::PMF::pmf(pmf_data,y1, y2, pmf_vals);

    state(i,j,k,DEF_Temp) = pmf_vals[0];;

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);


    state(i,j,k,Xvel) = 0;
#if ( AMREX_SPACEDIM == 2 ) 
    state(i,j,k,Yvel) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Yvel) = 0;
    state(i,j,k,Zvel) = pmf_vals[1]*1e-2;
#endif

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
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();

  if (sgn == 1) {
    pele::physics::PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

    if (!ac_parm.ctrl_active || time <= 0.0) {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = pmf_vals[1]*1e-2;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = pmf_vals[1]*1e-2;
#endif
    } else {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#endif
    }

    s_ext[DEF_Temp] = pmf_vals[0];

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

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
