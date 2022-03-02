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
void pelelm_initdata(int i, int j, int k,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5*Lx;,
                 const amrex::Real y_c = prob_lo[1] + 0.5*Ly;,
                 const amrex::Real z_c = prob_lo[2] + 0.5*Lz);

    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;


    AMREX_D_TERM(state(i,j,k,Xvel) = 0.0;,
                 state(i,j,k,Yvel) = 0.0;,
                 state(i,j,k,Zvel) = 0.0);
    state(i,j,k,Xvel+prob_parm.meanFlowDir) =  prob_parm.meanFlowMag;

    state(i,j,k,DEF_Temp) = prob_parm.T_mean;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,DEF_Temp), massfrac, h_cgs);
    state(i,j,k,DEF_RhoH) = h_cgs * 1.0e-4 * state(i,j,k,Density);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
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
  ProbParm const& prob_parm,
  ACParm const& /*ac_parm*/,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5*Lx;,
                 const amrex::Real y_c = prob_lo[1] + 0.5*Ly;,
                 const amrex::Real z_c = prob_lo[2] + 0.5*Lz);

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    auto eos = pele::physics::PhysicsType::eos();

    // When using turbulnce injection, s_ext comes in with
    // turbulence data. Only add the mean flow here.
    s_ext[prob_parm.meanFlowDir] +=  prob_parm.meanFlowMag;

    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    amrex::Real rad = 0.0;
    if ( prob_parm.meanFlowDir == 0 ) {
        rad = std::sqrt( (x[1]-y_c)*(x[1]-y_c) +
                         (x[2]-z_c)*(x[2]-z_c) );
    } else if ( prob_parm.meanFlowDir == 1 ) {
        rad = std::sqrt( (x[0]-x_c)*(x[0]-x_c) +
                         (x[2]-z_c)*(x[2]-z_c) );
    } else if ( prob_parm.meanFlowDir == 2 ) {
        rad = std::sqrt( (x[0]-x_c)*(x[0]-x_c) +
                         (x[1]-y_c)*(x[1]-y_c) );
    }

    if (rad < 0.15 * Lx) {
        massfrac[O2_ID] = 0.5;
        massfrac[N2_ID] = 0.5;
    }

    s_ext[DEF_Temp] = prob_parm.T_mean;

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
