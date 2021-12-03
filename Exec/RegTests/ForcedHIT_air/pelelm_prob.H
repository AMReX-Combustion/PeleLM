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
#include <Utilities.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    amrex::Real x[3] = {
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real y1, y2;

    // Fill in the velocities
    amrex::Real uinterp[3] = {0.0};
    // amrex::Print() << "turbulence_from_file? "<< prob_parm.turbulence_from_file << "\n";
    // if(prob_parm.turbulence_from_file){
      // Interpolation factors
      amrex::Real mod[3] = {0.0};
      int idx[3] = {0};
      int idxp1[3] = {0};
      amrex::Real slp[3] = {0.0};
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        mod[idim] = std::fmod(x[idim], prob_parm.Linput);
        locate(prob_parm.d_xarray, prob_parm.input_resolution, mod[idim], idx[idim]);
        idxp1[idim] = (idx[idim] + 1) % prob_parm.input_resolution;
        slp[idim] = (mod[idim] - prob_parm.d_xarray[idx[idim]]) / prob_parm.d_xdiff[idx[idim]];
      }

      int inSize = prob_parm.input_resolution;

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
      const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
      const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
      const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
      const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
      const amrex::Real f7 = slp[0] * slp[1] * slp[2];

      uinterp[0] =   prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0
                   + prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                   + prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                   + prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                   + prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                   + prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                   + prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                   + prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

      uinterp[1] =   prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 
                   + prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                   + prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                   + prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                   + prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                   + prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                   + prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                   + prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

      uinterp[2] =   prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0
                   + prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                   + prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                   + prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                   + prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                   + prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                   + prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                   + prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

      //
      // Fill Velocity
      //
      // state(i,j,k,Xvel) = uinterp[0];
      // state(i,j,k,Yvel) = uinterp[1];
      // state(i,j,k,Zvel) = uinterp[2];
      // }
      // else{
      //   // amrex::Print() << "Initializing velocity field with sine function \n";
      //   const amrex::Real Lx    = prob_hi[0] - prob_lo[0];
      //   const amrex::Real Ly    = prob_hi[1] - prob_lo[1];
      //   const amrex::Real Lz    = prob_hi[2] - prob_lo[1];

      //   uinterp[0] =  prob_parm.turb_scale * std::cos(prob_parm.TwoPi*x[1]/Ly) * std::cos(prob_parm.TwoPi*x[2]/Lz);
      //   uinterp[1] =  prob_parm.turb_scale * std::cos(prob_parm.TwoPi*x[0]/Lx) * std::cos(prob_parm.TwoPi*x[2]/Lz);
      //   uinterp[2] =  prob_parm.turb_scale * std::cos(prob_parm.TwoPi*x[0]/Lx) * std::cos(prob_parm.TwoPi*x[1]/Ly);
      // }


    state(i,j,k,DEF_Temp) = prob_parm.T_mean;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;
    

    state(i,j,k,Xvel) = 0.1*uinterp[0];
#if ( AMREX_SPACEDIM == 2 ) 
    state(i,j,k,Yvel) = 0.1*uinterp[1]; //[m/s]
#elif ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Yvel) = 0;
    state(i,j,k,Zvel) = 0.1*uinterp[2] + prob_parm.V_in; //[m/s]
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

    // amrex::Real rho_cgs, P_cgs;
    // P_cgs = prob_parm.P_mean * 10.0;

    // auto eos = pele::physics::PhysicsType::eos();
    // eos.PYT2R(P_cgs, Yt, state(i,j,k,DEF_Temp), rho_cgs);
    // state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    // eos.TY2H(state(i,j,k,DEF_Temp), Yt, state(i,j,k,DEF_RhoH));
    // state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    // for (int n = 0; n < NUM_SPECIES; n++) {
    //   state(i,j,k,DEF_first_spec+n) = Yt[n] * state(i,j,k,Density);
    // }
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
  ACParm const& ac_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real velx, vely, velz;
  amrex::Real rho, T, pres;

  amrex::Real Lx = prob_hi[0]-prob_lo[0];
  amrex::Real l_fluct = Lx/10.0; //corresponds to lt for the turbulence level we are aiming at 
  amrex::Real pert = prob_parm.V_in * 0.1; //velocity perturbation at the inlet. For now it is 10 % the mean value

  if (sgn == 1) {
    // pele::physics::PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);
    if (!ac_parm.ctrl_active || time <= 0.0) {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = 0.0; //[m/s]
#elif (AMREX_SPACEDIM == 3)
       // s_ext[Yvel] = prob_parm.V_in;
       s_ext[Yvel]   = sin(x[2]/l_fluct*3.14156)*pert + prob_parm.V_in;
       velz = pmf_vals[2]; //[m/s]
#endif
    } 
    else {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = prob_parm.V_in;
       velz = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#endif
    }

    s_ext[DEF_Temp] = prob_parm.T_mean;

    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;
    // auto eos = pele::physics::PhysicsType::eos();
    // eos.X2Y(molefrac, massfrac);

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
