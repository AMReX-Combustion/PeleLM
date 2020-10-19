#ifndef PeleLM_EF_K_H_
#define PeleLM_EF_K_H_

#include <AMReX_REAL.H>
#include <PeleLM_EF_Constant.H>
#include <mechanism.h>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getKappaE(int i, int j, int k,
          amrex::Array4<      amrex::Real> const& Ke) noexcept
{
   using namespace amrex::literals;

   Ke(i,j,k) = 0.4_rt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getKappaSp(int i, int j, int k,
           amrex::Real const mwt[NUM_SPECIES],
           amrex::GpuArray<amrex::Real, NUM_SPECIES> a_zk,
           amrex::Array4<      amrex::Real> const& rhoY,
           amrex::Array4<      amrex::Real> const& rhoD,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Ksp) noexcept
{
   using namespace amrex::literals;
   amrex::Real rho = 0.0;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   for (int n = 0; n < NUM_SPECIES; n++) {
      Ksp(i,j,k,n) = rhoD(i,j,k,n) * mwt[n] * 0.001 * a_zk[n] / ( rho * PP_RU_MKS * T(i,j,k) );
      if ( n == E_ID ) Ksp(i,j,k,n) = 0.0;         // Set the E species mobility to zero. Electron are treated separately.
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getDiffE(int i, int j, int k,
         const amrex::Real &factor,
         amrex::Array4<const amrex::Real> const& T,
         amrex::Array4<const amrex::Real> const& Ke,
         amrex::Array4<      amrex::Real> const& De) noexcept
{
   using namespace amrex::literals;

   De(i,j,k) = Ke(i,j,k) * T(i,j,k) * factor;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ef_edge_state(int i,int j,int k,
              int dir,
              amrex::Array4<amrex::Real const> const& state,
              amrex::Array4<amrex::Real const> const& effvel ) noexcept
{
   using namespace amrex::literals;
   int idx_hi[3] = {i,j,k};
   int idx_lo[3] = {i,j,k};
   idx_lo[dir] -= 1;
   amrex::Real state_ed;
   if ( effvel(i,j,k) > small_vel) {
      state_ed = state(idx_lo[0],idx_lo[1],idx_lo[2]);
   } else if ( effvel(i,j,k) < -small_vel) {
      state_ed = state(idx_hi[0],idx_hi[1],idx_hi[2]);
   } else {
      state_ed = 0.5_rt * ( state(idx_lo[0],idx_lo[1],idx_lo[2]) +
                            state(idx_hi[0],idx_hi[1],idx_hi[2]) );
   }
   return state_ed;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
cen2edg_upwind_cpp(int i, int j, int k,
                   int dir,
                   int ncomp,
                   const bool                       on_lo,
                   const bool                       on_hi,
                   amrex::Array4<const amrex::Real> const& ueff,
                   amrex::Array4<const amrex::Real> const& cfab,
                   amrex::Array4<      amrex::Real> const& efab) noexcept
{
   using namespace amrex::literals;

   // Default behavior: id_l -1 in dir
   // Handle the BCs
   // Need -1 in id_l and id_h in dir on low Dirichlet BC
   // Need nothing on high Dirichlet BC
   int  id_l[3] = {i,j,k};
   int  id_h[3] = {i,j,k};
   if (!on_lo ) {
      if ( !on_hi ) {
         id_l[dir] -= 1;
      }
   } else {
      id_l[dir] -= 1;
      id_h[dir] -= 1;
   }

   for (int n = 0; n < ncomp; n++) {
      if ( ueff(i,j,k) > small_vel ) {
         efab(i,j,k,n) = cfab(id_l[0],id_l[1],id_l[2],n);
      } else if ( ueff(i,j,k) < -small_vel ) {
         efab(i,j,k,n) = cfab(id_h[0],id_h[1],id_h[2],n);
      } else {
         efab(i,j,k,n) = 0.5_rt*(cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n));
      }
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
addLorentzForce(int i, int j, int k,
                int scomp,
                int ncomp,
                amrex::GpuArray<int,3> const bxlo,
                amrex::GpuArray<int,3> const bxhi,
                const amrex::Real &time,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const dx,
                amrex::GpuArray<amrex::Real, NUM_SPECIES> const a_zk,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& ne,
                amrex::Array4<const amrex::Real> const& phiV,
                amrex::Array4<      amrex::Real> const& force) noexcept
{
   using namespace amrex::literals;
   // Electric field. Use 1-sided stencil in the force ghost cells.
   amrex::Real EFx;
   if ( i == bxlo[0] ) {
      EFx = - 1.0_rt / dx[0] * ( phiV(i+1,j,k) - phiV(i,j,k) );
   } else if ( i == bxhi[0] ) {
      EFx = - 1.0_rt / dx[0] * ( phiV(i,j,k) - phiV(i-1,j,k) );
   } else {
      EFx = - 0.5_rt / dx[0] * ( phiV(i+1,j,k) - phiV(i-1,j,k) );
   }
   amrex::Real EFy;
   if ( j == bxlo[1] ) {
      EFy = - 1.0_rt / dx[1] * ( phiV(i,j+1,k) - phiV(i,j,k) );
   } else if ( j == bxhi[1] ) {
      EFy = - 1.0_rt / dx[1] * ( phiV(i,j,k) - phiV(i,j-1,k) );
   } else {
      EFy = - 0.5_rt / dx[1] * ( phiV(i,j+1,k) - phiV(i,j-1,k) );
   }
#if ( AMREX_SPACEDIM == 3 )
   amrex::Real EFz;
   if ( k == bxlo[2] ) {
      EFz = - 1.0_rt / dx[2] * ( phiV(i,j,k+1) - phiV(i,j,k) );
   } else if ( k == bxhi[2] ) {
      EFz = - 1.0_rt / dx[2] * ( phiV(i,j,k) - phiV(i,j,k-1) );
   } else {
      EFz = - 0.5_rt / dx[2] * ( phiV(i,j,k+1) - phiV(i,j,k-1) );
   }
#endif
   for (int n = 0; n < NUM_SPECIES; n++) {
      force(i,j,k,0) += a_zk[n] * rhoY(i,j,k,n) * EFx;
      force(i,j,k,1) += a_zk[n] * rhoY(i,j,k,n) * EFy;
#if ( AMREX_SPACEDIM == 3 )
      force(i,j,k,2) += a_zk[n] * rhoY(i,j,k,n) * EFz;
#endif
   }
   force(i,j,k,0) -= ne(i,j,k) * EFConst::elemCharge * EFx;
   force(i,j,k,1) -= ne(i,j,k) * EFConst::elemCharge * EFy;
#if ( AMREX_SPACEDIM == 3 )
   force(i,j,k,2) -= ne(i,j,k) * EFConst::elemCharge * EFz;
#endif
}

AMREX_GPU_HOST_DEVICE
inline
void
reactionRateRhoY_EF(int i, int j, int k,
                    amrex::Array4<const amrex::Real> const& rhoY,
                    amrex::Array4<const amrex::Real> const& rhoH,
                    amrex::Array4<const amrex::Real> const& T,
                    amrex::Array4<const amrex::Real> const& nE,
                    amrex::Array4<int>         const& mask,
                    amrex::Array4<      amrex::Real> const& rhoYdot,
                    amrex::Array4<      amrex::Real> const& nEdot) noexcept
{
   using namespace amrex::literals;

   // Get rho from rhoY. Is there a better way ?
   amrex::Real rho = 0.0_rt;
   amrex::Real rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   // Get Y[e] from nE
   amrex::Real mw[NUM_SPECIES];
   EOS::molecular_weight(mw);
   y[E_ID] = nE(i,j,k) / EFConst::Na * mw[E_ID] * 0.001 * rhoinv;

   // Get T from Y/H. Was done like that before.
   amrex::Real Tloc = T(i,j,k);
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;      // Include MKS -> CGS conversion
   EOS::HY2T(H, y, Tloc);

   // Get wdot. Don't like the temporary wdot there ...
   rho = rho * 0.001_rt;                                 // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES];
   EOS::RTY2WDOT(rho, Tloc, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;            // CGS -> MKS conversion
   }
   nEdot(i,j,k) = wdot[E_ID] * EFConst::Na * 1.0e6;
   rhoYdot(i,j,k,E_ID) = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata_ef (int i, int j, int k,
                    amrex::Array4<amrex::Real> const& state,
                    amrex::GeometryData const& geomdata)
{
    // Got the mass fractions from the pltfile
    // Get rid of E species and get the charge
    state(i,j,k,DEF_first_spec+N2_ID) += state(i,j,k,DEF_first_spec+E_ID);
    state(i,j,k,DEF_first_spec+E_ID) = 0.0;

    amrex::Real zk[NUM_SPECIES];
    EOS::charge_mass(zk);
    amrex::Real chargeDist = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
       chargeDist += state(i,j,k,DEF_first_spec+n) * zk[n];
    }

    // Clean up places where charge dist < 0
    if ( chargeDist < 0.0 ) {
       amrex::Real sum_pos_ion = 0.0;
       amrex::Real mean_zk_pos_ion = 0.0;
       for (int n = 0; n < NUM_SPECIES; n++) {
          if ( zk[n] > 0.0 ) {
             sum_pos_ion += state(i,j,k,DEF_first_spec+n);
             mean_zk_pos_ion += state(i,j,k,DEF_first_spec+n) * zk[n];
          }
       }
       for (int n = 0; n < NUM_SPECIES; n++) {
          if ( zk[n] > 0.0 ) {
             state(i,j,k,DEF_first_spec+N2_ID) += state(i,j,k,DEF_first_spec+n) 
                                                / sum_pos_ion * chargeDist / mean_zk_pos_ion * sum_pos_ion;
             state(i,j,k,DEF_first_spec+n) -= state(i,j,k,DEF_first_spec+n) 
                                             / sum_pos_ion * chargeDist / mean_zk_pos_ion * sum_pos_ion;
          }
       }
    }

    // Get the electron number density such that electro-neutral.
    chargeDist = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
        chargeDist += state(i,j,k,DEF_first_spec+n) * zk[n];
    }
    state(i,j,k,DEF_nE) = std::max(1.0e-24,chargeDist/EFConst::elemCharge);
    state(i,j,k,DEF_PhiV) = 0.0;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = state(i,j,k,DEF_first_spec+n);
    }
    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

    state(i,j,k,DEF_nE) *= state(i,j,k,Density);
}

#endif
