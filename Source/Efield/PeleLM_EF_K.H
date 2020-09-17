#ifndef PeleLM_EF_K_H_
#define PeleLM_EF_K_H_

#include <AMReX_REAL.H>
#include <PeleLM_EF_Constant.H>

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
   amrex::Real EFx = - 0.5_rt / dx[0] * ( phiV(i+1,j,k) - phiV(i-1,j,k) );
   if ( i == bxlo[0] ) EFx = - 1.0_rt / dx[0] * ( phiV(i+1,j,k) - phiV(i,j,k) );
   if ( i == bxhi[0] ) EFx = - 1.0_rt / dx[0] * ( phiV(i,j,k) - phiV(i-1,j,k) );
   amrex::Real EFy = - 0.5_rt / dx[1] * ( phiV(i,j+1,k) - phiV(i,j-1,k) );
   if ( j == bxlo[1] ) EFy = - 1.0_rt / dx[1] * ( phiV(i,j+1,k) - phiV(i,j,k) );
   if ( j == bxhi[1] ) EFy = - 1.0_rt / dx[1] * ( phiV(i,j,k) - phiV(i,j-1,k) );
#if ( AMREX_SPACEDIM == 3 )
   amrex::Real EFz = - 0.5_rt / dx[2] * ( phiV(i,j,k+1) - phiV(i,j,k-1) );
   if ( k == bxlo[2] ) EFz = - 1.0_rt / dx[2] * ( phiV(i,j,k+1) - phiV(i,j,k) );
   if ( k == bxhi[2] ) EFz = - 1.0_rt / dx[2] * ( phiV(i,j,k) - phiV(i,j,k-1) );
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
#endif
