#ifndef PeleLM_EF_K_H_
#define PeleLM_EF_K_H_

#include <AMReX_REAL.H>

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
#endif
