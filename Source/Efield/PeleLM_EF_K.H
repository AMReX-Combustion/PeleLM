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
ef_edge_state(int i,int j,int k, int n,
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

#endif
