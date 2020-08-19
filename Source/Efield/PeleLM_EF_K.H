#ifndef PeleLM_EF_K_H_
#define PeleLM_EF_K_H_

#include <AMReX_REAL.H>

AMREX_GPU_HOST_DEVICE
inline
void
getKappaE(int i, int j, int k,
          amrex::Array4<      amrex::Real> const& Ke) noexcept
{
   using namespace amrex::literals;

   Ke(i,j,k) = 0.4_rt;
}

AMREX_GPU_HOST_DEVICE
inline
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
#endif
