#ifndef PeleLM_K_H_
#define PeleLM_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <mechanism.h>
#include <EOS.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
reactionRateRhoY(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& rhoH,
                 amrex::Array4<int>         const& mask,
                 amrex::Array4<amrex::Real> const& rhoYdot) noexcept
{
   using namespace amrex::literals;

   // Get rho from rhoY. Is there a better way ?
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   // Get T from Y/H. Was done like that before.
   amrex::Real T;
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;  // Include MKS -> CGS convertion
   EOS::HY2T(H, y, T);

   // Get wdot. Don't like the temporary wdot there ...
   rho = rho * 0.001_rt;                             // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES];
   EOS::RTY2WDOT(rho, T, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {   
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;        // MKS -> CGS convertion
   }

//   printf("rho: %6.4f, T: %8.6f, wdot: %10.8f \n",rho,T,rhoYdot(i,j,k,0));  
}
#endif
