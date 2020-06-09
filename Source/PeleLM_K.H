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
                 amrex::Array4<      amrex::Real> const& rhoYdot) noexcept
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
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;      // Include MKS -> CGS conversion
   EOS::HY2T(H, y, T);

   // Get wdot. Don't like the temporary wdot there ...
   rho = rho * 0.001_rt;                                 // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES];
   EOS::RTY2WDOT(rho, T, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {   
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;            // CGS -> MKS conversion
   }
}

AMREX_GPU_HOST_DEVICE
inline
void
getPGivenRTY(int i, int j, int k,
             amrex::Array4<const amrex::Real> const& rho,
             amrex::Array4<const amrex::Real> const& rhoY,
             amrex::Array4<const amrex::Real> const& T,
             amrex::Array4<      amrex::Real> const& P) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real rho_cgs = rho(i,j,k) * 0.001_rt;

   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   EOS::TRY2P(rho_cgs, T(i,j,k), y, P(i,j,k));
   P(i,j,k) = P(i,j,k) * 0.1_rt;                         // CGS -> MKS conversion
}

// This kernel returns directly rhoHmix
AMREX_GPU_HOST_DEVICE
inline
void
getRHmixGivenTY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& Hmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   EOS::TY2H(T(i,j,k), y, Hmix(i,j,k));
   Hmix(i,j,k) = Hmix(i,j,k) * 1.0e-4_rt * rho(i,j,k);   // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
inline
void
getMwmixGivenRY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& Mwmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   EOS::Y2WBAR(y, Mwmix(i,j,k));
   Mwmix(i,j,k) = Mwmix(i,j,k) * 0.001_rt;            // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
inline
void
getHGivenT(int i, int j, int k,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Hi) noexcept
{
   using namespace amrex::literals;

   amrex::Real hi_spec[NUM_SPECIES];
   EOS::T2Hi(T(i,j,k), hi_spec);
   for (int n = 0; n < NUM_SPECIES; n++) {   
      Hi(i,j,k,n) = hi_spec[n] * 0.0001_rt;            // CGS -> MKS conversion
   }
   
}

AMREX_GPU_HOST_DEVICE
inline
void
buildAdvectionForcing(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& rho,
                      amrex::Array4<const amrex::Real> const& rhoY,
                      amrex::Array4<const amrex::Real> const& T,
                      amrex::Array4<const amrex::Real> const& dn,
                      amrex::Array4<const amrex::Real> const& ddn,
                      amrex::Array4<const amrex::Real> const& r,
                      amrex::Real                      const& dp0dt,
                      int                              const& closed_chamber,
                      amrex::Array4<      amrex::Real> const& forceY,
                      amrex::Array4<      amrex::Real> const& forceT) noexcept
{
   using namespace amrex::literals;

   // Get species enthalpy
   amrex::Real hi_spec[NUM_SPECIES];
   EOS::T2Hi(T(i,j,k), hi_spec);

   // For species mass fractions
   amrex::Real rhoinv = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES];

   forceT(i,j,k) = dn(i,j,k,NUM_SPECIES+1) + ddn(i,j,k);

   for (int n = 0; n < NUM_SPECIES; n++) {   
      y[n] = rhoY(i,j,k,n) * rhoinv;                                 // get y
      forceY(i,j,k,n) = r(i,j,k,n) + dn(i,j,k,n);                    // forceY = Rk + Dnk
      forceT(i,j,k) -= forceY(i,j,k,n) * hi_spec[n] * 0.0001_rt;     // forceT -= Sum{ hk. (Rk + Dnk) } with CGS -> MKS conversion of hk
   }

   for (int n = 0; n < NUM_SPECIES; n++) {   
   }

   if ( closed_chamber == 1 ) forceT(i,j,k) += dp0dt;                // Add dp0/dt

   amrex::Real cpmix, cpmixinv;
   EOS::TY2Cp(T(i,j,k), y, cpmix);
   cpmixinv = 1.0_rt / cpmix * 1.0e4_rt;                             // CGS -> MKS conversion
   forceT(i,j,k) *= rhoinv * cpmixinv;

}

AMREX_GPU_HOST_DEVICE
inline
void
getCpmixGivenRYT(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rho,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& T,
                 amrex::Array4<      amrex::Real> const& cpmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   EOS::TY2Cp(T(i,j,k), y, cpmix(i,j,k));
   cpmix(i,j,k) *= 0.0001_rt;                         // CGS -> MKS conversion
}
#endif
