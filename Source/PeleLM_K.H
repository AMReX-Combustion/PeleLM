#ifndef PeleLM_K_H_
#define PeleLM_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <mechanism.h>
#include <EOS.H>
#include <Transport.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
reactionRateRhoY(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& rhoH,
                 amrex::Array4<const amrex::Real> const& T,
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

   EOS::RTY2P(rho_cgs, T(i,j,k), y, P(i,j,k));
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

AMREX_GPU_HOST_DEVICE
inline
void
getTGivenHY(int i, int j, int k,
            amrex::Array4<const amrex::Real> const& rho,
            amrex::Array4<const amrex::Real> const& rhoY,
            amrex::Array4<const amrex::Real> const& rhoH,
            amrex::Array4<      amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;   // MKS -> CGS conversion

   EOS::HY2T(H, y, T(i,j,k));
}

AMREX_GPU_HOST_DEVICE
inline
void
compute_divu(int i, int j, int k,
             amrex::Array4<const amrex::Real> const& rhoY,
             amrex::Array4<const amrex::Real> const& T,
             amrex::Array4<const amrex::Real> const& viscTermT,
             amrex::Array4<const amrex::Real> const& viscTermY,
             amrex::Array4<const amrex::Real> const& rhoYdot,
             amrex::Array4<      amrex::Real> const& divu) noexcept
{
   using namespace amrex::literals;

   amrex::Real mwtinv[NUM_SPECIES];
   EOS::inv_molecular_weight(mwtinv);

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real cpmix, Wbar, hi[NUM_SPECIES];;
   EOS::TY2Cp(T(i,j,k), y, cpmix);
   EOS::Y2WBAR(y, Wbar);
   EOS::T2Hi(T(i,j,k), hi);
   cpmix *= 1.0e-4_rt;                                // CGS -> MKS conversion
   for (int n = 0; n < NUM_SPECIES; n++) {
      hi[n] *= 1.0e-4_rt;                             // CGS -> MKS conversion
   }

   amrex::Real denominv = 1.0_rt / ( rho * cpmix * T(i,j,k) );
   divu(i,j,k) = ( divu(i,j,k) + viscTermT(i,j,k) ) * denominv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      divu(i,j,k) +=  ( viscTermY(i,j,k,n) + rhoYdot(i,j,k,n) )
                    * ( mwtinv[n] * Wbar * rhoinv - hi[n] * denominv );
   }
}

AMREX_GPU_HOST_DEVICE
inline
void
fabMinMax(int i, int j, int k,
          int           ncomp,  
          amrex::Real   fmin,
          amrex::Real   fmax,
          amrex::Array4<amrex::Real> const& fab) noexcept
{
   using namespace amrex::literals;

   for (int n = 0; n < ncomp; n++) {
      fab(i,j,k,n) = amrex::max( fmin, amrex::min( fmax, fab(i,j,k,n) ) );
   }
}

AMREX_GPU_HOST_DEVICE
inline
amrex::Real
getGammaInv(int i, int j, int k,
         amrex::Array4<const amrex::Real> const& rhoY,
         amrex::Array4<const amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real cpmix, cvmix;
   EOS::TY2Cp(T(i,j,k), y, cpmix);
   EOS::TY2Cv(T(i,j,k), y, cvmix);
   
   amrex::Real gammainv = cvmix / cpmix;

   return gammainv;
}

AMREX_GPU_HOST_DEVICE
inline
void
getRhoGivenTYP(int i, int j, int k,
               amrex::Array4<const amrex::Real> const& T,
               amrex::Array4<const amrex::Real> const& Y,
               amrex::Real                      const& Patm,
               amrex::Array4<      amrex::Real> const& rho) noexcept
{
   using namespace amrex::literals;

   amrex::Real rho_cgs, P_cgs, y[NUM_SPECIES];
   P_cgs = Patm * EOS::PATM;
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = Y(i,j,k,n);                       // TODO: need to find a way not to do that
   }
   EOS::PYT2R(P_cgs, y, T(i,j,k), rho_cgs);
   rho(i,j,k) = rho_cgs * 1.0e3_rt;            // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
inline
void
getTransportCoeff(int i, int j, int k,
                  amrex::Array4<const amrex::Real> const& rhoY,
                  amrex::Array4<const amrex::Real> const& T,
                  amrex::Array4<      amrex::Real> const& rhoDi,
                  amrex::Array4<      amrex::Real> const& lambda,
                  amrex::Array4<      amrex::Real> const& mu) noexcept
{
   using namespace amrex::literals;

   amrex::Real mwtinv[NUM_SPECIES];
   EOS::inv_molecular_weight(mwtinv);

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real Wbar;
   EOS::Y2WBAR(y, Wbar);

   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real rhoDi_cgs[NUM_SPECIES], lambda_cgs, mu_cgs, dummy_xi;
   amrex::Real Tloc = T(i,j,k);

   bool get_xi = false, get_mu = true, get_lam = true, get_Ddiag = true;
   transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
             rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs);

   // Do CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = rhoDi_cgs[n] * Wbar * mwtinv[n] * 1.0e-1_rt;
   }
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_HOST_DEVICE
inline
void
getTransportCoeffUnityLe(int i, int j, int k,
                         amrex::Real                      const& ScInv,
                         amrex::Real                      const& PrInv,
                         amrex::Array4<const amrex::Real> const& rhoY,
                         amrex::Array4<      amrex::Real> const& T,
                         amrex::Array4<      amrex::Real> const& rhoDi,
                         amrex::Array4<      amrex::Real> const& lambda,
                         amrex::Array4<      amrex::Real> const& mu) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   rho *= 1.0e-3_rt;                                     // MKS -> CGS conversion
   amrex::Real rhoDi_cgs[NUM_SPECIES], lambda_cgs, mu_cgs, dummy_xi;

   bool get_xi = false, get_mu = true, get_lam = false, get_Ddiag = false;
   transport(get_xi, get_mu, get_lam, get_Ddiag, T(i,j,k),
             rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs);

   mu(i,j,k) = mu_cgs * 1.0e-1_rt;                       // CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = mu_cgs * 1.0e-1_rt * ScInv;       // Constant Schmidt number
   }

   amrex::Real cpmix;
   EOS::TY2Cp(T(i,j,k), y, cpmix);
   lambda(i,j,k) = mu_cgs * 1.0e-1_rt * PrInv / cpmix;   // Constant Prandtl number
}

#ifdef USE_WBAR
AMREX_GPU_HOST_DEVICE
inline
void
getBetaWbar(int i, int j, int k,
            amrex::Array4<const amrex::Real> const& rhoY,
            amrex::Array4<      amrex::Real> const& rhoDi,
            amrex::Array4<      amrex::Real> const& DWbar) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real Wbar;
   EOS::Y2WBAR(y, Wbar);

   for (int n = 0; n < NUM_SPECIES; n++) {
      DWbar(i,j,k,n) = rhoDi(i,j,k,n) * y[n] / Wbar;
   }
}
#endif

AMREX_GPU_HOST_DEVICE
inline
void
getVelViscosity(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& mu) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real dummy_rhoDi[NUM_SPECIES], dummy_lambda, mu_cgs, dummy_xi;

   bool get_xi = false, get_mu = true, get_lam = false, get_Ddiag = false;
   transport(get_xi, get_mu, get_lam, get_Ddiag, T(i,j,k),
             rho, y, dummy_rhoDi, mu_cgs, dummy_xi, dummy_lambda);

   // Do CGS -> MKS conversions
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_HOST_DEVICE
inline
void
getConductivity(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& lambda) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt, rhoinv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES];
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real dummy_rhoDi[NUM_SPECIES], dummy_mu, lambda_cgs, dummy_xi;

   bool get_xi = false, get_mu = false, get_lam = true, get_Ddiag = false;
   transport(get_xi, get_mu, get_lam, get_Ddiag, T(i,j,k),
             rho, y, dummy_rhoDi, dummy_mu, dummy_xi, lambda_cgs);

   // Do CGS -> MKS conversions
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
}

AMREX_GPU_HOST_DEVICE
inline
void
fluxDivergence(int i, int j, int k,
               int ncomp,
               D_DECL(amrex::Array4<const amrex::Real>     const& fluxX,
                      amrex::Array4<const amrex::Real>     const& fluxY,
                      amrex::Array4<const amrex::Real>     const& fluxZ),
               amrex::Array4<const        amrex::Real>     const& vol,
               const amrex::Real                                & scaling,
               amrex::Array4<      amrex::Real>            const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = D_TERM (  ( fluxX(i+1,j,k,n) - fluxX(i,j,k,n) ),
                             + ( fluxY(i,j+1,k,n) - fluxY(i,j,k,n) ),
                             + ( fluxZ(i,j,k+1,n) - fluxZ(i,j,k,n) ));
      div(i,j,k,n) *= factor;
   }
}


AMREX_GPU_HOST_DEVICE
inline
void
calcMixtFrac(int i, int j, int k,
                amrex::Real                      const& ZO,
                amrex::Real                      const& denom_inv,
                amrex::Real                      const  fact_Bilger[],
                amrex::Array4<const amrex::Real> const& density,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& mixt_frac) noexcept
{
        amrex::Real rho_inv = 1.0 / density(i,j,k);
        mixt_frac(i,j,k) = 0.0;
        for (int n=0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) += ( rhoY(i,j,k,n) * fact_Bilger[n] ) * rho_inv;
        }
         mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - ZO ) * denom_inv;
}

#endif
