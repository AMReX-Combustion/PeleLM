#ifndef PeleLM_K_H_
#define PeleLM_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <mechanism.H>
#include <PelePhysics.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
gpu_flatten(int icell, int i, int j, int k,  
            amrex::Array4<amrex::Real> const& rhoY,
            amrex::Array4<amrex::Real> const& frcExt,
            amrex::Array4<int>         const& EBmask,
            amrex::Real * tmp_vect, amrex::Real * tmp_src_vect,
            amrex::Real * tmp_vect_energy, amrex::Real * tmp_src_vect_energy,
            int * tmp_FCount, int * tmp_EBmask) noexcept
{
    amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
    amrex::Real * TVS_cell = tmp_src_vect + icell * (NUM_SPECIES);
    for (int n = 0; n < NUM_SPECIES; n++) {
       TV_cell[n]   = rhoY(i,j,k,n) * 1.e-3;
       TVS_cell[n]  = frcExt(i,j,k,n) * 1.e-3;
    }
    TV_cell[NUM_SPECIES]       = rhoY(i,j,k,NUM_SPECIES+1);
    tmp_vect_energy[icell]     = rhoY(i,j,k,NUM_SPECIES) * 10.0;
    tmp_src_vect_energy[icell] = frcExt(i,j,k,NUM_SPECIES) * 10.0;
    tmp_FCount[icell]          = -1;
    tmp_EBmask[icell]          = EBmask(i,j,k);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
gpu_unflatten(int icell, int i, int j, int k,
              amrex::Array4<amrex::Real> const& rhoY,
              amrex::Array4<amrex::Real> const& fcl,
              amrex::Real *tmp_vect,
              amrex::Real *tmp_vect_energy,
              int * tmp_FCount) noexcept
{
   amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoY(i,j,k,n) = TV_cell[n] * 1.e+3;
   }
   rhoY(i,j,k,NUM_SPECIES+1) = TV_cell[NUM_SPECIES];
   rhoY(i,j,k,NUM_SPECIES) = tmp_vect_energy[icell] * 1.e-01;
   fcl(i,j,k) = tmp_FCount[icell];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
reactionRateRhoY(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& rhoH,
                 amrex::Array4<const amrex::Real> const& T,
                 amrex::Array4<const int>         const& mask,
                 amrex::Array4<      amrex::Real> const& rhoYdot) noexcept
{
   using namespace amrex::literals;

   // mask < 0 is EB covered cells: no reactions.
   if ( mask(i,j,k) < 0 ) {
      for (int n = 0; n < NUM_SPECIES; n++) {
         rhoYdot(i,j,k,n) = 0.0;
      }
      return;
   }

   auto eos = pele::physics::PhysicsType::eos();
   // Get rho & Ys from rhoYs.
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   // Get T from Y/H.
   amrex::Real Tloc = T(i,j,k);
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;      // Include MKS -> CGS conversion
   eos.HY2T(H, y, Tloc);

   // Get wdot.
   rho = rho * 0.001_rt;                                 // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES] = {0.0_rt};
   eos.RTY2WDOT(rho, Tloc, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;            // CGS -> MKS conversion
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
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

   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   auto eos = pele::physics::PhysicsType::eos();
   eos.RTY2P(rho_cgs, T(i,j,k), y, P(i,j,k));
   P(i,j,k) = P(i,j,k) * 0.1_rt;                         // CGS -> MKS conversion
}

// This kernel returns directly rhoHmix
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getRHmixGivenTY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& Hmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   auto eos = pele::physics::PhysicsType::eos();
   eos.TY2H(T(i,j,k), y, Hmix(i,j,k));
   Hmix(i,j,k) = Hmix(i,j,k) * 1.0e-4_rt * rho(i,j,k);   // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getMwmixGivenRY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& Mwmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   auto eos = pele::physics::PhysicsType::eos();
   eos.Y2WBAR(y, Mwmix(i,j,k));
   Mwmix(i,j,k) = Mwmix(i,j,k) * 0.001_rt;            // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getHGivenT(int i, int j, int k,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Hi) noexcept
{
   using namespace amrex::literals;

   amrex::Real hi_spec[NUM_SPECIES] = {0.0_rt};
   auto eos = pele::physics::PhysicsType::eos();
   eos.T2Hi(T(i,j,k), hi_spec);
   for (int n = 0; n < NUM_SPECIES; n++) {
      Hi(i,j,k,n) = hi_spec[n] * 0.0001_rt;            // CGS -> MKS conversion
   }
   
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
buildAdvectionForcing(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& rho,
                      amrex::Array4<const amrex::Real> const& rhoY,
                      amrex::Array4<const amrex::Real> const& T,
                      amrex::Array4<const amrex::Real> const& dn,
                      amrex::Array4<const amrex::Real> const& ddn,
                      amrex::Array4<const amrex::Real> const& r,
                      amrex::Array4<const amrex::Real> const& extY,
                      amrex::Array4<const amrex::Real> const& extRhoH,
                      amrex::Real                      const& dp0dt,
                      int                              const& closed_chamber,
                      amrex::Array4<      amrex::Real> const& forceY,
                      amrex::Array4<      amrex::Real> const& forceT) noexcept
{
   using namespace amrex::literals;

   // Get species enthalpy
   amrex::Real hi_spec[NUM_SPECIES] = {0.0};
   auto eos = pele::physics::PhysicsType::eos();
   eos.T2Hi(T(i,j,k), hi_spec);

   // For species mass fractions
   amrex::Real rhoinv = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};

   forceT(i,j,k) = dn(i,j,k,NUM_SPECIES+1) + ddn(i,j,k);

   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;                                 // get y
      forceY(i,j,k,n) = r(i,j,k,n) + dn(i,j,k,n);                    // forceY = Rk + Dnk
      forceT(i,j,k) -= forceY(i,j,k,n) * hi_spec[n] * 0.0001_rt;     // forceT -= Sum{ hk. (Rk + Dnk) } with CGS -> MKS conversion of hk
   }
   for (int n = 0; n < NUM_SPECIES; n++) {
      forceY(i,j,k,n) += extY(i,j,k,n);
   }
   forceT(i,j,k) += extRhoH(i,j,k);

   if ( closed_chamber == 1 ) forceT(i,j,k) += dp0dt;                // Add dp0/dt

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real cpmixinv = 1.0_rt / cpmix * 1.0e4_rt;                 // CGS -> MKS conversion
   forceT(i,j,k) *= rhoinv * cpmixinv;

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getCpmixGivenRYT(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rho,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& T,
                 amrex::Array4<      amrex::Real> const& cpmix) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   auto eos = pele::physics::PhysicsType::eos();
   eos.TY2Cp(T(i,j,k), y, cpmix(i,j,k));
   cpmix(i,j,k) *= 0.0001_rt;                         // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getTGivenHY(int i, int j, int k,
            amrex::Array4<const amrex::Real> const& rho,
            amrex::Array4<const amrex::Real> const& rhoY,
            amrex::Array4<const amrex::Real> const& rhoH,
            amrex::Array4<      amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;   // MKS -> CGS conversion
   auto eos = pele::physics::PhysicsType::eos();
   eos.HY2T(H, y, T(i,j,k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
compute_divu(int i, int j, int k,
             amrex::Array4<const amrex::Real> const& rhoY,
             amrex::Array4<const amrex::Real> const& T,
             amrex::Array4<const amrex::Real> const& viscTermT,
             amrex::Array4<const amrex::Real> const& viscTermY,
             amrex::Array4<const amrex::Real> const& rhoYdot,
             amrex::Array4<const amrex::Real> const& extRho,
             amrex::Array4<const amrex::Real> const& extY,
             amrex::Array4<const amrex::Real> const& extRhoH,
             amrex::Array4<      amrex::Real> const& divu) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real mwtinv[NUM_SPECIES] = {0.0};
   eos.inv_molecular_weight(mwtinv);

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real Wbar = 0.0_rt;
   eos.Y2WBAR(y, Wbar);
   amrex::Real hi[NUM_SPECIES] = {0.0};
   eos.T2Hi(T(i,j,k), hi);
   cpmix *= 1.0e-4_rt;                                // CGS -> MKS conversion
   for (int n = 0; n < NUM_SPECIES; n++) {
      hi[n] *= 1.0e-4_rt;                             // CGS -> MKS conversion
   }

   amrex::Real denominv = 1.0_rt / ( rho * cpmix * T(i,j,k) );
   divu(i,j,k) = ( divu(i,j,k) + viscTermT(i,j,k) ) * denominv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      divu(i,j,k) +=  ( viscTermY(i,j,k,n) + rhoYdot(i,j,k,n) )
                    * ( mwtinv[n] * Wbar * rhoinv - hi[n] * denominv );
      divu(i,j,k) += extY(i,j,k,n) * mwtinv[n] * Wbar * rhoinv;
   }
   divu(i,j,k) += extRho(i,j,k) * rhoinv;
   divu(i,j,k) += extRhoH(i,j,k) * denominv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
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
AMREX_FORCE_INLINE
amrex::Real
getGammaInv(int i, int j, int k,
         amrex::Array4<const amrex::Real> const& rhoY,
         amrex::Array4<const amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real cvmix = 0.0_rt;
   eos.TY2Cv(T(i,j,k), y, cvmix);
   
   amrex::Real gammainv = cvmix / cpmix;

   return gammainv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getRhoGivenTYP(int i, int j, int k,
               amrex::Array4<const amrex::Real> const& T,
               amrex::Array4<const amrex::Real> const& Y,
               amrex::Real                      const& Patm,
               amrex::Array4<      amrex::Real> const& rho) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real P_cgs = Patm * pele::physics::Constants::PATM;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = Y(i,j,k,n);                       // TODO: need to find a way not to do that
   }
   amrex::Real rho_cgs = 0.0_rt;
   eos.PYT2R(P_cgs, y, T(i,j,k), rho_cgs);
   rho(i,j,k) = rho_cgs * 1.0e3_rt;            // CGS -> MKS conversion
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getTransportCoeff(int i, int j, int k,
                  amrex::Array4<const amrex::Real> const& rhoY,
                  amrex::Array4<const amrex::Real> const& T,
                  amrex::Array4<      amrex::Real> const& rhoDi,
                  amrex::Array4<      amrex::Real> const& lambda,
                  amrex::Array4<      amrex::Real> const& mu,
                  pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                      pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real mwtinv[NUM_SPECIES] = {0.0};
   eos.inv_molecular_weight(mwtinv);

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real Tloc = T(i,j,k);               // So that we can use const_array

   amrex::Real Wbar = 0.0_rt;
   eos.Y2WBAR(y, Wbar);

   amrex::Real rhoDi_cgs[NUM_SPECIES] = {0.0};
   amrex::Real lambda_cgs = 0.0_rt;
   amrex::Real mu_cgs = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;

   bool get_xi = false;
   bool get_mu = true;
   bool get_lam = true;
   bool get_Ddiag = true;
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
                   rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs, trans_parm);

   // Do CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = rhoDi_cgs[n] * Wbar * mwtinv[n] * 1.0e-1_rt;
   }
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getTransportCoeffUnityLe(int i, int j, int k,
                         amrex::Real                      const& ScInv,
                         amrex::Real                      const& PrInv,
                         amrex::Array4<const amrex::Real> const& rhoY,
                         amrex::Array4<const amrex::Real> const& T,
                         amrex::Array4<      amrex::Real> const& rhoDi,
                         amrex::Array4<      amrex::Real> const& lambda,
                         amrex::Array4<      amrex::Real> const& mu,
                         pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                             pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real Tloc = T(i,j,k);               // So that we can use const_array

   amrex::Real rhoDi_cgs[NUM_SPECIES] = {0.0};
   amrex::Real lambda_cgs = 0.0_rt;
   amrex::Real mu_cgs = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;

   bool get_xi = false;
   bool get_mu = true;
   bool get_lam = false;
   bool get_Ddiag = false;
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
                   rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs, trans_parm);

   mu(i,j,k) = mu_cgs * 1.0e-1_rt;                       // CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = mu_cgs * 1.0e-1_rt * ScInv;       // Constant Schmidt number
   }

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   lambda(i,j,k) = mu_cgs * 1.0e-1_rt * PrInv * cpmix;   // Constant Prandtl number
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getTfromHY(int i, int j, int k,
           amrex::Array4<const amrex::Real> const& rho,
           amrex::Array4<const amrex::Real> const& rhoY,
           amrex::Array4<const amrex::Real> const& rhoH,
           amrex::Array4<      amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   amrex::Real h_cgs = rhoH(i,j,k) * rhoinv * 1.0e4_rt;
   amrex::Real T_loc = T(i,j,k);
   eos.HY2T(h_cgs, y, T_loc);
   T(i,j,k) = T_loc;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getVelViscosity(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& mu,
                pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                    pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real Tloc = T(i,j,k);               // So that we can use const_array

   amrex::Real dummy_rhoDi[NUM_SPECIES] = {0.0};
   amrex::Real dummy_lambda = 0.0_rt;
   amrex::Real mu_cgs = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;

   bool get_xi = false;
   bool get_mu = true;
   bool get_lam = false;
   bool get_Ddiag = false;
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
                   rho, y, dummy_rhoDi, mu_cgs, dummy_xi, dummy_lambda, trans_parm);

   // Do CGS -> MKS conversions
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getConductivity(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& lambda,
                pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                    pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real Tloc = T(i,j,k);               // So that we can use const_array

   amrex::Real dummy_rhoDi[NUM_SPECIES] = {0.0};
   amrex::Real lambda_cgs = 0.0_rt;
   amrex::Real dummy_mu = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;

   bool get_xi = false;
   bool get_mu = false;
   bool get_lam = true;
   bool get_Ddiag = false;
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
                   rho, y, dummy_rhoDi, dummy_mu, dummy_xi, lambda_cgs, trans_parm);

   // Do CGS -> MKS conversions
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
fluxDivergence(int i, int j, int k,
               int ncomp,
               AMREX_D_DECL(amrex::Array4<const amrex::Real>     const& fluxX,
                            amrex::Array4<const amrex::Real>     const& fluxY,
                            amrex::Array4<const amrex::Real>     const& fluxZ),
               amrex::Array4<const        amrex::Real>     const& vol,
               const amrex::Real                                & scaling,
               amrex::Array4<      amrex::Real>            const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = AMREX_D_TERM (  ( fluxX(i+1,j,k,n) - fluxX(i,j,k,n) ),
                                   + ( fluxY(i,j+1,k,n) - fluxY(i,j,k,n) ),
                                   + ( fluxZ(i,j,k+1,n) - fluxZ(i,j,k,n) ));
      div(i,j,k,n) *= factor;
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
fluxDivergenceEB(int i, int j, int k,
                 int ncomp,
                 AMREX_D_DECL(amrex::Array4<const amrex::Real>     const& fluxX,
                              amrex::Array4<const amrex::Real>     const& fluxY,
                              amrex::Array4<const amrex::Real>     const& fluxZ),
               int                                                useEBFlux,
               amrex::Array4<const        amrex::Real>     const& ebflux,
               amrex::Array4<const        amrex::Real>     const& vol,
               const amrex::Real                                & volFracInv,
               const amrex::Real                                & scaling,
               amrex::Array4<      amrex::Real>            const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = AMREX_D_TERM (  ( fluxX(i+1,j,k,n) - fluxX(i,j,k,n) ),
                                   + ( fluxY(i,j+1,k,n) - fluxY(i,j,k,n) ),
                                   + ( fluxZ(i,j,k+1,n) - fluxZ(i,j,k,n) ));
      if (useEBFlux) div(i,j,k,n) += ebflux(i,j,k,n);
      div(i,j,k,n) *= factor * volFracInv;
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
buildDiffusionForcing(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& dn,
                      amrex::Array4<const amrex::Real> const& ddn,
                      amrex::Array4<const amrex::Real> const& dnp1k,
                      amrex::Array4<const amrex::Real> const& ddnp1k,
                      amrex::Array4<const amrex::Real> const& r,
                      amrex::Array4<const amrex::Real> const& a,
                      amrex::Real                      const& dp0dt,
                      int                              const& closed_chamber,
                      amrex::Array4<      amrex::Real> const& forceY,
                      amrex::Array4<      amrex::Real> const& forceT) noexcept
{
   using namespace amrex::literals;

   for (int n = 0; n < NUM_SPECIES; n++) {
      forceY(i,j,k,n) = a(i,j,k,n) + 0.5_rt * ( dn(i,j,k,n) - dnp1k(i,j,k,n) ) + r(i,j,k,n);
   }
   forceT(i,j,k) = a(i,j,k,NUM_SPECIES) + 0.5_rt * (  dn(i,j,k,NUM_SPECIES+1) + ddn(i,j,k) 
                             - dnp1k(i,j,k,NUM_SPECIES+1) - ddnp1k(i,j,k) );

   if ( closed_chamber == 1 ) forceT(i,j,k) += dp0dt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
calcMixtFrac(int i, int j, int k,
                amrex::Real                      const& ZO,
                amrex::Real                      const& denom_inv,
                amrex::Real                      const  fact_Bilger[],
                amrex::Array4<const amrex::Real> const& density,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& mixt_frac) noexcept
{
   using namespace amrex::literals;

   amrex::Real rho_inv = 1.0_rt / density(i,j,k);
   mixt_frac(i,j,k) = 0.0_rt;
   for (int n=0; n<NUM_SPECIES; ++n) {
       mixt_frac(i,j,k) += ( rhoY(i,j,k,n) * fact_Bilger[n] ) * rho_inv;
   }
   mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - ZO ) * denom_inv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
check_divu_dt(int i, int j, int k,
              int                                          const& check_type,
              amrex::Real                                  const& dtfactor,
              amrex::Real                                  const& rhomin,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const dxinv,
              amrex::Array4<const amrex::Real>             const& rho,
              amrex::Array4<const amrex::Real>             const& vel,
              amrex::Array4<const amrex::Real>             const& divu,
              amrex::Array4<const amrex::Real>             const& vol,
              AMREX_D_DECL(amrex::Array4<const amrex::Real>      const& areaX,
                           amrex::Array4<const amrex::Real>      const& areaY,
                           amrex::Array4<const amrex::Real>      const& areaZ),
              amrex::Real                                  const& dt) noexcept
{
   using namespace amrex::literals;

   amrex::Real dtcell = 1.0e12_rt;
   amrex::Real denom  = 0.0_rt;
   bool rho_is_okay = (rho(i,j,k) > rhomin);

   if ( check_type == 1 ) {                           // Check based on divU of DivU_type
      if ( divu(i,j,k) > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( 1.0_rt - rhomin / rho(i,j,k) ) / divu(i,j,k) : 1.0_rt / divu(i,j,k);
      }
   } else if ( check_type == 2 ) {                    // Check based rho * divU + u \cdot grad rho
      denom = rho(i,j,k) * divu(i,j,k) +
              AMREX_D_TERM(  vel(i,j,k,0) * ( rho(i+1,j,k) - rho(i-1,j,k) ) * dxinv[0],
                           + vel(i,j,k,1) * ( rho(i,j+1,k) - rho(i,j-1,k) ) * dxinv[1],
                           + vel(i,j,k,2) * ( rho(i,j,k+1) - rho(i,j,k-1) ) * dxinv[2] );
      if ( denom > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( rho(i,j,k) - rhomin ) / denom : std::abs( rho(i,j,k) ) / denom;
      }
   } else if ( check_type == 3 ) {                    // Check based on div (rho U)
      AMREX_D_TERM( amrex::Real fxlo = 0.25_rt * ( (rho(i,j,k)+rho(i-1,j,k)) * (vel(i,j,k,0)+vel(i-1,j,k,0)) );,
                    amrex::Real fylo = 0.25_rt * ( (rho(i,j,k)+rho(i,j-1,k)) * (vel(i,j,k,1)+vel(i,j-1,k,1)) );,
                    amrex::Real fzlo = 0.25_rt * ( (rho(i,j,k)+rho(i,j,k-1)) * (vel(i,j,k,2)+vel(i,j,k-1,2)) ););
      AMREX_D_TERM( amrex::Real fxhi = 0.25_rt * ( (rho(i,j,k)+rho(i+1,j,k)) * (vel(i,j,k,0)+vel(i+1,j,k,0)) );,
                    amrex::Real fyhi = 0.25_rt * ( (rho(i,j,k)+rho(i,j+1,k)) * (vel(i,j,k,1)+vel(i,j+1,k,1)) );,
                    amrex::Real fzhi = 0.25_rt * ( (rho(i,j,k)+rho(i,j,k+1)) * (vel(i,j,k,2)+vel(i,j,k+1,2)) ););
      denom = AMREX_D_TERM (  areaX(i+1,j,k) * fxhi - areaX(i,j,k) * fxlo,
                            + areaY(i,j+1,k) * fyhi - areaY(i,j,k) * fylo,
                            + areaZ(i,j,k+1) * fzhi - areaZ(i,j,k) * fzlo );
      denom /= vol(i,j,k);
      if ( denom > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( rho(i,j,k) - rhomin ) / denom : std::abs( rho(i,j,k) ) / denom;
      }
   }
   if (dt>dtcell) {
      printf(" ERROR: check_divu_dt - [%i,%i,%i] dtcell %f8.6 > dt %f8.6", i, j, k, dtcell, dt);
   } else if (dt>dtcell*dtfactor) {
      printf(" WARNING: check_divu_dt - [%i,%i,%i] dtcell*dtfactor %f8.6 > dt %f8.6", i, j, k, dtcell*dtfactor, dt);
   }
}

// Full est_divu_dt function - TODO: not available due to limitation of amrex::ReduceMin
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
est_divu_dt(int i, int j, int k,
            int                                          const& step_type,
            amrex::Real                                  const& dtfactor,
            amrex::Real                                  const& rhomin,
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const dxinv,
            amrex::Array4<const amrex::Real>             const& rho,
            amrex::Array4<const amrex::Real>             const& vel,
            amrex::Array4<const amrex::Real>             const& divu,
            amrex::Array4<const amrex::Real>             const& vol,
            AMREX_D_DECL(amrex::Array4<const amrex::Real>      const& areaX,
                         amrex::Array4<const amrex::Real>      const& areaY,
                         amrex::Array4<const amrex::Real>      const& areaZ)) noexcept
{
   using namespace amrex::literals;

   amrex::Real dtcell = 1.0e20_rt;
   amrex::Real denom  = 0.0_rt;
   bool rho_is_okay = (rho(i,j,k) > rhomin);

   if ( step_type == 1 ) {                           // Check based on divU of DivU_type
      if ( divu(i,j,k) > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( 1.0_rt - rhomin / rho(i,j,k) ) / divu(i,j,k) : 1.0_rt / divu(i,j,k);
         dtcell *= dtfactor;
      }
   } else if ( step_type == 2 ) {                    // Check based rho * divU + u \cdot grad rho
      denom = rho(i,j,k) * divu(i,j,k) +
              AMREX_D_TERM(  vel(i,j,k,0) * ( rho(i+1,j,k) - rho(i-1,j,k) ) * dxinv[0],
                           + vel(i,j,k,1) * ( rho(i,j+1,k) - rho(i,j-1,k) ) * dxinv[1],
                           + vel(i,j,k,2) * ( rho(i,j,k+1) - rho(i,j,k-1) ) * dxinv[2] );
      if ( denom > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( rho(i,j,k) - rhomin ) / denom : std::abs( rho(i,j,k) ) / denom;
         dtcell *= dtfactor;
      }
   } else if ( step_type == 3 ) {                    // Check based on div (rho U)
      AMREX_D_TERM( amrex::Real fxlo = 0.25_rt * ( (rho(i,j,k)+rho(i-1,j,k)) * (vel(i,j,k,0)+vel(i-1,j,k,0)) );,
                    amrex::Real fylo = 0.25_rt * ( (rho(i,j,k)+rho(i,j-1,k)) * (vel(i,j,k,1)+vel(i,j-1,k,1)) );,
                    amrex::Real fzlo = 0.25_rt * ( (rho(i,j,k)+rho(i,j,k-1)) * (vel(i,j,k,2)+vel(i,j,k-1,2)) ););
      AMREX_D_TERM( amrex::Real fxhi = 0.25_rt * ( (rho(i,j,k)+rho(i+1,j,k)) * (vel(i,j,k,0)+vel(i+1,j,k,0)) );,
                    amrex::Real fyhi = 0.25_rt * ( (rho(i,j,k)+rho(i,j+1,k)) * (vel(i,j,k,1)+vel(i,j+1,k,1)) );,
                    amrex::Real fzhi = 0.25_rt * ( (rho(i,j,k)+rho(i,j,k+1)) * (vel(i,j,k,2)+vel(i,j,k+1,2)) ););
      denom = AMREX_D_TERM (  areaX(i+1,j,k) * fxhi - areaX(i,j,k) * fxlo,
                            + areaY(i,j+1,k) * fyhi - areaY(i,j,k) * fylo,
                            + areaZ(i,j,k+1) * fzhi - areaZ(i,j,k) * fzlo );
      denom /= vol(i,j,k);
      if ( denom > 0.0_rt ) {
         dtcell = (rho_is_okay) ? ( rho(i,j,k) - rhomin ) / denom : std::abs( rho(i,j,k) ) / denom;
         dtcell *= dtfactor;
      }
   }
   return dtcell;
}

// DivuCeiling type 1 est_divu_dt function
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
est_divu_dt_1(int i, int j, int k,
              int                                          const& step_type,
              amrex::Real                                  const& dtfactor,
              amrex::Real                                  const& rhomin,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const /*dxinv*/,
              amrex::Array4<const amrex::Real>             const& rho,
              amrex::Array4<const amrex::Real>             const& divu) noexcept
{
   using namespace amrex::literals;
   AMREX_ASSERT(step_type==1); 

   amrex::Real dtcell = 1.0e20_rt;
   bool rho_is_okay = (rho(i,j,k) > rhomin);

   if ( divu(i,j,k) > 0.0_rt ) {
      dtcell = (rho_is_okay) ? ( 1.0_rt - rhomin / rho(i,j,k) ) / divu(i,j,k) : 1.0_rt / divu(i,j,k);
      dtcell *= dtfactor;
   }
   return dtcell;
}

// DivuCeiling type 2 est_divu_dt function
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
est_divu_dt_2(int i, int j, int k,
              int                                          const& step_type,
              amrex::Real                                  const& dtfactor,
              amrex::Real                                  const& rhomin,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const dxinv,
              amrex::Array4<const amrex::Real>             const& rho,
              amrex::Array4<const amrex::Real>             const& vel,
              amrex::Array4<const amrex::Real>             const& divu) noexcept
{
   using namespace amrex::literals;
   AMREX_ASSERT(step_type==2); 

   amrex::Real dtcell = 1.0e20_rt;
   amrex::Real denom  = 0.0_rt;
   bool rho_is_okay = (rho(i,j,k) > rhomin);

   denom = rho(i,j,k) * divu(i,j,k) +
           AMREX_D_TERM(  vel(i,j,k,0) * ( rho(i+1,j,k) - rho(i-1,j,k) ) * dxinv[0],
                        + vel(i,j,k,1) * ( rho(i,j+1,k) - rho(i,j-1,k) ) * dxinv[1],
                        + vel(i,j,k,2) * ( rho(i,j,k+1) - rho(i,j,k-1) ) * dxinv[2] );
   if ( denom > 0.0_rt ) {
      dtcell = (rho_is_okay) ? ( rho(i,j,k) - rhomin ) / denom : std::abs( rho(i,j,k) ) / denom;
      dtcell *= dtfactor;
   }
   return dtcell;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
repair_flux(int i, int j, int k,
            int dir,
            const bool                       on_lo,
            const bool                       on_hi,
            amrex::Array4<const amrex::Real> const& rhoY,
            amrex::Array4<      amrex::Real> const& flux ) noexcept
{
   using namespace amrex::literals;

   // Handle the BCs : need a -1 on id_l in dir away from the BC
   // Need -1 in id_l and id_h in dir on low Dirichlet BC
   // Need nothing on high Dirichlet BC
   int  id_l[3] = {i,j,k};
   int  id_h[3] = {i,j,k};
   if ( !on_lo ) {
      if ( !on_hi ) {
        id_l[dir] -= 1;
      }
   } else {
      id_l[dir] -= 1;
      id_h[dir] -= 1;
   }

   amrex::Real sumFlux = 0.0_rt;
   amrex::Real sumRhoYe = 0.0_rt;
   amrex::Real RhoYe[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      sumFlux += flux(i,j,k,n);
      RhoYe[n] = 0.5_rt * (rhoY(id_l[0],id_l[1],id_l[2],n) + rhoY(id_h[0],id_h[1],id_h[2],n));
      sumRhoYe += RhoYe[n];
   }
   sumRhoYe = 1.0_rt / sumRhoYe;
   for (int n = 0; n < NUM_SPECIES; n++) {
      flux(i,j,k,n) -= sumFlux * RhoYe[n] * sumRhoYe;
   }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
repair_flux_eb(int i, int j, int k,
               int dir,
               const bool                       on_lo,
               const bool                       on_hi,
               amrex::Array4<const amrex::Real> const& rhoY,
               amrex::Array4<const amrex::Real> const& rhoYed,
               amrex::Array4<const amrex::Real> const& areafrac,
               amrex::Array4<      amrex::Real> const& flux ) noexcept
{
   using namespace amrex::literals;

   // Handle the BCs : need a -1 on id_l in dir on Dirichlet BC
   int  id_l[3] = {i,j,k};
   id_l[dir] -= 1;

   if ( areafrac(i,j,k) > 0.0_rt ) {
      amrex::Real sumFlux = 0.0_rt;
      amrex::Real sumRhoYe = 0.0_rt;
      amrex::Real RhoYe[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
         sumFlux += flux(i,j,k,n);
         RhoYe[n] = rhoYed(i,j,k,n);
         if ( on_lo ) RhoYe[n] = rhoY(id_l[0],id_l[1],id_l[2],n);
         if ( on_hi ) RhoYe[n] = rhoY(i,j,k,n);
         sumRhoYe += RhoYe[n];
      }
      sumRhoYe = 1.0_rt / sumRhoYe;
      for (int n = 0; n < NUM_SPECIES; n++) {
         flux(i,j,k,n) -= sumFlux * RhoYe[n] * sumRhoYe;
      }
   } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
         flux(i,j,k,n) = 0.0_rt;
      }
   }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
cen2edg_cpp(int i, int j, int k,
            int dir,
            int ncomp,
            const bool                       do_harmonic,
            const bool                       on_lo,
            const bool                       on_hi,
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

   if (do_harmonic) {
      for (int n = 0; n < ncomp; n++) {
         if ( ( cfab(id_l[0],id_l[1],id_l[2],n) * cfab(id_h[0],id_h[1],id_h[2],n) ) > 0.0_rt ) {
            efab(i,j,k,n) = 2.0_rt * ( cfab(id_l[0],id_l[1],id_l[2],n) * cfab(id_h[0],id_h[1],id_h[2],n) )
                                   / ( cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n) );
         } else {
            efab(i,j,k,n) = 0.0_rt;
         }
      }
   } else {
      for (int n = 0; n < ncomp; n++) {
         efab(i,j,k,n) = 0.5_rt*(cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n));
      }
   }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
makeForce(int i, int j, int k,
          int                                     scomp,
          int                                     ncomp,
          int                                     pseudo_gravity,
          const amrex::Real                     & /*time*/,
          const amrex::Real                     & gravity,
          const amrex::Real                     & dV_control,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const /*dx*/,
          amrex::Array4<const amrex::Real> const& /*vel*/,
          amrex::Array4<const amrex::Real> const& scal,
          amrex::Array4<      amrex::Real> const& force) noexcept
{
   using namespace amrex::literals;
   // TODO: component handling here is awfull. Will be updated.

   if ( scomp == 0 ) {
      if ( std::abs(gravity) > 0.0001_rt ) {
         force(i,j,k,0) = 0.0_rt;
#if ( AMREX_SPACEDIM == 2 )
         force(i,j,k,1) = gravity*scal(i,j,k,0);
#elif ( AMREX_SPACEDIM == 3 )
         force(i,j,k,1) = 0.0_rt;
         force(i,j,k,2) = gravity*scal(i,j,k,0);
#endif
      } else {
         for (int n = 0; n < AMREX_SPACEDIM; n++) {
            force(i,j,k,n) = 0.0_rt;
         }
      }

      if (pseudo_gravity) {
         force(i,j,k,AMREX_SPACEDIM-1) += dV_control * scal(i,j,k,0);
      }
   }
   if ( (scomp+ncomp) > AMREX_SPACEDIM - 1) {
      for (int n = amrex::max(scomp,AMREX_SPACEDIM); n < scomp+ncomp; n++) {
         force(i,j,k,n) = 0.0_rt;
      }
   }
}

#endif
