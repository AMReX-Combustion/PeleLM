#include "PeleLM.H"
#include "PeleLM_K.H"
#include "PeleLM_derive.H"
#include <EOS.H>
#include <TransportParams.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#endif

#include <mechanism.h>

using namespace amrex;


//
// Multiply by Rho
//

void pelelm_dermprho (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 2);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,1)*in_dat(i,j,k,0);
    });
}

//
// Divide by Rho
//

void pelelm_derdvrho (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 2);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,1)/in_dat(i,j,k,0);
    });
}


//
// rhoY_n, identity ---------- FIXME?
//

void pelelm_derRhoY (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n);
    });
}

//
// Extract species mass fractions Y_n
//

void pelelm_dermassfrac (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = datfab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        der(i,j,k,n) = in_dat(i,j,k,n+1) * rhoinv;
    });
}

//
// Extract species mole fractions X_n
//

void pelelm_dermolefrac (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real Xt[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        for (int n = 0; n < NUM_SPECIES; n++) {
          Yt[n] = in_dat(i,j,k,n+1) * rhoinv;
        } 
        EOS::Y2X(Yt,Xt);
        for (int n = 0; n < NUM_SPECIES; n++) {
          der(i,j,k,n) = Xt[n];
        }
    });

}

//
// Compute the mixture mean molecular weight
//

void pelelm_dermolweight (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        for (int n = 0; n < NUM_SPECIES; n++) {
          Yt[n] = in_dat(i,j,k,n+1) * rhoinv;
        } 
        amrex::Real Wbar = 0.0;
        EOS::Y2WBAR(Yt,Wbar);
        der(i,j,k) = Wbar;
    });

}

//
// Compute the mixture mean heat capacity at cst pressure
//

void pelelm_dercpmix (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& datfab, const Geometry& /*geomdata*/,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+2);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        for (int n = 0; n < NUM_SPECIES; n++) {
          Yt[n] = in_dat(i,j,k,n+2) * rhoinv;
        }
        amrex::Real Temp = in_dat(i,j,k,1);
        amrex::Real Cpmix = 0.0;
        EOS::TY2Cp(Temp,Yt,Cpmix);
        der(i,j,k) = Cpmix * 1.0e-4; // CGS -> MKS ;
    });

}

//
// Compute rho - Sum_n(rhoY_n)
//

void pelelm_drhomry (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& datfab, const Geometry& /*geomdata*/,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    // we put the density in the component 0 of der array
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,0);
    });
    
    // we substract to rho with every rhoY(n)
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k) -= in_dat(i,j,k,n+1);
    });
}

//
// Compute sum of rhoY_dot
//

void pelelm_dsrhoydot (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                       const FArrayBox& datfab, const Geometry& /*geomdata*/,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    AMREX_ASSERT(in_dat.nComp()==NUM_SPECIES);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = 0.0;
    });

    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k) += in_dat(i,j,k,n);
    });
}


//
//  Compute cell-centered pressure as average of the 
//  surrounding nodal values 
//

void pelelm_deravgpres (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& datfab, const Geometry& /*geomdata*/,
                        Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.25;
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.125;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)
                                + in_dat(i+1,j+1,k)   + in_dat(i,j+1,k)
#if (AMREX_SPACEDIM == 3 )
                                + in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)
                                + in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1)
#endif
                                );
    });
}


//
//   Compute cell centered gradient of the nodal pressure in direction x
//

void pelelm_dergrdpx (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    const auto dxinv = geomdata.InvCellSizeArray();    
#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.5*dxinv[0];
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.25*dxinv[0];
#endif
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
#if (AMREX_SPACEDIM == 2 )
        der(i,j,k) = factor*(in_dat(i+1,j,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i,j+1,k));
#elif (AMREX_SPACEDIM == 3 )
        der(i,j,k) = factor*( in_dat(i+1,j,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i,j+1,k) +
                              in_dat(i+1,j,k+1)-in_dat(i,j,k+1)+in_dat(i+1,j+1,k+1)-in_dat(i,j+1,k+1));
#endif 
    });
}

//
//   Compute cell centered gradient of the nodal pressure in direction y
//

void pelelm_dergrdpy (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    const auto dxinv = geomdata.InvCellSizeArray();
#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.5*dxinv[1];
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.25*dxinv[1];
#endif
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
#if (AMREX_SPACEDIM == 2 )
        der(i,j,k) = factor*(in_dat(i,j+1,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i+1,j,k));
#elif (AMREX_SPACEDIM == 3 )
        der(i,j,k) = factor*( in_dat(i,j+1,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i+1,j,k) +
                              in_dat(i,j+1,k+1)-in_dat(i,j,k+1)+in_dat(i+1,j+1,k+1)-in_dat(i+1,j,k+1));
#endif
    });
}


//
//   Compute cell centered gradient of the nodal pressure in direction z
//

void pelelm_dergrdpz (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    const auto dxinv = geomdata.InvCellSizeArray();
    Real factor = 0.25*dxinv[2];
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = factor*( in_dat(i,  j,k+1)-in_dat(i,  j,k)+in_dat(i,  j+1,k+1)-in_dat(i,  j+1,k) +
                              in_dat(i+1,j,k+1)-in_dat(i+1,j,k)+in_dat(i+1,j+1,k+1)-in_dat(i+1,j+1,k));
    });
}

//
//  Compute transport coefficient: D_n, \lambda, \mus
//

void pelelm_dertransportcoeff (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+2);
    AMREX_ASSERT(ncomp == NUM_SPECIES+2);
    auto const T       = datfab.array(1);
    auto const rhoY    = datfab.array(2);
    auto       rhoD    = derfab.array(dcomp);
    auto       lambda  = derfab.array(dcomp+NUM_SPECIES);
    auto       mu      = derfab.array(dcomp+NUM_SPECIES+1);

    // Get the transport GPU data pointer
    TransParm const* ltransparm = trans_parm_g;

    amrex::ParallelFor(bx,
    [T, rhoY, rhoD, lambda, mu, ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getTransportCoeff( i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
    });

}

//
//  Compute both mixt. fraction and scalar diss. rate
//
void pelelm_dermixanddiss (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+2);
    AMREX_ASSERT(ncomp == 2);

    if (!PeleLM::mixture_fraction_ready) amrex::Abort("Mixture fraction not initialized");

    auto const bx_dat     = datfab.box();     //input box (grown to compute mixture fraction --grow_bow_by_one ?)
    auto const bx_der     = derfab.box();     //output box (ungrown)

    auto const density    = datfab.array(0);
    auto const temp       = datfab.array(1);
    auto const rhoY       = datfab.array(2);
    auto       mixt_frac  = derfab.array(dcomp);
    auto       scalar_dis = derfab.array(dcomp+1);

    AMREX_ASSERT(bx_dat.contains(grow(bx,1))); // check that data box is grown
    
    FArrayBox  tmp(bx_dat,1);    // def temporary array for mixt fraction
    //Elixir tmp_e = tmp.elixir();  // not sure we need this
    auto       mixt_frac_tmp  = tmp.array(0);

    const auto dxinv = geomdata.InvCellSizeArray();    
    amrex::Real factor = 0.5*dxinv[0];

    // TODO probably better way to do this ?
    amrex::Real Zox_lcl = PeleLM::Zox;
    amrex::Real Zfu_lcl = PeleLM::Zfu;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_lcl[n] = PeleLM::spec_Bilger_fact[n];
    }

    // Compute Z first -- on grown box
    amrex::ParallelFor(bx_dat,
    [density, rhoY, mixt_frac_tmp, fact_lcl, Zox_lcl, Zfu_lcl] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Get Y from rhoY
        amrex::Real rhoinv = 1.0 / density(i,j,k);
        amrex::Real y[NUM_SPECIES] = {0.0};
        for (int n = 0; n < NUM_SPECIES; n++) {
            y[n] = rhoY(i,j,k,n) * rhoinv;
        }

        mixt_frac_tmp(i,j,k) = 0.0;
        for (int n=0; n<NUM_SPECIES; ++n) {
            mixt_frac_tmp(i,j,k) += y[n] * fact_lcl[n];
        }
        mixt_frac_tmp(i,j,k) = ( mixt_frac_tmp(i,j,k) - Zox_lcl ) / ( Zfu_lcl - Zox_lcl ) ;
    });

    // Get the transport GPU data pointer
    TransParm const* ltransparm = trans_parm_g;

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Copy mixt frac on ungrown box
        mixt_frac(i,j,k) = mixt_frac_tmp(i,j,k);

        // Get Y from rhoY
        amrex::Real rhoinv = 1.0 / density(i,j,k);
        amrex::Real y[NUM_SPECIES] = {0.0};
        for (int n = 0; n < NUM_SPECIES; n++) {
            y[n] = rhoY(i,j,k,n) * rhoinv;
        }

        // get only conductivity
        amrex::Real Tloc = temp(i,j,k);
        amrex::Real rho  = density(i,j,k);

        amrex::Real dummy_rhoDi[NUM_SPECIES] = {0.0};
        amrex::Real dummy_mu = 0.0;
        amrex::Real lambda_cgs = 0.0;
        amrex::Real dummy_xi = 0.0;

        transport(false, false, true, false, 
                  Tloc, rho, y, dummy_rhoDi, dummy_mu, dummy_xi, lambda_cgs, ltransparm);
        amrex::Real lambda = lambda_cgs * 1.0e-5;  // CGS -> MKS 
        
        amrex::Real cpmix = 0.0;
        EOS::TY2Cp(Tloc, y, cpmix);
        cpmix *= 0.0001;                         // CGS -> MKS conversion

        //grad mixt. fraction
        amrex::Real grad[3] = {0.0};
        grad[0] = factor * ( mixt_frac_tmp(i+1,j,k)-mixt_frac_tmp(i-1,j,k) );
        grad[1] = factor * ( mixt_frac_tmp(i,j+1,k)-mixt_frac_tmp(i,j-1,k) );
#if (AMREX_SPACEDIM == 3 )
        grad[2] = factor * ( mixt_frac_tmp(i,j,k+1)-mixt_frac_tmp(i,j,k-1) );
#endif
        scalar_dis(i,j,k) = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
        scalar_dis(i,j,k) *= 2.0 * lambda * rhoinv / cpmix;
    });

}

//
//  Compute Bilger's element based mixture fraction
//

void pelelm_dermixfrac (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == 1);

    if (!PeleLM::mixture_fraction_ready) amrex::Abort("Mixture fraction not initialized");

    auto const density   = datfab.array(0);
    auto const rhoY      = datfab.array(1);
    auto       mixt_frac = derfab.array(0);

    // TODO probably better way to do this ?
    amrex::Real Zox_lcl   = PeleLM::Zox;
    amrex::Real Zfu_lcl   = PeleLM::Zfu;
    amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_lcl[n] = PeleLM::spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, rhoY, mixt_frac, fact_lcl, Zox_lcl, denom_inv] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        calcMixtFrac(i,j,k,
                     Zox_lcl, denom_inv, fact_lcl.data(),
                     density, rhoY, mixt_frac);
    });
}
        

//
// Compute species concentrations C_n
//

void pelelm_derconcentration (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                              const FArrayBox& datfab, const Geometry& /*geomdata*/,
                              Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+2);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real Ct[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        for (int n = 0; n < NUM_SPECIES; n++) {
          Yt[n] = in_dat(i,j,k,n+2) * rhoinv;
        }
        amrex::Real Temp = in_dat(i,j,k,1);
        amrex::Real Rho = in_dat(i,j,k,0) * 1.0e-3; // ! kg/m^3 -> g/cm^3
        EOS::RTY2C(Rho,Temp,Yt,Ct);
        for (int n = 0; n < NUM_SPECIES; n++) {
          der(i,j,k,n) = Ct[n] * 1.0e6; // cm^(-3) -> m^(-3)
        }
    });
}

//
//  Compute the heat release rate -Sum_n(rhoYn_dot * Hn(T))
//
void pelelm_dhrr (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                              const FArrayBox& datfab, const Geometry& /*geomdata*/,
                              Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);

    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        amrex::Real hi_spec[NUM_SPECIES] = {0.0};
        EOS::T2Hi(in_dat(i,j,k,0),hi_spec);
        der(i,j,k) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
            hi_spec[n] *= 0.0001;   // erg/g -> J/kg
            der(i,j,k) -= in_dat(i,j,k,n+1)*hi_spec[n];
        }
    });
}


//
//  Compute the heat release rate -Sum_n(rhoYn_dot * Hn(T)) + the mixture fraction (dcma)
//
void pelelm_dcma (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                              const FArrayBox& datfab, const Geometry& /*geomdata*/,
                              Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 2*NUM_SPECIES+2);
    AMREX_ASSERT(ncomp == 4);
    auto const density   = datfab.array(0);
    auto const rhoY      = datfab.array(1);
    auto const Temp      = derfab.array(NUM_SPECIES+1);
    auto const rhoY_dot  = derfab.array(NUM_SPECIES+2);
    auto       mixt_frac = derfab.array(dcomp);
    auto       hr        = derfab.array(dcomp+1);

#ifdef C12H25O2_ID
    auto       Y1        = derfab.array(dcomp+2);
    auto       Y2        = derfab.array(dcomp+3);
    int        OH        = OH_ID;
    int        RO2       = C12H25O2_ID;
#else
    amrex::Abort("C12H25O2 is not present in your chemistry: do not use dcma.");
#endif

    // TODO probably better way to do this ?
    amrex::Real Zox_lcl   = PeleLM::Zox;
    amrex::Real Zfu_lcl   = PeleLM::Zfu;
    amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_lcl[n] = PeleLM::spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        // Z
        calcMixtFrac(i,j,k,
                     Zox_lcl, denom_inv, fact_lcl.data(),
                     density, rhoY, mixt_frac);

        // HR
        amrex::Real hi_spec[NUM_SPECIES] = {0.0};
        EOS::T2Hi(Temp(i,j,k),hi_spec);
        hr(i,j,k) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
            hi_spec[n] *= 0.0001;   // erg/g -> J/kg
            hr(i,j,k) -= rhoY_dot(i,j,k,n)*hi_spec[n];
        }

        // Y
#ifdef C12H25O2_ID
        amrex::Real rho_inv = 1.0 / density(i,j,k);
        Y1(i,j,k) = rhoY(i,j,k,OH)  * rho_inv; 
        Y2(i,j,k) = rhoY(i,j,k,RO2) * rho_inv; 
#endif
    });
}

//
//  Compute vorticity
//
void pelelm_mgvort (const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    AMREX_D_TERM(const amrex::Real idx = geomdata.InvCellSize(0);,
                 const amrex::Real idy = geomdata.InvCellSize(1);,
                 const amrex::Real idz = geomdata.InvCellSize(2););

    amrex::Array4<amrex::Real const> const& dat_arr = datfab.const_array();
    amrex::Array4<amrex::Real>       const&vort_arr = derfab.array();

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(datfab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();
    auto typ = flags.getType(bx);
    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vort_arr(i,j,k) = 0.0;
        });
    } else if (typ == FabType::singlevalued)
    {
       const auto& flag_fab = flags.const_array();
       amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
          constexpr amrex::Real c0 = -1.5;
          constexpr amrex::Real c1 = 2.0;
          constexpr amrex::Real c2 = -0.5;
          if (flag_fab(i,j,k).isCovered()) {
             vort_arr(i,j,k) = 0.0;
          } else {
             amrex::Real vx = 0.0;
             amrex::Real uy = 0.0;
#if ( AMREX_SPACEDIM == 2 )
             // Need to check if there are covered cells in neighbours --
             // -- if so, use one-sided difference computation (but still quadratic)
             if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                vx = - (c0 * dat_arr(i  ,j,k,1)
                      + c1 * dat_arr(i-1,j,k,1)
                      + c2 * dat_arr(i-2,j,k,1)) * idx;
             } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                vx = (c0 * dat_arr(i  ,j,k,1)
                    + c1 * dat_arr(i+1,j,k,1)
                    + c2 * dat_arr(i+2,j,k,1)) * idx;
             } else {
                vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
             }
             // Do the same in y-direction
             if (!flag_fab(i,j,k).isConnected( 0,1,0)) {
                uy = - (c0 * dat_arr(i,j  ,k,0)
                      + c1 * dat_arr(i,j-1,k,0)
                      + c2 * dat_arr(i,j-2,k,0)) * idy;
             } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                uy = (c0 * dat_arr(i,j  ,k,0)
                    + c1 * dat_arr(i,j+1,k,0)
                    + c2 * dat_arr(i,j+2,k,0)) * idy;
             } else {
                uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
             }
             vort_arr(i,j,k) = vx-uy;
#elif ( AMREX_SPACEDIM == 3 )
             amrex::Real wx = 0.0;
             amrex::Real wy = 0.0;
             amrex::Real uz = 0.0;
             amrex::Real vz = 0.0;
             // Need to check if there are covered cells in neighbours --
             // -- if so, use one-sided difference computation (but still quadratic)
             if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                // Covered cell to the right, go fish left
                vx = - (c0 * dat_arr(i  ,j,k,1)
                      + c1 * dat_arr(i-1,j,k,1)
                      + c2 * dat_arr(i-2,j,k,1)) * idx;
                wx = - (c0 * dat_arr(i  ,j,k,2)
                      + c1 * dat_arr(i-1,j,k,2)
                      + c2 * dat_arr(i-2,j,k,2)) * idx;
             } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                // Covered cell to the left, go fish right
                vx = (c0 * dat_arr(i  ,j,k,1)
                    + c1 * dat_arr(i+1,j,k,1)
                    + c2 * dat_arr(i+2,j,k,1)) * idx;
                wx = (c0 * dat_arr(i  ,j,k,2)
                    + c1 * dat_arr(i+1,j,k,2)
                    + c2 * dat_arr(i+2,j,k,2)) * idx;
             } else {
                // No covered cells right or left, use standard stencil
                vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
             }
             // Do the same in y-direction
             if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                 uy = - (c0 * dat_arr(i,j  ,k,0)
                       + c1 * dat_arr(i,j-1,k,0)
                       + c2 * dat_arr(i,j-2,k,0)) * idy;
                 wy = - (c0 * dat_arr(i,j  ,k,2)
                       + c1 * dat_arr(i,j-1,k,2)
                       + c2 * dat_arr(i,j-2,k,2)) * idy;
             } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                 uy = (c0 * dat_arr(i,j  ,k,0)
                     + c1 * dat_arr(i,j+1,k,0)
                     + c2 * dat_arr(i,j+2,k,0)) * idy;
                 wy = (c0 * dat_arr(i,j  ,k,2)
                     + c1 * dat_arr(i,j+1,k,2)
                     + c2 * dat_arr(i,j+2,k,2)) * idy;
             } else {
                 uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                 wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
             }
             // Do the same in z-direction
             if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                 uz = - (c0 * dat_arr(i,j,k  ,0)
                       + c1 * dat_arr(i,j,k-1,0)
                       + c2 * dat_arr(i,j,k-2,0)) * idz;
                 vz = - (c0 * dat_arr(i,j,k  ,1)
                       + c1 * dat_arr(i,j,k-1,1)
                       + c2 * dat_arr(i,j,k-2,1)) * idz;
             } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                 uz = (c0 * dat_arr(i,j,k  ,0)
                     + c1 * dat_arr(i,j,k+1,0)
                     + c2 * dat_arr(i,j,k+2,0)) * idz;
                 vz = (c0 * dat_arr(i,j,k  ,1)
                     + c1 * dat_arr(i,j,k+1,1)
                     + c2 * dat_arr(i,j,k+2,1)) * idz;
             } else {
                 uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                 vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
             }
             vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
          }
       });
    } else
#endif
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            vort_arr(i,j,k) = vx-uy;

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }
}
