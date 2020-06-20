#include "PeleLM.H"
#include "PeleLM_K.H"
#include "PeleLM_derive.H"

#include <mechanism.h>

using namespace amrex;


//
// Multiply by Rho
//

void pelelm_dermprho (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,dcomp) = in_dat(i,j,k,1)*in_dat(i,j,k,0);
    });
}

//
// Divide by Rho
//

void pelelm_derdvrho (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,dcomp) = in_dat(i,j,k,1)/in_dat(i,j,k,0);
    });
}


//
// Extract species mass rhoY_n
//

void pelelm_derRhoY (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n);
    });
}

//
// Extract species mass rhoY_n
//

void pelelm_dermassfrac (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,0);
        der(i,j,k,n) = in_dat(i,j,k,n+1) * rhoinv;
    });
}


//
// Compute rho - Sum_n(rhoY_n)
//

void pelelm_drhomry (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

    // we put the density in the component 0 of der array
    amrex::ParallelFor(bx, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n);
    });
    
    int nspec_comp = in_dat.nComp() - 1;  //here we get back the correct number of species

    // we substract to rho with every rhoY(n)
    amrex::ParallelFor(bx, nspec_comp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,0) -= in_dat(i,j,k,n+1);
    });
}

//
// Compute sum of rhoY_dot
//

void pelelm_dsrhoydot (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

    amrex::ParallelFor(bx, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = 0.;
    });

    int nspec_comp = in_dat.nComp();  //here we get the correct number of species

    amrex::ParallelFor(bx, nspec_comp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,0) += in_dat(i,j,k,n);
    });
}


//
//  Compute cell-centered pressure as average of the 
//  surrounding nodal values 
//

void pelelm_deravgpres (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.25;
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.125;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,dcomp) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)  
                                      + in_dat(i+1,j+1,k)   + in_dat(i,j+1,k) 
#if (AMREX_SPACEDIM == 3 )
                                      + in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)  
                                      + in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1) 
#endif
                                 );
    });

}


//
//   Compute node centered pressure gradient in direction x
//

void pelelm_dergrdpx (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{   
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

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
        der(i,j,k,dcomp) = factor*(in_dat(i+1,j,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i,j+1,k));
#elif (AMREX_SPACEDIM == 3 )
        der(i,j,k,dcomp) = factor*( in_dat(i+1,j,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i,j+1,k) + 
                                    in_dat(i+1,j,k+1)-in_dat(i,j,k+1)+in_dat(i+1,j+1,k+1)-in_dat(i,j+1,k+1));
#endif 

    });

}

//
//   Compute node centered pressure gradient in direction y
//

void pelelm_dergrdpy (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

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
        der(i,j,k,dcomp) = factor*(in_dat(i,j+1,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i+1,j,k));
#elif (AMREX_SPACEDIM == 3 )
        der(i,j,k,dcomp) = factor*( in_dat(i,j+1,k)-in_dat(i,j,k)+in_dat(i+1,j+1,k)-in_dat(i+1,j,k) + 
                                    in_dat(i,j+1,k+1)-in_dat(i,j,k+1)+in_dat(i+1,j+1,k+1)-in_dat(i+1,j,k+1));
#endif

    });

}


//
//   Compute node centered pressure gradient in direction z
//

void pelelm_dergrdpz (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const in_dat = datfab.array();
    auto       der = derfab.array();

    const auto dxinv = geomdata.InvCellSizeArray();
    Real factor = 0.25*dxinv[2];

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

        der(i,j,k,dcomp) = factor*( in_dat(i,  j,k+1)-in_dat(i,  j,k)+in_dat(i,  j+1,k+1)-in_dat(i,  j+1,k) + 
                                    in_dat(i+1,j,k+1)-in_dat(i+1,j,k)+in_dat(i+1,j+1,k+1)-in_dat(i+1,j+1,k));

    });

}

//
//  Compute transport coefficient: D_n, \lambda, \mus
//

void pelelm_dertransportcoeff (const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    auto const density = datfab.array(0);
    auto const T       = datfab.array(1);
    auto const rhoY    = datfab.array(2);
    auto       rhoD    = derfab.array(0);
    auto       lambda  = derfab.array(NUM_SPECIES);
    auto       mu      = derfab.array(NUM_SPECIES+1);

    amrex::ParallelFor(bx,
    [density, T, rhoY, rhoD, lambda, mu] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getTransportCoeff( i, j, k, rhoY, T, rhoD, lambda, mu);
    });

}

//
//  Compute both mixt. fraction and scalar diss. rate
//
void pelelm_dermixanddiss (const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    if (!PeleLM::mixture_fraction_ready) amrex::Abort("Mixture fraction not initialized");

    auto const density    = datfab.array(0);
    auto const temp       = datfab.array(1);
    auto const rhoY       = datfab.array(2);
    auto       mixt_frac  = derfab.array(0);
    auto       scalar_dis = derfab.array(1);

    const auto dxinv = geomdata.InvCellSizeArray();    
#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.5*dxinv[0];
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.25*dxinv[0];
#endif

    // TODO probably better way to do this ?
    amrex::Real Zox_lcl = PeleLM::Zox;
    amrex::Real Zfu_lcl = PeleLM::Zfu;
    amrex::Real fact_lcl[NUM_SPECIES];
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_lcl[n] = PeleLM::spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, temp, rhoY, mixt_frac, fact_lcl, Zox_lcl, Zfu_lcl] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

        // Get Y from rhoY
        amrex::Real rhoinv;
        rhoinv = 1.0 / density(i,j,k);
        amrex::Real y[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
            y[n] = rhoY(i,j,k,n) * rhoinv;
        }

        mixt_frac(i,j,k) = 0.0;
        for (int n=0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) = mixt_frac(i,j,k) + y[n] * fact_lcl[n];
        }
        mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - Zox_lcl ) / ( Zfu_lcl - Zox_lcl ) ;

        amrex::Real lambda;
        getConductivity(i, j, k, rhoY, temp, lambda);

        amrex::Real cpmix;
        EOS::TY2Cp(temp(i,j,k), y, cpmix);

        cpmix *= 0.0001;                         // CGS -> MKS conversion

        //grad mixt. fraction
        amrex::Real grad[3];
        grad[1] = factor * ( mixt_frac(i+1,j,k)-mixt_frac(i-1,j,k) );
        grad[2] = factor * ( mixt_frac(i,j+1,k)-mixt_frac(i,j-1,k) );
#if (AMREX_SPACEDIM == 3 )
        grad(3) = factor * ( mixt_frac(i,j,k+1)-mixt_frac(i,j,k-1) );
#endif
        scalar_dis(i,j,k) = grad[1]*grad[1] + grad[2]*grad[2] + grad[3]*grad[3];
        scalar_dis(i,j,k) = 2.0 * scalar_dis(i,j,k) * lambda * rhoinv / cpmix;
    });

}

//
//  Compute Bilger's element based mixture fraction
//

void pelelm_dermixfrac (const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)

{
    if (!PeleLM::mixture_fraction_ready) amrex::Abort("Mixture fraction not initialized");

    auto const density   = datfab.array(0);
    auto const rhoY      = datfab.array(1);
    auto       mixt_frac = derfab.array(0);

    // TODO probably better way to do this ?
    amrex::Real Zox_lcl = PeleLM::Zox;
    amrex::Real Zfu_lcl = PeleLM::Zfu;
    amrex::Real fact_lcl[NUM_SPECIES];
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_lcl[n] = PeleLM::spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, rhoY, mixt_frac, fact_lcl, Zox_lcl, Zfu_lcl] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        mixt_frac(i,j,k) = 0.0;
        for (int n=0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) = mixt_frac(i,j,k) + ( rhoY(i,j,k,n) * fact_lcl[n] ) / density(i,j,k);
        }
         mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - Zox_lcl ) / ( Zfu_lcl - Zox_lcl ) ;
        
    });
}
