#include "PeleLM.H"
#include "PeleLM_derive.H"

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
        der(i,j,k,0) = der(i,j,k,0) - in_dat(i,j,k,n+1);
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
        der(i,j,k,0) = der(i,j,k,0) + in_dat(i,j,k,n);
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

#if (AMREX_SPACEDIM >= 2 )
     Real factor = 0.25;
#elif (AMREX_SPACEDIM == 3 )
     Real factor = 0.125;
#endif

    amrex::ParallelFor(bx, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) =  factor * (  in_dat(i+1,j,k,n)     + in_dat(i,j,k,n)  
                                  + in_dat(i+1,j+1,k,n)   + in_dat(i,j+1,k,n) 
#if (AMREX_SPACEDIM == 3 )
                                  + in_dat(i+1,j,k+1,n)   + in_dat(i,j,k+1,n)  
                                  + in_dat(i+1,j+1,k+1,n) + in_dat(i,j+1,k+1,n) 
#endif
                                 );
    });

}





