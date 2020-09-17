#include <PeleLM_EF_Constant.H>

//
//  Compute electric efield in X direction
//
void pelelm_derefx (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der     = derfab.box();     //output box (ungrown)
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[0];

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       bool on_lo = ( ( bcrec[0] == EXT_DIR ) && i <= domain.smallEnd(0) );
       bool on_hi = ( ( bcrec[0+AMREX_SPACEDIM] == EXT_DIR ) && i >= domain.bigEnd(0) );
       der(i,j,k) = factor * ( in_dat(i+1,j,k) - in_dat(i-1,j,k) );
       if ( on_lo ) der(i,j,k) = factor * ( in_dat(i+1,j,k) + in_dat(i,j,k) - 2.0 * in_dat(i-1,j,k) ) ;
       if ( on_hi ) der(i,j,k) = factor * ( 2.0 * in_dat(i+1,j,k) - in_dat(i,j,k) - in_dat(i-1,j,k) ) ;
    });
}

//
//  Compute electric efield in Y direction
//
void pelelm_derefy (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der     = derfab.box();     //output box (ungrown)
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[1];

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       bool on_lo = ( ( bcrec[1] == EXT_DIR ) && j <= domain.smallEnd(1) );
       bool on_hi = ( ( bcrec[1+AMREX_SPACEDIM] == EXT_DIR ) && j >= domain.bigEnd(1) );
       der(i,j,k) = factor * ( in_dat(i,j+1,k) - in_dat(i,j-1,k) );
       if ( on_lo ) der(i,j,k) = factor * ( in_dat(i,j+1,k) + in_dat(i,j,k) - 2.0 * in_dat(i,j-1,k) ) ;
       if ( on_hi ) der(i,j,k) = factor * ( 2.0 * in_dat(i,j+1,k) - in_dat(i,j,k) - in_dat(i,j-1,k) ) ;
    });
}

#if ( AMREX_SPACEDIM == 3 )
//
//  Compute electric field in Z direction
//
void pelelm_derefz (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der     = derfab.box();     //output box (ungrown)
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[2];

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       bool on_lo = ( ( bcrec[2] == EXT_DIR ) && k <= domain.smallEnd(2) );
       bool on_hi = ( ( bcrec[2+AMREX_SPACEDIM] == EXT_DIR ) && k >= domain.bigEnd(2) );
       der(i,j,k) = factor * ( in_dat(i,j,k+1) - in_dat(i,j,k-1) );
       if ( on_lo ) der(i,j,k) = factor * ( in_dat(i,j,k+1) + in_dat(i,j,k) - 2.0 * in_dat(i,j,k-1) ) ;
       if ( on_hi ) der(i,j,k) = factor * ( 2.0 * in_dat(i,j,k+1) - in_dat(i,j,k) - in_dat(i,j,k-1) ) ;
    });
}
#endif

//
//  Compute charge distribution
//
void pelelm_dercd (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der   = derfab.box();     //output box (ungrown)
    auto const ne_arr   = datfab.array(0);
    auto const rhoY_arr = datfab.array(1);
    auto          der   = derfab.array(dcomp);

    amrex::GpuArray<amrex::Real,NUM_SPECIES> zk_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
       zk_lcl[n] = PeleLM::zk[n];
    }
    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       der(i,j,k) = - ne_arr(i,j,k) * EFConst::elemCharge;
       for (int n = 0; n < NUM_SPECIES; ++n ) {
          der(i,j,k) += zk_lcl[n] * rhoY_arr(i,j,k,n);
       }
    });
}

//
//  Compute Lorentz forces in X direction
//
void pelelm_derLorentzx (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der   = derfab.box();     //output box (ungrown)
    auto const phi_arr  = datfab.array(0);
    auto const ne_arr   = datfab.array(1);
    auto const rhoY_arr = datfab.array(2);
    auto          der   = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[0];

    amrex::GpuArray<amrex::Real,NUM_SPECIES> zk_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
       zk_lcl[n] = PeleLM::zk[n];
    }

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       // Get gradient of PhiV
       amrex::Real EFx;
       bool on_lo = ( ( bcrec[0] == EXT_DIR ) && i <= domain.smallEnd(0) );
       bool on_hi = ( ( bcrec[0+AMREX_SPACEDIM] == EXT_DIR ) && i >= domain.bigEnd(0) );
       EFx = factor * ( phi_arr(i+1,j,k) - phi_arr(i-1,j,k) );
       if ( on_lo ) EFx = factor * ( phi_arr(i+1,j,k) + phi_arr(i,j,k) - 2.0 * phi_arr(i-1,j,k) ) ;
       if ( on_hi ) EFx = factor * ( 2.0 * phi_arr(i+1,j,k) - phi_arr(i,j,k) - phi_arr(i-1,j,k) ) ;

       // Assemble Lorentz force in X
       der(i,j,k) = - ne_arr(i,j,k) * EFConst::elemCharge * EFx;
       for (int n = 0; n < NUM_SPECIES; ++n ) {
          der(i,j,k) += zk_lcl[n] * rhoY_arr(i,j,k,n) * EFx;
       }
    });
}

//
//  Compute Lorentz forces in Y direction
//
void pelelm_derLorentzy (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der   = derfab.box();     //output box (ungrown)
    auto const phi_arr  = datfab.array(0);
    auto const ne_arr   = datfab.array(1);
    auto const rhoY_arr = datfab.array(2);
    auto          der   = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[0];

    amrex::GpuArray<amrex::Real,NUM_SPECIES> zk_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
       zk_lcl[n] = PeleLM::zk[n];
    }

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       // Get gradient of PhiV
       amrex::Real EFy;
       bool on_lo = ( ( bcrec[1] == EXT_DIR ) && j <= domain.smallEnd(1) );
       bool on_hi = ( ( bcrec[1+AMREX_SPACEDIM] == EXT_DIR ) && j >= domain.bigEnd(1) );
       EFy = factor * ( phi_arr(i,j+1,k) - phi_arr(i,j-1,k) );
       if ( on_lo ) EFy = factor * ( phi_arr(i,j+1,k) + phi_arr(i,j,k) - 2.0 * phi_arr(i,j-1,k) ) ;
       if ( on_hi ) EFy = factor * ( 2.0 * phi_arr(i,j+1,k) - phi_arr(i,j,k) - phi_arr(i,j-1,k) ) ;

       // Assemble Lorentz force in Y
       der(i,j,k) = - ne_arr(i,j,k) * EFConst::elemCharge * EFy;
       for (int n = 0; n < NUM_SPECIES; ++n ) {
          der(i,j,k) += zk_lcl[n] * rhoY_arr(i,j,k,n) * EFy;
       }
    });
}

#if ( AMREX_SPACEDIM == 3 )
//
//  Compute Lorentz forces in Z direction
//
void pelelm_derLorentzz (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real time, const int* bcrec, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    auto const bx_der   = derfab.box();     //output box (ungrown)
    auto const phi_arr  = datfab.array(0);
    auto const ne_arr   = datfab.array(1);
    auto const rhoY_arr = datfab.array(2);
    auto          der   = derfab.array(dcomp);

    const auto dxinv = geomdata.InvCellSizeArray();
    const auto domain = geomdata.Domain();
    amrex::Real factor = -0.5*dxinv[0];

    amrex::GpuArray<amrex::Real,NUM_SPECIES> zk_lcl;
    for (int n=0; n<NUM_SPECIES; ++n) {
       zk_lcl[n] = PeleLM::zk[n];
    }

    amrex::ParallelFor(bx_der,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       // Get gradient of PhiV
       amrex::Real EFz;
       bool on_lo = ( ( bcrec[2] == EXT_DIR ) && k <= domain.smallEnd(2) );
       bool on_hi = ( ( bcrec[2+AMREX_SPACEDIM] == EXT_DIR ) && k >= domain.bigEnd(2) );
       EFz = factor * ( phi_arr(i,j,k+1) - phi_arr(i,j,k-1) );
       if ( on_lo ) EFz = factor * ( phi_arr(i,j,k+1) + phi_arr(i,j,k) - 2.0 * phi_arr(i,j,k-1) ) ;
       if ( on_hi ) EFz = factor * ( 2.0 * phi_arr(i,j,k+1) - phi_arr(i,j,k) - phi_arr(i,j,k-1) ) ;

       // Assemble Lorentz force in Z
       der(i,j,k) = - ne_arr(i,j,k) * EFConst::elemCharge * EFz;
       for (int n = 0; n < NUM_SPECIES; ++n ) {
          der(i,j,k) += zk_lcl[n] * rhoY_arr(i,j,k,n) * EFz;
       }
    });
}
#endif
