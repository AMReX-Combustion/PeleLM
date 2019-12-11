#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_BCRec.H>

#include <PeleLM.H>

using namespace amrex;

namespace cfinterp
{
void
InterpCCtoFcentEB (const Box& a_bx,
                   D_DECL( FArrayBox& edgestate_x,
                           FArrayBox& edgestate_y,
                           FArrayBox& edgestate_z),
                   const FArrayBox& a_state,
                   const int a_comp,
                   const int a_ncomp,
                   const Box& a_domain,
                   const Vector<BCRec>& a_bcs,
                   D_DECL( const FArrayBox& a_afracx,
                           const FArrayBox& a_afracy,
                           const FArrayBox& a_afracz),
                   D_DECL( const FArrayBox& a_face_centx,
                           const FArrayBox& a_face_centy,
                           const FArrayBox& a_face_centz),
                   const IArrayBox& a_cc_mask,
                   const EBCellFlagFab& a_flags)
{
    const Dim3 domlo = amrex::lbound(a_domain);
    const Dim3 domhi = amrex::ubound(a_domain);

    const auto& state = a_state.array();

    D_TERM( const Box ubx = amrex::surroundingNodes(a_bx,0);,
            const Box vbx = amrex::surroundingNodes(a_bx,1);,
            const Box wbx = amrex::surroundingNodes(a_bx,2););

    D_TERM( const Box ubx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),0);,
            const Box vbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),1);,
            const Box wbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),2););

    D_TERM( FArrayBox s_on_x_face(ubx_grown, a_ncomp);,
            FArrayBox s_on_y_face(vbx_grown, a_ncomp);,
            FArrayBox s_on_z_face(wbx_grown, a_ncomp););

    // These lines ensure that the temporary Fabs above aren't destroyed
    //   before we're done with them when running with GPUs
    D_TERM( Elixir eli_x = s_on_x_face.elixir();,
            Elixir eli_y = s_on_y_face.elixir();,
            Elixir eli_z = s_on_z_face.elixir(););

    D_TERM( const auto& areafrac_x = a_afracx.array();,
            const auto& areafrac_y = a_afracy.array();,
            const auto& areafrac_z = a_afracz.array(););

    D_TERM( const auto& fcx_fab = a_face_centx.array();,
            const auto& fcy_fab = a_face_centy.array();,
            const auto& fcz_fab = a_face_centz.array(););

    D_TERM( const auto& sx = s_on_x_face.array();,
            const auto& sy = s_on_y_face.array();,
            const auto& sz = s_on_z_face.array(););

    D_TERM( const auto& edgs_x = edgestate_x.array();,
            const auto& edgs_y = edgestate_y.array();,
            const auto& edgs_z = edgestate_z.array(););
    
    const auto& ccm_fab = a_cc_mask.const_array();

    const auto bc = a_bcs.dataPtr();

    //
    // ===================== X =====================
    //
    AMREX_FOR_4D(ubx_grown, a_ncomp, i, j, k, n,
    {
      if( areafrac_x(i,j,k) > 0 )
      {
        if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
        {
          sx(i,j,k,n) = state(domlo.x-1,j,k,a_comp+n);
        }
        else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
        {
          sx(i,j,k,n) = state(domhi.x+1,j,k,a_comp+n);
        }
        else
        {
          sx(i,j,k,n) = 0.5 * ( state(i,j,k,a_comp+n) + state(i-1,j,k,a_comp+n) );
        }
      }
      else
      {
        sx(i,j,k,n) = COVERED_VAL;
      }
    });

    // Interpolate to face centroid
    AMREX_FOR_4D(ubx, a_ncomp, i, j, k, n,
    {
      if( areafrac_x(i,j,k) > 0 )
      {
        int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
#if (AMREX_SPACEDIM==3)
        int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));
#else
        int kk = k;
#endif

        Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM==3)
        Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;
#else
        Real fracz(0.0);
#endif

        edgs_x(i,j,k,n) = (1.0-fracy)*(1.0-fracz)*sx(i, j,k ,n)+
            fracy *(1.0-fracz)*sx(i,jj,k ,n)+
            fracz *(1.0-fracy)*sx(i, j,kk,n)+
            fracy *     fracz *sx(i,jj,kk,n);
      }
      else
      {
        edgs_x(i,j,k,n) = COVERED_VAL;
      }
    });

    //
    // ===================== Y =====================
    //
    AMREX_FOR_4D(vbx_grown, a_ncomp, i, j, k, n,
    {
      if( areafrac_y(i,j,k) > 0 )
      {
        if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
          sy(i,j,k,n) = state(i,domlo.y-1,k,a_comp+n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
          sy(i,j,k,n) = state(i,domhi.y+1,k,a_comp+n);
        }
        else
        {
          sy(i,j,k,n) = 0.5 * (state(i,j  ,k,a_comp+n) + state(i,j-1,k,a_comp+n));
        }
      }
      else
      {
        sy(i,j,k,n) = COVERED_VAL;
      }
    });

    AMREX_FOR_4D(vbx, a_ncomp, i, j, k, n,
    {
      if ( areafrac_y(i,j,k) > 0 )
      {
        int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
#if (AMREX_SPACEDIM==3)
        int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));
#else
        int kk = k;
#endif

        Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM==3)
        Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;
#else
        Real fracz(0.0);
#endif

        edgs_y(i,j,k,n) = (1.0-fracx)*(1.0-fracz)*sy(i ,j,k ,n)+
          fracx *(1.0-fracz)*sy(ii,j,k ,n)+
          fracz *(1.0-fracx)*sy(i ,j,kk,n)+
          fracx *     fracz *sy(ii,j,kk,n);
        }
      else
      {
        edgs_y(i,j,k,n) = COVERED_VAL;
      }
    });


    //
    // ===================== Z =====================
    //
#if ( AMREX_SPACEDIM == 3 )
    AMREX_FOR_4D(wbx_grown, a_ncomp, i, j, k, n,
    {
      Real wpls(0.0);
      Real wmns(0.0);

      if( areafrac_z(i,j,k) > 0 )
      {
        if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
        {
          sz(i,j,k,n) = state(i,j,domlo.z-1,a_comp+n);
        }
        else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
        {
          sz(i,j,k,n) = state(i,j,domhi.z+1,a_comp+n);
        }
        else
        {
          sz(i,j,k,n) = 0.5 * ( state(i,j,k,a_comp+n) + state(i,j,k-1,a_comp+n) );
        }
      }
      else
      {
        sz(i,j,k,n) = COVERED_VAL;
      }
    });

    AMREX_FOR_4D(wbx, a_ncomp, i, j, k, n,
    {
      if( areafrac_z(i,j,k) > 0 )
      {
        int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
        int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

        Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
        Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

        edgs_z(i,j,k,n) = (1.0-fracx)*(1.0-fracy)*sz(i ,j ,k,n)+
          fracx *(1.0-fracy)*sz(ii,j ,k,n)+
          fracy *(1.0-fracx)*sz(i ,jj,k,n)+
          fracx *     fracy *sz(ii,jj,k,n);
      }
      else
      {
        edgs_z(i,j,k,n) = COVERED_VAL;
      }
    });
#endif
}

void
InterpCCtoFcent (const Box& a_bx,
                 D_DECL( FArrayBox& edgestate_x,
                         FArrayBox& edgestate_y,
                         FArrayBox& edgestate_z),
                 const FArrayBox& a_state,
                 const int a_comp,
                 const int a_ncomp,
                 const Box& a_domain,
                 const Vector<BCRec>& a_bcs)
{
    const Dim3 domlo = amrex::lbound(a_domain);
    const Dim3 domhi = amrex::ubound(a_domain);

    const auto& state = a_state.array();

    D_TERM( const Box ubx = amrex::surroundingNodes(a_bx,0);,
            const Box vbx = amrex::surroundingNodes(a_bx,1);,
            const Box wbx = amrex::surroundingNodes(a_bx,2););

    D_TERM( const Box ubx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),0);,
            const Box vbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),1);,
            const Box wbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),2););

    D_TERM( const auto& edgs_x = edgestate_x.array();,
            const auto& edgs_y = edgestate_y.array();,
            const auto& edgs_z = edgestate_z.array(););
    
    const auto bc = a_bcs.dataPtr();

    //
    // ===================== X =====================
    //
    AMREX_FOR_4D(ubx_grown, a_ncomp, i, j, k, n,
    {
      if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
      {
        edgs_x(i,j,k,n) = state(domlo.x-1,j,k,a_comp+n);
      }
      else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
      {
        edgs_x(i,j,k,n) = state(domhi.x+1,j,k,a_comp+n);
      }
      else
      {
        edgs_x(i,j,k,n) = 0.5 * ( state(i,j,k,a_comp+n) + state(i-1,j,k,a_comp+n) );
      }
    });

    //
    // ===================== Y =====================
    //
    AMREX_FOR_4D(vbx_grown, a_ncomp, i, j, k, n,
    {
      if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
      {
        edgs_y(i,j,k,n) = state(i,domlo.y-1,k,a_comp+n);
      }
      else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
      {
        edgs_y(i,j,k,n) = state(i,domhi.y+1,k,a_comp+n);
      }
      else
      {
        edgs_y(i,j,k,n) = 0.5 * (state(i,j  ,k,a_comp+n) + state(i,j-1,k,a_comp+n));
      }
    });

    //
    // ===================== Z =====================
    //
#if ( AMREX_SPACEDIM == 3 )
    AMREX_FOR_4D(wbx_grown, a_ncomp, i, j, k, n,
    {
      if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
      {
        edgs_z(i,j,k,n) = state(i,j,domlo.z-1,a_comp+n);
      }
      else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
      {
        edgs_z(i,j,k,n) = state(i,j,domhi.z+1,a_comp+n);
      }
      else
      {
        edgs_z(i,j,k,n) = 0.5 * ( state(i,j,k,a_comp+n) + state(i,j,k-1,a_comp+n) );
      }
    });
#endif
}
} // end of namespace cfinterp

void
PeleLM::InterpCCtoFcent (D_DECL( MultiFab& a_edgestate_x,
                                 MultiFab& a_edgestate_y,
                                 MultiFab& a_edgestate_z),
                         const MultiFab& a_state,
                         const int a_comp,
                         const int a_ncomp,
                         const Box& a_domain,
                         const Geometry& a_geom,
                         const Vector<BCRec>& a_bcs)
{
    AMREX_ALWAYS_ASSERT(a_state.hasEBFabFactory());
    AMREX_ALWAYS_ASSERT(a_state.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(a_bcs.size() == a_ncomp );

    // FIXME: For now use 4 ghost nodes
    const int nghost(4);

    Box domain(a_geom.Domain());

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_state.Factory());

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

    areafrac  =   ebfactory.getAreaFrac();
    facecent  =   ebfactory.getFaceCent();

    // Create cc_mask
    iMultiFab cc_mask(a_state.boxArray(), a_state.DistributionMap(), 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = a_geom.periodicity().shiftIntVect();
        const BoxArray& ba = cc_mask.boxArray();
        for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
        {
            Array4<int> const& fab = cc_mask.array(mfi);

            const Box& bx = mfi.fabbox();
            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);
                for (const auto& is : isects)
                {
                    const Box& b = is.second-iv;
                    AMREX_FOR_3D ( b, i, j, k,
                    {
                        fab(i,j,k) = 1;
                    });
                }
            }
            // NOTE: here we do not need host-device synchronization since it
            // is already included in the MFIter destructor
        }
    }

    // Initialize edge state
    D_TERM(a_edgestate_x.setVal(COVERED_VAL);,
           a_edgestate_y.setVal(COVERED_VAL);,
           a_edgestate_z.setVal(COVERED_VAL););

    for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox ();

        const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
        const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
        {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
            {
              cfinterp::InterpCCtoFcent(bx, D_DECL(a_edgestate_x[mfi], a_edgestate_y[mfi], a_edgestate_z[mfi]),
                                        a_state[mfi], a_comp, a_ncomp, domain, a_bcs);              
            }
            else
            {
              cfinterp::InterpCCtoFcentEB(bx, D_DECL(a_edgestate_x[mfi], a_edgestate_y[mfi], a_edgestate_z[mfi]),
                                          a_state[mfi], a_comp, a_ncomp, domain, a_bcs,
                                          D_DECL((*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi]),
                                          D_DECL((*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi]),
                                          cc_mask[mfi], flags);
            }
        }

    }

    // MR: incflo does not have this: should it be added?
    a_edgestate_x.FillBoundary(a_geom.periodicity());
    a_edgestate_y.FillBoundary(a_geom.periodicity());
#if ( AMREX_SPACEDIM == 3 )
    a_edgestate_z.FillBoundary(a_geom.periodicity());
#endif
}
