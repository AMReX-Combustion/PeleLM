#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <Prob_F.H>
#include <AMReX_ArrayLim.H>
#include "mechanism.h"

module derive_PLM_nd

  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use fuego_chemistry

  implicit none

  private

  public :: dermgvort, dermgdivu, & 
            dhrr, dcma

  REAL_T, dimension(NUM_SPECIES,NUM_ELEMENTS) :: coeff_mix
  REAL_T, dimension(NUM_ELEMENTS) :: beta_mix
  REAL_T, dimension(NUM_SPECIES)  :: fact
  REAL_T :: Zfu, Zox
  logical :: init_mixture = .FALSE.

contains
 
!=========================================================
!  Compute the amagnitude of the vorticity from the 
!  velocity field
!=========================================================

   subroutine dermgvort (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgvort")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T :: uy, uz, vx, vz, wx, wy, dx, dy, dz
      REAL_T :: uycen, uzcen, uylo, uyhi, uzlo, uzhi
      REAL_T :: vxcen, vzcen, vxlo, vxhi, vzlo, vzhi
      REAL_T :: wxcen, wycen, wxlo, wxhi, wylo, wyhi
      REAL_T :: vorfun

      logical :: fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical :: fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical :: fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z

      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)

!
!     ::::: statement functions that implement stencil
!
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

#if ( AMREX_SPACEDIM == 3 )
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)
#endif

#if ( AMREX_SPACEDIM == 2 )
      vorfun(uy,uz,vx,vz,wx,wy) = vx - uy
#elif ( AMREX_SPACEDIM == 3 )
      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! Init all logical tests on BC to false
      fixvlo_x = .FALSE. ; fixwlo_x = .FALSE. ; fixvhi_x = .FALSE. ; fixwhi_x = .FALSE.
      fixulo_y = .FALSE. ; fixwlo_y = .FALSE. ; fixuhi_y = .FALSE. ; fixwhi_y = .FALSE.
      fixulo_z = .FALSE. ; fixvlo_z = .FALSE. ; fixuhi_z = .FALSE. ; fixvhi_z = .FALSE.

      ! Init all vorticity comp. In 2d uz, vz, wx, wy will alway be zero
      uy = 0.0d0
      uz = 0.0d0
      vx = 0.0d0
      vz = 0.0d0
      wx = 0.0d0
      wy = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )
#endif

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
#endif

!
!     First do all the faces
!
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
#endif

!
!     Next do all the edges
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

!
!     Finally do all the corners
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

#     undef U
#     undef V
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

   end subroutine dermgvort

!=========================================================
!  Compute the magnitude of the velocity divergence
!=========================================================

   subroutine dermgdivu (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgdivu")

      use bc_fill_nd_module, only : xvel_fill, yvel_fill, zvel_fill

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: ux, vy, wz, dx, dy, dz
      REAL_T  :: uxcen, uxlo, uxhi
      REAL_T  :: vycen, vylo, vyhi
#if ( AMREX_SPACEDIM == 3 )
      REAL_T  :: wzcen, wzlo, wzhi
#endif
      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)

!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k) = (eight*U(i,j,k)-six*U(i+1,j,k)+U(i+2,j,k))/(three*dx)
      uxhi(i,j,k) = (eight*U(i,j,k)-six*U(i-1,j,k)+U(i-2,j,k))/(three*dx)

#if ( AMREX_SPACEDIM >= 2 )
      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k) = (eight*V(i,j,k)-six*V(i,j+1,k)+V(i,j+2,k))/(three*dy)
      vyhi(i,j,k) = (eight*V(i,j,k)-six*V(i,j-1,k)+V(i,j-2,k))/(three*dy)

#if ( AMREX_SPACEDIM == 3 )
      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k+1)+W(i,j,k+2))/(three*dz)
      wzhi(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k-1)+W(i,j,k-2))/(three*dz)
#endif
#endif

      call xvel_fill (dat(:,:,:,1), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,1))
#if ( AMREX_SPACEDIM >= 2 )
      call yvel_fill (dat(:,:,:,2), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,2))
#if ( AMREX_SPACEDIM == 3 )
      call zvel_fill (dat(:,:,:,3), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

!
!     :: at physical bndries where an edge value is prescribed,
!     :: set the value in the outside cell so that a central
!     :: difference formula is equivalent to the higher order
!     :: one sided formula
!
      if (lo(1) == domlo(1)) then
         i = lo(1)
         if (ULOX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = two*U(i-1,j,k) - U(i,j,k)
               end do
            end do
         else if (ULOX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = uxlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(1) == domhi(1)) then
         i = hi(1)
         if (UHIX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = two*U(i+1,j,k) - U(i,j,k)
               end do
            end do
         else if (UHIX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = uxhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM >= 2 )
      if (lo(2) == domlo(2)) then
         j = lo(2)
         if (VLOY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = two*V(i,j-1,k) - V(i,j,k)
               end do
            end do
         else if (VLOY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = vylo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(2) == domhi(2)) then
         j = hi(2)
         if (VHIY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = two*V(i,j+1,k) - V(i,j,k)
               end do
            end do
         else if (VHIY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = vyhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (lo(3) == domlo(3)) then
         k = lo(3)
         if (WLOZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = two*W(i,j,k-1) - W(i,j,k)
               end do
            end do
         else if (WLOZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = wzlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(3) == domhi(3)) then
         k = hi(3)
         if (WHIZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = two*W(i,j,k+1) - W(i,j,k)
               end do
            end do
         else if (WHIZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = wzhi(i,j,k)
               end do
            end do
         end if
      end if
#endif
#endif

      ux = 0.0d0
      vy = 0.0d0
      wz = 0.0d0
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
#if ( AMREX_SPACEDIM >= 2 )
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wz = wzcen(i,j,k)
#endif
#endif
               e(i,j,k,1) = ux + vy + wz
            end do
         end do
      end do

!
! we overwrote the ghost cells above, so set them back below
!
      call xvel_fill (dat(:,:,:,1), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,1))
#if ( AMREX_SPACEDIM >= 2 )
      call yvel_fill (dat(:,:,:,2), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,2))
#if ( AMREX_SPACEDIM == 3 )
      call zvel_fill (dat(:,:,:,3), d_lo, d_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif
#endif

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY
#     undef WLOZ
#     undef WHIZ

   end subroutine dermgdivu


!=========================================================
!  Compute CMA
!=========================================================

   subroutine dcma (e,   e_lo, e_hi, nv, &
                    dat, d_lo, d_hi, ncomp, &
                    lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                    level, grid_no) &
                    bind(C, name="dcma")

      use fuego_chemistry, only : get_species_index

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

! Local
      REAL_T  :: rhoinv
      integer :: rho, T, fS, OH, RO2

      integer :: i, j, k

      rho = 1
      fS  = 2
      T   = 2 + NUM_SPECIES

! TODO : this should be specified somewhere in the probin file
      OH = get_species_index('OH')
      RO2 = get_species_index('C12H25O2')

      !CALL  dermixfrac(e(:,:,:,1),  e_lo, e_hi, 1, &
      !                 dat(:,:,:,rho:NUM_SPECIES+1), d_lo, d_hi, NUM_SPECIES+1, &
      !                 lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
      !                 level, grid_no)

      CALL  dhrr(e(:,:,:,2), e_lo, e_hi, 1, &
                 dat(:,:,:,T:ncomp), d_lo, d_hi, NUM_SPECIES+1, &
                 lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                 level, grid_no)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoinv = dat(i,j,k,rho)
               e(i,j,k,3) = dat(i,j,k,fS+OH-1) * rhoinv
               e(i,j,k,4) = dat(i,j,k,fS+RO2-1) * rhoinv
            enddo
         enddo
      enddo

   end subroutine dcma

!=========================================================
!  Compute forcing term
!  Reference : ??   
!=========================================================

   subroutine FORT_DERFORCING (e,   e_lo, e_hi, nv, &
                               dat, d_lo, d_hi, ncomp, &
                               lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                               level, grid_no) &
                               bind(C, name="FORT_DERFORCING")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

! Local
#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>

      integer :: ilo, jlo, klo
      integer :: ihi, jhi, khi
      REAL_T  :: x, y, z
      REAL_T  :: hx, hy, hz
      REAL_T  :: sga, cga
      REAL_T  :: f1, f2, f3
      REAL_T  :: twicePi
      REAL_T  :: force_time
      REAL_T  :: kxd, kyd, kzd
      REAL_T  :: xt, yt, zt
      REAL_T  :: HLx_inv, HLy_inv, HLz_inv
      REAL_T  :: Lx, Ly, Lz, Lmin, kappa, kappaMax
      REAL_T  :: rho, u, v, w
      integer :: kx, ky, kz, mode_count, xstep, ystep, zstep
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho

      integer :: i, j, k, n

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi = two*Pi

!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif

      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz==1) then
         Lz = Lz/two
      endif

      Lmin = min(Lx,Ly,Lz)
      kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(0.5d0+Lx/Lmin)
      nymodes = nmodes*int(0.5d0+Ly/Lmin)
      nzmodes = nmodes*int(0.5d0+Lz/Lmin)

      xstep = int(Lx/Lmin+0.5d0)
      ystep = int(Ly/Lmin+0.5d0)
      zstep = int(Lz/Lmin+0.5d0)

      if (forcing_twice_wavelength==1) then
         HLx_inv = two/Lx
         HLy_inv = two/Ly
         HLz_inv = two/Lz
      else
         HLx_inv = 1.0d0/Lx
         HLy_inv = 1.0d0/Ly
         HLz_inv = 1.0d0/Lz
      endif

      do k = klo, khi
         z = xlo(3) + hz*(float(k-klo) + half)
         do j = jlo, jhi
            y = xlo(2) + hy*(float(j-jlo) + half)
            do i = ilo, ihi
               x = xlo(1) + hx*(float(i-ilo) + half)
               f1 = zero
               f2 = zero
               f3 = zero
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dfloat(kz)
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dfloat(ky)
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f1 = f1 + xT * (  FAZ(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) &
                                              - FAY(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * (  FAX(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) &
                                              - FAZ(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * (  FAY(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) &
                                              - FAX(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT * FAX(kx,ky,kz) * cos(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                              f2 = f2 + xT * FAY(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             cos(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                              f3 = f3 + xT * FAZ(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             cos(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               do kz = 1, zstep - 1
                  kzd = dfloat(kz)
                  do ky = mode_start, nymodes
                     kyd = dfloat(ky)
                     do kx = mode_start, nxmodes
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f1 = f1 + xT * (  FAZ(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) &
                                              - FAY(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * (  FAX(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) &
                                              - FAZ(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * (  FAY(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) &
                                              - FAX(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT * FAX(kx,ky,kz) * cos(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                            sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                            sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                              f2 = f2 + xT * FAY(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             cos(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                              f3 = f3 + xT * FAZ(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             cos(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing==1) then
                  rho = dat(i,j,k,1)
                  u   = dat(i,j,k,2)
                  v   = dat(i,j,k,3)
                  w   = dat(i,j,k,4)
                  e(i,j,k,1) = f1*rho*u
                  e(i,j,k,2) = f2*rho*v
                  e(i,j,k,3) = f3*rho*w
               else
                  e(i,j,k,1) = f1*u
                  e(i,j,k,2) = f2*v
                  e(i,j,k,3) = f3*w
               endif
            enddo
         enddo
      enddo
#endif

   end subroutine FORT_DERFORCING 

!=========================================================
!  Compute forcing term in x-dir
!  Reference : ??   
!=========================================================

   subroutine FORT_DERFORCEX (e,   e_lo, e_hi, nv, &
                              dat, d_lo, d_hi, ncomp, &
                              lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                              level, grid_no) &
                              bind(C, name="FORT_DERFORCEX")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

! Local
#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer :: ilo, jlo, klo
      integer :: ihi, jhi, khi
      integer :: a2, a3, a4, a5
      REAL_T  :: x, y, z
      REAL_T  :: hx, hy, hz
      REAL_T  :: sga, cga
      REAL_T  :: f1
      REAL_T  :: twicePi
      REAL_T  :: force_time
      REAL_T  :: kxd, kyd, kzd
      REAL_T  :: xt, yt, zt
      REAL_T  :: HLx_inv, HLy_inv, HLz_inv
      REAL_T  :: Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer :: kx, ky, kz, mode_count, xstep, ystep, zstep
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho

      integer :: i, j, k, n

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi = two*Pi
      
!     Adjust z offset for probtype 15
      if (time_offset>(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif
      
      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz==1) then
         Lz = Lz/two
      endif

      Lmin = min(Lx,Ly,Lz)
      kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(0.5+Lx/Lmin)
      nymodes = nmodes*int(0.5+Ly/Lmin)
      nzmodes = nmodes*int(0.5+Lz/Lmin)

      xstep = int(Lx/Lmin+0.5)
      ystep = int(Ly/Lmin+0.5)
      zstep = int(Lz/Lmin+0.5)

      if (forcing_twice_wavelength==1) then
         HLx_inv = two/Lx
         HLy_inv = two/Ly
         HLz_inv = two/Lz
      else
         HLx_inv = 1.0d0/Lx
         HLy_inv = 1.0d0/Ly
         HLz_inv = 1.0d0/Lz
      endif

      do k = klo, khi
         z = xlo(3) + hz*(float(k-klo) + half)
         do j = jlo, jhi
            y = xlo(2) + hy*(float(j-jlo) + half)
            do i = ilo, ihi
               x = xlo(1) + hx*(float(i-ilo) + half)
               f1 = zero
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dfloat(kz)
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dfloat(ky)
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f1 = f1 + xT * (  FAZ(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) &
                                              - FAY(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT * FAX(kx,ky,kz) * cos(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               do kz = 1, zstep - 1
                  kzd = dfloat(kz)
                  do ky = mode_start, nymodes
                     kyd = dfloat(ky)
                     do kx = mode_start, nxmodes
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f1 = f1 + xT * (  FAZ(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) &
                                              - FAY(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT * FAX(kx,ky,kz) * cos(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing==1) then
                  e(i,j,k,1) = f1*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f1
               endif
            enddo
         enddo
      enddo
#endif

  end subroutine FORT_DERFORCEX

!=========================================================
!  Compute forcing term in y-dir
!  Reference : ??   
!=========================================================

   subroutine FORT_DERFORCEY (e,   e_lo, e_hi, nv, &
                              dat, d_lo, d_hi, ncomp, &
                              lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                              level, grid_no) &
                              bind(C, name="FORT_DERFORCEY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

! Local
#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer :: ilo, jlo, klo
      integer :: ihi, jhi, khi
      integer :: a2, a3, a4, a5
      REAL_T  :: x, y, z
      REAL_T  :: hx, hy, hz
      REAL_T  :: sga, cga
      REAL_T  :: f2
      REAL_T  :: twicePi
      REAL_T  :: force_time
      REAL_T  :: kxd, kyd, kzd
      REAL_T  :: xt, yt, zt
      REAL_T  :: HLx_inv, HLy_inv, HLz_inv
      REAL_T  :: Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer :: kx, ky, kz, mode_count, xstep, ystep, zstep
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho

      integer :: i, j, k, n

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi = two*Pi

!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif

      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz==1) then
         Lz = Lz/two
      endif

      Lmin = min(Lx,Ly,Lz)
      kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(0.5+Lx/Lmin)
      nymodes = nmodes*int(0.5+Ly/Lmin)
      nzmodes = nmodes*int(0.5+Lz/Lmin)

      xstep = int(Lx/Lmin+0.5)
      ystep = int(Ly/Lmin+0.5)
      zstep = int(Lz/Lmin+0.5)

      if (forcing_twice_wavelength==1) then
         HLx_inv = two/Lx
         HLy_inv = two/Ly
         HLz_inv = two/Lz
      else
         HLx_inv = 1.0d0/Lx
         HLy_inv = 1.0d0/Ly
         HLz_inv = 1.0d0/Lz
      endif

      do k = klo, khi
         z = xlo(3) + hz*(float(k-klo) + half)
         do j = jlo, jhi
            y = xlo(2) + hy*(float(j-jlo) + half)
            do i = ilo, ihi
               x = xlo(1) + hx*(float(i-ilo) + half)
               f2 = zero
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dfloat(kz)
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dfloat(ky)
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f2 = f2 + xT * (  FAX(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) &
                                              - FAZ(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT * FAY(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             cos(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               do kz = 1, zstep - 1
                  kzd = dfloat(kz)
                  do ky = mode_start, nymodes
                     kyd = dfloat(ky)
                     do kx = mode_start, nxmodes
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f2 = f2 + xT * (  FAX(kx,ky,kz) * twicePi * (kzd*HLz_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                cos(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) &
                                              - FAZ(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPZX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPZY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT * FAY(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             cos(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             sin(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing==1) then
                  e(i,j,k,1) = f2*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f2
               endif
            enddo
         enddo
      enddo
#endif

   end subroutine FORT_DERFORCEY

!=========================================================
!  Compute forcing term in z-dir
!  Reference : ??   
!=========================================================

   subroutine FORT_DERFORCEZ (e,   e_lo, e_hi, nv, &
                              dat, d_lo, d_hi, ncomp, &
                              lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                              level, grid_no) &
                              bind(C, name="FORT_DERFORCEZ")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

! Local
#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer :: ilo, jlo, klo
      integer :: ihi, jhi, khi
      integer :: a2, a3, a4, a5
      REAL_T  :: x, y, z
      REAL_T  :: hx, hy, hz
      REAL_T  :: sga, cga
      REAL_T  :: f3
      REAL_T  :: twicePi
      REAL_T  :: force_time
      REAL_T  :: kxd, kyd, kzd
      REAL_T  :: xt, yt, zt
      REAL_T  :: HLx_inv, HLy_inv, HLz_inv
      REAL_T  :: Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer :: kx, ky, kz, mode_count, xstep, ystep, zstep
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho

      integer :: i, j, k, n

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi = two * Pi

!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif

      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz==1) then
         Lz = Lz/two
      endif

      Lmin = min(Lx,Ly,Lz)
      kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(0.5+Lx/Lmin)
      nymodes = nmodes*int(0.5+Ly/Lmin)
      nzmodes = nmodes*int(0.5+Lz/Lmin)

      xstep = int(Lx/Lmin+0.5)
      ystep = int(Ly/Lmin+0.5)
      zstep = int(Lz/Lmin+0.5)

      if (forcing_twice_wavelength==1) then
         HLx_inv = two/Lx
         HLy_inv = two/Ly
         HLz_inv = two/Lz
      else
         HLx_inv = 1.0d0/Lx
         HLy_inv = 1.0d0/Ly
         HLz_inv = 1.0d0/Lz
      endif

      do k = klo, khi
         z = xlo(3) + hz*(float(k-klo) + half)
         do j = jlo, jhi
            y = xlo(2) + hy*(float(j-jlo) + half)
            do i = ilo, ihi
               x = xlo(1) + hx*(float(i-ilo) + half)
               f3 = zero
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dfloat(kz)
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dfloat(ky)
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f3 = f3 + xT * (  FAY(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) &
                                              - FAX(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT * FAZ(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             cos(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               do kz = 1, zstep - 1
                  kzd = dfloat(kz)
                  do ky = mode_start, nymodes
                     kyd = dfloat(ky)
                     do kx = mode_start, nxmodes
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa<=kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force==1) then
                              f3 = f3 + xT * (  FAY(kx,ky,kz) * twicePi * (kxd*HLx_inv) * &
                                                cos(twicePi*kxd*x*HLx_inv+FPYX(kx,ky,kz)) * &
                                                sin(twicePi*kyd*y*HLy_inv+FPYY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPYZ(kx,ky,kz)) &
                                              - FAX(kx,ky,kz) * twicePi * (kyd*HLy_inv) * &
                                                sin(twicePi*kxd*x*HLx_inv+FPXX(kx,ky,kz)) * &
                                                cos(twicePi*kyd*y*HLy_inv+FPXY(kx,ky,kz)) * &
                                                sin(twicePi*kzd*z*HLz_inv+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT * FAZ(kx,ky,kz) * sin(twicePi*kxd*x*HLx_inv+FPX(kx,ky,kz)) * &
                                                             sin(twicePi*kyd*y*HLy_inv+FPY(kx,ky,kz)) * &
                                                             cos(twicePi*kzd*z*HLz_inv+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing==1) then
                  e(i,j,k,1) = f3*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f3
               endif
            enddo
         enddo
      enddo
#endif

   end subroutine FORT_DERFORCEZ

end module derive_PLM_nd
