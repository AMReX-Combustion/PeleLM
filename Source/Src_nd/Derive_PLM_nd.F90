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

  public :: dermgvort, dermgdivu, derdvrho, dermprho, &
            deravgpres, dergrdpx, dergrdpy, dergrdpz, &
            drhomry, dsrhoydot, dermassfrac, drhort, &
            dermolefrac, derconcentration, dertransportcoeff, dermolweight, &
            dhrr, dermixanddiss, dcma, init_mixture_fraction

  REAL_T, dimension(NUM_SPECIES,NUM_ELEMENTS) :: coeff_mix
  REAL_T, dimension(NUM_ELEMENTS) :: beta_mix
  REAL_T, dimension(NUM_SPECIES)  :: fact
  REAL_T :: Zfu, Zox
  REAL_T :: Zstoic = -1.0d0
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
!  Compute C / rho
!=========================================================

   subroutine derdvrho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="derdvrho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)/dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine derdvrho

!=========================================================
!  Compute rho * C 
!=========================================================

   subroutine dermprho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dermprho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)*dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine dermprho

!=========================================================
!  Compute cell-centered pressure as average of the 
!  surrounding nodal values 
!=========================================================

   subroutine deravgpres (e,   e_lo, e_hi, nv, &
                          dat, d_lo, d_hi, ncomp, &
                          lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                          level, grid_no) &
                          bind(C, name="deravgpres")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      integer :: i, j, k

      factor = 0.5d0
#if (AMREX_SPACEDIM >= 2 )
      factor = 0.25d0
#if (AMREX_SPACEDIM == 3 )
      factor = 0.125d0
#endif
#endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) =  factor * (  dat(i+1,j,k,1)     + dat(i,j,k,1)  &
#if (AMREX_SPACEDIM >= 2 )
                                       + dat(i+1,j+1,k,1)   + dat(i,j+1,k,1) &
#if (AMREX_SPACEDIM == 3 )
                                       + dat(i+1,j,k+1,1)   + dat(i,j,k+1,1)  &
                                       + dat(i+1,j+1,k+1,1) + dat(i,j+1,k+1,1) &
#endif
#endif
                                      )
            end do
         end do
      end do

   end subroutine deravgpres

!=========================================================
!  Compute node centered pressure gradient in direction dir
!=========================================================

   subroutine gradp_dir (p, p_lo, p_hi, &
                         gp, g_lo, g_hi, &
                         lo, hi, dir, dx) &
                         bind(C, name="gradp_dir")

      implicit none

! In/Out
      integer :: lo(3),  hi(3)
      integer :: p_lo(3), p_hi(3)
      integer :: g_lo(3), g_hi(3)
      integer :: dir
      REAL_T  :: dx
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p
      REAL_T, dimension(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3)) :: gp

! Local
      integer :: i, j, k
      REAL_T  :: d

      d = one/dx
#if (AMREX_SPACEDIM >= 2 )
      d = half/dx
#if (AMREX_SPACEDIM == 3 )
      d = fourth/dx
#endif
#endif

!
!     ::::: compute gradient on interior
!
      if (dir == 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 1 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k))
#elif (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k)+p(i+1,j+1,k)-p(i,j+1,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  ) + &
                                  p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
#endif 

               end do
            end do
         end do
      else if (dir == 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i,j+1,k)-p(i,j,k)+p(i+1,j+1,k)-p(i+1,j,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  ) + &
                                  p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
#endif
               end do
            end do
         end do
      else if (dir == 2) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k) + &
                                  p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
               end do
            end do
         end do
      else
         call amrex_abort("gradp_dir: invalid dir = ")
      end if

   end subroutine gradp_dir

!=========================================================
!  Compute node centered pressure gradient in X-dir
!=========================================================

   subroutine dergrdpx (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpx")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 0, delta(1))

   end subroutine dergrdpx

!=========================================================
!  Compute node centered pressure gradient in Y-dir
!=========================================================

   subroutine dergrdpy (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpy")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

#if (AMREX_SPACEDIM < 2 )
      call amrex_abort("dergrdpy called but AMREX_SPACEDIM<2 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 1, delta(2))

   end subroutine dergrdpy

!=========================================================
!  Compute node centered pressure gradient in Z-dir
!=========================================================

   subroutine dergrdpz (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpz")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no


#if (AMREX_SPACEDIM < 3 )
      call amrex_abort("dergrdpz called but AMREX_SPACEDIM<3 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 2, delta(3))

   end subroutine dergrdpz

!=========================================================
!  Compute rho - Sum_n(rhoY_n)
!=========================================================

   subroutine drhomry (e,   e_lo, e_hi, nv, &
                       dat, d_lo, d_hi, ncomp, &
                       lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                       level, grid_no) &
                       bind(C, name="drhomry")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: nxlo, nxhi, nylo, nyhi, nzlo,nzhi

      integer :: i, j, k, n

      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         CALL amrex_abort("drhomry: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1)
            enddo
         enddo
      enddo

! TODO: weird loop nesting ?
      do n=2,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)-dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

   end subroutine drhomry

!=========================================================
!  Compute sum of rhoY_dot
!=========================================================

   subroutine dsrhoydot (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dsrhoydot")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: nxlo, nxhi, nylo, nyhi, nzlo,nzhi

      integer :: i, j, k, n

      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         CALL amrex_abort("dsrhoydot: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = 0.0d0
            enddo
         enddo
      enddo

! TODO: weird loop nesting ?
      do n=1,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)+dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

   end subroutine dsrhoydot

!=========================================================
!  Compute rho*R*T
!=========================================================

   subroutine drhort (e,   e_lo, e_hi, nv, &
                      dat, d_lo, d_hi, ncomp, &
                      lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                      level, grid_no) &
                      bind(C, name="drhort")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: fS, rho, T
      integer :: nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T, dimension(NUM_SPECIES) :: Yt
      REAL_T  :: rhoinv, rhoCGS, Press

      integer :: i, j, k, n

      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         CALL amrex_abort("drhort: outside domain")
      endif

!
!     Set pointers into state (these must agree with setup for this derived quant).
!
      rho = 1
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               do n=1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1) * rhoinv
               end do
               rhoCGS = dat(i,j,k,rho) * 1.0d-3 ! kg/m^3 -> g/cm^3
               CALL CKPY(rhoCGS, dat(i,j,k,T), Yt, Press)
               e(i,j,k,1) = Press * 1.0d-1
            end do
         end do
      end do

   end subroutine drhort

!=========================================================
!  Compute species mass fractions Y_n
!=========================================================

   subroutine dermassfrac (e,   e_lo, e_hi, nv, &
                           dat, d_lo, d_hi, ncomp, &
                           lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                           level, grid_no) &
                           bind(C, name="dermassfrac")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: fS, rho
      REAL_T  :: rhoinv

      integer :: i, j, k, n

      rho = 1
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               do n = 1,NUM_SPECIES
                  e(i,j,k,n) = dat(i,j,k,fS+n-1) * rhoinv
               enddo
            enddo
         enddo
      enddo

   end subroutine dermassfrac

!=========================================================
!  Compute species mole fractions X_n
!=========================================================

   subroutine dermolefrac (e,   e_lo, e_hi, nv, &
                           dat, d_lo, d_hi, ncomp, &
                           lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                           level, grid_no) &
                           bind(C, name="dermolefrac")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: fS, rho
      REAL_T, dimension(NUM_SPECIES) :: Yt, Xt
      REAL_T  :: rhoinv

      integer :: i, j, k, n

      rho = 1
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               do n = 1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1) * rhoinv
               end do
               CALL CKYTX(Yt,Xt)
               do n = 1,NUM_SPECIES
                  e(i,j,k,n) = Xt(n)
               end do
            enddo
         enddo
      enddo

   end subroutine dermolefrac

!=========================================================
!  Compute species concentrations C_n
!=========================================================

   subroutine derconcentration (e,   e_lo, e_hi, nv, &
                                dat, d_lo, d_hi, ncomp, &
                                lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                                level, grid_no) &
                                bind(C, name="derconcentration")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: fS, rho, T
      REAL_T, dimension(NUM_SPECIES) :: Yt, Ct
      REAL_T  :: rhoinv, rhoCGS

      integer :: i, j, k, n

      rho = 1
      T = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               rhoCGS = dat(i,j,k,rho) * 1.0d-3 ! kg/m^3 -> g/cm^3
               do n = 1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1) * rhoinv
               end do
               CALL CKYTCR(rhoCGS,dat(i,j,k,T),Yt,Ct)
               do n = 1,NUM_SPECIES
                  e(i,j,k,n) = Ct(n) * 1.0d6 ! cm^(-3) -> m^(-3)
               end do
            enddo
         enddo
      enddo

   end subroutine derconcentration

!=========================================================
!  Extract species mass rhoY_n
!=========================================================

   subroutine derRhoY (e,   e_lo, e_hi, nv, &
                       dat, d_lo, d_hi, ncomp, &
                       lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                       level, grid_no) &
                       bind(C, name="derRhoY")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: fS

      integer :: i, j, k, n

      fS  = 1

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,NUM_SPECIES
                  e(i,j,k,n) = dat(i,j,k,fS+n-1)
               enddo
            enddo
         enddo
      enddo

   end subroutine derRhoY

!=========================================================
!  Compute transport coefficients: D_n, \lambda, \mu
!=========================================================

   subroutine dertransportcoeff (e,   e_lo, e_hi, nv, &
                                 dat, d_lo, d_hi, ncomp, &
                                 lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                                 level, grid_no) &
                                 bind(C, name="dertransportcoeff")

      use transport_module, only : get_transport_coeffs_F
      use PeleLM_nd,        only : vel_visc
      use mod_Fvar_def,     only : Pr, Sc, LeEQ1, thickFac

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T, dimension(NUM_SPECIES) :: Yt, D, invmwt
      REAL_T, dimension(1)        :: rho_dummy, MU, XI, LAM, Wavg
      REAL_T                      :: rhoinv, cpmix, Tfac, Yfac
      REAL_T, dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) :: mu_le1
      integer :: fS, rho, T
      integer :: lo_chem(3), hi_chem(3)

      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      integer :: i, j, k, n

      rho = 1
      T   = 2
      fS  = 3

      call get_imw(invmwt)
      if (.not. LeEQ1) then

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoinv = 1.0d0/dat(i,j,k,rho)
               do n = 1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1)*rhoinv
               enddo
               rho_dummy(1) = dat(i,j,k,rho) * 1.d-3

                  CALL CKMMWY(Yt,Wavg(1))

                  CALL get_transport_coeffs_F(lo_chem, hi_chem, &
                                            Yt,           lo_chem,hi_chem,  &
                                            dat(i,j,k,T), lo_chem,hi_chem,  &
                                            rho_dummy(1), lo_chem,hi_chem,  &
                                            D,            lo_chem,hi_chem,  &
                                            MU(1),        lo_chem,hi_chem,  &
                                            XI(1),        lo_chem,hi_chem,  &
                                            LAM(1),       lo_chem,hi_chem)

                  do n = 1,NUM_SPECIES
                     e(i,j,k,n) = Wavg(1) * invmwt(n) * D(n) * 0.1d0
                  enddo

                  e(i,j,k,NUM_SPECIES+1) = LAM(1) * 1.0d-05
                  e(i,j,k,NUM_SPECIES+2) = MU(1)  * 0.1d0

               enddo
            enddo
         enddo

      else

         Tfac = thickFac / Pr
         Yfac = thickFac / Sc

         call vel_visc(lo, hi,&
                 dat(lo(1),lo(2),lo(3),T), d_lo, d_hi, & 
                 dat(lo(1),lo(2),lo(3),fS), d_lo, d_hi, &
                 mu_le1, lo, hi)

         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  rhoinv = 1.0d0/dat(i,j,k,rho)
                  do n=1,NUM_SPECIES
                     Yt(n) = dat(i,j,k,fS+n-1)*rhoinv
                  end do
                  CALL CKCPBS(dat(i,j,k,T),Yt,cpmix)
                  cpmix = cpmix * 1.0d-4 
                  do n = 1,NUM_SPECIES
                     e(i,j,k,n) = mu_le1(i,j,k) * Yfac
                  enddo
                  e(i,j,k,NUM_SPECIES+1)    = mu_le1(i,j,k)*cpmix*Tfac
                  e(i,j,k,NUM_SPECIES+2)    = mu_le1(i,j,k)
               enddo
            enddo
         enddo


      endif

   end subroutine dertransportcoeff

!=========================================================
!  Compute the mixture mean molecular weight
!=========================================================

   subroutine dermolweight (e,   e_lo, e_hi, nv, &
                            dat, d_lo, d_hi, ncomp, &
                            lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                            level, grid_no) &
                            bind(C, name="dermolweight")

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
      REAL_T, dimension(NUM_SPECIES) :: Yt
      integer :: fS,rho

      integer :: i, j, k, n

      rho = 1
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho)
               enddo
               CALL CKMMWY(Yt,e(i,j,k,1))
            enddo
         enddo
      enddo

   end subroutine dermolweight

!=========================================================
!  Compute the mixture mean heat capacity at cst pressure
!=========================================================

   subroutine dercpmix (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dercpmix")

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
      REAL_T, dimension(NUM_SPECIES) :: Yt
      integer :: fS,rho, T

      integer :: i, j, k, n

      rho = 1
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho)
               enddo
               CALL CKCPBS(dat(i,j,k,T),Yt,e(i,j,k,1))
               e(i,j,k,1) = e(i,j,k,1) * 1.0d-4 ! CGS -> MKS
            enddo
         enddo
      enddo

   end subroutine dercpmix

!=========================================================
!  Init Bilger's element based mixture fraction
!=========================================================

   subroutine init_mixture_fraction(Yfu, Yox) bind(C,name='init_mixture_fraction')

      use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor
      use fuego_chemistry, only : elem_names, spec_names

      implicit none

! In/Out
      REAL_T, intent(in), dimension(NUM_SPECIES) :: Yfu, Yox

! Local
      REAL_T, dimension(NUM_SPECIES)  :: WtS
      REAL_T, dimension(NUM_ELEMENTS) :: WtE
      integer, dimension(NUM_ELEMENTS,NUM_SPECIES) :: ELTinSp
      integer :: i, j
      logical :: is_ioproc
      REAL_T, parameter :: tol = ten*tiny(zero)

      is_ioproc = amrex_pd_ioprocessor()

      ! Print stream composition
      if (is_ioproc) then
        write(6,'(2x,a)') 'Initialise mixture fraction'
        write(6,'(4x,a)') 'Fuel-stream mass fractions:'
        do i = 1, nspecies
          if (Yfu(i) .gt. 1e-14) then
            write(6,'(4x,a22,1x,f12.7)') adjustl(spec_names(i)), Yfu(i)
          endif
        enddo
        write(6,'(4x,a)') 'Oxidizer-stream mass fractions:'
        do i = 1, nspecies
          if (Yox(i) .gt. 1e-14) then
            write(6,'(4x,a22,1x,f12.7)') adjustl(spec_names(i)), Yox(i)
          endif
        enddo
      endif

      ! Sanity checks
      if      (abs(sum(Yfu) - one) .gt. tol) then
        call amrex_abort('sum(Yfu) != 1')
      else if (abs(sum(Yox) - one) .gt. tol) then
        call amrex_abort('sum(Yox) != 1')
      endif

      CALL ckwt(WtS)
      CALL ckawt(WtE)
      CALL ckncf(NUM_ELEMENTS,ELTinSP)

      do i = 1, NUM_SPECIES
         do j = 1, NUM_ELEMENTS
            coeff_mix(i,j) = ELTinSP(j,i)*WtE(j)/WtS(i)
         enddo
      enddo

      ! Bilger coeffs
      do i = 1,NUM_ELEMENTS
         beta_mix(i) = 0
         if(elem_names(i).eq.'C') beta_mix(i) = 2.0/WtE(i)
         if(elem_names(i).eq.'H') beta_mix(i) = 1.0/(2.0*WtE(i))
         if(elem_names(i).eq.'O') beta_mix(i) = -1.0/WtE(i)
      enddo

      Zfu = 0.0d0
      Zox = 0.0d0
      do i=1,NUM_SPECIES
         fact(i) = 0.0d0
         do j=1,NUM_ELEMENTS
            fact(i) = fact(i) + beta_mix(j)*coeff_mix(i,j)
         enddo
         Zfu = Zfu+fact(i)*Yfu(i)
         Zox = Zox+fact(i)*Yox(i)
      enddo

      Zstoic = (zero - Zox)/(Zfu - Zox)

      if (is_ioproc) then
        ! Print stream composition
        write(6,'(4x,a)') 'Stoichiometric mixture fraction:'
        write(6,'(8x,a,1x,f12.7)') 'Zstoic = ', Zstoic
      endif

      init_mixture = .TRUE.

   end subroutine init_mixture_fraction

!=========================================================
!  Compute Bilger's element based mixture fraction
!=========================================================

   subroutine dermixfrac (e,   e_lo, e_hi, nv, &
                          dat, d_lo, d_hi, ncomp, &
                          lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                          level, grid_no) &
                          bind(C, name="dermixfrac")

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
      REAL_T :: rhoinv
      integer :: fS,rho

      integer :: i, j, k, n

      rho = 1
      fS  = 2

      if (.not.init_mixture) call amrex_abort("mixture fraction not initialized")

      do k=lo(3), hi(3)
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               e(i,j,k,1) = 0.0d0
               do n=1,NUM_SPECIES
                  e(i,j,k,1) = e(i,j,k,1) + dat(i,j,k,fs+n-1) * rhoinv * fact(n)
               enddo
               e(i,j,k,1) = ( e(i,j,k,1) - Zox ) / ( Zfu - Zox )
            enddo
         enddo
      enddo

   end subroutine dermixfrac

!=========================================================
!  Compute the heat release rate -Sum_n(rhoYn_dot * Hn(T))
!=========================================================

   subroutine dhrr (e,   e_lo, e_hi, nv, &
                    dat, d_lo, d_hi, ncomp, &
                    lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                    level, grid_no) &
                    bind(C, name="dhrr")

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
      integer :: T, rhoYd
      REAL_T, dimension(NUM_SPECIES) :: Ht
      integer :: i, j, k, n

      T = 1
      rhoYd = 2

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = 0.0d0
               CALL CKHMS(dat(i,j,k,T),Ht)
               do n = 1, NUM_SPECIES
                  Ht(n) = Ht(n) * 1.0d-4     ! erg/g -> J/kg
                  e(i,j,k,1) = e(i,j,k,1) - dat(i,j,k,rhoYd+n-1)*Ht(n)
               enddo
            enddo
         enddo
      enddo

   end subroutine dhrr

!=========================================================
!  Compute both mixt. fraction and scalar diss. rate
!=========================================================

   subroutine dermixanddiss (e,   e_lo, e_hi, nv, &
                             dat, d_lo, d_hi, ncomp, &
                             lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                             level, grid_no) &
                             bind(C, name="dermixanddiss")

      use transport_module, only : get_transport_coeffs_F

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
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: mixfrac
      REAL_T, dimension(NUM_SPECIES) :: Yt, D
      REAL_T, dimension(3) :: grad
      REAL_T, dimension(1) :: rhoCGS, MU, XI, LAM
      REAL_T :: cpmix, rhoinv
      integer :: rho, T, fS

      integer :: lo_chem(3), hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      integer :: i, j, k, n

      rho = 1
      T   = 2
      fS  = 3

      if (.not.init_mixture) call amrex_abort("mixture fraction not initialized")

      grad(1) = 0.0d0; grad(2) = 0.0d0; grad(3) = 0.0d0

!     Grown box will be given for mixture fraction because derivative
!     needs to be calculated
      do k = d_lo(3), d_hi(3)
         do j = d_lo(2), d_hi(2)
            do i = d_lo(1), d_hi(1)
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               mixfrac(i,j,k) = 0.0d0
               do n = 1,NUM_SPECIES
                  mixfrac(i,j,k) = mixfrac(i,j,k) + dat(i,j,k,fS+n-1) * rhoinv * fact(n)
               enddo
               mixfrac(i,j,k) = ( mixfrac(i,j,k) - Zox ) / ( Zfu - Zox )
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ! Fill mixture fraction array
               e(i,j,k,1) = mixfrac(i,j,k)

               ! Compute scalar dissipation rate
               ! cpmix 
               rhoinv = 1.0d0 / dat(i,j,k,rho)
               do n = 1, NUM_SPECIES
                  Yt(n) = dat(i,j,k,fS+n-1) * rhoinv
               enddo
               CALL CKCPBS(dat(i,j,k,T), Yt, cpmix)
               cpmix = cpmix * 1.0d-4 ! erg/g.K -> J/kg.K

               ! lambda mix 
               rhoCGS(1) = dat(i,j,k,rho) * 1.0d-3 ! kg/m^3 -> g/cm^3
               CALL get_transport_coeffs_F(lo_chem, hi_chem, &
                                         Yt,           lo_chem,hi_chem, &
                                         dat(i,j,k,T), lo_chem,hi_chem, &
                                         rhoCGS(1),    lo_chem,hi_chem, &
                                         D,            lo_chem,hi_chem, &
                                         MU(1),        lo_chem,hi_chem, &
                                         XI(1),        lo_chem,hi_chem, &
                                         LAM(1),       lo_chem,hi_chem)
               LAM(1) = LAM(1) * 1.0d-05 ! erg/(s.cm.K) -> J/(s.m.K)

               ! grad mixt. fraction
               grad(1) = 0.5*(mixfrac(i+1,j,k)-mixfrac(i-1,j,k))/delta(1)
               grad(2) = 0.5*(mixfrac(i,j+1,k)-mixfrac(i,j-1,k))/delta(2)
               if ( delta(3) > 0.0d0 ) grad(3) = 0.5*(mixfrac(i,j,k+1)-mixfrac(i,j,k-1))/delta(3)

               ! Gather
               e(i,j,k,2) = grad(1)**2.0d0 + grad(2)**2.0d0 + grad(3)**2.0d0
               e(i,j,k,2) = 2.0d0 * e(i,j,k,2) * LAM(1) * rhoinv / cpmix
            enddo
         enddo
      enddo

   end subroutine dermixanddiss

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

      CALL  dermixfrac(e(:,:,:,1),  e_lo, e_hi, 1, &
                       dat(:,:,:,rho:NUM_SPECIES+1), d_lo, d_hi, NUM_SPECIES+1, &
                       lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                       level, grid_no)

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
