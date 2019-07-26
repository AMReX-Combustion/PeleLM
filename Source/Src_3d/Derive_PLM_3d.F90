#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <Prob_F.H>
#include <AMReX_ArrayLim.H>

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

module derive_PLM_3D

  use mod_Fvar_def, only : dim
  
  implicit none

  private
 
  public :: derdvrho, dermprho, dermgvort, dermgdivu, &
            deravgpres, dergrdpx, dergrdpy, dergrdpz, &
            drhomry, dsrhoydot, drhort, dermassfrac, & 
            dermolefrac, derconcentration, dertransportcoeff, &
            dermolweight

contains


 subroutine derdvrho (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                           lo,hi,domlo,domhi,delta,xlo,time,dt, &
                           bc,level, grid_no) bind(C, name="derdvrho")
      implicit none
!c
!c ::: This routine will derive C/RHO
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)/dat(i,j,k,1)
            end do
         end do
      end do

      end subroutine derdvrho

      subroutine dermprho (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt, &
                               bc,level, grid_no) bind(C, name="dermprho")
      implicit none
!c
!c ::: This routine will derive RHO*C
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)*dat(i,j,k,1)
            end do
         end do
      end do

      end subroutine dermprho
      
      

   subroutine dermgvort (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                bc,level,grid_no) bind(C,name="dermgvort")
      implicit none
!c
!c ::: This routine will derive magnitude of vorticity from
!c ::: the velocity field
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    uy, uz, vx, vz, wx, wy, dx, dy, dz
      REAL_T    uycen, uzcen, uylo, uyhi, uzlo, uzhi
      REAL_T    vxcen, vzcen, vxlo, vxhi, vzlo, vzhi
      REAL_T    wxcen, wycen, wxlo, wxhi, wylo, wyhi
      REAL_T    vorfun

      logical   fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z
!c
!c     ::::: some useful macro definitions
!c
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
!c
!c     ::::: statement functions that implement stencil
!c
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               uz = uzcen(i,j,k)
               vx = vxcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
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
!c
!c     First do all the faces
!c
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

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
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
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
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
!c
!c     Next do all the edges
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if
!c
!c     Finally do all the corners
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
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


  subroutine dermgdivu (divu,DIMS(divu),nv,dat,DIMS(dat),ncomp, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                bc,level,grid_no) bind(C, name="dermgdivu")
      
      use bc_fill_3d_module, only : xvel_fill, yvel_fill, zvel_fill
      
      implicit none
!c
!c ::: This routine will derive magnitude of the divergence of velocity
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(divu)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     divu(DIMV(divu),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    ux, vy, wz, dx, dy, dz
      REAL_T    uxcen, uxlo, uxhi
      REAL_T    vycen, vylo, vyhi
      REAL_T    wzcen, wzlo, wzhi
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define WLOZ bc(3,1,2)
#     define WHIZ bc(3,2,2)
!c
!c     ::::: statement functions that implement stencil
!c
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k) = (eight*U(i,j,k)-six*U(i+1,j,k)+U(i+2,j,k))/(three*dx)
      uxhi(i,j,k) = (eight*U(i,j,k)-six*U(i-1,j,k)+U(i-2,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k) = (eight*V(i,j,k)-six*V(i,j+1,k)+V(i,j+2,k))/(three*dy)
      vyhi(i,j,k) = (eight*V(i,j,k)-six*V(i,j-1,k)+V(i,j-2,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k+1)+W(i,j,k+2))/(three*dz)
      wzhi(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k-1)+W(i,j,k-2))/(three*dz)

      call xvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,1))
      call yvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),2),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,2))
      call zvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),3),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,3))

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
!c
!c     :: at physical bndries where an edge value is prescribed,
!c     :: set the value in the outside cell so that a central
!c     :: difference formula is equivalent to the higher order
!c     :: one sided formula
!c
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (ULOX.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = two*U(i-1,j,k) - U(i,j,k)
               end do
            end do
         else if (ULOX.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = uxlo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (UHIX.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = two*U(i+1,j,k) - U(i,j,k)
               end do
            end do
         else if (UHIX.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = uxhi(i,j,k)
               end do
            end do
	 end if
      end if
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
	 if (VLOY.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = two*V(i,j-1,k) - V(i,j,k)
               end do
            end do
         else if (VLOY.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = vylo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(2) .eq. domhi(2)) then
         j = hi(2)
	 if (VHIY.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = two*V(i,j+1,k) - V(i,j,k)
               end do
            end do
	 else if (VHIY.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = vyhi(i,j,k)
               end do
            end do
	 end if
      end if
      if (lo(3) .eq. domlo(3)) then
         k = lo(3)
	 if (WLOZ.eq.EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = two*W(i,j,k-1) - W(i,j,k)
               end do
            end do
	 else if (WLOZ.eq.HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = wzlo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(3) .eq. domhi(3)) then
         k = hi(3)
	 if (WHIZ.eq.EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = two*W(i,j,k+1) - W(i,j,k)
               end do
            end do
	 else if (WHIZ.eq.HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = wzhi(i,j,k)
               end do
            end do
	 end if
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               vy = vycen(i,j,k)
               wz = wzcen(i,j,k)
               divu(i,j,k,1) = ux + vy + wz
            end do
         end do
      end do

!c
!c we overwrote the ghost cells above, so set them back below
!c
      call xvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,1))
      call yvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,2))
      call zvel_fill(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat), &
                        domlo,domhi,delta,xlo,time,bc(1,1,3))

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

      subroutine gradp_dir  ( &
          p,DIMS(p), &
          gp,DIMS(gp), &
          lo,hi,dir,dx) &
          bind(C, name="gradp_dir")

      implicit none
!c
!c     compute a node centered pressure gradient in direction (dir)
!c
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer     lo(dim),  hi(dim)
      integer    dir
      REAL_T     dx
      REAL_T   p(DIMV(p))
      REAL_T  gp(DIMV(gp))

      integer    i,j,k
      REAL_T     d

      d = fourth/dx
!c
!c     ::::: compute gradient on interior
!c
      if (dir .eq. 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( &
                      p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+ &
                      p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
               end do
            end do
         end do
      else if (dir .eq. 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( &
                      p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+ &
                      p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
               end do
            end do
         end do
      else if (dir .eq. 2) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( &
                      p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+ &
                      p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
               end do
            end do
         end do
      else
	 call bl_abort("FORT_GRADP_DIR: invalid dir = ")
      end if

      end subroutine gradp_dir

      subroutine dergrdpx (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt, &
                               bc,level,grid_no) bind(C, name="dergrdpx")
      implicit none
!c
!c     This routine computes pressure gradient in x direciton
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     grdpx(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
!c
      call gradp_dir ( &
          dat,DIMS(dat),grdpx,DIMS(gp), &
          lo,hi,0,delta(1))

      end subroutine dergrdpx

      subroutine dergrdpy (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt,&
                               bc,level,grid_no) bind(C, name="dergrdpy")
      implicit none
!c
!c     This routine computes pressure gradient in Y direciton
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     grdpy(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
!c
      call gradp_dir ( &
          dat,DIMS(dat),grdpy,DIMS(gp), &
          lo,hi,1,delta(2))

      end subroutine dergrdpy

      subroutine dergrdpz (grdpz,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt, &
                               bc,level,grid_no) bind(C, name="dergrdpz")
      implicit none
!c
!c     This routine computes pressure gradient in Z direciton
!c
      integer    lo(dim), hi(dim)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim)
      REAL_T     time, dt
      REAL_T     grdpz(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
!c
      call gradp_dir ( &
          dat,DIMS(dat),grdpz,DIMS(gp), &
          lo,hi,2,delta(3))

      end subroutine dergrdpz


      subroutine deravgpres (avgpres,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                                 lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                 bc,level,grid_no) bind(C, name="deravgpres")
      implicit none
!c
!c     This routine computes cell-centered pressure as average of the eight
!c       surrounding nodal values.
!c
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      REAL_T  avgpres(DIMV(gp))
      REAL_T  dat(DIMV(dat))
      integer nv, ncomp
      integer lo(dim), hi(dim)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim)
      REAL_T  xlo(dim)
      REAL_T  time, dt
      integer bc(dim,2,ncomp)
      integer level
      integer grid_no

      integer i,j,k

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            avgpres(i,j,k) = eighth*( &
                          dat(i+1,j,k)     + dat(i,j,k)  &
                        + dat(i+1,j+1,k)   + dat(i,j+1,k) &
                        + dat(i+1,j,k+1)   + dat(i,j,k+1)  &
                        + dat(i+1,j+1,k+1) + dat(i,j+1,k+1) )
          end do
        end do
      end do

      end subroutine deravgpres

  
  subroutine drhomry (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                      lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                      level,grid_no) &
                      bind(C, name="drhomry")


    implicit none
!
! ::: This routine will computes rho - sum (rho*Y)
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHOMINUSSUMRHOY: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1)
            enddo
         enddo
      enddo

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

  subroutine dsrhoydot (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                        lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                        level,grid_no) &
                        bind(C, name="dsrhoydot")

    implicit none
!
! ::: This routine will computes sum (rhoYdot or Ydot)
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo, nxhi, nylo, nyhi, nzlo,nzhi
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERSUMYDOT: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = 0.0d0
            enddo
         enddo
      enddo
      
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

  subroutine drhort (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                     lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                     level,grid_no) &
                     bind(C, name="drhort")

    use network,        only : nspec
    use PeleLM_3D, only: pphys_PfromRTY
    
    implicit none

!     
! ::: This routine will derive rho*R*T
!
      
      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i, j, k, n, rho, T, fS
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      REAL_T     Yt(nspec)
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHORT: outside domain")
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
               do n=1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1) / dat(i,j,k,rho)
               end do
               call pphys_PfromRTY(lo_chem, hi_chem, &
                   e(i,j,k,1),     ARLIM(lo_chem),ARLIM(hi_chem), &
                   dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
                   dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
                   Yt,             ARLIM(lo_chem),ARLIM(hi_chem))
            end do
         end do
      end do
  end subroutine drhort

!=========================================================

  subroutine dermassfrac (x,DIMS(x),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                          level,grid_no) &
                          bind(C, name="dermassfrac")

    use network,        only : nspec

    implicit none

      integer    lo(dim), hi(dim)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      integer fS,rho
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  x(i,j,k,n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho)
               enddo
            enddo
         enddo
      enddo

  end subroutine dermassfrac

!=========================================================

  subroutine derRhoY (x,DIMS(x),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                          level,grid_no) &
                          bind(C, name="derRhoY")

    use network,        only : nspec

    implicit none

      integer    lo(dim), hi(dim)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      integer fS,rho
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      fS  = 1

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  x(i,j,k,n) = dat(i,j,k,fS+n-1)
               enddo
            enddo
         enddo
      enddo

  end subroutine derRhoY

!=========================================================

  subroutine dermolefrac (x,DIMS(x),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                          level,grid_no) &
                          bind(C, name="dermolefrac")

    use network,        only : nspec
    use PeleLM_3D, only : pphys_mass_to_mole

    implicit none
    
      integer    lo(dim), hi(dim)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(nspec),Xt(nspec)
      integer fS,rho
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
               call pphys_mass_to_mole(lo_chem, hi_chem, &
                          Yt, ARLIM(lo_chem),ARLIM(hi_chem), &
                          Xt, ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  x(i,j,k,n) = Xt(n)
               enddo
            enddo
         enddo
      enddo

  end subroutine dermolefrac

!=========================================================
  
  subroutine derconcentration (C,DIMS(C),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                               level,grid_no) &
                               bind(C, name="derconcentration")

    use network,        only : nspec
    use PeleLM_3D, only: pphys_massr_to_conc
                               
    implicit none

      integer    lo(dim), hi(dim)
      integer    DIMDEC(C)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     C(DIMV(C),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(nspec),Ct(nspec)
      integer fS,rho,T
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
               call pphys_massr_to_conc(lo_chem,hi_chem, &
                  Yt,           ARLIM(lo_chem),ARLIM(hi_chem), &
                  dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
                  dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
                  Ct,           ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)
               enddo
            enddo
         enddo
      enddo

  end subroutine derconcentration

  subroutine dertransportcoeff (C,DIMS(C),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                               level,grid_no) &
                               bind(C, name="dertransportcoeff")

    use network,        only : nspec
    use transport_module, only : get_transport_coeffs

    implicit none

    integer    lo(dim), hi(dim)
    integer    DIMDEC(C)
    integer    DIMDEC(dat)
    integer    domlo(dim), domhi(dim)
    integer    nv, ncomp
    integer    bc(dim,2,ncomp)
    REAL_T     delta(dim), xlo(dim), time, dt
    REAL_T     C(DIMV(C),nv)
    REAL_T     dat(DIMV(dat),ncomp)
    integer    level, grid_no

    integer i,j,k,n
    REAL_T Yt(nspec), rho_dummy(1), D(Nspec), MU(1), XI(1), LAM(1)
    integer fS,rho,T
    integer lo_chem(3),hi_chem(3)
    data lo_chem /1,1,1/
    data hi_chem /1,1,1/

    rho = 1
    T   = 2
    fS  = 3

    do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

              do n = 1,Nspec
                Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho)
              enddo
              rho_dummy(1) = dat(i,j,k,rho) * 1.d-3

              call get_transport_coeffs(lo_chem, hi_chem, &
                                        Yt,           lo_chem,hi_chem,  &
                                        dat(i,j,k,T),   lo_chem,hi_chem,  &
                                        rho_dummy(1), lo_chem,hi_chem,  &
                                        D,            lo_chem,hi_chem,  &
                                        MU(1),        lo_chem,hi_chem,  &
                                        XI(1),        lo_chem,hi_chem,  &
                                        LAM(1),      lo_chem,hi_chem)

              do n = 1,Nspec
                C(i,j,k,n) = D(n) * 0.1d0
              enddo

              C(i,j,k,Nspec+1) = LAM(1) * 1.0d-05
              C(i,j,k,Nspec+2) = MU(1)  * 0.1d0

            enddo
         enddo
      enddo

  end subroutine dertransportcoeff
  
!=========================================================

  subroutine dermolweight (x,DIMS(x),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                          level,grid_no) &
                          bind(C, name="dermolweight")

    use network,        only : nspec
    use fuego_chemistry

    implicit none

      integer    lo(dim), hi(dim)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      REAL_T Yt(nspec)
      integer    level, grid_no

      integer i,j,k,n
      integer fS,rho
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho)
               enddo
               
               CALL CKMMWY(Yt,x(i,j,k,1))
               
            enddo
         enddo
      enddo

  end subroutine dermolweight

!=========================================================

  subroutine FORT_DERFORCING (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                              level,grid_no) &
                              bind(C, name="FORT_DERFORCING")
     
      implicit none

!
! ::: This routine will computes the forcing term
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>

      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1, f2, f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      REAL_T  rho, u, v, w
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

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
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif
      
      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz.eq.1) then
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

      if (forcing_twice_wavelength.eq.1) then
         HLx = Lx/two
         HLy = Ly/two
         HLz = Lz/two
      else
         HLx = Lx
         HLy = Ly
         HLz = Lz
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
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
  
  subroutine FORT_DERFORCEX (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEX")
                             
      implicit none

!
! ::: This routine will computes the forcing term
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

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
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif
      
      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz.eq.1) then
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

      if (forcing_twice_wavelength.eq.1) then
         HLx = Lx/two
         HLy = Ly/two
         HLz = Lz/two
      else
         HLx = Lx
         HLy = Ly
         HLz = Lz
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
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

  subroutine FORT_DERFORCEY (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEY")
                             
      implicit none
!
! ::: This routine will computes the forcing term
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f2
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

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
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif
      
      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz.eq.1) then
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

      if (forcing_twice_wavelength.eq.1) then
         HLx = Lx/two
         HLy = Ly/two
         HLz = Lz/two
      else
         HLx = Lx
         HLy = Ly
         HLz = Lz
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
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

  subroutine FORT_DERFORCEZ (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEZ")
                             
      implicit none

!
! ::: This routine will computes the forcing term
!

      integer    lo(dim), hi(dim)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(dim), domhi(dim)
      integer    nv, ncomp
      integer    bc(dim,2,ncomp)
      REAL_T     delta(dim), xlo(dim), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

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
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
      if (time_offset.gt.(-half)) then
         force_time = time + time_offset
      else
         force_time = time
      endif
      
      Lx = domnhi(1)-domnlo(1)
      Ly = domnhi(2)-domnlo(2)
      Lz = domnhi(3)-domnlo(3)

      if (hack_lz.eq.1) then
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

      if (forcing_twice_wavelength.eq.1) then
         HLx = Lx/two
         HLy = Ly/two
         HLz = Lz/two
      else
         HLx = Lx
         HLy = Ly
         HLz = Lz
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                           if (div_free_force.eq.1) then
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
                  e(i,j,k,1) = f3*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f3
               endif
            enddo
         enddo
      enddo
#endif
  end subroutine FORT_DERFORCEZ

end module derive_PLM_3D
