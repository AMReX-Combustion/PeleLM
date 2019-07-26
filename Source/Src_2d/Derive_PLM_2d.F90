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

module derive_PLM_2D

  use mod_Fvar_def, only : dim
  
  implicit none

  private
 
  public :: derdvrho, dermprho, dermgvort, dermgdivu, deravgpres, dergrdpx, dergrdpy, &
            drhomry, dsrhoydot, drhort, dermassfrac, dermolefrac, derconcentration, &
            dertransportcoeff, dermolweight

contains
 
!
! ::: This routine will derive C/RHO
!

 subroutine derdvrho (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no)bind(C,name="derdvrho")

   implicit none

   integer    lo(2), hi(2)
   integer    DIMDEC(e)
   integer    DIMDEC(dat)
   integer    domlo(2), domhi(2)
   integer    nv, ncomp
   integer    bc(2,2,ncomp)
   REAL_T     delta(2), xlo(2), time, dt
   REAL_T     e(DIMV(e),nv)
   REAL_T     dat(DIMV(dat),ncomp)
   integer    level, grid_no
 
   integer    i,j
 
   do j = lo(2), hi(2)
     do i = lo(1), hi(1)
	     e(i,j,1) = dat(i,j,2)/dat(i,j,1)
	   end do
   end do

 end subroutine derdvrho
 
!
! ::: This routine will derive RHO*C
!
 
 subroutine dermprho (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="dermprho")
   implicit none

   integer    lo(2), hi(2)
   integer    DIMDEC(e)
   integer    DIMDEC(dat)
   integer    domlo(2), domhi(2)
   integer    nv, ncomp
   integer    bc(2,2,ncomp)
   REAL_T     delta(2), xlo(2), time, dt
   REAL_T     e(DIMV(e),nv)
   REAL_T     dat(DIMV(dat),ncomp)
   integer    level, grid_no

   integer    i,j

   do j = lo(2), hi(2)
     do i = lo(1), hi(1)
	     e(i,j,1) = dat(i,j,2)*dat(i,j,1)
	   end do
   end do

 end subroutine dermprho
 
! 
! ::: This routine will derive magnitude of vorticity from
! ::: the velocity field
!
    
 subroutine dermgvort (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp,&
                                lo,hi,domlo,domhi,delta,xlo,time,dt,&
                                bc,level,grid_no) bind(C,name="dermgvort")
   implicit none

   integer    lo(2), hi(2)
   integer    DIMDEC(vort)
   integer    DIMDEC(dat)
   integer    domlo(2), domhi(2)
   integer    nv, ncomp
   integer    bc(2,2,ncomp)
   REAL_T     delta(2), xlo(2), time, dt
   REAL_T     vort(DIMV(vort),nv)
   REAL_T     dat(DIMV(dat),ncomp)
   integer    level, grid_no

   integer   i,j
   REAL_T    vx, uy, dx, dy
   logical   fixlft, fixrgt, fixbot, fixtop
   REAL_T    vxcen, vxlft, vxrgt, uycen, uybot, uytop, vorfun
   
!
!     ::::: some useful macro definitions
!
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)
#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
!
!     ::::: statement functions that implement stencil
!
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)
      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)
      vorfun(vx,uy) = vx - uy

      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	         vx  = vxcen(i,j)
	         uy  = uycen(i,j)
	         vort(i,j,1) = vorfun(vx,uy)
	       end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. &
                (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. &
                (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. &
                (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. &
                (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
!
!     ::::: handle special bndry conditions at the left edge
!
      if (fixlft) then
        i = lo(1)
        do j = lo(2), hi(2)
	        vx  = vxlft(i,j)
	        uy  = uycen(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	      end do
      end if
!
!     ::::: handle special bndry conditions on the right
!
      if (fixrgt) then
        i = hi(1)
        do j = lo(2), hi(2)
	        vx  = vxrgt(i,j)
	        uy  = uycen(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
        end do
      end if
!
!     ::::: handle special bndry conditions on bottom
!
      if (fixbot) then
        j = lo(2)
        do i = lo(1), hi(1)
	        vx  = vxcen(i,j)
	        uy  = uybot(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	      end do
      end if
!
!     ::::: handle special bndry conditions on top
!
      if (fixtop) then
        j = hi(2)
        do i = lo(1), hi(1)
	        vx  = vxcen(i,j)
	        uy  = uytop(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	      end do
      end if
!
!     ::::: check corners
!
      if (fixlft .and. fixbot) then
        i = lo(1)
        j = lo(2)
        vx = vxlft(i,j)
	      uy = uybot(i,j)
	      vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixlft .and. fixtop) then
        i = lo(1)
	      j = hi(2)
        vx = vxlft(i,j)
        uy = uytop(i,j)
	      vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixrgt .and. fixtop) then
        i = hi(1)
	      j = hi(2)
	      vx = vxrgt(i,j)
        uy = uytop(i,j)
	      vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixrgt .and. fixbot) then
        i = hi(1)
	      j = lo(2)
	      vx = vxrgt(i,j)
	      uy = uybot(i,j)
	      vort(i,j,1) = vorfun(vx,uy)
      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

    end subroutine dermgvort   
 
!
! ::: This routine will derive magnitude of the divergence of velocity
!
 
   subroutine dermgdivu (divu,DIMS(divu),nv,dat,DIMS(dat),ncomp,&
                                lo,hi,domlo,domhi,delta,xlo,time,dt,&
                                bc,level,grid_no) bind(C,name="dermgdivu")
      
      use bc_fill_2d_module, only : xvel_fill, yvel_fill
      
      implicit none

      integer    lo(2), hi(2)
      integer    DIMDEC(divu)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     divu(DIMV(divu),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j
      REAL_T    ux, vy, dx, dy
      REAL_T    uxcen, uxlo, uxhi
      REAL_T    vycen, vylo, vyhi
!
!     ::::: some useful macro definitions
!
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)

#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlo(i,j) = (eight*U(i,j)-six*U(i+1,j)+U(i+2,j))/(three*dx)
      uxhi(i,j) = (eight*U(i,j)-six*U(i-1,j)+U(i-2,j))/(three*dx)

      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vylo(i,j) = (eight*V(i,j)-six*V(i,j+1)+V(i,j+2))/(three*dy)
      vyhi(i,j) = (eight*V(i,j)-six*V(i,j-1)+V(i,j-2))/(three*dy)

      call xvel_fill(dat(ARG_L1(dat),ARG_L2(dat),1),DIMS(dat),&
                        domlo,domhi,delta,xlo,time,bc(1,1,1))
      call yvel_fill(dat(ARG_L1(dat),ARG_L2(dat),2),DIMS(dat),&
                        domlo,domhi,delta,xlo,time,bc(1,1,2))

      dx = delta(1)
      dy = delta(2)
!
!     :: at physical bndries where an edge value is prescribed,
!     :: set the value in the outside cell so that a central
!     :: difference formula is equivalent to the higher order
!     :: one sided formula
!
!     boundaries handled by fill routines
!
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (ULOX.eq.EXT_DIR) then
           do j = lo(2), hi(2)
              U(i-1,j) = two*U(i-1,j)-U(i,j)
           end do
         else if (ULOX.eq.HOEXTRAP) then
           do j = lo(2), hi(2)
              U(i-1,j) = uxlo(i,j)
           end do
         end if
      end if
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (UHIX.eq.EXT_DIR) then
           do j = lo(2), hi(2)
              U(i+1,j) = two*U(i+1,j)-U(i,j)
           end do
         else if (UHIX.eq.HOEXTRAP) then
           do j = lo(2), hi(2)
              U(i+1,j) = uxhi(i,j)
           end do
         end if
      end if
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
         if (VLOY.eq.EXT_DIR) then
           do i = lo(1), hi(1)
              V(i,j-1) = two*V(i,j-1)-V(i,j)
           end do
         else if (VLOY.eq.HOEXTRAP) then
           do i = lo(1), hi(1)
              V(i,j-1) = vylo(i,j)
           end do
         end if
      end if
      if (hi(2) .eq. domhi(2)) then
         j = hi(2)
         if (VHIY.eq.EXT_DIR) then
           do i = lo(1), hi(1)
              V(i,j+1) = two*V(i,j+1)-V(i,j)
           end do
         else if (VHIY.eq.HOEXTRAP) then
           do i = lo(1), hi(1)
              V(i,j+1) = vyhi(i,j)
           end do
         end if
      end if

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux = uxcen(i,j)
            vy = vycen(i,j)
            divu(i,j,1) = ux + vy
         end do
      end do
!
! we overwrote the ghost cells above, so set them back below
!
      call xvel_fill(dat(ARG_L1(dat),ARG_L2(dat),1),DIMS(dat),&
                        domlo,domhi,delta,xlo,time,bc(1,1,1))
      call yvel_fill(dat(ARG_L1(dat),ARG_L2(dat),2),DIMS(dat),&
                        domlo,domhi,delta,xlo,time,bc(1,1,2))

#     undef U
#     undef V      
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY

    end subroutine dermgdivu
 
!
!     This routine computes cell-centered pressure as average of the four
!       surrounding nodal values.
!
 
   subroutine deravgpres (avgpres,DIMS(gp),nv,dat,DIMS(dat),ncomp,&
                                 lo,hi,domlo,domhi,delta,xlo,time,dt,&
                                 bc,level,grid_no) bind(C,name="deravgpres")
      implicit none

      integer DIMDEC(gp)
      integer DIMDEC(dat)
      REAL_T  avgpres(DIMV(gp))
      REAL_T  dat(DIMV(dat))
      integer nv, ncomp
      integer lo(2), hi(2)
      integer domlo(2), domhi(2)
      REAL_T  delta(2)
      REAL_T  xlo(2)
      REAL_T  time, dt
      integer bc(2,2,ncomp)
      integer level
      integer grid_no

      integer i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            avgpres(i,j) = fourth*( dat(i+1,j  )+dat(i,j  ) +&
                                   dat(i+1,j+1)+dat(i,j+1) )
         end do
      end do

    end subroutine deravgpres

!
!     compute a node centered pressure gradient in direction (dir)
!

    subroutine gradp_dir (p,DIMS(p), gp,DIMS(gp),&
           lo,hi,dir,dx)

      implicit none

      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer    lo(dim),  hi(dim)
      integer    dir
      REAL_T     dx
      REAL_T     p(DIMV(p))
      REAL_T     gp(DIMV(gp))
      integer    i,j
      REAL_T     d

#if (BL_PRVERSION == 9)
      d = half
#else
      d = half/dx
#endif

      if (dir .eq. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j) = d*(p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
            end do
         end do
      else if (dir .eq. 1) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j) = d*(p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
            end do
         end do
      else
         call bl_abort("FORT_GRADP_DIR: invalid dir")
      end if

    end subroutine gradp_dir

!
!     This routine computes pressure gradient in X direction
!
    
    subroutine dergrdpx (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,&
                               lo,hi,domlo,domhi,delta,xlo,time,dt,&
                               bc,level,grid_no) bind(C,name="dergrdpx")
      implicit none

      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpx(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no

      call gradp_dir (&
          dat,DIMS(dat),grdpx,DIMS(gp),&
          lo,hi,0,delta(1))

    end subroutine dergrdpx

!
!     This routine computes pressure gradient in Y direction
!
    
    subroutine dergrdpy (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp,&
                               lo,hi,domlo,domhi,delta,xlo,time,dt,&
                               bc,level,grid_no) bind(C,name="dergrdpy")
      implicit none

      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpy(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no

      call gradp_dir (&
          dat,DIMS(dat),grdpy,DIMS(gp),&
          lo,hi,1,delta(2))

    end subroutine dergrdpy
 
!
! ::: This routine will computes rho - sum (rho*Y)
!
 

  subroutine drhomry (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                      lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                      level,grid_no) &
                      bind(C, name="drhomry")
                                  
    implicit none

    integer    lo(2), hi(2)
    integer    DIMDEC(e)
    integer    DIMDEC(dat)
    integer    domlo(2), domhi(2)
    integer    nv, ncomp
    integer    bc(2,2,ncomp)
    REAL_T     delta(2), xlo(2), time, dt
    REAL_T     e(DIMV(e),nv)
    REAL_T     dat(DIMV(dat),ncomp)
    integer    level, grid_no

    integer    i,j,k
    integer    nlft, nrgt, nbot, ntop
      
!     ::::: lets punt if not in domain interior
    nlft = max(0,domlo(1)-lo(1))
    nrgt = max(0,hi(1)-domhi(1))
    nbot = max(0,domlo(2)-lo(2))
    ntop = max(0,hi(2)-domhi(2))

    if (nlft+nrgt+nbot+ntop .gt. 0) then
      call bl_abort("FORT_DERRHOMINUSSUMRHOY: outside domain")
    endif

    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        e(i,j,1) = dat(i,j,1)
      enddo
    enddo

    do k=2,ncomp
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          e(i,j,1) = e(i,j,1)-dat(i,j,k)
        enddo
       enddo
    enddo
    
    return

  end subroutine drhomry

!
! ::: This routine will computes sum (rhoYdot or Ydot)
!

  subroutine dsrhoydot (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                        lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                        level,grid_no) &
                        bind(C, name="dsrhoydot")

    implicit none

    integer    lo(2), hi(2)
    integer    DIMDEC(e)
    integer    DIMDEC(dat)
    integer    domlo(2), domhi(2)
    integer    nv, ncomp
    integer    bc(2,2,ncomp)
    REAL_T     delta(2), xlo(2), time, dt
    REAL_T     e(DIMV(e),nv)
    REAL_T     dat(DIMV(dat),ncomp)
    integer    level, grid_no

    integer    i,j,k
    integer    nlft, nrgt, nbot, ntop
      
!     ::::: lets punt if not in domain interior
    nlft = max(0,domlo(1)-lo(1))
    nrgt = max(0,hi(1)-domhi(1))
    nbot = max(0,domlo(2)-lo(2))
    ntop = max(0,hi(2)-domhi(2))

    if (nlft+nrgt+nbot+ntop .gt. 0) then
       call bl_abort("FORT_DERSUMRHOYDOT: outside domain")
    endif

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
         e(i,j,1) = 0.0
       enddo
    enddo

    do k=1,ncomp
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          e(i,j,1) = e(i,j,1)+dat(i,j,k)
        enddo
      enddo
    enddo

    return

  end subroutine dsrhoydot
  
!     
! ::: This routine will derive rho*R*T
!

  subroutine drhort (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                     lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                     level,grid_no) &
                     bind(C, name="drhort")

    use network,        only : nspec
    use PeleLM_2d, only: pphys_PfromRTY
    
    implicit none
      
    integer    lo(2), hi(2)
    integer    DIMDEC(e)
    integer    DIMDEC(dat)
    integer    domlo(2), domhi(2)
    integer    nv, ncomp
    integer    bc(2,2,ncomp)
    REAL_T     delta(2), xlo(2), time, dt
    REAL_T     e(DIMV(e),nv)
    REAL_T     dat(DIMV(dat),ncomp)
    integer    level, grid_no

    integer    i, j, n, rho, T, fS
    integer    nlft, nrgt, nbot, ntop
    REAL_T     Yt(nspec)
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/
      
!     ::::: lets punt if not in domain interior
    nlft = max(0,domlo(1)-lo(1))
    nrgt = max(0,hi(1)-domhi(1))
    nbot = max(0,domlo(2)-lo(2))
    ntop = max(0,hi(2)-domhi(2))

    if (nlft+nrgt+nbot+ntop .gt. 0) then
      call bl_abort("FORT_DERRHORT: outside domain")
    endif
      
!     Set pointers into state (these must agree with setup for this derived quant)
    rho = 1
    T = 2
    fS = 3

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n=1,Nspec
          Yt(n) = dat(i,j,fS+n-1) / dat(i,j,rho)
         end do

         call pphys_PfromRTY(lo_chem, hi_chem, &
              e(i,j,1),     ARLIM(lo_chem),ARLIM(hi_chem), &
              dat(i,j,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
              dat(i,j,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
              Yt,           ARLIM(lo_chem),ARLIM(hi_chem))
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

    integer i,j,n
    integer fS,rho
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/

    rho = 1
    fS = 2

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          x(i,j,n) = dat(i,j,fS+n-1)/dat(i,j,rho)
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

    integer i,j,n
    integer fS,rho
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/

    fS = 1

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          x(i,j,n) = dat(i,j,fS+n-1)
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
    use PeleLM_2D, only : pphys_mass_to_mole

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

    integer i,j,n
    REAL_T Yt(nspec),Xt(nspec)
    integer fS,rho
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/

    rho = 1 
    fS = 2

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          Yt(n) = dat(i,j,fS+n-1)/dat(i,j,rho) 
         enddo

         call pphys_mass_to_mole(lo_chem, hi_chem, &
                          Yt, ARLIM(lo_chem),ARLIM(hi_chem), &
                          Xt, ARLIM(lo_chem),ARLIM(hi_chem))
          do n = 1,Nspec
            x(i,j,n) = Xt(n)
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
    use PeleLM_2D, only: pphys_massr_to_conc
                               
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

    integer i,j,n
    REAL_T Yt(nspec),Ct(nspec)
    integer fS,rho,T
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/

    rho = 1 
    T   = 2
    fS  = 3

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          Yt(n) = dat(i,j,fS+n-1)/dat(i,j,rho) 
        enddo

        call pphys_massr_to_conc(lo_chem,hi_chem, &
                  Yt,           ARLIM(lo_chem),ARLIM(hi_chem), &
                  dat(i,j,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
                  dat(i,j,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
                  Ct,           ARLIM(lo_chem),ARLIM(hi_chem))
        do n = 1,Nspec
          C(i,j,n) = Ct(n)
        enddo
      enddo
    enddo

  end subroutine derconcentration
  
!=========================================================
  
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

    integer i,j,n
    REAL_T Yt(nspec), rho_dummy(1), D(Nspec), MU(1), XI(1), LAM(1)
    integer fS,rho,T
    integer lo_chem(3),hi_chem(3)
    data lo_chem /1,1,1/
    data hi_chem /1,1,1/

    rho = 1 
    T   = 2
    fS  = 3

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          Yt(n) = dat(i,j,fS+n-1)/dat(i,j,rho) 
        enddo
        rho_dummy(1) = dat(i,j,rho) * 1.d-3
        
        call get_transport_coeffs(lo_chem,hi_chem, &
                                  Yt,    lo_chem,hi_chem,  &
                                  dat(i,j,T),  lo_chem,hi_chem,  &
                                  rho_dummy(1), lo_chem,hi_chem,  &
                                  D,      lo_chem,hi_chem,  &
                                  MU(1),     lo_chem,hi_chem,  &
                                  XI(1),     lo_chem,hi_chem,  &
                                  LAM(1),    lo_chem,hi_chem)

        do n = 1,Nspec
          C(i,j,n) = D(n) * 0.1d0
        enddo

       C(i,j,Nspec+1) = LAM(1) * 1.0d-05
       C(i,j,Nspec+2) = MU(1) * 0.1d0

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

    integer i,j,n
    integer fS,rho
    integer lo_chem(2),hi_chem(2)
    data lo_chem /1,1/
    data hi_chem /1,1/

    rho = 1
    fS = 2

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n = 1,Nspec
          Yt(n) = dat(i,j,fS+n-1)/dat(i,j,rho)
         enddo
         
         CALL CKMMWY(Yt,x(i,j,1))
         
      enddo
    enddo

  end subroutine dermolweight
  
end module derive_PLM_2D
