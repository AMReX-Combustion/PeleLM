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
 
  public :: drhomry, dsrhoydot, drhort, dermolefrac, derconcentration

contains
 
  subroutine drhomry (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                      lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                      level,grid_no) &
                      bind(C, name="drhomry")
                                     

    implicit none

!
! ::: This routine will computes rho - sum (rho*Y)
!

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

!     
!     end of routine
!     
    return

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
!
!     end of routine
!
    return

  end subroutine dsrhoydot

!=========================================================

  subroutine drhort (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                     lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                     level,grid_no) &
                     bind(C, name="drhort")

    use network,        only : nspec
    use PeleLM_2d, only: pphys_PfromRTY
    
    implicit none

!     
! ::: This routine will derive rho*R*T
!
      
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

end module derive_PLM_2D
