#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module MakeForce_2d_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none
  
  private
  
  public :: FORT_MAKEFORCE

contains
  

!
!
! ::: -----------------------------------------------------------
!
!     This routine add the forcing terms to the momentum equation
!

      subroutine FORT_MAKEFORCE(time,force, &
                               vel, &
                               scal, &
                               DIMS(force), &
                               DIMS(vel), &
                               DIMS(scal), &
                               dx,xlo,xhi,gravity,scomp,ncomp, &
                               nscal,getForceVerbose &
      )bind(C, name="FORT_MAKEFORCE")

      use mod_Fvar_def, only : dv_control, pseudo_gravity

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     time, dx(dim)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     gravity
      integer    DIMDEC(vel)
      integer    getForceVerbose, nscal
      REAL_T     vel    (DIMV(vel),0:dim-1)
      REAL_T     scal   (DIMV(scal),0:nscal-1)

      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      REAL_T  hx, hy
      integer isioproc
      integer nXvel, nYvel, nRho, nTrac
      integer nRhoScal

      REAL_T  forcemin(scomp:scomp+ncomp-1)
      REAL_T  forcemax(scomp:scomp+ncomp-1)

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)

      ilo = force_l1
      jlo = force_l2
      ihi = force_h1
      jhi = force_h2
 
!     Assumes components are in the following order
      nXvel = 0
      nYvel = 1
      nRho  = 2
      nTrac = 3

      nRhoScal   = nRho-dim

      if (scomp.eq.0) then
         if (abs(gravity).gt.0.0001) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = gravity*scal(i,j,nRhoScal)
               enddo
            enddo
!     else to zero
         else
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = zero
               enddo
            enddo
         endif
!     Add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nYvel) = force(i,j,nYvel) + dV_control*scal(i,j,nRhoScal)
               enddo
            enddo
         endif
!     End of velocity forcing
      endif

      if ((scomp+ncomp).gt.BL_SPACEDIM) then
!     Scalar forcing
         do n = max(scomp,nRho), scomp+ncomp-1
            if (n.eq.nRho) then
!     Density
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
!     Tracer
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
!     Other scalar
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif

      if (getForceVerbose.gt.0 .and. isioproc .eq. 1) then
         do n = scomp,scomp+ncomp-1
            forcemin(n) = 1.d234
            forcemax(n) = -1.d234
         enddo
         do j = jlo, jhi
            do i = ilo, ihi
               do n = scomp,ncomp+scomp-1
                  forcemin(n) = min(forcemin(n),force(i,j,n))
                  forcemax(n) = max(forcemax(n),force(i,j,n))
               enddo
            enddo
         enddo
         do n = scomp,ncomp+scomp-1
            write (6,*) "forcemin (",n,") = ",forcemin(n)
            write (6,*) "forcemax (",n,") = ",forcemax(n)
         enddo
      endif

  end subroutine FORT_MAKEFORCE




end module MakeForce_2d_module

