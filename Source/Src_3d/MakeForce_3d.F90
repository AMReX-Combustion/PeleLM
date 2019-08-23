#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module MakeForce_3d_module

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
      integer    nscal, getForceVerbose
      REAL_T     vel    (DIMV(vel),0:dim-1)
      REAL_T     scal   (DIMV(scal),0:nscal-1)

!
!     ::::: local variables
!
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  hx, hy, hz
      integer isioproc, count
      integer nXvel, nYvel, nZvel, nRho, nTrac, nRhoScal

      REAL_T  velmin(0:dim-1)
      REAL_T  velmax(0:dim-1)
      REAL_T  scalmin(0:nscal-1)
      REAL_T  scalmax(0:nscal-1)
      REAL_T  forcemin(scomp:scomp+ncomp-1)
      REAL_T  forcemax(scomp:scomp+ncomp-1)


      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      ilo = force_l1
      jlo = force_l2
      klo = force_l3
      ihi = force_h1
      jhi = force_h2
      khi = force_h3

!     Assumes components are in the following order
      nXvel  = 0
      nYvel  = 1
      nZvel  = 2
      nRho   = 3
      nTrac  = 4

      nRhoScal   = nRho-dim

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
 
            write (6,*) "In MAKEFORCE"

            write (6,*) "gravity = ",gravity
            write (6,*) "scomp = ",scomp
            write (6,*) "ncomp = ",ncomp
            write (6,*) "nscal = ",nscal

            do n = 0, dim-1
               velmin(n) = 1.d234
               velmax(n) = -1.d234
            enddo
            do n = 0, nscal-1
               scalmin(n) = 1.d234
               scalmax(n) = -1.d234
            enddo

            count = 0
!c     Get min/max
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
!c     Velocities
                     do n = 0, dim-1
                        if (vel(i,j,k,n).gt.velmax(n)) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if (vel(i,j,k,n).lt.velmin(n)) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo
!c     Scalars
                     if (scal(i,j,k,0).lt.0.001) then
                        count=count+1
!c                        write (*,*) i,j,k,scal(i,j,k,n)
                     endif
                     do n = 0, nscal-1
                        if (scal(i,j,k,n).gt.scalmax(n)) then
                           scalmax(n)=scal(i,j,k,n)
                        endif
                        if (scal(i,j,k,n).lt.scalmin(n)) then
                           scalmin(n)=scal(i,j,k,n)
                        endif
                     enddo

                  enddo
               enddo
            enddo

            write(*,*) "DODGY DENSITY COUNT = ",count

            do n = 0, dim-1
               write (6,*) "velmin (",n,") = ",velmin(n)
               write (6,*) "velmax (",n,") = ",velmax(n)
            enddo
            do n = 0, nscal-1
               write (6,*) "scalmin(",n,") = ",scalmin(n)
               write (6,*) "scalmax(",n,") = ",scalmax(n)
            enddo
         endif
      endif

!c
!c     Here's where the forcing actually gets done
!c


      if (scomp.eq.0) then
        if (abs(gravity).gt.0.0001) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
                  enddo
               enddo
            enddo
!     else to zero
         else
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = zero
                  enddo
               enddo
            enddo
         endif
!     Let's add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*scal(i,j,k,nRhoScal)
                  enddo
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
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac) then
!     Tracer
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else
!     Other scalar
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

      if (getForceVerbose.gt.0 .and. isioproc.eq.1) then
         do n = scomp,scomp+ncomp-1
            forcemin(n) = 1.d234
            forcemax(n) = -1.d234
         enddo
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi
                  do n = scomp,ncomp+scomp-1
                     forcemin(n) = min(forcemin(n),force(i,j,k,n))
                     forcemax(n) = max(forcemax(n),force(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
         do n = scomp,ncomp+scomp-1
            write (6,*) "forcemin (",n,") = ",forcemin(n)
            write (6,*) "forcemax (",n,") = ",forcemax(n)
         enddo
      endif

  end subroutine FORT_MAKEFORCE


end module MakeForce_3d_module

