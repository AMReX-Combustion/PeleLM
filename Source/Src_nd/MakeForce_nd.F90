#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module MakeForce_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public :: FORT_MAKEFORCE

contains

!
! ::: -----------------------------------------------------------
!     This routine add the forcing terms to the momentum equation

   subroutine FORT_MAKEFORCE( time, &
                              force, f_lo, f_hi,&
                              vel, v_lo, v_hi,&
                              scal, s_lo, s_hi,&
                              dx, xlo, xhi, gravity, scomp, ncomp, &
                              nscal, getForceVerbose ) &
                              bind(C, name="FORT_MAKEFORCE")

      use mod_Fvar_def, only : dv_control, pseudo_gravity

      implicit none

! In/Out
      integer :: f_lo(3), f_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: s_lo(3), s_hi(3)
      integer :: scomp, ncomp
      integer :: nscal, getForceVerbose
      REAL_T  :: time, dx(3)
      REAL_T  :: xlo(3), xhi(3)
      REAL_T  :: gravity
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),scomp:scomp+ncomp-1) :: force
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),0:dim-1) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nscal-1) :: scal

! Local
      REAL_T  :: velmin(0:dim-1)
      REAL_T  :: velmax(0:dim-1)
      REAL_T  :: scalmin(0:nscal-1)
      REAL_T  :: scalmax(0:nscal-1)
      REAL_T  :: forcemin(scomp:scomp+ncomp-1)
      REAL_T  :: forcemax(scomp:scomp+ncomp-1)
      REAL_T  :: hx, hy, hz
      integer :: isioproc, count
      integer :: nXvel, nYvel, nZvel, nRho, nTrac, nRhoScal

      integer :: i, j, k, n

      call bl_pd_is_ioproc(isioproc)

      if (isioproc==1 .and. pseudo_gravity==1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

!     Assumes components are in the following order
      nXvel  = 0
      nYvel  = 1
      nZvel  = dim-1
      nRho   = dim
      nTrac  = dim+1

      nRhoScal   = nRho-dim

      if ( getForceVerbose > 0 ) then
         call bl_pd_is_ioproc(isioproc)
         if ( isioproc == 1 ) then

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

            ! Get min/max
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)

                     ! Velocities
                     do n = 0, dim-1
                        if ( vel(i,j,k,n) > velmax(n) ) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if ( vel(i,j,k,n) < velmin(n) ) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo

                     ! Scalars
                     if ( scal(i,j,k,0) < 0.001 ) then
                        count = count + 1
!                        write (*,*) i,j,k,scal(i,j,k,n)
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

!
! Here's where the forcing actually gets done
!

      if ( scomp == 0 ) then
        if ( abs(gravity) > 0.0001d0 ) then
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nXvel) = zero
#if ( AMREX_SPACEDIM == 2 )
                     force(i,j,k,nYvel) = gravity*scal(i,j,k,nRhoScal)
#elif ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
#endif
                  enddo
               enddo
            enddo
         ! else to zero
         else
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
#if ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nZvel) = zero
#endif
                  enddo
               enddo
            enddo
         endif

         ! Let's add the pseudo gravity afterwards...
         if ( pseudo_gravity == 1 ) then
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
#if ( AMREX_SPACEDIM == 2 )
                     force(i,j,k,nYvel) = force(i,j,k,nYvel)  + dV_control*scal(i,j,k,nRhoScal)
#elif ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*scal(i,j,k,nRhoScal)
#endif
                  enddo
               enddo
            enddo
         endif
         ! End of velocity forcing
      endif

      if ( (scomp+ncomp) > AMREX_SPACEDIM) then
         ! Scalar forcing
         do n = max(scomp,nRho), scomp+ncomp-1
            if ( n == nRho) then
               ! Density
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if ( n == nTrac ) then
               ! Tracer
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else
               ! Other scalar
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

      if ( getForceVerbose>0 .and. isioproc==1) then
         do n = scomp,scomp+ncomp-1
            forcemin(n) = 1.d234
            forcemax(n) = -1.d234
         enddo
         do k = f_lo(3), f_hi(3)
            do j = f_lo(2), f_hi(2)
               do i = f_lo(1), f_hi(1)
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

end module MakeForce_nd_module

