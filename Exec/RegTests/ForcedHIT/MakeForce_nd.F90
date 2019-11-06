#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module MakeForce_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort

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

      use mod_Fvar_def, only : dv_control, pseudo_gravity, domnlo, domnhi
      use probdata_module, only : FTX, TAT, FPX, FPY, FPZ, &
                                  FAX, FAY, FAZ, FPXX, FPYX, FPXY, FPYY, &
                                  FPZY, FPYZ, FPZX, FPXZ, FPZZ, &
                                  blrandseed, div_free_force, &
                                  mode_start, moderate_zero_modes, nmodes, time_offset, &
                                  use_rho_in_forcing

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
      integer :: i, j, k, n
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho
      integer :: xstep, ystep, zstep
      integer :: nxmodes, nymodes, nzmodes
      integer :: kx, ky, kz
      integer :: count, nRhoScal   

      double precision :: f1, f2, f3
      double precision :: kappa, kappaMax
      double precision :: kxd, kyd, kzd
      double precision :: Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      double precision :: HLx_inv, HLy_inv, HLz_inv
      double precision :: force_time
      double precision :: twicepi, xT, x, y, z
      double precision :: hx, hy, hz

      double precision ::   velmin(0:dim-1)
      double precision ::   velmax(0:dim-1)
      double precision ::   scalmin(0:nscal-1)
      double precision ::   scalmax(0:nscal-1)
      double precision ::   forcemin(scomp:scomp+ncomp-1)
      double precision ::   forcemax(scomp:scomp+ncomp-1)

      call bl_pd_is_ioproc(isioproc)

#if ( AMREX_SPACEDIM == 2 )
      call amrex_abort("Turbulent Forcing is only available in 3D")
#endif

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
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
!     Here's where the forcing actually gets done
!

      if (scomp.eq.0) then
         ! Homogeneous Isotropic Turbulence
         twicePi = two * Pi

         ! Adjust time offset if controlled:  FIXME: What is this really doing??
         if (time_offset.gt.(-half)) then
            force_time = time + time_offset
         else
            force_time = time
         endif

         Lx = domnhi(1)-domnlo(1)
         Ly = domnhi(2)-domnlo(2)
         Lz = domnhi(3)-domnlo(3)

         Lmin = min(Lx,Ly,Lz)
         kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5+Lx/Lmin)
         nymodes = nmodes*int(0.5+Ly/Lmin)
         nzmodes = nmodes*int(0.5+Lz/Lmin)

         xstep = int(Lx/Lmin+0.5)
         ystep = int(Ly/Lmin+0.5)
         zstep = int(Lz/Lmin+0.5)

         HLx_inv = 1.0d0/Lx
         HLy_inv = 1.0d0/Ly
         HLz_inv = 1.0d0/Lz

         do k = f_lo(3), f_hi(3)
            z = xlo(3) + dx(3)*(float(k-f_lo(3)) + half)
            do j = f_lo(2), f_hi(2)
               y = xlo(2) + dx(2)*(float(j-f_lo(2)) + half)
               do i = f_lo(1), f_hi(1)
                  x = xlo(1) + dx(1)*(float(i-f_lo(1)) + half)

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

                  !          write(*,*) 'DEBUG F1, F2, F3',kx,ky,kz,kappa,xT,f1,f2,f3
                  !          write(*,*) 'DEBUG x,y,z',x,y,z,kxd,kyd,kzd
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

                  if (use_rho_in_forcing.eq.1) then
                       force(i,j,k,nXvel) = f1*scal(i,j,k,nRhoScal)
                       force(i,j,k,nYvel) = f2*scal(i,j,k,nRhoScal)
                       force(i,j,k,nZvel) = f3*scal(i,j,k,nRhoScal)
                  else
                       force(i,j,k,nXvel) = f1
                       force(i,j,k,nYvel) = f2
                       force(i,j,k,nZvel) = f3
                  endif
                  !force(i,j,k,nXvel) = 0.0d0
                  !force(i,j,k,nYvel) = 0.0d0
                  !force(i,j,k,nZvel) = 0.0d0
               enddo
            enddo
         enddo

         ! Let's add the pseudo gravity afterwards...
         if ( pseudo_gravity == 1 ) then
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*scal(i,j,k,nRhoScal)
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
