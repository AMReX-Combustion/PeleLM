#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module MakeForce_3d_module

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

      use mod_Fvar_def, only : dv_control, pseudo_gravity, dim, domnlo, domnhi
      use probdata_module, only : FTX, TAT, FPX, FPY, FPZ, &
                                  FAX, FAY, FAZ, FPXX, FPYX, FPXY, FPYY, &
                                  FPZY, FPYZ, FPZX, FPXZ, FPZZ, &
                                  blrandseed, div_free_force, &
                                  mode_start, moderate_zero_modes, nmodes, time_offset, &
                                  use_rho_in_forcing

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
      integer :: i, j, k, n
      integer :: ilo, jlo, klo
      integer :: ihi, jhi, khi
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho, nTrac
      integer :: xstep, ystep, zstep
      integer :: nxmodes, nymodes, nzmodes
      integer :: kx, ky, kz
      integer :: count, nRhoScal   
 
      double precision :: f1, f2, f3
      double precision :: kappa, kappaMax
      double precision :: kxd, kyd, kzd
      double precision :: Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      double precision :: force_time
      double precision :: twicepi, xT, x, y, z
  
      double precision ::   velmin(0:dim-1)
      double precision ::   velmax(0:dim-1)
      double precision ::   scalmin(0:nscal-1)
      double precision ::   scalmax(0:nscal-1)
      double precision ::   forcemin(scomp:scomp+ncomp-1)
      double precision ::   forcemax(scomp:scomp+ncomp-1)


      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

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
!     Get min/max
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
!     Velocities
                     do n = 0, dim-1
                        if (vel(i,j,k,n).gt.velmax(n)) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if (vel(i,j,k,n).lt.velmin(n)) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo
!     Scalars
                     if (scal(i,j,k,0).lt.0.001) then
                        count=count+1
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
!     Homogeneous Isotropic Turbulence
       twicePi=two*Pi

!     Adjust time offset if controlled:  FIXME: What is this really doing??
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

       HLx = Lx
       HLy = Ly
       HLz = Lz


       do k = klo, khi
         z = xlo(3) + dx(3)*(float(k-klo) + half)
         do j = jlo, jhi
           y = xlo(2) + dx(2)*(float(j-jlo) + half)
           do i = ilo, ihi
             x = xlo(1) + dx(1)*(float(i-ilo) + half)
             
             
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
                             f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                         -    FAY(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                             f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                         -    FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                             f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                         -    FAX(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )

                         else
                             f1 = f1 + xT*FAX(kx,ky,kz)* cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                             f2 = f2 + xT*FAY(kx,ky,kz)* sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                             f3 = f3 + xT*FAZ(kx,ky,kz)* sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

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
                      if (kappa.le.kappaMax) then
                       write(*,*) 'COUCOU'
                         xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                         if (div_free_force.eq.1) then
                             f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                         -    FAY(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                             f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                         -    FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                             f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                         -    FAX(kx,ky,kz)*twicePi*(kyd/HLy) & 
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )

                         else
                             f1 = f1 + xT*FAX(kx,ky,kz)* cos(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                             f2 = f2 + xT*FAY(kx,ky,kz)* sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * cos(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                             f3 = f3 + xT*FAZ(kx,ky,kz)* sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * cos(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
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
           enddo
         enddo
       enddo


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

