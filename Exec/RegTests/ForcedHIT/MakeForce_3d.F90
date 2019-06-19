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

  subroutine FORT_MAKEFORCE(time,force,rho, &
                            DIMS(istate),DIMS(state), &
                            dx,xlo,xhi,gravity,scomp,ncomp) &
                            bind(C,name="FORT_MAKEFORCE")

      use mod_Fvar_def, only : dv_control, pseudo_gravity, dim

      implicit none

      integer    DIMDEC(state)
      integer    DIMDEC(istate)
      integer    scomp, ncomp
      REAL_T     time, dx(dim)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
      REAL_T     rho    (DIMV(state))
      REAL_T     gravity


!
!     ::::: local variables
!
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  hx, hy, hz
      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      ilo = istate_l1
      jlo = istate_l2
      klo = istate_l3
      ihi = istate_h1
      jhi = istate_h2
      khi = istate_h3

!     Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nZvel = 3
      nRho  = 4
      nTrac = 5

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

                       f1 = f1 + xT *   &
                            ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                            *   sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                            *   cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                            *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                            - FAY(kx,ky,kz)*twicePi*(kzd/HLz)         &
                            *   sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                            *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                            *   cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                       f2 = f2 + xT * &
                            ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                            *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                            *   sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                            *   cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                            - FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                            *   cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                            *   sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                            *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                       f3 = f3 + xT *  &
                            ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                            *   cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                            *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                            *   sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                            - FAX(kx,ky,kz)*twicePi*(kyd/HLy)  &
                            *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                            *   cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                            *   sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )

                     else

                       f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                               *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                               *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                       f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                               *                     cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                               *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                       f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                               *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                               *                     cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
    
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
                        f1 = f1 + xT *
                             ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)
                             *   sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz))
                             *   cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz))
                             *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) 
                             - FAY(kx,ky,kz)*twicePi*(kzd/HLz)
                             *   sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz))
                             *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz))
                             *   cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                        f2 = f2 + xT *
                             ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)
                             *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz))
                             *   sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz))
                             *   cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) 
                             - FAZ(kx,ky,kz)*twicePi*(kxd/HLx)
                             *   cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz))
                             *   sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz))
                             *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                        f3 = f3 + xT *
                             ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)
                             *   cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz))
                             *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz))
                             *   sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) 
                             - FAX(kx,ky,kz)*twicePi*(kyd/HLy)
                             *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz))
                             *   cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz))
                             *   sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                      else

                        f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz))
                             *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz))
                             *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                        f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz))
                             *                     cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz))
                             *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                        f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz))
                             *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz))
                             *                     cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                      endif
                    endif
                 enddo
               enddo
             enddo

             if (use_rho_in_forcing.eq.1) then
                  force(i,j,k,nXvel) = f1*rho(i,j,k)
                  force(i,j,k,nYvel) = f2*rho(i,j,k)
                  force(i,j,k,nZvel) = f3*rho(i,j,k)
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
               force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*rho(i,j,k)
             enddo
           enddo
         enddo
       endif

      if ((scomp+ncomp).gt.BL_SPACEDIM) then
!     Scalar forcing
         do n = max(scomp+1,nRho), scomp+ncomp
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

  end subroutine FORT_MAKEFORCE


end module MakeForce_3d_module

