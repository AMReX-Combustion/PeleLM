
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3

module PeleLM_3d

  implicit none

  use fuego_chemistry

  private

  public ::  calc_divu_fortran, calc_gamma_pinv, floor_spec, enth_diff_terms, &
             compute_rho_dgrad_hdot_grad_Y, vel_visc, spec_temp_visc, &
             est_divu_dt, check_divu_dt, dqrad_fill, divu_fill, &
             dsdt_fill, ydot_fill, rhoYdot_fill, fab_minmax, repair_flux, &
             incrwext_flx_div, flux_div, compute_ugradp, conservative_T_floor, &
             part_cnt_err, mcurve, smooth, grad_wbar, recomp_update, &
             pphys_PfromRTY, pphys_mass_to_mole, pphys_massr_to_conc, pphys_HfromT, &
             pphys_HMIXfromTY, pphys_RHOfromPTY 

contains
             
  subroutine calc_divu_fortran(lo, hi, &
                              divu, DIMS(divu), rYdot, DIMS(rYdot), &
                              vtY,  DIMS(vtY),    vtT, DIMS(vtT), &
                              rhoY, DIMS(rhoY),     T, DIMS(T)) &
                              bind(C, name="calc_divu_fortran")

      use network,        only : nspec

      implicit none

#include <visc.H>
#include <htdata.H>

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rYdot)
      integer DIMDEC(vtY)
      integer DIMDEC(vtT)
      integer DIMDEC(rhoY)
      integer DIMDEC(T)
      REAL_T  divu(DIMV(divu))
      REAL_T  rYdot(DIMV(rYdot),1:Nspec)
      REAL_T  vtY(DIMV(vtY),1:Nspec)
      REAL_T  vtT(DIMV(vtT))
      REAL_T  rhoY(DIMV(rhoY),1:Nspec)
      REAL_T  T(DIMV(T))
      
      integer i, j, k, n
      REAL_T Y(nspec), H(nspec), cpmix, rhoInv, tmp, mmw, invmtw(nspec)

      call CKWT(invmtw)
      do n=1,Nspec
         invmtw(n) = one / invmtw(n)
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoInv = 0.d0
               do n=1,Nspec
                  rhoInv = rhoInv + rhoY(i,j,k,n)
               enddo
               rhoInv = 1.d0 / rhoInv
               do n=1,Nspec
                  Y(n) = rhoInv*rhoY(i,j,k,n)
               enddo
               CALL CKCPBS(T(i,j,k),Y,cpmix)
               CALL CKHMS(T(i,j,k),H)
               CALL CKMMWY(Y,mmw)

               cpmix = cpmix*1.d-4
               do n=1,Nspec
                  H(n) = H(n)*1.d-4
               enddo

               tmp = 0.d0
               divu(i,j,k) = vtT(i,j,k)
               do n=1,Nspec
                  tmp = tmp + (rYdot(i,j,k,n)+vtY(i,j,k,n))*invmtw(n)
                  divu(i,j,k) = divu(i,j,k) - rYdot(i,j,k,n)*H(n)
               enddo
               divu(i,j,k) = ( divu(i,j,k)/(cpmix*T(i,j,k)) + tmp*mmw ) * rhoInv
            enddo
         enddo
      enddo

  end subroutine calc_divu_fortran
  
!----------------------------------------- 

  subroutine calc_gamma_pinv(lo, hi, &
                             theta, DIMS(theta), &
                             rhoY, DIMS(rhoY), &
                             T, DIMS(T), &
                             Pamb_in)&
                             bind(C, name="calc_gamma_pinv")
                                                        
      use network,        only : nspec

      implicit none
      
#include <visc.H>
#include <htdata.H>
      
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(theta)
      integer DIMDEC(rhoY)
      integer DIMDEC(T)
      REAL_T  theta(DIMV(theta))
      REAL_T  rhoY(DIMV(rhoY),1:Nspec)
      REAL_T  T(DIMV(T))
      REAL_T  Pamb_in
      
      integer i, j, k, n
      REAL_T Y(nspec), cpmix, cvmix, rhoInv

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoInv = 0.d0
               do n=1,Nspec
                  rhoInv = rhoInv + rhoY(i,j,k,n)
               enddo
               rhoInv = 1.d0 / rhoInv
               do n=1,Nspec
                  Y(n) = rhoInv*rhoY(i,j,k,n)
               enddo
               CALL CKCPBS(T(i,j,k),Y,cpmix)
               cpmix = cpmix*1.d-4
               CALL CKCVBS(T(i,j,k),Y,cvmix)
               cvmix = cvmix*1.d-4
               
               theta(i,j,k) = cvmix / (cpmix*Pamb_in)
               
            enddo
         enddo
      enddo
  end subroutine calc_gamma_pinv

!-------------------------------------  
  
  subroutine floor_spec(lo, hi,spec,  DIMS(spec))bind(C, name="floor_spec")

    use network,        only : nspec

    implicit none
      
#include <visc.H>
#include <htdata.H>
      
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(spec)
      REAL_T  spec(DIMV(spec),1:Nspec)
      
      integer i, j, k, n

      do n=1,Nspec
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  spec(i,j,k,n) = max(0.d0,spec(i,j,k,n))
               end do
            enddo
         enddo
      enddo

  end subroutine floor_spec
  
!-----------------------------------  
  
  subroutine enth_diff_terms (lo, hi, dlo, dhi, dx, &
                              T, DIMS(T), RhoY, DIMS(RhoY), &
                              rhoDx, DIMS(rhoDx), Fx, DIMS(Fx), Ax, DIMS(Ax), &
                              rhoDy, DIMS(rhoDy), Fy, DIMS(Fy), Ay, DIMS(Ay), &
                              rhoDz, DIMS(rhoDz), Fz, DIMS(Fz), Az, DIMS(Az), &
                              FiGHi, DIMS(FiGHi), Tbc ) &
                              bind(C, name="enth_diff_terms")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM), dlo(SDIM), dhi(SDIM), Tbc(SDIM,2)
      REAL_T  dx(SDIM)
      integer DIMDEC(T)
      REAL_T  T(DIMV(T))
      integer DIMDEC(RhoY)
      REAL_T  RhoY(DIMV(RhoY),Nspec)

      integer DIMDEC(rhoDx)
      REAL_T  rhoDx(DIMV(rhoDx),Nspec+2)
      integer DIMDEC(Fx)
      REAL_T  Fx(DIMV(Fx),Nspec+3)
      integer DIMDEC(Ax)
      REAL_T  Ax(DIMV(Ax))

      integer DIMDEC(rhoDy)
      REAL_T  rhoDy(DIMV(rhoDy),Nspec+2)
      integer DIMDEC(Fy)
      REAL_T  Fy(DIMV(Fy),Nspec+3)
      integer DIMDEC(Ay)
      REAL_T  Ay(DIMV(Ay))

      integer DIMDEC(rhoDz)
      REAL_T  rhoDz(DIMV(rhoDz),Nspec+2)
      integer DIMDEC(Fz)
      REAL_T  Fz(DIMV(Fz),Nspec+3)
      integer DIMDEC(Az)
      REAL_T  Az(DIMV(Az))

      integer DIMDEC(FiGHi)
      REAL_T  FiGHi(DIMV(FiGHi))

      REAL_T, allocatable :: H(:,:,:,:), AD(:,:,:,:)

      integer i, j, k, d, n
      integer lob(SDIM), hib(SDIM)
      REAL_T AxDxInv_lo, AxDxInv_hi, dxInv
      REAL_T AyDyInv_lo, AyDyInv_hi, dyInv
      REAL_T AzDzInv_lo, AzDzInv_hi, dzInv
      logical fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi

      REAL_T, allocatable :: rhoInv(:,:,:)
      REAL_T gradY

      fix_xlo = .false.
      fix_xhi = .false.
      fix_ylo = .false.
      fix_yhi = .false.
      fix_zlo = .false.
      fix_zhi = .false.

!     Compute species enthalpies on box grown by one
      do d=1,SDIM
         lob(d) = lo(d)-1
         hib(d) = hi(d)+1
      enddo

!     Make space for Hi, use T box, since this better be big enough as well.
!     Note that any cells on a physical boundary with Dirichlet conditions will 
!     actually be centered on the edge, so the stencils below must reflect this

      allocate( H(DIMV(T),1:Nspec) )

      call pphys_HfromT(lob, hib, H, DIMS(T), T, DIMS(T))

!     On entry, Fx(1:Nspec) = spec flux
!     On exit:
!     Fx(Nspec+1) = untouched
!     Fx(Nspec+2) = sum[ (species flux).(species enthalpy) ]
!     Fx(Nspec+3) = extensive heat conduction

!     Compute lambda.Grad(T)
      dxInv = 1.d0 / dx(1)
      dyInv = 1.d0 / dx(2)
      dzInv = 1.d0 / dx(3)

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               Fx(i,j,k,Nspec+3) = - rhoDx(i,j,k,Nspec+2)*(T(i,j,k) - T(i-1,j,k))* dxInv * Ax(i,j,k)
            enddo
         enddo
      enddo

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = - rhoDy(i,j,k,Nspec+2)*(T(i,j,k) - T(i,j-1,k)) * dyInv * Ay(i,j,k)
            enddo
         enddo
      enddo

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = - rhoDz(i,j,k,Nspec+2)*(T(i,j,k) - T(i,j,k-1)) * dzInv * Az(i,j,k)
            enddo
         enddo
      enddo

!     xlo
      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               Fx(i,j,k,Nspec+3) = 2*Fx(i,j,k,Nspec+3)
            enddo
         enddo
      endif
!     xhi
      if (hi(1)+1.eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               Fx(i,j,k,Nspec+3) = 2*Fx(i,j,k,Nspec+3)
            enddo
         enddo
      endif
!     ylo
      if (lo(2).eq.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
         j=lo(2)
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = 2*Fy(i,j,k,Nspec+3)
            enddo
         enddo
      endif
!     yhi
      if (hi(2)+1.eq.dhi(2)+1 .and. Tbc(2,2).eq.EXT_DIR) then
         j=hi(2)+1
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = 2*Fy(i,j,k,Nspec+3)
            enddo
         enddo
      endif
!     zlo
      if (lo(3).eq.dlo(3) .and. Tbc(3,1).eq.EXT_DIR) then
         k=lo(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = 2*Fz(i,j,k,Nspec+3)
            enddo
         enddo
      endif
!     zhi
      if (hi(3)+1.eq.dhi(3)+1 .and. Tbc(3,2).eq.EXT_DIR) then
         k=hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = 2*Fz(i,j,k,Nspec+3)
            enddo
         enddo
      endif

!     Compute enthalpy flux as hi*(Fi+(lambda/cp).Grad(Yi))

      Fx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),Nspec+2) = 0.d0
      Fy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),Nspec+2) = 0.d0
      Fz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,Nspec+2) = 0.d0

      allocate(rhoInv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

      rhoInv = 0.0d0
      do n=1,Nspec
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
               rhoInv(i,j,k) = rhoInv(i,j,k) + RhoY(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      rhoInv(:,:,:) = 1.0D0/rhoInv(:,:,:)

      do n=1,Nspec
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
         Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2) &
             + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*(H(i,j,k,n)+H(i-1,j,k,n))*0.5d0
      enddo
      enddo
      enddo
      enddo

      do n=1,Nspec
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
         Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2) &
             + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*(H(i,j,k,n)+H(i,j-1,k,n))*0.5d0
      enddo
      enddo
      enddo
      enddo

      do n=1,Nspec
      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
         Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2) &
             + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*(H(i,j,k,n)+H(i,j,k-1,n))*0.5d0
      enddo
      enddo
      enddo
      enddo

!     xlo
      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         Fx(i,lo(2):hi(2),lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
                  Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2) &
                      + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*H(i-1,j,k,n)
               enddo
            enddo
         enddo
      endif
!     xhi
      if (hi(1)+1.eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         Fx(i,lo(2):hi(2),lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
                  Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2) &
                      + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
!     ylo
      if (lo(2).eq.dlo(2)  .and.  Tbc(2,1).eq.EXT_DIR) then
         j = dlo(2)
         Fy(lo(1):hi(1),j,lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
                  Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2) &
                      + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*H(i,j-1,k,n)
               enddo
            enddo
         enddo
      endif
!     yhi
      if (hi(2)+1.eq.dhi(2)+1  .and.  Tbc(2,2).eq.EXT_DIR) then
         j = dhi(2)+1
         Fy(lo(1):hi(1),j,lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
                  Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2) &
                      + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
!     zlo
      if (lo(3).eq.dlo(3)  .and.  Tbc(3,1).eq.EXT_DIR) then
         k = dlo(3)
         Fz(lo(1):hi(1),lo(2):hi(2),k,Nspec+2) = 0.d0
         do n=1,Nspec
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
                  Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2) &
                      + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*H(i,j,k-1,n)
               enddo
            enddo
         enddo
      endif
!     zhi
      if (hi(3)+1.eq.dhi(3)+1  .and.  Tbc(3,2).eq.EXT_DIR) then
         k = dhi(3)+1
         Fz(lo(1):hi(1),lo(2):hi(2),k,Nspec+2) = 0.d0
         do n=1,Nspec
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
                  Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2) &
                      + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
      deallocate(rhoInv)


!     Set FiGHi = (species flux) dot Grad(species enthalpy)
!        compute Grad(H) on each face, and average across faces in each coordinate
!        Fi is extensive here, so need to remove area.  Also, assume that we are 
!        away from domain boundary, fix afterward
!
!     FIXME: This will fail for r-z since Ax(dlo(1),:)=0
!

      !
      ! Use AD to cut down on some divides.
      ! Unfortunately it ups memory usage a bit.
      !
      allocate(AD(6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         AD(1,i,j,k) = 1.d0/(Ax(i,j,k)*dx(1))
         AD(2,i,j,k) = 1.d0/(Ax(i+1,j,k)*dx(1))
         AD(3,i,j,k) = 1.d0/(Ay(i,j,k)*dx(2))
         AD(4,i,j,k) = 1.d0/(Ay(i,j+1,k)*dx(2))
         AD(5,i,j,k) = 1.d0/(Az(i,j,k)*dx(3))
         AD(6,i,j,k) = 1.d0/(Az(i,j,k+1)*dx(3))
      enddo
      enddo
      enddo

      FiGHi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.0D0


      do n=1,Nspec
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
             ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AD(2,i,j,k) &
             + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AD(1,i,j,k) &
             + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AD(4,i,j,k) &
             + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AD(3,i,j,k) &
             + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AD(6,i,j,k) &
             + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AD(5,i,j,k) )
      enddo
      enddo
      enddo
      enddo

      deallocate(AD)

!     xlo
      if (fix_xlo) then
         i = lo(1)
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo &
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

!     xhi
      if (fix_xhi) then
         i = hi(1)
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo & 
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

!     ylo
      if (fix_ylo) then
         j = lo(2)
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo &
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

!     yhi
      if (fix_yhi) then
         j = hi(2)
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo &
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

!     zlo
      if (fix_zlo) then
         k = lo(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo & 
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

!     zhi
      if (fix_zhi) then
         k = hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            FiGHi(i,j,k) = 0.d0
         enddo
         enddo
         do n=1,Nspec
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            AxDxInv_lo = 2.d0/(Ax(i,j,k)*dx(1))
            AxDxInv_hi = 1.d0/(Ax(i+1,j,k)*dx(1))
            AyDyInv_lo = 1.d0/(Ay(i,j,k)*dx(2))
            AyDyInv_hi = 1.d0/(Ay(i,j+1,k)*dx(2))
            AzDzInv_lo = 1.d0/(Az(i,j,k)*dx(3))
            AzDzInv_hi = 1.d0/(Az(i,j,k+1)*dx(3))
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0* &
                ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi &
                + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo &
                + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi &
                + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo &
                + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi &
                + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

      deallocate(H)
 
  end subroutine enth_diff_terms

!-------------------------------------

  subroutine pphys_RRATERHOY(lo,hi,RhoY,DIMS(RhoY),RhoH,DIMS(RhoH),T,DIMS(T), &
                       RhoYdot,DIMS(RhoYdot))&
                       bind(C, name="pphys_RRATERHOY")
      
      use network,        only : nspec
      use PeleLM_F,       only : pphys_calc_src_sdc

      implicit none

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(RhoY)
      integer DIMDEC(RhoH)
      integer DIMDEC(T)
      integer DIMDEC(RhoYdot)
      REAL_T RhoY(DIMV(RhoY),nspec)
      REAL_T RhoH(DIMV(RhoH))
      REAL_T T(DIMV(T))
      REAL_T RhoYdot(DIMV(RhoYdot),nspec)

      REAL_T Zt(nspec+1),Zdott(nspec+1)
      REAL_T Temperature, TIME
      integer i,j,k,n

      TIME = 0.

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
 
               Zt(Nspec+1) = RhoH(i,j,k)
               do n=1,Nspec
                  Zt(n) = RhoY(i,j,k,n)
               end do
               Temperature = T(i,j,k)

               call pphys_calc_src_sdc(nspec,TIME,Temperature,Zt,Zdott)

               do n=1,Nspec
                  RhoYdot(i,j,k,n) = Zdott(n)
               end do
            end do
         end do
      end do
      
  end subroutine pphys_RRATERHOY

!-------------------------------------

  subroutine pphys_PfromRTY(lo, hi, P, DIMS(P), RHO, DIMS(RHO), &
                      T, DIMS(T), Y, DIMS(Y))&
                      bind(C, name="pphys_PfromRTY")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(P)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T P(DIMV(P))
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspec), RHOt, SCAL, SCAL1
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
!           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
      SCAL = 1.d-1
      SCAL1 = SCAL**3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n=1,nspec
                   Yt(n) = Y(i,j,k,n)
                end do

                RHOt = RHO(i,j,k) * SCAL1
                CALL CKPY(RHOt,T(i,j,k),Yt,P(i,j,k))

                P(i,j,k) = P(i,j,k) * SCAL

            end do
         end do
      end do

  end subroutine pphys_PfromRTY

!-------------------------------------

  subroutine pphys_mass_to_mole(lo, hi, Y, DIMS(Y), X, DIMS(X)) &
                          bind(C, name="pphys_mass_to_mole")

      use network,        only : nspec

      implicit none

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(X)
      REAL_T Y(DIMV(Y),*)
      REAL_T X(DIMV(X),*)

      REAL_T Xt(nspec), Yt(nspec)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      
  end subroutine pphys_mass_to_mole

!-------------------------------------

  subroutine pphys_massr_to_conc(lo, hi, Y, DIMS(Y), &
                           T, DIMS(T), RHO, DIMS(RHO), C, DIMS(C))&
                           bind(C, name="pphys_massr_to_conc")
                                   
      use network,        only : nspec

      implicit none

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      integer DIMDEC(RHO)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      REAL_T RHO(DIMV(RHO))

      REAL_T Yt(nspec), Ct(nspec), rhoScl
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               rhoScl = RHO(i,j,k)*1.e-3
               CALL CKYTCR(rhoScl,T(i,j,k),Yt,Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*1.e6
               end do
            end do
         end do
      end do

  end subroutine pphys_massr_to_conc

!-------------------------------------


  subroutine pphys_HMIXfromTY(lo, hi, HMIX, DIMS(HMIX), T, DIMS(T), &
                        Y, DIMS(Y))&
                        bind(C, name="pphys_HMIXfromTY")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(HMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T HMIX(DIMV(HMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspec), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n=1,Nspec
                   Yt(n) = Y(i,j,k,n)
                end do

                CALL CKHBMS(T(i,j,k),Yt,HMIX(i,j,k))

                HMIX(i,j,k) = HMIX(i,j,k) * SCAL

            end do
         end do
      end do

  end subroutine pphys_HMIXfromTY

!-------------------------------------

  subroutine pphys_RHOfromPTY(lo, hi, RHO, DIMS(RHO), T, DIMS(T), &
                        Y, DIMS(Y), Patm) &
                        bind(C, name="pphys_RHOfromPTY")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(nspec), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      CALL CKRP(RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                do n=1,Nspec
                   Yt(n) = Y(i,j,k,n)
                end do
                CALL CKRHOY(Ptmp,T(i,j,k),Yt,RHO(i,j,k))
                RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
  end subroutine pphys_RHOfromPTY

!-------------------------------------

  subroutine pphys_HfromT(lo, hi, H, DIMS(H), T, DIMS(T))&
                    bind(C, name="pphys_HfromT")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(H)
      integer DIMDEC(T)
      REAL_T H(DIMV(H),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, Ht(nspec)
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),Ht)
               do n=1,Nspec
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do

  end subroutine pphys_HfromT

!-------------------------------------

  subroutine pphys_MWMIXfromY(lo, hi, MWMIX, DIMS(MWMIX), Y, DIMS(Y))&
                        bind(C, name="pphys_MWMIXfromY")

      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(MWMIX)
      integer DIMDEC(Y)
      REAL_T MWMIX(DIMV(MWMIX))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspec)

!     Returns mean molecular weight in kg/kmole

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,nspec
                  Yt(n) = Y(i,j,k,n)
               end do

               CALL CKMMWY(Yt,MWMIX(i,j,k))

            end do
         end do
      end do

  end subroutine pphys_MWMIXfromY

!-------------------------------------

  subroutine pphys_CPMIXfromTY(lo, hi, &
                               CPMIX, cmix_lo, cmix_hi, &
                               T, T_lo, T_hi, &
                               Y, Y_lo, Y_hi )&
                         bind(C,name="pphys_CPMIXfromTY")

      use network,        only : nspec
                         
      implicit none

      integer         , intent(in   ) ::     lo(3),      hi(3)
      integer         , intent(in   ) ::   T_lo(3),    T_hi(3)
      integer         , intent(in   ) ::   Y_lo(3),    Y_hi(3)
      integer         , intent(in   ) ::cmix_lo(3),cmix_hi(3)
      REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T, intent(in   ) :: Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspec)
      REAL_T, intent(out  ) :: CPMIX(cmix_lo(1):cmix_hi(1),cmix_lo(2):cmix_hi(2),cmix_lo(3):cmix_hi(3))
      
      integer i, j, k, n
      REAL_T Yt(nspec), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            do n=1,nspec
               Yt(n) = Y(i,j,k,n)
            end do
            CALL CKCPBS(T(i,j,k),Yt,CPMIX(i,j,k))
            CPMIX(i,j,k) = CPMIX(i,j,k) * SCAL
          end do
        end do
      end do

  end subroutine pphys_CPMIXfromTY

!-------------------------------------

  integer function pphys_TfromHY(lo, hi, T, DIMS(T), &
                           HMIX, DIMS(HMIX), Y, DIMS(Y), &
                           errMax, NiterMAX, res) &
                           bind(C, name="pphys_TfromHY")
                           
      use network,        only : nspec
      use PeleLM_F,       only : pphys_TfromHYpt
      
      implicit none

      integer lo(SDIM), hi(SDIM)
      integer NiterMAX
      integer DIMDEC(T)
      integer DIMDEC(HMIX)
      integer DIMDEC(Y)
      REAL_T T(DIMV(T))
      REAL_T HMIX(DIMV(HMIX))
      REAL_T Y(DIMV(Y),*)
      REAL_T errMAX
      REAL_T res(0:NiterMAX-1)

      REAL_T Yt(nspec), lres(0:NiterMAX-1)
      integer i, j, k, n, lierr, Niter,MAXiters
      REAL_T HMIX_CGS, Tguess

      MAXiters = 0
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,nspec
                  Yt(n) = Y(i,j,k,n)
               end do

               call pphys_TfromHYpt(T(i,j,k),HMIX(i,j,k),Yt,errMax,NiterMAX,lres,Niter)

               if (Niter .lt. 0) then
                       call bl_abort(" Something went wrong in pphys_TfromHYpt ")
               end if

               if (Niter .gt. MAXiters) then
                       MAXiters = Niter
               end if

            end do
         end do
      end do

!     Set max iters taken during this solve, and exit
      pphys_TfromHY = MAXiters 

      return

  end function pphys_TfromHY

!-------------------------------------

  subroutine compute_rho_dgrad_hdot_grad_Y(dx, &
              lo, hi, DIMS(species), species, &
              DIMS(h), h, DIMS(betax), betax, &
              DIMS(betay), betay, DIMS(betaz), betaz,  &
              DIMS(rdghdgy), rdghdgy) &
               bind(C, name="compute_rho_dgrad_hdot_grad_y")

      implicit none

      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(species)
      integer DIMDEC(h)
      REAL_T  species(DIMV(species))
      REAL_T  h(DIMV(h))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
      integer DIMDEC(betaz)
      REAL_T betax(DIMV(betax))
      REAL_T betay(DIMV(betay))
      REAL_T betaz(DIMV(betaz))
      integer DIMDEC(rdghdgy)

      REAL_T rdghdgy(DIMV(rdghdgy))

      integer i,j,k
      REAL_T  dxsqr, dysqr, dzsqr
      REAL_T  bdotxlo, bdotxhi, bdotylo, bdotyhi
      REAL_T  bdotzlo, bdotzhi

      dxsqr = 1.0D0/dx(1)**2
      dysqr = 1.0D0/dx(2)**2
      dzsqr = 1.0D0/dx(3)**2

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            bdotxlo = betax(i,j,k  )*(h(i,j,k)-h(i-1,j,k)) &
                 *(species(i,j,k)-species(i-1,j,k))
            bdotxhi = betax(i+1,j,k)*(h(i+1,j,k)-h(i,j,k)) &
                *(species(i+1,j,k)-species(i,j,k))
            bdotylo = betay(i,j,k  )*(h(i,j,k)-h(i,j-1,k)) &
                *(species(i,j,k)-species(i,j-1,k))
            bdotyhi = betay(i,j+1,k)*(h(i,j+1,k)-h(i,j,k)) &
                *(species(i,j+1,k)-species(i,j,k))
            bdotzlo = betaz(i,j,k  )*(h(i,j,k)-h(i,j,k-1)) &
                *(species(i,j,k)-species(i,j,k-1))
            bdotzhi = betaz(i,j,k+1)*(h(i,j,k+1)-h(i,j,k)) &
                *(species(i,j,k+1)-species(i,j,k))
            rdghdgy(i,j,k) =  half*((bdotxlo + bdotxhi)*dxsqr + &
                                   (bdotylo + bdotyhi)*dysqr + &
                                   (bdotzlo + bdotzhi)*dysqr)
          enddo
        enddo
      enddo

  end subroutine compute_rho_dgrad_hdot_grad_y

!--------------------------------------------------------------

  subroutine vel_visc (lo,hi, &
                       T, T_lo, T_hi, &
                       Y, Y_lo, Y_hi, &
                       mu, mu_lo, mu_hi) &
                       bind(C, name="vel_visc")
                                 
      use network,          only : nspec
      use transport_module, only: get_visco_coeffs

      implicit none
      
#include <visc.H>

      integer         , intent(in   ) ::     lo(3),    hi(3)
      integer         , intent(in   ) ::   T_lo(3),  T_hi(3)
      integer         , intent(in   ) ::   Y_lo(3),  Y_hi(3)
      integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
      REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T, intent(in   ) :: Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspec)
      REAL_T, intent(out  ) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))

      integer i, j, k

      if (.not.use_constant_mu) then
         if (.not. LeEQ1) then
            call get_visco_coeffs(lo, hi, &
                              Y,  Y_lo,  Y_hi, &
                              T,  T_lo,  T_hi, &
                              mu, mu_lo, mu_hi)
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     mu(i,j,k) = mu(i,j,k) * 1.0d-1
                  end do
               end do
            end do
         else 
            do k=lo(3), hi(3)
              do j=lo(2), hi(2)
                do i=lo(1), hi(1)
                  mu(i,j,k) = 1.85e-5*(T(i,j,k)/298.0)**.7
                end do
              end do
            end do
         end if
      else
         do k=lo(3), hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  mu(i,j,k) = constant_mu_val
               end do
            end do
         end do
      end if
  end subroutine vel_visc 

!-------------------------------------------------

  subroutine spec_temp_visc(lo,hi, &
                            T, T_lo, T_hi, &
                            Y, Y_lo, Y_hi, &
                            rhoD, rhoD_lo, rhoD_hi, &
                            ncompd, P1ATM_MKS, do_temp, do_VelVisc, &
                            Pamb_in) &
                            bind(C, name="spec_temp_visc")
       
      use network,          only : nspec
      use transport_module, only : get_transport_coeffs
    
      implicit none
      
#include <visc.H>
#include <htdata.H>

      integer         , intent(in   ) ::     lo(3),    hi(3)
      integer         , intent(in   ) ::   T_lo(3),  T_hi(3)
      integer         , intent(in   ) ::   Y_lo(3),  Y_hi(3)
      integer         , intent(in   ) ::rhoD_lo(3),  rhoD_hi(3)
      REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T, intent(in   ) :: Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspec)
      REAL_T, intent(out  ) :: rhoD(rhoD_lo(1):rhoD_hi(1),rhoD_lo(2):rhoD_hi(2),rhoD_lo(3):rhoD_hi(3),nspec+2)
      integer ncompd, do_temp, do_VelVisc
      REAL_T  Pamb_in, P1ATM_MKS

      integer i, j, k, n, nc, ncs
      REAL_T Patm, Wavg
      REAL_T Yt(nspec), invmwt(nspec)
      REAL_T cpmix(1,1,1), Tfac, Yfac
      REAL_T  Y_real(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspec)
      REAL_T  RHO(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T  D(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3),nspec)
      REAL_T  MU(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T  XI(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T  LAM(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      integer lo_chem(3),hi_chem(3)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

!     Nspec+1 is Temp stuff, if requested
!     Nspec+2 is Velocity stuff, if requested
      nc = Nspec
      ncs = nc
      if (do_temp.eq.1) then
         nc = nc + 1
         ncs = nc
         if (do_VelVisc.eq.1) then
            nc = nc + 1
         endif
      endif
      if (ncompd.lt.nc) then
         call bl_pd_abort("not enough components in rhoD")
      endif
      if (do_temp.eq.0 .and. do_VelVisc.eq.1) then
         call bl_pd_abort("if do_VelVisc, rhoD requires do_temp")
      endif

!     Warning, FORT_VELVISC is called separately from this routine, so if there's
!     any hacking to be done on viscosity, be sure to catch it there as well.
      Tfac = thickFacTR / Pr
      Yfac = thickFacTR / Sc

      call CKWT(invmwt)
      do n=1,Nspec
          invmwt(n) = one / invmwt(n)
      end do

      if (.not.use_constant_rhoD) then
          do k=lo(3),hi(3) 
            do j=lo(2),hi(2)
              do i=lo(1),hi(1)
                RHO(i,j,k) = 0.d0
                do n=1,nspec
                    RHO(i,j,k) = RHO(i,j,k) + Y(i,j,k,n)
                end do
                do n=1,nspec
                    Y_real(i,j,k,n) = Y(i,j,k,n) / RHO(i,j,k)
                end do
                RHO(i,j,k) = RHO(i,j,k) * 1.d-3
              end do
            end do
          end do
          if (.not. LeEQ1) then
            Patm = Pamb_in / P1ATM_MKS
            call get_transport_coeffs(lo,   hi, &
                                  Y_real, Y_lo,  Y_hi,  &
                                  T,      T_lo,  T_hi,  &
                                  RHO,    T_lo,  T_hi,  &
                                  D,      T_lo,  T_hi,  &
                                  MU,     T_lo,  T_hi,  &
                                  XI,     T_lo,  T_hi,  &
                                  LAM,    T_lo,  T_hi)
            do k=lo(3),hi(3) 
              do j=lo(2), hi(2)
                do i=lo(1), hi(1)
                  CALL CKMMWY(Y_real(i,j,k,:),Wavg)
                  do n=1,nspec
                    rhoD(i,j,k,n) = Wavg * invmwt(n) * D(i,j,k,n)  * 1.0d-1 
                    !rhoD(i,j,k,n) = rhoD(i,j,k,n) * thickFacTR
                  end do
                  if (do_temp .ne. 0) then 
                    rhoD(i,j,k,nspec+1) = LAM(i,j,k) * (one / 100000.0D0)
                  end if
                  if (do_VelVisc .ne. 0) then 
                    rhoD(i,j,k,nspec+2) = MU(i,j,k) * 1.0d-1
                  end if
                end do
              end do
            end do
          else
            call vel_visc(lo,      hi,&
                      T,       T_lo,    T_hi, &
                      Y_real,  Y_lo,    Y_hi, &
                      rhoD,    rhoD_lo, rhoD_hi)
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     do n=Nspec+1,nc
                        rhoD(i,j,k,n) = rhoD(i,j,k,1)
                     end do
                     rhoD(i,j,k,1) = rhoD(i,j,k,1) * Yfac
                  end do
               end do
            end do
            do n=2,Nspec
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        rhoD(i,j,k,n) = rhoD(i,j,k,1)
                     end do
                  end do
               end do
            end do

            if (do_temp .ne. 0) then
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        do n=1,Nspec
                           Yt(n) = Y_real(i,j,k,n)
                        end do
                        CALL pphys_CPMIXfromTY(lo,       hi, & 
                                       cpmix,    lo_chem, hi_chem, &
                                       T(i,j,k), lo_chem, hi_chem, &
                                       Yt,       lo_chem, hi_chem)
                        rhoD(i,j,k,Nspec+1) = rhoD(i,j,k,Nspec+1)*cpmix(1,1,1)*Tfac
                     end do
                  end do
               end do
            end if
          end if
      else
         do n=1,ncs 
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     rhoD(i,j,k,n) = constant_rhoD_val
                  end do
               end do
            end do
         end do
      end if
      
  end subroutine spec_temp_visc

!-------------------------------------------

  subroutine est_divu_dt(flag, dtfactor, delta, divu, DIMS(divu), &
                         dsdt, DIMS(dsdt), rho, DIMS(rho),  &
                         u, DIMS(u), &
                         volume, DIMS(volume), &
                         areax,  DIMS(areax), &
                         areay,  DIMS(areay), &
                         areaz,  DIMS(areaz), &
                         lo, hi, dt, rhomin) &
                         bind(C, name="est_divu_dt")

      implicit none

      integer flag
      integer lo(SDIM), hi(SDIM)
      REAL_T  delta(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      integer DIMDEC(dsdt)
      REAL_T  rho(DIMV(rho))      
      REAL_T  u(DIMV(u),BL_SPACEDIM)      
      REAL_T  rhomin, dtfactor
      REAL_T  divu(DIMV(divu))
      REAL_T  dsdt(DIMV(dsdt))

      integer DIMDEC(volume)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer DIMDEC(areaz)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))
      REAL_T  areaz(DIMV(areaz))

      REAL_T dt

      integer i,j,k
      REAL_T  dtcell, dtcell2, denom, rhominij, rhoij
      REAL_T  fluxxlo, fluxxhi, fluxylo, fluxyhi, fluxzlo, fluxzhi
      REAL_T  a,b,c

      dt = 1.0D20

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dtcell = dt
               if (flag.eq.1) then
                  if(divu(i,j,k).gt.zero) then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(one-rhomin/rho(i,j,k))/divu(i,j,k)
                     else
                        dtcell = dtfactor*0.5D0/divu(i,j,k)
                     endif
                     if (dsdt(i,j,k).gt.1.0D-20) then
                        if (abs(rho(i,j,k)).gt.rhomin) then
                           rhominij = rhomin
                        else
                           rhominij = 0.9D0*abs(rho(i,j,k)) 
                        endif
                        rhoij = abs(rho(i,j,k))
!
!     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!     We use the latter because it is more robust
!
                        a = rhoij*dsdt(i,j,k)*half
                        b = rhoij*divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))

                        dtcell2 = dtfactor*dtcell2
                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
                  if(dtcell.le.zero)then
                     write(6,*)'aha'
                  endif
               else if (flag.eq.2) then
                  denom = rho(i,j,k)*divu(i,j,k)+ &
                      u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) + &
                      u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) + &
                      u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = dtfactor*abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k)) &
                                   *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k)) &
                                   *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k)) &
                                   *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k)) &
                                   *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k-1))  &
                                   *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k+1)) &
                                   *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+ &
                          (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+ &
                          (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo)) &
                           /volume(i,j,k)
                  
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = dtfactor*abs(rho(i,j,k))/denom
                     endif
                  endif
               endif
#if 0 
               write(6,*)'i,j,k,dtcell=',i,j,k,dtcell
#endif
               dt = min(dtcell,dt)
            enddo
         enddo
      enddo

  end subroutine est_divu_dt

!---------------------------------------------------------------

  subroutine check_divu_dt(flag, dtfactor, delta, divu, DIMS(divu), &
                           dsdt, rho, DIMS(rho),  &
                           u, DIMS(u),  &
                           volume, DIMS(volume), &
                           areax,  DIMS(areax), &
                           areay,  DIMS(areay), &
                           areaz,  DIMS(areaz), &
                           lo, hi, &
                           dt, rhomin) &
                           bind(C, name="check_divu_dt")

      implicit none

      integer flag
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      REAL_T  delta(SDIM)
      REAL_T  rho(DIMV(rho))      
      REAL_T  u(DIMV(u),BL_SPACEDIM)      
      REAL_T  rhomin
      REAL_T  divu(DIMV(divu))
      REAL_T  dsdt(DIMV(divu))
      REAL_T  dt, dtfactor

      integer DIMDEC(volume)
      integer DIMDEC(areax)
      integer DIMDEC(areay)
      integer DIMDEC(areaz)
      REAL_T  volume(DIMV(volume))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))
      REAL_T  areaz(DIMV(areaz))

      integer i,j,k
      REAL_T  dtcell, denom
      REAL_T  fluxxlo,fluxxhi,fluxylo,fluxyhi,fluxzlo,fluxzhi
      REAL_T  a,b,c,dtcell2,rhominij,rhoij

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dtcell = bigreal
               if (flag.eq.1) then
                  if(divu(i,j,k).gt.zero) then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (one-rhomin/rho(i,j,k))/divu(i,j,k)
                     else
                        dtcell = one/divu(i,j,k)
                     endif
                     if (dsdt(i,j,k).gt.1.0D-20) then
                        if (abs(rho(i,j,k)).gt.rhomin) then
                           rhominij = rhomin
                        else
                           rhominij = 0.9D0*abs(rho(i,j,k)) 
                        endif
                        rhoij = abs(rho(i,j,k))
!
!     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!     We use the latter because it is more robust
!
                        a = rhoij*dsdt(i,j,k)*half
                        b = rhoij*divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))

                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
               else if (flag.eq.2) then
                  denom = rho(i,j,k)*divu(i,j,k)+ &
                      u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) + &
                      u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) + &
                      u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k)) &
                                   *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k)) &
                                   *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k)) &
                                   *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k)) &
                                   *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k-1)) &
                                   *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k+1)) &
                                   *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+ &
                          (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+ &
                          (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo)) &
                           /volume(i,j,k)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               endif
               if (dt.gt.dtcell) then
                  write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell = ', &
                      i,j,k,dt,dtcell
               else if (dt.gt.dtcell*dtfactor) then
                  write(6,*)'WARNING: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell*dtfactor = ', &
                      i,j,k,dt,dtcell*dtfactor
               endif
            enddo
         enddo
      enddo

  end subroutine check_divu_dt

!===================================================================
!
! ... The routines YDOTFILL, DIVUFILL, DQRADFILL, and DSDTFILL
!     are here instead of in the problem dependent code because
!     we always fill the quantitities ydot, divu, dqrad, and dsdt
!     the same way
!
!===================================================================

  subroutine dqrad_fill  (dqrad,DIMS(dqrad),domlo,domhi,delta, &
                          xlo,time,bc )bind(C, name="dqrad_fill")

      integer    DIMDEC(dqrad)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dqrad(DIMV(dqrad))

      call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)
      call fillEdges(dqrad,DIMS(dqrad),domlo,domhi,bc)
      
  end subroutine dqrad_fill

!--------------------------------------------------

  subroutine divu_fill  (divu,DIMS(divu),domlo,domhi,delta, &
                         xlo,time,bc )bind(C, name="divu_fill")

      integer    DIMDEC(divu)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      call fillEdges(divu,DIMS(divu),domlo,domhi,bc)

  end subroutine divu_fill

!--------------------------------------------

  subroutine dsdt_fill (dsdt,DIMS(dsdt),domlo,domhi,delta, &
                        xlo,time,bc )bind(C, name="dsdt_fill")

      integer    DIMDEC(dsdt)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(dsdt,DIMS(dsdt),domlo,domhi,bc)

  end subroutine dsdt_fill

!--------------------------------------

  subroutine ydot_fill (ydot,DIMS(ydot),domlo,domhi,delta, &
                        xlo,time,bc)bind(C, name="ydot_fill")

      integer    DIMDEC(ydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     ydot(DIMV(ydot))

      call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(ydot,DIMS(ydot),domlo,domhi,bc)

  end subroutine ydot_fill

!-------------------------------------------

  subroutine rhoYdot_fill (rhoydot,DIMS(rhoydot),domlo,domhi,delta, &
                           xlo,time,bc)bind(C, name="rhoYdot_fill")

      integer    DIMDEC(rhoydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     rhoydot(DIMV(rhoydot))

      call filcc (rhoydot,DIMS(rhoydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(rhoydot,DIMS(rhoydot),domlo,domhi,bc)

  end subroutine rhoYdot_fill

!------------------------------------------------

  subroutine fab_minmax(lo, hi, &
                        fab, DIMS(fab), &
                        fmin, fmax, nc) &
                        bind(C, name="fab_minmax")
                        
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fab)
      REAL_T  fab(DIMV(fab),nc)
      REAL_T  fmin, fmax

      integer i,j,k,n


      do n = 1,nc
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  fab(i,j,k,n) = MAX( fmin, MIN( fmax, fab(i,j,k,n) ) )
               end do
            end do
         end do
      end do

  end subroutine fab_minmax

!--------------------------------------------------

  subroutine repair_flux (lo, hi, dlo, dhi, &
                          flux, DIMS(flux), &
                          RhoY, DIMS(RhoY), dir, Ybc)&
                          bind(C, name="repair_flux")
                          
      use network,        only : nspec

      implicit none

      integer lo(SDIM), hi(SDIM), dlo(SDIM), dhi(SDIM), dir, Ybc(SDIM,2)
      integer DIMDEC(flux)
      integer DIMDEC(RhoY)
      REAL_T flux(DIMV(flux),Nspec)
      REAL_T RhoY(DIMV(RhoY),Nspec)
      
      integer i, j, k, n
      REAL_T sumFlux, RhoYe(Nspec), sumRhoYe

      if (dir.eq.0) then

!     First, assume away from physical boundaries, then use boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,Nspec
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i-1,j,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  end do
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,Nspec
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                  end do
               end do
            end do
         end do
!     xlo
         if (Ybc(1,1).eq.EXT_DIR.and.lo(1).le.dlo(1)) then
            do i = lo(1),dlo(1)
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i-1,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i-1,j,k,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif
!     xhi
         if (Ybc(1,2).eq.EXT_DIR.and.hi(1).ge.dhi(1)) then
            do i = dhi(1),hi(1)
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i,j,k,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif

      else if (dir.eq.1) then

!     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,Nspec
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j-1,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,Nspec
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                  end do
               end do
            end do
         end do
!     ylo
         if (Ybc(2,1).eq.EXT_DIR.and.lo(2).le.dlo(2)) then
            do j = lo(2),dlo(2)
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j-1,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i,j-1,k,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif
!     yhi
         if (Ybc(2,2).eq.EXT_DIR.and.hi(2).ge.dhi(2)) then
            do j = dhi(2),hi(2)
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i,j,k,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif

      else if (dir.eq.2) then

!     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,Nspec
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j,k-1,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,Nspec
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                  end do
               end do
            end do
         end do
!     zlo
         if (Ybc(3,1).eq.EXT_DIR.and.lo(3).le.dlo(3)) then
            do k = lo(3),dlo(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k-1,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i,j,k-1,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif
!     yzi
         if (Ybc(3,2).eq.EXT_DIR.and.hi(3).ge.dhi(3)) then
            do k = dhi(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     sumFlux = 0.d0
                     sumRhoYe = 0.d0
                     do n=1,Nspec
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,Nspec
                        flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoY(i,j,k,n)*sumRhoYe
                     enddo
                  enddo
               enddo
            enddo
         endif

      endif
      
  end subroutine repair_flux

!-----------------------------------------

  subroutine incrwext_flx_div(lo, hi, &
                              xflux,  DIMS(xflux), &
                              yflux,  DIMS(yflux), &
                              zflux,  DIMS(zflux), &
                              stateo, DIMS(stateo), &
                              staten, DIMS(staten), &
                              vol,    DIMS(vol), &
                              nc, dt) &
                              bind(C, name="incrwext_flx_div")
                              
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(stateo)
      integer DIMDEC(staten)
      integer DIMDEC(vol)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T stateo(DIMV(stateo))
      REAL_T staten(DIMV(staten))
      REAL_T vol(DIMV(vol))
      REAL_T dt

      integer i, j, k, n
      REAL_T dF

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dF = zero
               do n=1,nc
                  dF = dF + ( xflux(i+1,j,k,n) - xflux(i,j,k,n) ) &
                      +    ( yflux(i,j+1,k,n) - yflux(i,j,k,n) ) &
                      +    ( zflux(i,j,k+1,n) - zflux(i,j,k,n) )
               end do
               staten(i,j,k) = stateo(i,j,k) + dF*dt/vol(i,j,k)
            end do
         end do
      end do

  end subroutine incrwext_flx_div

!----------------------------------------------

  subroutine flux_div (lo, hi, &
                       update, DIMS(update), &
                       xflux,  DIMS(xflux), &
                       yflux,  DIMS(yflux), &
                       zflux,  DIMS(zflux), &
                       vol,    DIMS(vol), &
                       nc, scal) &
                       bind(C, name="flux_div")
     
     
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T vol(DIMV(vol))
      REAL_T scal
      REAL_T, allocatable :: ivol(:,:,:)

      integer i, j, k, n

      allocate(ivol(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ivol(i,j,k) = scal / vol(i,j,k)
            end do
         end do
      end do

      do n=1,nc
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  update(i,j,k,n) = &
                      ( (xflux(i+1,j,k,n)-xflux(i,j,k,n)) &
                      + (yflux(i,j+1,k,n)-yflux(i,j,k,n))  &
                      + (zflux(i,j,k+1,n)-zflux(i,j,k,n)) ) * ivol(i,j,k)
               end do
            end do
         end do
      end do

  end subroutine flux_div

!-----------------------------------------  
  
  subroutine fillEdges(dat,DIMS(dat),domlo,domhi,bc )bind(C, name="fillEdges")
      
!  This routine fills ghost cells with the value from the nearest
!  interior cell.  
      integer    DIMDEC(dat)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dat(DIMV(dat))

      integer    lo(SDIM), hi(SDIM)
      integer    i, j, k
      integer    ilo, ihi, jlo, jhi, klo, khi
      logical    do_xlo,do_xhi,do_ylo,do_yhi,do_zlo,do_zhi

      lo(1) = ARG_L1(dat)
      hi(1) = ARG_H1(dat)
      lo(2) = ARG_L2(dat)
      hi(2) = ARG_H2(dat)
      lo(3) = ARG_L3(dat)
      hi(3) = ARG_H3(dat)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      klo = max(lo(3),domlo(3))
      khi = min(hi(3),domhi(3))

      do_xlo = lo(1).lt.domlo(1)
      do_xhi = hi(1).gt.domhi(1)
      do_ylo = lo(2).lt.domlo(2)
      do_yhi = hi(2).gt.domhi(2)
      do_zlo = lo(3).lt.domlo(3)
      do_zhi = hi(3).gt.domhi(3)

      if (bc(1,1).eq.EXT_DIR.and. do_xlo) then
         do k = klo,khi
            do j = jlo, jhi
               do i = lo(1), domlo(1)-1
                  dat(i,j,k) = dat(domlo(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            


      if (bc(1,2).eq.EXT_DIR.and. do_xhi) then
         do k = klo, khi
            do j = jlo, jhi
               do i = domhi(1)+1,hi(1)
                  dat(i,j,k) = dat(domhi(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1,hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(2,1).eq.EXT_DIR.and.do_ylo) then
         do k = klo, khi
            do j = lo(2), domlo(2)-1
               do i = ilo,ihi
                  dat(i,j,k) = dat(i,domlo(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = lo(2), domlo(2)-1
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = lo(2),domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = lo(2),domlo(2)-1
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            
      
      if (bc(2,2).eq.EXT_DIR.and.do_yhi) then
         do k = klo, khi
            do j = domhi(2)+1, hi(2)
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,domhi(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1, hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1, hi(2)
                     do i = lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = domhi(2)+1, hi(2)
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
            if (do_zlo) then
               do k = lo(3),domlo(3)-1
                  do j = domhi(2)+1, hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, hi(3)
                  do j = domhi(2)+1, hi(2)
                     do i = domhi(1)+1,hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_zlo) then
            do k = lo(3), domlo(3)-1
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, hi(3)
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(3,1).eq.EXT_DIR.and. do_zlo) then
         do k = lo(3), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,j,domlo(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = lo(3), domlo(3)-1
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = lo(3), domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = lo(3), domlo(3)-1
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = lo(3), domlo(3)-1
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = lo(3), domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = lo(3), domlo(3)-1
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(3,2).eq.EXT_DIR.and. do_zhi) then
         do k = domhi(3)+1,hi(3)
            do j = jlo, jhi
               do i = ilo,ihi
                  dat(i,j,k) = dat(i,j,domhi(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = domhi(3)+1,hi(3)
               do j = lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1,hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1,hi(3)
                  do j = lo(2), domlo(2)-1
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = domhi(3)+1,hi(3)
               do j = domhi(2)+1, hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1,hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1,hi(3)
                  do j = domhi(2)+1,hi(2)
                     do i = domhi(1)+1, hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = domhi(3)+1,hi(3)
               do j = jlo, jhi
                  do i = lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = domhi(3)+1,hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1, hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif
      
  end subroutine fillEdges

!-----------------------------------------  
  
  subroutine fillWithZeros(dat,DIMS(dat),domlo,domhi,bc)&
                           bind(C, name="fillWithZeros")

      integer    DIMDEC(dat)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dat(DIMV(dat))

      integer    lo(SDIM), hi(SDIM)
      integer    i, j, k

      lo(1) = ARG_L1(dat)
      hi(1) = ARG_H1(dat)
      lo(2) = ARG_L2(dat)
      hi(2) = ARG_H2(dat)
      lo(3) = ARG_L3(dat)
      hi(3) = ARG_H3(dat)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), domlo(1)-1
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = domhi(1)+1,hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do k=lo(3),hi(3)
            do j = lo(2), domlo(2)-1
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do k=lo(3),hi(3)
            do j = domhi(2)+1, hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k=lo(3),domlo(3)-1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif            
      
      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k=domhi(3)+1,hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif
      
  end subroutine fillWithZeros

!-----------------------------------------

  subroutine compute_ugradp(p, DIMS(p), ugradp, DIMS(ugp), &
                            umac,  DIMS(umac), &
                            vmac,  DIMS(vmac), &
                            wmac,  DIMS(wmac), &
                            lo, hi, dx)&
                            bind(C, name="compute_ugradp")

      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(ugp)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T      p(DIMV(p))
      REAL_T ugradp(DIMV(ugp))
      REAL_T dx(SDIM)

      integer i, j, k
      REAL_T uadv, vadv, wadv
      REAL_T p_x_lo, p_x_hi
      REAL_T p_y_lo, p_y_hi
      REAL_T p_z_lo, p_z_hi

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               uadv = half*(umac(i,j,k) + umac(i+1,j,k))
               vadv = half*(vmac(i,j,k) + vmac(i,j+1,k))
               wadv = half*(wmac(i,j,k) + wmac(i,j,k+1))
               p_x_hi = merge(p(i  ,j,k),p(i+1,j,k),umac(i+1,j,k)>=zero)
               p_x_lo = merge(p(i-1,j,k),p(i  ,j,k),umac(i  ,j,k)>=zero)
               p_y_hi = merge(p(i,j  ,k),p(i,j+1,k),vmac(i,j+1,k)>=zero)
               p_y_lo = merge(p(i,j-1,k),p(i,j  ,k),vmac(i,j  ,k)>=zero)
               p_z_hi = merge(p(i,j,k  ),p(i,j,k+1),wmac(i,j,k+1)>=zero)
               p_z_lo = merge(p(i,j,k-1),p(i,j,k  ),wmac(i,j,k  )>=zero)
               ugradp(i,j,k) = uadv * (p_x_hi - p_x_lo) / dx(1) + &
                   vadv * (p_y_hi - p_y_lo) / dx(2) +    &
                   wadv * (p_z_hi - p_z_lo) / dx(3)
            end do
         end do
      end do

  end subroutine compute_ugradp

!------------------------------------------------

  integer function conservative_T_floor ( &
              loF,hiF,state,DIMS(state), &
              min_T, Tcomp, Rcomp, first_spec, last_spec, RhoH, &
              ratio, tmp, nt) &
              bind(C, name="conservative_T_floor")
              
      implicit none

      integer loF(SDIM), hiF(SDIM)
      integer DIMDEC(state)
      REAL_T  state(DIMV(state),0:*)      
      integer Tcomp, Rcomp, first_spec, last_spec, RhoH, ratio(SDIM), nt
      REAL_T  min_T, tmp(0:nt-1)
      integer n,i,j,k, loC(SDIM),hiC(SDIM),ii,jj,kk,iii,jjj,kkk,ncells
      Real ncellsInv
      logical bad_T

!     Returns the number of fine cells fixed up
      conservative_T_floor = 0

      ncells = 1
      do n=1,SDIM
         loC(n) = loF(n)/ratio(n)
         hiC(n) = (hiF(n)+1)/ratio(n) - 1
         ncells = ncells*ratio(n)
      enddo
      ncellsInv = 1.d0 / ncells

      do k=loC(3),hiC(3)
         do j=loC(2),hiC(2)
            do i=loC(1),hiC(1)

               bad_T = .false.
               do kk=0,ratio(3)-1
                  kkk = ratio(3)*k + kk
                  do jj=0,ratio(2)-1
                     jjj = ratio(2)*j + jj
                     do ii=0,ratio(1)-1
                        iii = ratio(1)*i + ii
                        if (state(iii,jjj,kkk,Tcomp).lt.min_T) then
                           bad_T = .true.
                        endif
                     enddo
                  enddo
               enddo

               if (bad_T .eqv. .true.) then

                  tmp(Rcomp) = 0.d0
                  do n=first_spec,last_spec
                     tmp(n) = 0.d0
                  enddo
                  tmp(RhoH) = 0.d0


                  do kk=0,ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj=0,ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii=0,ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           tmp(Rcomp) = tmp(Rcomp) + state(iii,jjj,kkk,Rcomp)
                           do n=first_spec,last_spec
                              tmp(n) = tmp(n) + state(iii,jjj,kkk,n)
                           enddo
                           tmp(RhoH) = tmp(RhoH) + state(iii,jjj,kkk,RhoH)
                           
                        enddo
                     enddo
                  enddo

                  conservative_T_floor = conservative_T_floor + ncells
                  tmp(Rcomp) = tmp(Rcomp) * ncellsInv
                  do n=first_spec,last_spec
                     tmp(n) = tmp(n) * ncellsInv
                  enddo
                  tmp(RhoH) = tmp(RhoH)* ncellsInv
                  
                  do kk=0,ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj=0,ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii=0,ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           state(iii,jjj,kkk,Rcomp) = tmp(Rcomp)
                           do n=first_spec,last_spec
                              state(iii,jjj,kkk,n) = tmp(n)
                           enddo
                           state(iii,jjj,kkk,RhoH) = tmp(RhoH)
                           
                        enddo
                     enddo
                  enddo


               endif

            enddo
         enddo
      enddo
      
  end function conservative_T_floor

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on whether or not
! ::: they contain any particles.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

  subroutine part_cnt_err(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                          set,clear,  &
                          var,varl1,varl2,varl3,varh1,varh2,varh3, &
                          lo,hi,nd,domlo,domhi,  &
                          delta,xlo,problo,time,level) &
                          bind(C, name="part_cnt_err")

      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision :: var(varl1:varh1,varl2:varh2,varl3:varh3,nd)
      double precision :: delta(3), xlo(3), problo(3), time
      integer          :: i,j,k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         if (var(i,j,k,1) .gt. zero) tag(i,j,k) = set
         !if (var(i,j,k,1) .gt. zero) print *,'TAGGING ',i,j,k
      end do
      end do
      end do

  end subroutine part_cnt_err

!------------------------

  subroutine mcurve(lo, hi, T, DIMS(T), curv, DIMS(curv), &
          wrk, DIMS(wrk), delta)bind(C, name="mcurve")

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(curv)
      integer DIMDEC(wrk)
      REAL_T    T(DIMV(T))
      REAL_T curv(DIMV(curv))
      REAL_T wrk(DIMV(wrk),SDIM)
      REAL_T delta(SDIM)

      integer i,j,k
      REAL_T mag

!     Fill normal on nodes (assumes 1 grow cell properly filled)

      REAL_T Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Txz,tdxI(3),dxSqI(3)

      do i=1,3
         tdxI(i) =  one / (two*delta(i))
         dxSqI(i) = one / (delta(i)*delta(i))
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               Tx = tdxI(1)*(T(i+1,j,k) - T(i-1,j,k))
               Ty = tdxI(2)*(T(i,j+1,k) - T(i,j-1,k))
               Tz = tdxI(3)*(T(i,j,k+1) - T(i,j,k-1))

               Txx = dxSqI(1)*(T(i+1,j,k) - two*T(i,j,k) + T(i-1,j,k))
               Tyy = dxSqI(2)*(T(i,j+1,k) - two*T(i,j,k) + T(i,j-1,k))
               Tzz = dxSqI(3)*(T(i,j,k+1) - two*T(i,j,k) + T(i,j,k-1))

               Txy = tdxI(1)*tdxI(2)*(T(i+1,j+1,k) - T(i-1,j+1,k) - T(i+1,j-1,k) + T(i-1,j-1,k))
               Txz = tdxI(1)*tdxI(3)*(T(i+1,j,k+1) - T(i-1,j,k+1) - T(i+1,j,k-1) + T(i-1,j,k-1))
               Tyz = tdxI(2)*tdxI(3)*(T(i,j+1,k+1) - T(i,j-1,k+1) - T(i,j+1,k-1) + T(i,j-1,k-1))

               mag = max(1.0d-12, SQRT(Tx*Tx + Ty*Ty + Tz*Tz))

               curv(i,j,k) = -half*(Txx + Tyy + Tzz &
                   - ( Tx*(Tx*Txx + Ty*Txy + Tz*Txz) &
                     + Ty*(Tx*Txy + Ty*Tyy + Tz*Tyz) &
                     + Tz*(Tx*Txz + Ty*Tyz + Tz*Tzz) )/mag**2 )/mag
            end do
         end do
      end do

  end subroutine mcurve

!------------------------------------------

  subroutine smooth(lo, hi, Tin, DIMS(Tin), Tout, DIMS(Tout))bind(C, name="smooth")
  
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Tin)
      integer DIMDEC(Tout)
      REAL_T   Tin(DIMV(Tin))
      REAL_T  Tout(DIMV(Tout))

      integer i,j,k,ii,jj,kk

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Tout(i,j,k) = zero
               do kk=0,1
                  do jj=0,1
                     do ii=0,1
                        Tout(i,j,k) = Tout(i,j,k) &
                            + Tin(i+ii,j+jj,  k+kk-1) + Tin(i+ii-1,j+jj,  k+kk-1) &
                            + Tin(i+ii,j+jj-1,k+kk-1) + Tin(i+ii-1,j+jj-1,k+kk-1) &
                            + Tin(i+ii,j+jj,  k+kk)   + Tin(i+ii-1,j+jj,  k+kk  ) &
                            + Tin(i+ii,j+jj-1,k+kk  ) + Tin(i+ii-1,j+jj-1,k+kk  )
                     end do
                  end do
               end do
               Tout(i,j,k) = Tout(i,j,k) / 64.d0
            end do
         end do
      end do

  end subroutine smooth

!-----------------------------------------

  subroutine grad_wbar(lo, hi, Wbar, DIMS(Wbar), &
                       rDe, DIMS(rDe), flux, DIMS(flux), &
                       area, DIMS(area), dx, dir, mult, inc) &
                       bind(C, name="grad_wbar")
                       
      implicit none
      
      integer lo(SDIM), hi(SDIM)
      integer dir
      integer DIMDEC(Wbar)
      integer DIMDEC(rDe)
      integer DIMDEC(flux)
      integer DIMDEC(area)
      REAL_T  Wbar(DIMV(Wbar))
      REAL_T  rDe(DIMV(rDe))
      REAL_T  flux(DIMV(flux))
      REAL_T  area(DIMV(area))
      REAL_T  dx, mult, inc

      REAL_T Wgr, fac
      integer i,j,k

      fac = mult / dx

      if (inc .eq. 0) then

!     compute grad wbar fluxes

         if (dir.eq.0) then
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir.eq.1) then

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir.eq.2) then

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else
            call bl_pd_abort('Bad dir in FORT_GRADWBAR')
         endif

      else

!     increment grad wbar fluxes by a factor of inc (can be negative)

         if (dir.eq.0) then

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir.eq.1) then

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir.eq.2) then

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else
            call bl_pd_abort('Bad dir in FORT_GRADWBAR')
         endif

      end if

  end subroutine grad_wbar

!-----------------------------------

  subroutine recomp_update(lo, hi, &
                           update, DIMS(update), &
                           xflux,  DIMS(xflux), &
                           yflux,  DIMS(yflux), &
                           zflux,  DIMS(zflux), &
                           vol,    DIMS(vol), &
                           nc) &
                           bind(C, name="recomp_update")
                           
      implicit none
      
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T vol(DIMV(vol))

      integer i, j, k, n

      do n=1,nc
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  update(i,j,k,n)=-((xflux(i+1,j,k,n)-xflux(i,j,k,n)) &
                      +            (yflux(i,j+1,k,n)-yflux(i,j,k,n)) &
                      +            (zflux(i,j,k+1,n)-zflux(i,j,k,n))) &
                      /vol(i,j,k)
               end do
            end do
         end do
      end do

  end subroutine recomp_update

end module PeleLM_3d
