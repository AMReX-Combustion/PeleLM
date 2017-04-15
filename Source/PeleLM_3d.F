
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <ChemDriver_F.H>
#include <AMReX_ArrayLim.H>

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3

      subroutine FORT_CALCDIVU(lo, hi,
     &                         divu, DIMS(divu), rYdot, DIMS(rYdot),
     &                         vtY,  DIMS(vtY),    vtT, DIMS(vtT),
     &                         rhoY, DIMS(rhoY),     T, DIMS(T))
      implicit none
#include <cdwrk.H>
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
      REAL_T Y(maxspec), H(maxspec), cpmix, rhoInv, tmp, mmw, invmtw(maxspec)

      call CKWT(IWRK(ckbi),RWRK(ckbr),invmtw)
      do n=1,Nspec
         invmtw(n) = one / invmtw(n)
      end do

!$omp parallel do private(i,j,k,n,rhoInv,Y,cpmix,H,mmw,tmp)
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
               CALL CKCPBS(T(i,j,k),Y,IWRK(ckbi),RWRK(ckbr),cpmix)
               CALL CKHMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),H)
               CALL CKMMWY(Y,IWRK(ckbi),RWRK(ckbr),mmw)

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
!$omp end parallel do
      end

      subroutine FORT_CALCGAMMAPINV(lo, hi,
     &                              theta, DIMS(theta),
     &                              rhoY, DIMS(rhoY),
     &                              T, DIMS(T),
     &                              Pamb_in)
      implicit none
#include <cdwrk.H>
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
      REAL_T Y(maxspec), cpmix, cvmix, rhoInv

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
               CALL CKCPBS(T(i,j,k),Y,IWRK(ckbi),RWRK(ckbr),cpmix)
               cpmix = cpmix*1.d-4
               CALL CKCVBS(T(i,j,k),Y,IWRK(ckbi),RWRK(ckbr),cvmix)
               cvmix = cvmix*1.d-4
               
               theta(i,j,k) = cvmix / (cpmix*Pamb_in)
               
            enddo
         enddo
      enddo
      end

      subroutine FORT_FLOORSPEC(lo, hi,
     &                          spec,  DIMS(spec))
      implicit none
#include <cdwrk.H>
#include <visc.H>
#include <htdata.H>
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(spec)
      REAL_T  spec(DIMV(spec),1:Nspec)
      
      integer i, j, k, n

!$omp parallel do private(i,j,k,n)
      do n=1,Nspec
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  spec(i,j,k,n) = max(0.d0,spec(i,j,k,n))
               end do
            enddo
         enddo
      enddo
!$omp end parallel do
      end

      subroutine FORT_ENTH_DIFF_TERMS(lo, hi, dlo, dhi, dx, T, DIMS(T), RhoY, DIMS(RhoY),
     &                                rhoDx, DIMS(rhoDx), Fx, DIMS(Fx), Ax, DIMS(Ax), 
     &                                rhoDy, DIMS(rhoDy), Fy, DIMS(Fy), Ay, DIMS(Ay), 
     &                                rhoDz, DIMS(rhoDz), Fz, DIMS(Fz), Az, DIMS(Az), 
     &                                FiGHi, DIMS(FiGHi), Tbc )

      implicit none

#include <cdwrk.H>      

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

c     Compute species enthalpies on box grown by one
      do d=1,SDIM
         lob(d) = lo(d)-1
         hib(d) = hi(d)+1
      enddo

c     Make space for Hi, use T box, since this better be big enough as well.
c     Note that any cells on a physical boundary with Dirichlet conditions will 
c     actually be centered on the edge, so the stencils below must reflect this

      allocate( H(DIMV(T),1:Nspec) )

      call FORT_HfromT(lob, hib, H, DIMS(T), T, DIMS(T))

c     On entry, Fx(1:Nspec) = spec flux
c     On exit:
c     Fx(Nspec+1) = untouched
c     Fx(Nspec+2) = sum[ (species flux).(species enthalpy) ]
c     Fx(Nspec+3) = extensive heat conduction

c     Compute lambda.Grad(T)
      dxInv = 1.d0 / dx(1)
      dyInv = 1.d0 / dx(2)
      dzInv = 1.d0 / dx(3)

!$omp parallel do private(i,j,k)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               Fx(i,j,k,Nspec+3) = - rhoDx(i,j,k,Nspec+2)*(T(i,j,k) - T(i-1,j,k))* dxInv * Ax(i,j,k)
            enddo
         enddo
      enddo
!$omp end parallel do

!$omp parallel do private(i,j,k)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = - rhoDy(i,j,k,Nspec+2)*(T(i,j,k) - T(i,j-1,k)) * dyInv * Ay(i,j,k)
            enddo
         enddo
      enddo
!$omp end parallel do

!$omp parallel do private(i,j,k)
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = - rhoDz(i,j,k,Nspec+2)*(T(i,j,k) - T(i,j,k-1)) * dzInv * Az(i,j,k)
            enddo
         enddo
      enddo
!$omp end parallel do

c     xlo
      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               Fx(i,j,k,Nspec+3) = 2*Fx(i,j,k,Nspec+3)
            enddo
         enddo
      endif
c     xhi
      if (hi(1).eq.dhi(1)  .and.  Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               Fx(i,j,k,Nspec+3) = 2*Fx(i,j,k,Nspec+3)
            enddo
         enddo
      endif
c     ylo
      if (lo(2).eq.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
         j=lo(2)
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = 2*Fy(i,j,k,Nspec+3)
            enddo
         enddo
      endif
c     yhi
      if (hi(2).eq.dhi(2) .and. Tbc(2,2).eq.EXT_DIR) then
         j=hi(2)+1
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               Fy(i,j,k,Nspec+3) = 2*Fy(i,j,k,Nspec+3)
            enddo
         enddo
      endif
c     zlo
      if (lo(3).eq.dlo(3) .and. Tbc(3,1).eq.EXT_DIR) then
         k=lo(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = 2*Fz(i,j,k,Nspec+3)
            enddo
         enddo
      endif
c     zhi
      if (hi(3).eq.dhi(3) .and. Tbc(3,2).eq.EXT_DIR) then
         k=hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Fz(i,j,k,Nspec+3) = 2*Fz(i,j,k,Nspec+3)
            enddo
         enddo
      endif

c     Compute enthalpy flux as hi*(Fi+(lambda/cp).Grad(Yi))

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

!$omp parallel do private(i,j,k)
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               rhoInv(i,j,k) = 1.0D0/rhoInv(i,j,k)
            enddo
         enddo
      enddo
!$omp end parallel do

!$omp parallel private(i,j,k,n,gradY)
      do n=1,Nspec
!$omp do
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
         Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2)
     &        + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*(H(i,j,k,n)+H(i-1,j,k,n))*0.5d0
      enddo
      enddo
      enddo
!$omp end do
      enddo
!$omp end parallel

!$omp parallel private(i,j,k,n,gradY)
      do n=1,Nspec
!$omp do
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
         Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2)
     &        + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*(H(i,j,k,n)+H(i,j-1,k,n))*0.5d0
      enddo
      enddo
      enddo
!$omp end do
      enddo
!$omp end parallel

!$omp parallel private(i,j,k,n,gradY)
      do n=1,Nspec
!$omp do
      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         gradY = (RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
         Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2)
     &        + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*(H(i,j,k,n)+H(i,j,k-1,n))*0.5d0
      enddo
      enddo
      enddo
!$omp end do
      enddo
!$omp end parallel

c     xlo
      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         Fx(i,lo(2):hi(2),lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
                  Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2)
     &                 + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*H(i-1,j,k,n)
               enddo
            enddo
         enddo
      endif
c     xhi
      if (hi(1).eq.dhi(1)  .and.  Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         Fx(i,lo(2):hi(2),lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i-1,j,k,n)*rhoInv(i-1,j,k))*dxInv
                  Fx(i,j,k,Nspec+2) = Fx(i,j,k,Nspec+2)
     &                 + (Fx(i,j,k,n) + rhoDx(i,j,k,Nspec+1)*gradY*Ax(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
c     ylo
      if (lo(2).eq.dlo(2)  .and.  Tbc(2,1).eq.EXT_DIR) then
         j = dlo(2)
         Fy(lo(1):hi(1),j,lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
                  Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2)
     &                 + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*H(i,j-1,k,n)
               enddo
            enddo
         enddo
      endif
c     yhi
      if (hi(2).eq.dhi(2)  .and.  Tbc(2,2).eq.EXT_DIR) then
         j = dhi(2)+1
         Fy(lo(1):hi(1),j,lo(3):hi(3),Nspec+2) = 0.d0
         do n=1,Nspec
            do k=lo(3),hi(3)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j-1,k,n)*rhoInv(i,j-1,k))*dyInv
                  Fy(i,j,k,Nspec+2) = Fy(i,j,k,Nspec+2)
     &                 + (Fy(i,j,k,n) + rhoDy(i,j,k,Nspec+1)*gradY*Ay(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
c     zlo
      if (lo(3).eq.dlo(3)  .and.  Tbc(3,1).eq.EXT_DIR) then
         k = dlo(3)
         Fz(lo(1):hi(1),lo(2):hi(2),k,Nspec+2) = 0.d0
         do n=1,Nspec
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
                  Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2)
     &                 + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*H(i,j,k-1,n)
               enddo
            enddo
         enddo
      endif
c     zhi
      if (hi(3).eq.dhi(3)  .and.  Tbc(3,2).eq.EXT_DIR) then
         k = dhi(3)+1
         Fz(lo(1):hi(1),lo(2):hi(2),k,Nspec+2) = 0.d0
         do n=1,Nspec
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  gradY = 2*(RhoY(i,j,k,n)*rhoInv(i,j,k) - RhoY(i,j,k-1,n)*rhoInv(i,j,k-1))*dzInv
                  Fz(i,j,k,Nspec+2) = Fz(i,j,k,Nspec+2)
     &                 + (Fz(i,j,k,n) + rhoDz(i,j,k,Nspec+1)*gradY*Az(i,j,k))*H(i,j,k,n)
               enddo
            enddo
         enddo
      endif
      deallocate(rhoInv)


c     Set FiGHi = (species flux) dot Grad(species enthalpy)
c        compute Grad(H) on each face, and average across faces in each coordinate
c        Fi is extensive here, so need to remove area.  Also, assume that we are 
c        away from domain boundary, fix afterward
c
c     FIXME: This will fail for r-z since Ax(dlo(1),:)=0
c
      !
      ! Use AD to cut down on some divides.
      ! Unfortunately it ups memory usage a bit.
      !
      allocate(AD(6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

!$omp parallel do private(i,j,k)
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
!$omp end parallel do

      FiGHi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.0D0

!$omp parallel private(i,j,k,n)
      do n=1,Nspec
!$omp do
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &        ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AD(2,i,j,k)
     &        + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AD(1,i,j,k)
     &        + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AD(4,i,j,k)
     &        + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AD(3,i,j,k)
     &        + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AD(6,i,j,k)
     &        + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AD(5,i,j,k) )
      enddo
      enddo
      enddo
!$omp end do
      enddo
!$omp end parallel

      deallocate(AD)

c     xlo
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

c     xhi
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

c     ylo
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

c     yhi
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

c     zlo
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

c     zhi
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
            FiGHi(i,j,k) = FiGHi(i,j,k) - 0.5d0*
     &           ( ( H(i+1,j,k,n) - H(i,j,k,n) )*Fx(i+1,j  ,k  ,n)*AxDxInv_hi
     &           + ( H(i,j,k,n) - H(i-1,j,k,n) )*Fx(i  ,j  ,k  ,n)*AxDxInv_lo
     &           + ( H(i,j+1,k,n) - H(i,j,k,n) )*Fy(i  ,j+1,k  ,n)*AyDyInv_hi
     &           + ( H(i,j,k,n) - H(i,j-1,k,n) )*Fy(i  ,j  ,k  ,n)*AyDyInv_lo 
     &           + ( H(i,j,k+1,n) - H(i,j,k,n) )*Fz(i  ,j  ,k+1,n)*AzDzInv_hi
     &           + ( H(i,j,k,n) - H(i,j,k-1,n) )*Fz(i  ,j  ,k  ,n)*AzDzInv_lo )
         enddo
         enddo
         enddo
      endif

      deallocate(H)
      end

      subroutine FORT_COMPUTE_RHODGRADHDOTGRADY(dx,
     &         lo, hi, DIMS(species), species,
     &         DIMS(h), h, DIMS(betax), betax,
     &         DIMS(betay), betay, DIMS(betaz), betaz, 
     &         DIMS(rdghdgy), rdghdgy)

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

!$omp parallel do private(i,j,k,bdotxlo,bdotxhi,bdotylo)
!$omp&private(bdotyhi,bdotzlo,bdotzhi)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            bdotxlo = betax(i,j,k  )*(h(i,j,k)-h(i-1,j,k))
     &            *(species(i,j,k)-species(i-1,j,k))
            bdotxhi = betax(i+1,j,k)*(h(i+1,j,k)-h(i,j,k))
     &           *(species(i+1,j,k)-species(i,j,k))
            bdotylo = betay(i,j,k  )*(h(i,j,k)-h(i,j-1,k))
     &           *(species(i,j,k)-species(i,j-1,k))
            bdotyhi = betay(i,j+1,k)*(h(i,j+1,k)-h(i,j,k))
     &           *(species(i,j+1,k)-species(i,j,k))
            bdotzlo = betaz(i,j,k  )*(h(i,j,k)-h(i,j,k-1))
     &           *(species(i,j,k)-species(i,j,k-1))
            bdotzhi = betaz(i,j,k+1)*(h(i,j,k+1)-h(i,j,k))
     &           *(species(i,j,k+1)-species(i,j,k))
            rdghdgy(i,j,k) =  half*((bdotxlo + bdotxhi)*dxsqr +
     &                              (bdotylo + bdotyhi)*dysqr +
     &                              (bdotzlo + bdotzhi)*dysqr)
          enddo
        enddo
      enddo
!$omp end parallel do

      end

      subroutine FORT_VELVISC(lo,hi,
     &                        DIMS(T), T,
     &                        DIMS(Y), Y,
     &                        DIMS(mu), mu)
      implicit none
#include <visc.H>
#include <cdwrk.H>
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T  T(DIMV(T))      
      REAL_T  Y(DIMV(Y),*)      
      integer DIMDEC(mu)
      REAL_T  mu(DIMV(mu))
      
      integer i, j, k

      if (.not.use_constant_mu) then
         if (.not. LeEQ1) then
            call FORT_MIX_SHEAR_VISC(lo, hi, mu, DIMS(mu),
     &                               T, DIMS(T), Y, DIMS(Y))       
         else
!$omp parallel do private(i,j,k)
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     mu(i,j,k) = 1.85D-5*(T(i,j,k)/298.D0)**0.7D0
                  end do
               end do
            end do
!$omp end parallel do
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
      end

      subroutine FORT_SPECTEMPVISC(lo,hi,
     &     DIMS(T), T,
     &     DIMS(Y), Y,
     &     DIMS(rhoD), rhoD,
     &     ncompd, P1ATM_MKS, do_temp, do_VelVisc,
     &     Pamb_in)
      implicit none
#include <cdwrk.H>
#include <visc.H>
#include <htdata.H>
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(rhoD)
      integer ncompd, do_temp, do_VelVisc
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),*)
      REAL_T  rhoD(DIMV(rhoD),ncompd)
      REAL_T  Pamb_in
      
      integer i, j, k, n, nc, ncs
      REAL_T P1ATM_MKS, Patm
      REAL_T Yt(maxspec)
      REAL_T cpmix, Tfac, Yfac
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

c     Nspec+1 is Temp stuff, if requested
c     Nspec+2 is Velocity stuff, if requested
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

c     Warning, FORT_VELVISC is called separately from this routine, so if there's
c     any hacking to be done on viscosity, be sure to catch it there as well.
      Tfac = thickFacTR / Pr
      Yfac = thickFacTR / Sc

      if (.not.use_constant_rhoD) then
         if (.not. LeEQ1) then
            Patm = Pamb_in / P1ATM_MKS
            call FORT_MIXAVG_RHODIFF_TEMP(lo, hi, rhoD, DIMS(rhoD),
     &           T, DIMS(T), Y, DIMS(Y), Patm, do_temp, do_VelVisc)
            if (thickFacTR.ne.1.d0) then
!$omp parallel do private(i,j,k,n)
               do n=1,ncs
                  do k=lo(3), hi(3)
                     do j=lo(2), hi(2)
                        do i=lo(1), hi(1)
                           rhoD(i,j,k,n) = rhoD(i,j,k,n)*thickFacTR
                        end do
                     end do
                  end do
               end do
!$omp end parallel do
            endif
         else
            call FORT_VELVISC(lo,hi,DIMS(T),T,DIMS(Y),Y,
     &           DIMS(rhoD),rhoD)
!$omp parallel do private(i,j,k,n)
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
!$omp end parallel do

!$omp parallel do private(i,j,k,n)
            do n=2,Nspec
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        rhoD(i,j,k,n) = rhoD(i,j,k,1)
                     end do
                  end do
               end do
            end do
!$omp end parallel do

            if (do_temp .ne. 0) then
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        do n=1,Nspec
                           Yt(n) = Y(i,j,k,n)
                        end do
                        CALL FORT_CPMIXfromTY(lo_chem, hi_chem, 
     &                       cpmix,  ARLIM(lo_chem), ARLIM(hi_chem),
     &                       T(i,j,k), ARLIM(lo_chem), ARLIM(hi_chem),
     &                       Yt,     ARLIM(lo_chem), ARLIM(hi_chem))
                        rhoD(i,j,k,Nspec+1) = rhoD(i,j,k,Nspec+1)*cpmix*Tfac
                     end do
                  end do
               end do
            end if
         end if
      else
         do n=1,Nspec
            do k=lo(3), hi(3)
               do j=lo(2), hi(2)
                  do i=lo(1), hi(1)
                     rhoD(i,j,k,n) = constant_rhoD_val
                  end do
               end do
            end do
         end do
      end if
      end

      subroutine FORT_EST_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           areaz,  DIMS(areaz),
     &                           lo, hi, dt, rhomin)

      implicit none

      integer flag
      integer lo(SDIM), hi(SDIM)
      REAL_T  delta(SDIM)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      REAL_T  rho(DIMV(rho))      
      REAL_T  u(DIMV(u),BL_SPACEDIM)      
      REAL_T  rhomin, dtfactor
      REAL_T  divu(DIMV(divu))
      REAL_T  dsdt(DIMV(divu))

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

!$omp parallel do reduction(min : dt) private(i,j,k,a,b,c)
!$omp&private(dtcell,dtcell2,denom,rhominij,rhoij)
!$omp&private(fluxxlo,fluxxhi,fluxylo,fluxyhi,fluxzlo,fluxzhi)
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
c
c     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c     We use the latter because it is more robust
c
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
                  denom = rho(i,j,k)*divu(i,j,k)+
     &                 u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) +
     &                 u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) +
     &                 u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = dtfactor*(rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = dtfactor*abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k))
     &                              *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k))
     &                              *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k))
     &                              *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k))
     &                              *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k-1))
     &                              *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k+1))
     &                              *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+
     &                     (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+
     &                     (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo))
     &                      /volume(i,j,k)
                  
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
!$omp end parallel do

      end

      subroutine FORT_CHECK_DIVU_DT(flag, dtfactor, delta, divu, DIMS(divu), 
     &                           dsdt, rho, DIMS(rho), 
     &                           u, DIMS(u), 
     &                           volume, DIMS(volume), 
     &                           areax,  DIMS(areax),
     &                           areay,  DIMS(areay),
     &                           areaz,  DIMS(areaz),
     &                           lo, hi,
     &                           dt, rhomin)

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

!$omp parallel do private(i,j,k,a,b,c)
!$omp&private(dtcell,dtcell2,denom,rhominij,rhoij)
!$omp&private(fluxxlo,fluxxhi,fluxylo,fluxyhi,fluxzlo,fluxzhi)
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
c
c     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
c     We use the latter because it is more robust
c
                        a = rhoij*dsdt(i,j,k)*half
                        b = rhoij*divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))

                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
               else if (flag.eq.2) then
                  denom = rho(i,j,k)*divu(i,j,k)+
     &                 u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1) +
     &                 u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2) +
     &                 u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               else if (flag.eq.3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k))
     &                              *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k))
     &                              *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k))
     &                              *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k))
     &                              *(u(i,j,k,2)+u(i,j+1,k,2))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k-1))
     &                              *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k+1))
     &                              *(u(i,j,k,3)+u(i,j,k+1,3))
                  denom = ((areax(i+1,j,k)*fluxxhi-areax(i,j,k)*fluxxlo)+
     &                     (areay(i,j+1,k)*fluxyhi-areay(i,j,k)*fluxylo)+
     &                     (areaz(i,j,k+1)*fluxzhi-areaz(i,j,k)*fluxzlo))
     &                      /volume(i,j,k)
                  if(denom.gt.zero)then
                     if(rho(i,j,k).gt.rhomin) then
                        dtcell = (rho(i,j,k)-rhomin)/denom
                     else
                        dtcell = abs(rho(i,j,k))/denom
                     endif
                  endif
               endif
               if (dt.gt.dtcell) then
                  write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell = ',
     &                 i,j,k,dt,dtcell
               else if (dt.gt.dtcell*dtfactor) then
                  write(6,*)'WARNING: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell*dtfactor = ',
     &                 i,j,k,dt,dtcell*dtfactor
               endif
            enddo
         enddo
      enddo
!$omp end parallel do

      end

c===================================================================

c ... The routines YDOTFILL, DIVUFILL, DQRADFILL, and DSDTFILL
c     are here instead of in the problem dependent code because
c     we always fill the quantitities ydot, divu, dqrad, and dsdt
c     the same way

c===================================================================

      subroutine FORT_DQRADFILL (dqrad,DIMS(dqrad),domlo,domhi,delta,
     &                         xlo,time,bc )

      integer    DIMDEC(dqrad)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dqrad(DIMV(dqrad))

      call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)
      call fillEdges(dqrad,DIMS(dqrad),domlo,domhi,bc)
      end

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,delta,
     &                        xlo,time,bc )

      integer    DIMDEC(divu)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      call fillEdges(divu,DIMS(divu),domlo,domhi,bc)

      end

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,delta,
     &                         xlo,time,bc )

      integer    DIMDEC(dsdt)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(dsdt,DIMS(dsdt),domlo,domhi,bc)

      end

      subroutine FORT_YDOTFILL (ydot,DIMS(ydot),domlo,domhi,delta,
     &                          xlo,time,bc)

      integer    DIMDEC(ydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     ydot(DIMV(ydot))

      call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(ydot,DIMS(ydot),domlo,domhi,bc)

      end

      subroutine FORT_RHOYDOTFILL (rhoydot,DIMS(rhoydot),domlo,domhi,delta,
     &                             xlo,time,bc)

      integer    DIMDEC(rhoydot), bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     rhoydot(DIMV(rhoydot))

      call filcc (rhoydot,DIMS(rhoydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(rhoydot,DIMS(rhoydot),domlo,domhi,bc)

      end

      subroutine FORT_FABMINMAX(lo, hi, 
     &                          fab, DIMS(fab),
     &                          fmin, fmax, nc)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fab)
      REAL_T  fab(DIMV(fab),nc)
      REAL_T  fmin, fmax

      integer i,j,k,n

!$omp parallel private(i,j,k,n)
      do n = 1,nc
!$omp do
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  fab(i,j,k,n) = MAX( fmin, MIN( fmax, fab(i,j,k,n) ) )
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel
      end

      subroutine FORT_REPAIR_FLUX(lo, hi, dlo, dhi,
     &                            flux, DIMS(flux),
     &                            RhoY, DIMS(RhoY), dir, Ybc)
      implicit none

#include <cdwrk.H>

      integer lo(SDIM), hi(SDIM), dlo(SDIM), dhi(SDIM), dir, Ybc(SDIM,2)
      integer DIMDEC(flux)
      integer DIMDEC(RhoY)
      REAL_T flux(DIMV(flux),Nspec)
      REAL_T RhoY(DIMV(RhoY),Nspec)
      
      integer i, j, k, n
      REAL_T sumFlux, RhoYe(Nspec), sumRhoYe

      if (dir.eq.0) then

c     First, assume away from physical boundaries, then use boundary-aware version below if applicable

!$omp parallel do private(i,j,k,n,sumFlux,RhoYe,sumRhoYe)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
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
!$omp end parallel do
c     xlo
         if (Ybc(1,1).eq.EXT_DIR) then
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
c     xhi
         if (Ybc(1,2).eq.EXT_DIR) then
            do i = dhi(1)+1,hi(1)+1
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

c     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable

!$omp parallel do private(i,j,k,n,sumFlux,RhoYe,sumRhoYe)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
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
!$omp end parallel do
c     ylo
         if (Ybc(2,1).eq.EXT_DIR) then
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
c     yhi
         if (Ybc(2,2).eq.EXT_DIR) then
            do j = dhi(2)+1,hi(2)+1
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

c     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable
!$omp parallel do private(i,j,k,n,sumFlux,RhoYe,sumRhoYe)
         do k = lo(3),hi(3)+1
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
!$omp end parallel do
c     zlo
         if (Ybc(3,1).eq.EXT_DIR) then
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
c     yzi
         if (Ybc(3,2).eq.EXT_DIR) then
            do k = dhi(3)+1,hi(3)+1
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
      end

      subroutine  FORT_INCRWEXTFLXDIV(lo, hi,
     &                                xflux,  DIMS(xflux),
     &                                yflux,  DIMS(yflux),
     &                                zflux,  DIMS(zflux),
     &                                stateo, DIMS(stateo),
     &                                staten, DIMS(staten),
     &                                vol,    DIMS(vol),
     &                                nc, dt)
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

!$omp parallel do private(i,j,k,n,dF)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               dF = zero
               do n=1,nc
                  dF = dF + ( xflux(i+1,j,k,n) - xflux(i,j,k,n) )
     &                 +    ( yflux(i,j+1,k,n) - yflux(i,j,k,n) )
     &                 +    ( zflux(i,j,k+1,n) - zflux(i,j,k,n) )
               end do
               staten(i,j,k) = stateo(i,j,k) + dF*dt/vol(i,j,k)
            end do
         end do
      end do
!$omp end parallel do
      end

      subroutine  FORT_FLUXDIV(lo, hi,
     &                         update, DIMS(update),
     &                         xflux,  DIMS(xflux),
     &                         yflux,  DIMS(yflux),
     &                         zflux,  DIMS(zflux),
     &                         vol,    DIMS(vol),
     &                         nc, scal)
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

!$omp parallel private(i,j,k,n)

!$omp do
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ivol(i,j,k) = scal / vol(i,j,k)
            end do
         end do
      end do
!$omp end do

      do n=1,nc
!$omp do
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  update(i,j,k,n) = 
     &                 ( (xflux(i+1,j,k,n)-xflux(i,j,k,n))
     &                 + (yflux(i,j+1,k,n)-yflux(i,j,k,n)) 
     &                 + (zflux(i,j,k+1,n)-zflux(i,j,k,n)) ) * ivol(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do

!$omp end parallel

      end

      subroutine fillEdges(dat,DIMS(dat),domlo,domhi,bc )
c  This routine fills ghost cells with the value from the nearest
c  interior cell.  
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
      end

      subroutine fillWithZeros(dat,DIMS(dat),domlo,domhi,bc)

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
      end

      subroutine  FORT_COMPUTE_UGRADP(p, DIMS(p), ugradp, DIMS(ugp),
     &                                umac,  DIMS(umac),
     &                                vmac,  DIMS(vmac),
     &                                wmac,  DIMS(wmac),
     &                                lo, hi, dx)

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

!$omp parallel do private(i,j,k,uadv,vadv,wadv)
!$omp&private(p_x_lo,p_x_hi,p_y_lo,p_y_hi,p_z_lo,p_z_hi)
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
               ugradp(i,j,k) = uadv * (p_x_hi - p_x_lo) / dx(1) +
     &              vadv * (p_y_hi - p_y_lo) / dx(2) +    
     &              wadv * (p_z_hi - p_z_lo) / dx(3)
            end do
         end do
      end do
!$omp end parallel do

      end

      integer function FORT_CONSERVATIVE_T_FLOOR(
     &     loF,hiF,state,DIMS(state),
     &     min_T, Tcomp, Rcomp, first_spec, last_spec, RhoH,
     &     ratio, tmp, nt)
      implicit none

      integer loF(SDIM), hiF(SDIM)
      integer DIMDEC(state)
      REAL_T  state(DIMV(state),0:*)      
      integer Tcomp, Rcomp, first_spec, last_spec, RhoH, ratio(SDIM), nt
      REAL_T  min_T, tmp(0:nt-1)
      integer n,i,j,k, loC(SDIM),hiC(SDIM),ii,jj,kk,iii,jjj,kkk,ncells
      Real ncellsInv
      logical bad_T

c     Returns the number of fine cells fixed up
      FORT_CONSERVATIVE_T_FLOOR = 0

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

                  FORT_CONSERVATIVE_T_FLOOR = FORT_CONSERVATIVE_T_FLOOR + ncells
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
      end

      subroutine FORT_HTDDRELAX(lo, hi, S, DIMS(S), yc, Tc, hc, rc,
     &     L, DIMS(L), a, DIMS(a), R, DIMS(R), thetaDt, fac, maxRes, maxCor,
     &     for_T0_H1, res_only, mult)
      implicit none
#include <cdwrk.H>      
      integer lo(SDIM), hi(SDIM), yc, Tc, hc, rc
      integer DIMDEC(S)
      integer DIMDEC(L)
      integer DIMDEC(a)
      integer DIMDEC(R)
      REAL_T S(DIMV(s),0:*)
      REAL_T L(DIMV(L),0:*)
      REAL_T a(DIMV(a),0:*)
      REAL_T R(DIMV(R),0:*)
      REAL_T thetaDt, fac(0:*), maxRes(0:*), maxCor(0:*)
      integer for_T0_H1, res_only
      REAL_T mult, cor
      integer i, j, k, n

#ifdef HT_SKIP_NITROGEN
      REAL_T, allocatable :: tfab(:,:,:)
#endif

c
c     Helper function to compute:
c     (a) dR = Rhs + (rho.phi - theta.dt.Lphi)
c     (b) Res = Rhs - A(S) = Rhs - (rho.phi - theta.dt.Lphi) for RhoH
c         or  Res = Rhs - (phi - theta.dt.Lphi) for Temp
c               and then relax: phi = phi + Res/alpha
c     -- note that if for Temp both Lphi and alpha have been scaled by (1/rho.Cp)^nph
c     
c     NOTE: Assumes maxRes and maxCor have been initialized properly
c

      do n=0,Nspec-1
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  L(i,j,k,n) = R(i,j,k,n) + mult*(S(i,j,k,rc)*S(i,j,k,yc+n) - thetaDt*L(i,j,k,n))
                  maxRes(n) = MAX(maxRes(n),ABS(L(i,j,k,n)))
               enddo
            enddo
         enddo
      enddo

#ifdef HT_SKIP_NITROGEN
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               L(i,j,k,Nspec-1) = 0.d0
            enddo
         enddo
      enddo
      maxRes(Nspec-1) = 0.d0
#endif

      if (for_T0_H1.eq.0) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  L(i,j,k,Nspec) = R(i,j,k,Nspec) + mult*(S(i,j,k,Tc) - thetaDt*L(i,j,k,Nspec))
                  maxRes(Nspec) = MAX(maxRes(Nspec),ABS(L(i,j,k,Nspec)))
               enddo
            enddo
         enddo
      else
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  L(i,j,k,Nspec) = R(i,j,k,Nspec) + mult*(S(i,j,k,rc)*S(i,j,k,hc) - thetaDt*L(i,j,k,Nspec))
                  maxRes(Nspec) = MAX(maxRes(Nspec),ABS(L(i,j,k,Nspec)))
               enddo
            enddo
         enddo
      endif

      if (res_only.ne.1) then
         do n=0,Nspec-1
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     cor = fac(n) * L(i,j,k,n) / (S(i,j,k,rc) - thetaDt*a(i,j,k,n))
                     maxCor(n) = MAX(maxCor(n),ABS(cor))
                     S(i,j,k,yc+n) = S(i,j,k,yc+n)  +  cor
                  enddo
               enddo
            enddo
         enddo
#ifdef HT_SKIP_NITROGEN
         allocate(tfab(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         tfab = zero
         do n=0,Nspec-2
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     tfab(i,j,k) = tfab(i,j,k) + S(i,j,k,yc+n)
                  enddo
               enddo
            enddo
         end do
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  S(i,j,k,Nspec-1) = 1.d0 - tfab(i,j,k)
               enddo
            enddo
         enddo
         maxCor(Nspec-1) = 0.d0
         deallocate(tfab)
#endif

         if (for_T0_H1.eq.0) then
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     cor = fac(Nspec) * L(i,j,k,Nspec) / (1.d0 - thetaDt*a(i,j,k,Nspec))
                     maxCor(Nspec) = MAX(maxCor(Nspec),ABS(cor))
                     S(i,j,k,Tc) = S(i,j,k,Tc)  +  cor
                  enddo
               enddo
            enddo
         else
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     cor = fac(Nspec) * L(i,j,k,Nspec) / (S(i,j,k,rc) - thetaDt*a(i,j,k,Nspec))
                     maxCor(Nspec) = MAX(maxCor(Nspec),ABS(cor))
                     S(i,j,k,hc) = S(i,j,k,hc)  +  cor
                  enddo
               enddo
            enddo
         endif
      endif
      end

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
      subroutine FORT_PART_CNT_ERR(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3,
     &     set,clear,
     &     var,varl1,varl2,varl3,varh1,varh2,varh3,
     &     lo,hi,nd,domlo,domhi,
     &     delta,xlo,problo,time,level)

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

      end subroutine

      subroutine FORT_MCURVE(lo, hi, T, DIMS(T), curv, DIMS(curv),
     &     wrk, DIMS(wrk), delta)

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

c     Fill normal on nodes (assumes 1 grow cell properly filled)

      REAL_T Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Txz,tdxI(3),dxSqI(3)

      do i=1,3
         tdxI(i) =  one / (two*delta(i))
         dxSqI(i) = one / (delta(i)*delta(i))
      enddo

!$omp parallel do private(i,j,k,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Txz,Tyz,mag)
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

               curv(i,j,k) = -half*(Txx + Tyy + Tzz
     &              - ( Tx*(Tx*Txx + Ty*Txy + Tz*Txz)
     &                + Ty*(Tx*Txy + Ty*Tyy + Tz*Tyz)
     &                + Tz*(Tx*Txz + Ty*Tyz + Tz*Tzz) )/mag**2 )/mag
            end do
         end do
      end do
!$omp end parallel do

      end

      subroutine FORT_SMOOTH(lo, hi, Tin, DIMS(Tin), Tout, DIMS(Tout))
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Tin)
      integer DIMDEC(Tout)
      REAL_T   Tin(DIMV(Tin))
      REAL_T  Tout(DIMV(Tout))

      integer i,j,k,ii,jj,kk

!$omp parallel do private(i,j,k,ii,jj,kk)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Tout(i,j,k) = zero
               do kk=0,1
                  do jj=0,1
                     do ii=0,1
                        Tout(i,j,k) = Tout(i,j,k)
     &                       + Tin(i+ii,j+jj,  k+kk-1) + Tin(i+ii-1,j+jj,  k+kk-1)
     &                       + Tin(i+ii,j+jj-1,k+kk-1) + Tin(i+ii-1,j+jj-1,k+kk-1)
     &                       + Tin(i+ii,j+jj,  k+kk)   + Tin(i+ii-1,j+jj,  k+kk  )
     &                       + Tin(i+ii,j+jj-1,k+kk  ) + Tin(i+ii-1,j+jj-1,k+kk  )
                     end do
                  end do
               end do
               Tout(i,j,k) = Tout(i,j,k) / 64.d0
            end do
         end do
      end do
!$omp end parallel do

      end

      subroutine FORT_GRADWBAR(lo, hi, Wbar, DIMS(Wbar),
     &     rDe, DIMS(rDe), flux, DIMS(flux),
     &     area, DIMS(area), dx, dir, mult, inc)
      implicit none
#include <cdwrk.H>
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
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else if (dir.eq.1) then
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else if (dir.eq.2) then
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else
            call bl_pd_abort('Bad dir in FORT_GRADWBAR')
         endif

      else

!     increment grad wbar fluxes by a factor of inc (can be negative)

         if (dir.eq.0) then
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else if (dir.eq.1) then
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else if (dir.eq.2) then
!     $omp parallel do private(i,j,k,Wgr)
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
!     $omp end parallel do
         else
            call bl_pd_abort('Bad dir in FORT_GRADWBAR')
         endif

      end if

      end

      subroutine  FORT_RECOMP_UPDATE(lo, hi,
     &                               update, DIMS(update),
     &                               xflux,  DIMS(xflux),
     &                               yflux,  DIMS(yflux),
     &                               zflux,  DIMS(zflux),
     &                               vol,    DIMS(vol),
     &                               nc)
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

!$omp parallel private(i,j,k,n)
      do n=1,nc
!$omp do
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  update(i,j,k,n)=-((xflux(i+1,j,k,n)-xflux(i,j,k,n))
     &                 +            (yflux(i,j+1,k,n)-yflux(i,j,k,n))
     &                 +            (zflux(i,j,k+1,n)-zflux(i,j,k,n)))
     &                 /vol(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel
      end
