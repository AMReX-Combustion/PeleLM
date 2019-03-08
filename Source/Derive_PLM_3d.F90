#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <Prob_F.H>
#include <AMReX_ArrayLim.H>
#include <ChemDriver_F.H>

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 3

module derive_PLM_3D

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

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHOMINUSSUMRHOY: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1)
            enddo
         enddo
      enddo

      do n=2,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)-dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

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

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n
      integer    nxlo, nxhi, nylo, nyhi, nzlo,nzhi
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERSUMYDOT: outside domain")
      endif

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = 0.0d0
            enddo
         enddo
      enddo
      
      do n=1,ncomp
         do k=lo(3),hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  e(i,j,k,1) = e(i,j,k,1)+dat(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

  end subroutine dsrhoydot

!=========================================================

  subroutine drhort (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                     lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                     level,grid_no) &
                     bind(C, name="drhort")
                     
      use chem_driver_3D, only: PfromRTY
      implicit none
!     
! ::: This routine will derive rho*R*T
!

#include <cdwrk.H>
#include <htdata.H>
      
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i, j, k, n, rho, T, fS
      integer    nxlo,nxhi,nylo,nyhi,nzlo,nzhi
      REAL_T     Yt(maxspec)
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/
      
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
	 call bl_abort("FORT_DERRHORT: outside domain")
      endif
!      
!     Set pointers into state (these must agree with setup for this derived quant).
!
      rho = 1
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1) / dat(i,j,k,rho)
               end do
               call PfromRTY(lo_chem, hi_chem, &
                   e(i,j,k,1),     ARLIM(lo_chem),ARLIM(hi_chem), &
                   dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
                   dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
                   Yt,             ARLIM(lo_chem),ARLIM(hi_chem))
            end do
         end do
      end do
  end subroutine drhort

!=========================================================

  subroutine dermolefrac(x,DIMS(x),nv,dat,DIMS(dat),ncomp, &
                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                         level,grid_no) &
                         bind(C, name="dermolefrac")
                          
      use chem_driver_3D, only : mass_to_mole
      implicit none

#include <cdwrk.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(x)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     x(DIMV(x),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(maxspec),Xt(maxspec)
      integer fS,rho
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      fS  = 2

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
               call mass_to_mole(lo_chem, hi_chem, &
                          Yt, ARLIM(lo_chem),ARLIM(hi_chem), &
                          Xt, ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  x(i,j,k,n) = Xt(n)
               enddo
            enddo
         enddo
      enddo


  end subroutine dermolefrac

!=========================================================
  
  subroutine derconcentration (C,DIMS(C),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                              level,grid_no) &
                               bind(C, name="derconcentration")
                               
      use chem_driver_3D, only: MASSR_TO_CONC
      implicit none

#include <cdwrk.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(C)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     C(DIMV(C),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer i,j,k,n
      REAL_T Yt(maxspec),Ct(maxspec)
      integer fS,rho,T
      integer lo_chem(SDIM),hi_chem(SDIM)
      data lo_chem /1,1,1/
      data hi_chem /1,1,1/

      rho = 1 
      T   = 2
      fS  = 3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = dat(i,j,k,fS+n-1)/dat(i,j,k,rho) 
               enddo
               call MASSR_TO_CONC(lo_chem,hi_chem, &
                Yt,             ARLIM(lo_chem),ARLIM(hi_chem), &
                dat(i,j,k,T),   ARLIM(lo_chem),ARLIM(hi_chem), &
                dat(i,j,k,rho), ARLIM(lo_chem),ARLIM(hi_chem), &
                Ct,             ARLIM(lo_chem),ARLIM(hi_chem))
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)
               enddo
            enddo
         enddo
      enddo


  end subroutine derconcentration

!=========================================================

  subroutine FORT_DERFORCING (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                              level,grid_no) &
                              bind(C, name="FORT_DERFORCING")
     
      implicit none

!
! ::: This routine will computes the forcing term
!

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>

      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1, f2, f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      REAL_T  rho, u, v, w
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
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
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
                  rho = dat(i,j,k,1)
                  u   = dat(i,j,k,2)
                  v   = dat(i,j,k,3)
                  w   = dat(i,j,k,4)
                  e(i,j,k,1) = f1*rho*u
                  e(i,j,k,2) = f2*rho*v
                  e(i,j,k,3) = f3*rho*w
               else
                  e(i,j,k,1) = f1*u
                  e(i,j,k,2) = f2*v
                  e(i,j,k,3) = f3*w
               endif
            enddo
         enddo
      enddo
#endif
  end subroutine FORT_DERFORCING 

!=========================================================  
  
  subroutine FORT_DERFORCEX (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEX")
                             
      implicit none

!
! ::: This routine will computes the forcing term
!

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
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
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                              f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                  -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                           else
                              f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
                  e(i,j,k,1) = f1*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f1
               endif
            enddo
         enddo
      enddo
#endif
  end subroutine FORT_DERFORCEX

!=========================================================

  subroutine FORT_DERFORCEY (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEY")
                             
      implicit none
!
! ::: This routine will computes the forcing term
!

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f2
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
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
               f2 = zero
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
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                              f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                  -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                           else
                              f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
                  e(i,j,k,1) = f2*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f2
               endif
            enddo
         enddo
      enddo
#endif
  end subroutine FORT_DERFORCEY

!=========================================================

  subroutine FORT_DERFORCEZ (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                             level,grid_no) &
                             bind(C, name="FORT_DERFORCEZ")
                             
      implicit none

!
! ::: This routine will computes the forcing term
!

#include <cdwrk.H>
#include <htdata.H>

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

#ifdef DO_LMC_FORCE
#include <probdata.H>
#include <forcedata.H>
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  HLx, HLy, HLz
      REAL_T  Lx, Ly, Lz, Lmin, kappa, kappaMax
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc
      integer nXvel, nYvel, nZvel, nRho, nTrac

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      ilo = lo(1)
      jlo = lo(2)
      klo = lo(3)
      ihi = hi(1)
      jhi = hi(2)
      khi = hi(3)

!     Homogeneous Isotropic Turbulence
      twicePi=two*Pi
      
!     Adjust z offset for probtype 15
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
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
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
                              f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                  -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                           else
                              f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                           endif
                        endif
                     enddo
                  enddo
               enddo
               if (use_rho_in_forcing.eq.1) then
                  e(i,j,k,1) = f3*dat(i,j,k,1)
               else
                  e(i,j,k,1) = f3
               endif
            enddo
         enddo
      enddo
#endif
  end subroutine FORT_DERFORCEZ

end module derive_PLM_3D