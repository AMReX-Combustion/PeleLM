
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

#define SDIM 2

module PeleLM_2d

  use fuego_chemistry

  implicit none

  private

  public ::  calc_divu_fortran, calc_gamma_pinv, floor_spec, enth_diff_terms, &
             compute_rho_dgrad_hdot_grad_Y, vel_visc, spec_temp_visc, &
             est_divu_dt, check_divu_dt, dqrad_fill, divu_fill, &
             dsdt_fill, ydot_fill, rhoYdot_fill, fab_minmax, repair_flux, &
             incrwext_flx_div, flux_div, compute_ugradp, conservative_T_floor, &
             part_cnt_err, mcurve, smooth, grad_wbar, recomp_update, &
             valgt_error, vallt_error, magvort_error, diffgt_error, &
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

! Note that this routine has been renamed calc_divu_fortran
! to not override the calc_divu already declared in PeleLM.cpp

    integer lo(SDIM),hi(SDIM)
    integer DIMDEC(divu)
    integer DIMDEC(rYdot)
    integer DIMDEC(vtY)
    integer DIMDEC(vtT)
    integer DIMDEC(rhoY)
    integer DIMDEC(T)
    REAL_T  divu(DIMV(divu))
    REAL_T  rYdot(DIMV(rYdot),1:nspec)
    REAL_T  vtY(DIMV(vtY),1:nspec)
    REAL_T  vtT(DIMV(vtT))
    REAL_T  rhoY(DIMV(rhoY),1:nspec)
    REAL_T  T(DIMV(T))
      
    integer i, j, n
    REAL_T Y(nspec), H(nspec), cpmix, rhoInv, tmp, mmw, invmtw(nspec)

    call CKWT(invmtw)

    do n=1,Nspec
       invmtw(n) = one / invmtw(n)
    end do

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        rhoInv = 0.d0
        do n=1,Nspec
          rhoInv = rhoInv + rhoY(i,j,n)
        enddo
        rhoInv = 1.d0 / rhoInv
        do n=1,Nspec
          Y(n) = rhoInv*rhoY(i,j,n)
        enddo
        CALL CKCPBS(T(i,j),Y,cpmix)
        CALL CKHMS(T(i,j),H)
        CALL CKMMWY(Y,mmw)

        cpmix = cpmix*1.d-4
        do n=1,Nspec
          H(n) = H(n)*1.d-4
        enddo

        tmp = 0.d0
        divu(i,j) = vtT(i,j)
        do n=1,Nspec
          tmp = tmp + (rYdot(i,j,n)+vtY(i,j,n))*invmtw(n)
          divu(i,j) = divu(i,j) - rYdot(i,j,n)*H(n)
         enddo
         divu(i,j) = ( divu(i,j)/(cpmix*T(i,j)) + tmp*mmw ) * rhoInv
        enddo
      enddo

  end subroutine calc_divu_fortran

!----------------------------------------- 

  subroutine calc_gamma_pinv(lo, hi, &
                             theta, DIMS(theta), &
                             rhoY, DIMS(rhoY), &
                             T, DIMS(T), &
                             Pamb_in) &
                             bind(C, name="calc_gamma_pinv") 

    use network,        only : nspec

    implicit none

    integer lo(SDIM),hi(SDIM)
    integer DIMDEC(theta)
    integer DIMDEC(rhoY)
    integer DIMDEC(T)
    REAL_T  theta(DIMV(theta))
    REAL_T  rhoY(DIMV(rhoY),1:Nspec)
    REAL_T  T(DIMV(T))
    REAL_T  Pamb_in
      
    integer i, j, n
    REAL_T Y(nspec), cpmix, cvmix, rhoInv

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        rhoInv = 0.d0
        do n=1,Nspec
          rhoInv = rhoInv + rhoY(i,j,n)
        enddo
        rhoInv = 1.d0 / rhoInv
        do n=1,Nspec
          Y(n) = rhoInv*rhoY(i,j,n)
        enddo

        CALL CKCPBS(T(i,j),Y,cpmix)
        cpmix = cpmix*1.d-4
        CALL CKCVBS(T(i,j),Y,cvmix)
        cvmix = cvmix*1.d-4

        theta(i,j) = cvmix / (cpmix*Pamb_in)
            
      enddo
    enddo

  end subroutine calc_gamma_pinv

!-------------------------------------

  subroutine floor_spec(lo, hi, spec,  DIMS(spec))bind(C, name="floor_spec")

    use network,        only : nspec

    implicit none

    integer lo(SDIM),hi(SDIM)
    integer DIMDEC(spec)
    REAL_T  spec(DIMV(spec),1:Nspec)
      
    integer i, j, n

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        do n=1,Nspec
          spec(i,j,n) = max(0.d0,spec(i,j,n))
        end do
      enddo
    enddo

  end subroutine floor_spec

!-----------------------------------

  subroutine enth_diff_terms (lo, hi, dlo, dhi, dx, &
                              T, DIMS(T), RhoY, DIMS(RhoY), &
                              rhoDx, DIMS(rhoDx), Fx, DIMS(Fx), Ax, DIMS(Ax), &
                              rhoDy, DIMS(rhoDy), Fy, DIMS(Fy), Ay, DIMS(Ay), &
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

    integer DIMDEC(FiGHi)
    REAL_T  FiGHi(DIMV(FiGHi))

    REAL_T, allocatable :: H(:,:,:)

    integer i, j, d, n
    integer lob(SDIM), hib(SDIM)
    REAL_T AxDxInv_lo, AxDxInv_hi, dxInv
    REAL_T AyDyInv_lo, AyDyInv_hi, dyInv
    logical fix_xlo, fix_xhi, fix_ylo, fix_yhi

    REAL_T, allocatable :: rhoInv(:,:)
    REAL_T gradY

    fix_xlo = .false.
    fix_xhi = .false.
    fix_ylo = .false.
    fix_yhi = .false.

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

!     Comute lambda.Grad(T)
    dxInv = 1.d0 / dx(1)
    dyInv = 1.d0 / dx(2)

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
        Fx(i,j,Nspec+3) = - rhoDx(i,j,Nspec+2)*(T(i,j) - T(i-1,j))* dxInv * Ax(i,j)
      enddo
    enddo
    do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
        Fy(i,j,Nspec+3) = - rhoDy(i,j,Nspec+2)*(T(i,j) - T(i,j-1)) * dyInv * Ay(i,j)
      enddo
    enddo

!     xlo
    if (lo(1).le.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
      i = dlo(1)
      do j=lo(2),hi(2)
        Fx(i,j,Nspec+3) = 2*Fx(i,j,Nspec+3)
      enddo
    endif
!     xhi
    if (hi(1)+1.ge.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
      i = dhi(1)+1
      do j=lo(2),hi(2)
        Fx(i,j,Nspec+3) = 2*Fx(i,j,Nspec+3)
      enddo
    endif
!     ylo
    if (lo(2).le.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
      j=lo(2)
      do i=lo(1),hi(1)
        Fy(i,j,Nspec+3) = 2*Fy(i,j,Nspec+3)
      enddo
    endif
!     yhi
    if (hi(2)+1.ge.dhi(2)+1 .and. Tbc(2,2).eq.EXT_DIR) then
      j=hi(2)+1
      do i=lo(1),hi(1)
        Fy(i,j,Nspec+3) = 2*Fy(i,j,Nspec+3)
      enddo
    endif

#if 0  
!     Compute enthalpy flux as Fi.hi.(Lei-1)

    Fx(lo(1):hi(1)+1,lo(2):hi(2),Nspec+2) = 0.d0
    Fy(lo(1):hi(1),lo(2):hi(2)+1,Nspec+2) = 0.d0

    do n=1,Nspec
      do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
          FiHi = + 0.5d0*(H(i,j,n)+H(i-1,j,n))*Fx(i,j,n)
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) + FiHi*(rhoDx(i,j,Nspec+1)/rhoDx(i,j,n) - 1.d0)
        enddo
      enddo
    enddo
    do n=1,Nspec
      do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
          FiHi = + 0.5d0*(H(i,j,n)+H(i,j-1,n))*Fy(i,j,n)
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) + FiHi*(rhoDy(i,j,Nspec+1)/rhoDy(i,j,n) - 1.d0)
        enddo
      enddo
    enddo
!     xlo
    if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
      i = dlo(1)
      Fx(i:i,lo(2):hi(2),Nspec+2) = 0.d0
      do n=1,Nspec
        do j=lo(2),hi(2)
          FiHi = H(i,j,n)*Fx(i,j,n)
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) + FiHi*(rhoDx(i,j,Nspec+1)/rhoDx(i,j,n) - 1.d0)
        enddo
      enddo
    endif
!     xhi
    if (hi(1)+1.eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
      i = dhi(1)+1
      Fx(i:i,lo(2):hi(2),Nspec+2) = 0.d0
      do n=1,Nspec
        do j=lo(2),hi(2)
          FiHi = H(i,j,n)*Fx(i,j,n)
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) + FiHi*(rhoDx(i,j,Nspec+1)/rhoDx(i,j,n) - 1.d0)
        enddo
      enddo
    endif
!     ylo
    if (lo(2).eq.dlo(2)  .and.  Tbc(2,1).eq.EXT_DIR) then
      j = dlo(2)
      Fy(lo(1):hi(1),j:j,Nspec+2) = 0.d0
      do n=1,Nspec
        do i=lo(1),hi(1)
          FiHi = H(i,j,n)*Fy(i,j,n)
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) + FiHi*(rhoDy(i,j,Nspec+1)/rhoDy(i,j,n) - 1.d0)
        enddo
      enddo
    endif
!     yhi
    if (hi(2)+1.eq.dhi(2)  .and.  Tbc(2,2).eq.EXT_DIR) then
      j = dhi(2)+1
      Fy(lo(1):hi(1),j:j,Nspec+2) = 0.d0
      do n=1,Nspec
        do i=lo(1),hi(1)
          FiHi = H(i,j,n)*Fy(i,j,n)
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) + FiHi*(rhoDy(i,j,Nspec+1)/rhoDy(i,j,n) - 1.d0)
        enddo
      enddo
    endif
#else

!     Compute enthalpy flux as hi*(Fi+(lambda/cp).Grad(Yi))

    Fx(lo(1):hi(1)+1,lo(2):hi(2),Nspec+2) = 0.d0
    Fy(lo(1):hi(1),lo(2):hi(2)+1,Nspec+2) = 0.d0

    allocate(rhoInv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    rhoInv = 0
    do n=1,Nspec
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
          rhoInv(i,j) = rhoInv(i,j) + RhoY(i,j,n)
        enddo
      enddo
    enddo
    rhoInv(:,:) = 1.d0/rhoInv(:,:)

    do n=1,Nspec
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
          gradY = (RhoY(i,j,n)*rhoInv(i,j) - RhoY(i-1,j,n)*rhoInv(i-1,j))*dxInv
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) &
                   + (Fx(i,j,n) + rhoDx(i,j,Nspec+1)*gradY*Ax(i,j))*(H(i,j,n)+H(i-1,j,n))*0.5d0
        enddo
      enddo
    enddo

    do n=1,Nspec
      do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
          gradY = (RhoY(i,j,n)*rhoInv(i,j) - RhoY(i,j-1,n)*rhoInv(i,j-1))*dyInv
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) &
                   + (Fy(i,j,n) + rhoDy(i,j,Nspec+1)*gradY*Ay(i,j))*(H(i,j,n)+H(i,j-1,n))*0.5d0
        enddo
      enddo
    enddo

!     xlo
    if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
      i = dlo(1)
      Fx(i:i,lo(2):hi(2),Nspec+2) = 0.d0
      do n=1,Nspec
        do j=lo(2),hi(2)
          gradY = 2*(RhoY(i,j,n)*rhoInv(i,j) - RhoY(i-1,j,n)*rhoInv(i-1,j))*dxInv
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) &
                   + (Fx(i,j,n) + rhoDx(i,j,Nspec+1)*gradY*Ax(i,j))*H(i-1,j,n)
        enddo
      enddo
    endif
!     xhi
    if (hi(1).eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
      i = dhi(1)+1
      Fx(i:i,lo(2):hi(2),Nspec+2) = 0.d0
      do n=1,Nspec
        do j=lo(2),hi(2)
          gradY = 2*(RhoY(i,j,n)*rhoInv(i,j) - RhoY(i-1,j,n)*rhoInv(i-1,j))*dxInv
          Fx(i,j,Nspec+2) = Fx(i,j,Nspec+2) &
                   + (Fx(i,j,n) + rhoDx(i,j,Nspec+1)*gradY*Ax(i,j))*H(i,j,n)
        enddo
      enddo
    endif
!     ylo
    if (lo(2).eq.dlo(2)  .and.  Tbc(2,1).eq.EXT_DIR) then
      j = dlo(2)
      Fy(lo(1):hi(1),j:j,Nspec+2) = 0.d0
      do n=1,Nspec
        do i=lo(1),hi(1)
          gradY = 2*(RhoY(i,j,n)*rhoInv(i,j) - RhoY(i,j-1,n)*rhoInv(i,j-1))*dyInv
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) &
                   + (Fy(i,j,n) + rhoDy(i,j,Nspec+1)*gradY*Ay(i,j))*H(i,j-1,n)
        enddo
      enddo
    endif
!     yhi
    if (hi(2)+1.eq.dhi(2)+1  .and.  Tbc(2,2).eq.EXT_DIR) then
      j = dhi(2)+1
      Fy(lo(1):hi(1),j:j,Nspec+2) = 0.d0
      do n=1,Nspec
        do i=lo(1),hi(1)
          gradY = 2*(RhoY(i,j,n)*rhoInv(i,j) - RhoY(i,j-1,n)*rhoInv(i,j-1))*dyInv
          Fy(i,j,Nspec+2) = Fy(i,j,Nspec+2) &
                   + (Fy(i,j,n) + rhoDy(i,j,Nspec+1)*gradY*Ay(i,j))*H(i,j,n)
        enddo
      enddo
    endif
    deallocate(rhoInv)
#endif

!     Set FiGHi = (species flux) dot Grad(species enthalpy)
!        compute Grad(H) on each face, and average across faces in each coordinate
!        Fi is extensive here, so need to remove area.  Also, assume that we are 
!        away from domain boundary, fix afterward
!
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        FiGHi(i,j) = 0.d0
      enddo
    enddo
    do n=1,Nspec
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          if (Ax(i,j).eq.0.d0) then
            AxDxInv_lo = 0.d0
          else
            AxDxInv_lo = 1.d0/(Ax(i,j)*dx(1))
          endif
          if (Ax(i+1,j).eq.0.d0) then
            AxDxInv_hi = 0.d0
          else
            AxDxInv_hi = 1.d0/(Ax(i+1,j)*dx(1))
          endif
          if (Ay(i,j).eq.0.d0) then
            AyDyInv_lo = 0.d0
          else
            AyDyInv_lo = 1.d0/(Ay(i,j)*dx(2))
          endif
          if (Ay(i,j+1).eq.0.d0) then
            AyDyInv_hi = 0.d0
          else
            AyDyInv_hi = 1.d0/(Ay(i,j+1)*dx(2))
          endif
               FiGHi(i,j) = FiGHi(i,j) - 0.5d0*  &
                  ( ( H(i+1,j,n) - H(i,j,n) )*Fx(i+1,j,n)*AxDxInv_hi &
                  + ( H(i,j,n) - H(i-1,j,n) )*Fx(i  ,j,n)*AxDxInv_lo &
                  + ( H(i,j+1,n) - H(i,j,n) )*Fy(i,j+1,n)*AyDyInv_hi &
                  + ( H(i,j,n) - H(i,j-1,n) )*Fy(i  ,j,n)*AyDyInv_lo )
        enddo
      enddo
    enddo

!     xlo
    if (fix_xlo) then
      i = lo(1)
      do j=lo(2),hi(2)
        FiGHi(i,j) = 0.d0
      enddo
      do n=1,Nspec
        do j=lo(2),hi(2)
          if (Ax(i,j).eq.0.d0) then
            AxDxInv_lo = 0.d0
          else
            AxDxInv_lo = 2.d0/(Ax(i,j)*dx(1))
          endif
          if (Ax(i+1,j).eq.0.d0) then
            AxDxInv_hi = 0.d0
          else
            AxDxInv_hi = 1.d0/(Ax(i+1,j)*dx(1))
          endif
          if (Ay(i,j).eq.0.d0) then
            AyDyInv_lo = 0.d0
          else
            AyDyInv_lo = 1.d0/(Ay(i,j)*dx(2))
          endif
          if (Ay(i,j+1).eq.0.d0) then
            AyDyInv_hi = 0.d0
          else
            AyDyInv_hi = 1.d0/(Ay(i,j+1)*dx(2))
          endif
               FiGHi(i,j) = FiGHi(i,j) - 0.5d0*  &
                   ( ( H(i+1,j,n) - H(i,j,n) )*Fx(i+1,j,n)*AxDxInv_hi  &
                   + ( H(i,j,n) - H(i-1,j,n) )*Fx(i  ,j,n)*AxDxInv_lo  &
                   + ( H(i,j+1,n) - H(i,j,n) )*Fy(i,j+1,n)*AyDyInv_hi  &
                   + ( H(i,j,n) - H(i,j-1,n) )*Fy(i  ,j,n)*AyDyInv_lo )
        enddo
      enddo
    endif

!     xhi
    if (fix_xhi) then
      i = hi(1)
      do j=lo(2),hi(2)
        FiGHi(hi(1),j) = 0.d0
      enddo
      do n=1,Nspec
        do j=lo(2),hi(2)
          if (Ax(i,j).eq.0.d0) then
            AxDxInv_lo = 0.d0
          else
            AxDxInv_lo = 1.d0/(Ax(i,j)*dx(1))
          endif
          if (Ax(i+1,j).eq.0.d0) then
            AxDxInv_hi = 0.d0
          else
            AxDxInv_hi = 2.d0/(Ax(i+1,j)*dx(1))
          endif
          if (Ay(i,j).eq.0.d0) then
            AyDyInv_lo = 0.d0
          else
            AyDyInv_lo = 1.d0/(Ay(i,j)*dx(2))
          endif
          if (Ay(i,j+1).eq.0.d0) then
            AyDyInv_hi = 0.d0
          else
            AyDyInv_hi = 1.d0/(Ay(i,j+1)*dx(2))
          endif
               FiGHi(i,j) = FiGHi(i,j) - 0.5d0*  &
                   ( ( H(i+1,j,n) - H(i,j,n) )*Fx(i+1,j,n)*AxDxInv_hi &
                   + ( H(i,j,n) - H(i-1,j,n) )*Fx(i  ,j,n)*AxDxInv_lo &
                   + ( H(i,j+1,n) - H(i,j,n) )*Fy(i,j+1,n)*AyDyInv_hi &
                   + ( H(i,j,n) - H(i,j-1,n) )*Fy(i  ,j,n)*AyDyInv_lo )
        enddo
      enddo
    endif

!     ylo
    if (fix_ylo) then
      j = lo(2)
        do i=lo(1),hi(1)
          FiGHi(i,j) = 0.d0
        enddo
        do n=1,Nspec
          do i=lo(1),hi(1)
            if (Ax(i,j).eq.0.d0) then
              AxDxInv_lo = 0.d0
            else
              AxDxInv_lo = 1.d0/(Ax(i,j)*dx(1))
            endif
            if (Ax(i+1,j).eq.0.d0) then
              AxDxInv_hi = 0.d0
            else
              AxDxInv_hi = 1.d0/(Ax(i+1,j)*dx(1))
            endif
            if (Ay(i,j).eq.0.d0) then
              AyDyInv_lo = 0.d0
            else
              AyDyInv_lo = 2.d0/(Ay(i,j)*dx(2))
            endif
            if (Ay(i,j+1).eq.0.d0) then
              AyDyInv_hi = 0.d0
            else
              AyDyInv_hi = 1.d0/(Ay(i,j+1)*dx(2))
            endif
            FiGHi(i,j) = FiGHi(i,j) - 0.5d0* &
                   ( ( H(i+1,j,n) - H(i,j,n) )*Fx(i+1,j,n)*AxDxInv_hi &
                   + ( H(i,j,n) - H(i-1,j,n) )*Fx(i  ,j,n)*AxDxInv_lo &
                   + ( H(i,j+1,n) - H(i,j,n) )*Fy(i,j+1,n)*AyDyInv_hi &
                   + ( H(i,j,n) - H(i,j-1,n) )*Fy(i  ,j,n)*AyDyInv_lo )
      enddo
    enddo
  endif

!     yhi
    if (fix_yhi) then
      j = hi(2)
      do i=lo(1),hi(1)
        FiGHi(i,j) = 0.d0
      enddo
      do n=1,Nspec
        do i=lo(1),hi(1)
          if (Ax(i,j).eq.0.d0) then
            AxDxInv_lo = 0.d0
          else
            AxDxInv_lo = 1.d0/(Ax(i,j)*dx(1))
          endif
          if (Ax(i+1,j).eq.0.d0) then
            AxDxInv_hi = 0.d0
          else
            AxDxInv_hi = 1.d0/(Ax(i+1,j)*dx(1))
          endif
          if (Ay(i,j).eq.0.d0) then
            AyDyInv_lo = 0.d0
          else
            AyDyInv_lo = 1.d0/(Ay(i,j)*dx(2))
          endif
          if (Ay(i,j+1).eq.0.d0) then
            AyDyInv_hi = 0.d0
          else
            AyDyInv_hi = 2.d0/(Ay(i,j+1)*dx(2))
          endif
          FiGHi(i,j) = FiGHi(i,j) - 0.5d0* &
                   ( ( H(i+1,j,n) - H(i,j,n) )*Fx(i+1,j,n)*AxDxInv_hi &
                   + ( H(i,j,n) - H(i-1,j,n) )*Fx(i  ,j,n)*AxDxInv_lo &
                   + ( H(i,j+1,n) - H(i,j,n) )*Fy(i,j+1,n)*AyDyInv_hi &
                   + ( H(i,j,n) - H(i,j-1,n) )*Fy(i  ,j,n)*AyDyInv_lo )
        enddo
      enddo
    endif

     deallocate(H)

  end subroutine enth_diff_terms

!-------------------------------------

  integer function pphys_CONPsolv_SDC(lo, hi, &
      rhoYnew,   DIMS(rhoYnew),  &
      rhoHnew,   DIMS(rhoHnew), &
      Tnew,      DIMS(Tnew), &
      rhoYold,   DIMS(rhoYold),  &
      rhoHold,   DIMS(rhoHold),  &
      Told,      DIMS(Told), &
      const_src, DIMS(const_src), &
      FuncCount, DIMS(FuncCount), &
      dt, &
      diag, do_diag, do_stiff)&
      bind(C, name="pphys_CONPsolv_SDC")

      use network,        only : nspec
      use reactor_module, only : react
          
      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(rhoYold)
      integer DIMDEC(rhoHold)
      integer DIMDEC(Told)
      integer DIMDEC(rhoYnew)
      integer DIMDEC(rhoHnew)
      integer DIMDEC(Tnew)
      integer DIMDEC(const_src)
      integer DIMDEC(FuncCount)
      integer do_diag, do_stiff
      REAL_T rhoYold(DIMV(rhoYold),*)
      REAL_T rhoHold(DIMV(rhoHold))
      REAL_T Told(DIMV(Told))
      REAL_T rhoYnew(DIMV(rhoYnew),*)
      REAL_T rhoHnew(DIMV(rhoHnew))
      REAL_T Tnew(DIMV(rhoHnew))
      REAL_T const_src(DIMV(const_src),1:Nspec+1)
      REAL_T FuncCount(DIMV(FuncCount))
      REAL_T dt
      REAL_T diag(DIMV(FuncCount),*)


      integer i, j, lierr
      REAL_T  TT1, TT2
      double precision  Z(nspec+1), Ysrc(nspec)
      REAL_T  Yt(nspec)
      double precision energy, energy_src
      REAL_T  rho, rhoInv, HMIX_CGS
      double precision  pressure

      !print *, "IN CONPsolv_SDC "
      call flush
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            !print *, "IJK ?", i,j
            TT1                = zero
            TT2                = dt
            ! CONVERTING INTO CGS
            Z(1:Nspec)         = rhoYold(i,j,1:Nspec) * 1.D-3
            Z(Nspec+1)         = Told(i,j)
            rho                = sum(Z(1:Nspec))
            energy             = rhoHold(i,j) * 10.0
            Ysrc(1:Nspec)      = const_src(i,j,1:Nspec) * 1.d-3
            energy_src         = const_src(i,j,Nspec+1) * 10.0


            ! dummy pressure not used in PeleLM
            pressure = 1013250.0
            FuncCount(i,j) = react(Z, Ysrc,&
                           energy, energy_src,&
                           pressure,TT2,TT1,0)
               
            ! CONVERTING BACK TO MKS
            rhoYnew(i,j,1:Nspec) = Z(1:Nspec) * 1.d+03
            rho                  = sum(rhoYnew(i,j,1:Nspec))
            rhoHnew(i,j)         = energy * 1.D-1 
            Tnew(i,j)            = Z(Nspec+1)

            ! Necessary ?
            HMIX_CGS = rhoHnew(i,j) * 1.0d4 / rho
            Yt(1:nspec) = rhoYnew(i,j,1:nspec) / rho
            call get_t_given_hy(HMIX_CGS, Yt, Tnew(i,j), lierr);
            !print *, "rho, T, rhoenergy, rhoenergy_src ", rho, Tnew(i,j), rhoHnew(i,j), const_src(i,j,Nspec+1)

         end do
      end do

      pphys_CONPsolv_SDC = 1

  end function pphys_CONPsolv_SDC

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
      integer i,j,n

      TIME = 0.

      !print *, "IN RRATERHOY ..."
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            !print *, "IJK ? ", i,j

            Zt(Nspec+1) = RhoH(i,j)
            do n=1,Nspec
               Zt(n) = RhoY(i,j,n)
            end do
            Temperature = T(i,j)

            call pphys_calc_src_sdc(nspec,TIME,Temperature,Zt,Zdott)

            do n=1,Nspec
               RhoYdot(i,j,n) = Zdott(n)
               !print *," RhoYdot ", RhoYdot(i,j,n)
            end do

         end do
      end do
      !call flush
      
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
      
      integer i, j, n
      REAL_T Yt(nspec), RHOt, SCAL, SCAL1
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
!           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
      SCAL = 1.d-1
      SCAL1 = SCAL**3

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            do n=1,nspec
               Yt(n) = Y(i,j,n)
            end do

            RHOt = RHO(i,j) * SCAL1
            CALL CKPY(RHOt,T(i,j),Yt,P(i,j))

            P(i,j) = P(i,j) * SCAL

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
      integer i,j,n

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            do n = 1,Nspec
               Yt(n) = Y(i,j,n)
            end do
            CALL CKYTX(Yt,Xt)
            do n = 1,Nspec
               X(i,j,n) = Xt(n)
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
      integer i,j,n

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            do n = 1,Nspec
               Yt(n) = Y(i,j,n)
            end do
            rhoScl = RHO(i,j)*1.e-3
            CALL CKYTCR(rhoScl,T(i,j),Yt,Ct)
            do n = 1,Nspec
               C(i,j,n) = Ct(n)*1.e6
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
      
      integer i, j, n
      REAL_T Yt(nspec), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            do n=1,Nspec
               Yt(n) = Y(i,j,n)
            end do

            CALL CKHBMS(T(i,j),Yt,HMIX(i,j))

            HMIX(i,j) = HMIX(i,j) * SCAL

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
      
      integer i, j, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(nspec), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      CALL CKRP(RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            do n=1,Nspec
               Yt(n) = Y(i,j,n)
            end do
            CALL CKRHOY(Ptmp,T(i,j),Yt,RHO(i,j))
            RHO(i,j) = RHO(i,j) * SCAL
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
      
      integer i, j, n
      REAL_T SCAL, Ht(nspec)
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            CALL CKHMS(T(i,j),Ht)

            do n=1,Nspec
               H(i,j,n) = Ht(n) * SCAL
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
      
      integer i, j, n
      REAL_T Yt(nspec)

!     Returns mean molecular weight in kg/kmole
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            do n=1,nspec
               Yt(n) = Y(i,j,n)
            end do

            CALL CKMMWY(Yt,MWMIX(i,j))

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
      !CALL FLUSH

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

      REAL_T Yt(nspec)
      integer i, j, n, lierr, Niter,MAXiters
      REAL_T HMIX_CGS, Tguess

      MAXiters = 0
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            do n=1,nspec
               Yt(n) = Y(i,j,n)
            end do

            !HMIX_CGS = HMIX(i,j) * 1.0d4
            !call get_t_given_hy(HMIX_CGS, Yt, T(i,j), lierr);

            call pphys_TfromHYpt(T(i,j),HMIX(i,j),Yt,errMax,NiterMAX,res,Niter)

            if (Niter .lt. 0) then
                    call bl_abort(" Something went wrong in pphys_TfromHYpt ")
            end if

            if (Niter .gt. MAXiters) then
                    MAXiters = Niter
            end if
            
         end do
      end do

!     Set max iters taken during this solve, and exit
      pphys_TfromHY = MAXiters !10
      return

  end function pphys_TfromHY

!-------------------------------------

  subroutine compute_rho_dgrad_hdot_grad_Y(dx, &
               lo, hi, DIMS(species), species, &
               DIMS(h), h, DIMS(betax), betax, &
               DIMS(betay), betay, DIMS(rdghdgy), rdghdgy) &
               bind(C, name="compute_rho_dgrad_hdot_grad_y")

    implicit none

! ... inputs

    integer lo(SDIM), hi(SDIM)
    REAL_T  dx(SDIM)
    integer DIMDEC(species)
    integer DIMDEC(h)
    REAL_T  species(DIMV(species))
    REAL_T  h(DIMV(h))
    integer DIMDEC(betax)
    integer DIMDEC(betay)
    REAL_T betax(DIMV(betax))
    REAL_T betay(DIMV(betay))
    integer DIMDEC(rdghdgy)

! ... outputs

    REAL_T rdghdgy(DIMV(rdghdgy))

! ... local

    integer i,j
    REAL_T  dxsqr, dysqr
    REAL_T  betadotleft, betadotright
    REAL_T  betadottop,  betadotbot

    dxsqr = dx(1)**2
    dysqr = dx(2)**2
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        betadotleft  = betax(i,j)*(h(i,j)-h(i-1,j))*(species(i,j)-species(i-1,j))
        betadotright = betax(i+1,j)*(h(i+1,j)-h(i,j))*(species(i+1,j)-species(i,j))
        betadotbot   = betay(i,j)*(h(i,j)-h(i,j-1))*(species(i,j)-species(i,j-1))
        betadottop   = betay(i,j+1)*(h(i,j+1)-h(i,j))*(species(i,j+1)-species(i,j))
        rdghdgy(i,j) =  half*((betadotleft + betadotright)/dxsqr + &
                               (betadottop  + betadotbot)/dysqr)
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
    use mod_Fvar_def, only : use_constant_mu, constant_mu_val, LeEQ1

    implicit none

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
        ! use_eg 
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
    use mod_Fvar_def, only : use_constant_rhoD, constant_rhoD_val, use_constant_lambda,  &
                             constant_lambda_val
    use mod_Fvar_def, only : Pr, Sc, LeEQ1, thickFac
    
    implicit none

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
    REAL_T  Patm, Wavg
    REAL_T  Yt(nspec), invmwt(nspec)
    REAL_T  Tfac, Yfac, cpmix(1,1,1)
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
    Tfac = thickFac / Pr
    Yfac = thickFac / Sc

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
              end do
              if (do_temp .ne. 0) then 
                rhoD(i,j,k,nspec+1) = LAM(i,j,k) * (one / 100000.0D0)
              end if
              if (thickFac.ne.1.d0) then
                do n=1,nspec+1
                  rhoD(i,j,k,n) = rhoD(i,j,k,n)*thickFac
                enddo
              endif
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

        do k=lo(3),hi(3) 
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
          do k=lo(3),hi(3) 
            do j=lo(2), hi(2)
              do i=lo(1), hi(1)
                rhoD(i,j,k,n) = rhoD(i,j,k,1)
              end do
            end do
          end do
        end do

        if (do_temp .ne. 0) then
          do k=lo(3),hi(3) 
            do j=lo(2), hi(2)
              do i=lo(1), hi(1)
                do n=1,nspec
                  Yt(n) = Y_real(i,j,k,n)
                end do
                CALL pphys_CPMIXfromTY(lo,       hi, & 
                                       cpmix,    lo_chem, hi_chem, &
                                       T(i,j,k), lo_chem, hi_chem, &
                                       Yt,       lo_chem, hi_chem)
               rhoD(i,j,k,nspec+1) = rhoD(i,j,k,nspec+1)*cpmix(1,1,1)*Tfac
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
         if (do_temp .ne. 0) then
           if (.NOT. use_constant_lambda) THEN  
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
           else
             do k=lo(3), hi(3)
                do j=lo(2), hi(2)
                   do i=lo(1), hi(1)
                      rhoD(i,j,k,Nspec+1) = constant_lambda_val
                   end do
                end do
             end do
           end if
         end if
      end if

  end subroutine spec_temp_visc

!-------------------------------------------

  subroutine est_divu_dt(flag, dtfactor, delta, divu, DIMS(divu), &
                         dsdt, DIMS(dsdt), rho, DIMS(rho), & 
                         u, DIMS(u), & 
                         volume, DIMS(volume), &
                         areax,  DIMS(areax), &
                         areay,  DIMS(areay), &
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
    REAL_T  volume(DIMV(volume))
    REAL_T  areax(DIMV(areax))
    REAL_T  areay(DIMV(areay))

    REAL_T dt

    integer i,j
    REAL_T  dtcell, dtcell2, denom, rhominij, rhoij
    REAL_T  fluxtop, fluxbot, fluxleft, fluxright
    REAL_T  a,b,c

    dt = 1.0d20
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        dtcell = dt
        if (flag.eq.1) then
          if(divu(i,j).gt.zero) then
            if(rho(i,j).gt.rhomin) then
              dtcell = dtfactor*(one-rhomin/rho(i,j))/divu(i,j)
            else
              dtcell = dtfactor*.5/divu(i,j)
            endif
            if (dsdt(i,j).gt.1.0e-20) then
              if (abs(rho(i,j)).gt.rhomin) then
                rhominij = rhomin
              else
                rhominij = .9*abs(rho(i,j)) 
              endif
              rhoij = abs(rho(i,j))
#if 0
              dtcell2 = (-rhoij*divu(i,j) +  &
                          sqrt((rhoij*divu(i,j))**2+ &
                               two*dsdt(i,j)*rhoij*abs(rhoij-rhominij)))/ &
                         (rhoij*dsdt(i,j))
#else
! ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!           We use the latter because it is more robust
              a = rhoij*dsdt(i,j)*half
              b = rhoij*divu(i,j)
              c = rhominij - rhoij
              dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))
#endif
              dtcell2 = dtfactor*dtcell2
              dtcell = min(dtcell,dtcell2)
            endif
          endif
          if(dtcell.le.0.0d0)then
            write(6,*)'aha'
          endif
        else if (flag.eq.2) then
          denom = rho(i,j)*divu(i,j)+ &
                   u(i,j,1)*(rho(i+1,j)-rho(i-1,j))/delta(1) + &
                   u(i,j,2)*(rho(i,j+1)-rho(i,j-1))/delta(2)
          if(denom.gt.zero)then
            if(rho(i,j).gt.rhomin) then
              dtcell = dtfactor*(rho(i,j)-rhomin)/denom
            else
              dtcell = dtfactor*abs(rho(i,j))/denom
            endif
          endif
        else if (flag.eq.3) then
          fluxtop   = fourth*(rho(i,j)+rho(i,j+1))*(u(i,j,2)+u(i,j+1,2))
          fluxbot   = fourth*(rho(i,j)+rho(i,j-1))*(u(i,j,2)+u(i,j-1,2))
          fluxleft  = fourth*(rho(i,j)+rho(i-1,j))*(u(i,j,1)+u(i-1,j,1))
          fluxright = fourth*(rho(i,j)+rho(i+1,j))*(u(i,j,1)+u(i+1,j,1))
          denom = ((areax(i+1,j)*fluxright-areax(i,j)*fluxleft)+ &
                    (areay(i,j+1)*fluxtop-areay(i,j)*fluxbot))/volume(i,j)

          if(denom.gt.zero)then
            if(rho(i,j).gt.rhomin) then
              dtcell = dtfactor*(rho(i,j)-rhomin)/denom
            else
              dtcell = dtfactor*abs(rho(i,j))/denom
            endif
          endif
        endif
#if 0
        write(6,*)'i,j,dtcell=',i,j,dtcell
#endif
        dt = min(dtcell,dt)
      enddo
    enddo

    !write(*,*) 'DEBUG in fortran 1 ',dt
    
  end subroutine est_divu_dt

!---------------------------------------------------------------

  subroutine check_divu_dt(flag, dtfactor, delta, divu, DIMS(divu), &
                           dsdt, rho, DIMS(rho), &
                           u, DIMS(u), &
                           volume, DIMS(volume), &
                           areax,  DIMS(areax), &
                           areay,  DIMS(areay), &
                           lo, hi, dt, rhomin) &
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
    REAL_T  volume(DIMV(volume))
    REAL_T  areax(DIMV(areax))
    REAL_T  areay(DIMV(areay))

    integer i,j
    REAL_T  dtcell, denom
    REAL_T  fluxtop, fluxbot, fluxleft, fluxright
    REAL_T a,b,c,dtcell2,rhominij,rhoij

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        dtcell = bigreal
        if (flag.eq.1) then
          if(divu(i,j).gt.zero) then
            if(rho(i,j).gt.rhomin) then
              dtcell = (one-rhomin/rho(i,j))/divu(i,j)
            else
              dtcell = one/divu(i,j)
            endif
            if (dsdt(i,j).gt.1.0e-20) then
              if (abs(rho(i,j)).gt.rhomin) then
                rhominij = rhomin
              else
                rhominij = .9*abs(rho(i,j)) 
              endif
              rhoij = abs(rho(i,j))
#if 0
              dtcell2 = (-rhoij*divu(i,j) + & 
                          sqrt((rhoij*divu(i,j))**2+ &
                               two*dsdt(i,j)*rhoij*abs(rhoij-rhominij)))/ &
                         (rhoij*dsdt(i,j))
#else
! ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!           We use the latter because it is more robust
              a = rhoij*dsdt(i,j)*half
              b = rhoij*divu(i,j)
              c = rhominij - rhoij
              dtcell2 = two*c/(-b-sqrt(b**2-four*a*c))
#endif
              dtcell = min(dtcell,dtcell2)
            endif
          endif
        else if (flag.eq.2) then
          denom = rho(i,j)*divu(i,j)+  &
                   u(i,j,1)*(rho(i+1,j)-rho(i-1,j))/delta(1) + &
                   u(i,j,2)*(rho(i,j+1)-rho(i,j-1))/delta(2)
          if(denom.gt.zero)then
            if(rho(i,j).gt.rhomin) then
              dtcell = (rho(i,j)-rhomin)/denom
            else
              dtcell = abs(rho(i,j))/denom
            endif
          endif
        else if (flag.eq.3) then
          fluxtop = fourth*(rho(i,j)+rho(i,j+1))*(u(i,j,2)+u(i,j+1,2))
          fluxbot = fourth*(rho(i,j)+rho(i,j-1))*(u(i,j,2)+u(i,j-1,2))
          fluxleft = fourth*(rho(i,j)+rho(i-1,j))*(u(i,j,1)+u(i-1,j,1))
          fluxright = fourth*(rho(i,j)+rho(i+1,j))*(u(i,j,1)+u(i+1,j,1))
          denom = ((areax(i+1,j)*fluxright-areax(i,j)*fluxleft)+ &
                    (areay(i,j+1)*fluxtop-areay(i,j)*fluxbot))/volume(i,j)
          if(denom.gt.zero)then
            if(rho(i,j).gt.rhomin) then
              dtcell = (rho(i,j)-rhomin)/denom
            else
              dtcell = abs(rho(i,j))/denom
            endif
          endif
        endif
        if (dt.gt.dtcell) then
          write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,dt>dtcell = ', &
                 i,j,dt,dtcell
        else if (dt.gt.dtcell*dtfactor) then
            write(6,*)'WARNING: ', &
                 'FORT_CHECK_DIVU_DT : i,j,dt>dtcell*dtfactor = ', &
                 i,j,dt,dtcell*dtfactor
        endif
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

  subroutine dqrad_fill (dqrad,DIMS(dqrad),domlo,domhi,delta, &
                              xlo,time,bc )bind(C, name="dqrad_fill")

    integer    DIMDEC(dqrad), bc(SDIM,2)
    integer    domlo(SDIM), domhi(SDIM)
    integer    lo(SDIM), hi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     dqrad(DIMV(dqrad))

    integer    i, j
    integer    ilo, ihi, jlo, jhi

    lo(1) = ARG_L1(dqrad)
    hi(1) = ARG_H1(dqrad)
    lo(2) = ARG_L2(dqrad)
    hi(2) = ARG_H2(dqrad)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
      if(jlo.le.jhi)then
        do j = jlo, jhi
          do i = lo(1), domlo(1)-1
            dqrad(i,j) = dqrad(domlo(1),j)
          enddo
        enddo
      endif
      if (lo(2).lt.domlo(2)) then
        do j = lo(2), domlo(2)-1
          do i = lo(1), domlo(1)-1
            dqrad(i,j) = dqrad(domlo(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(2).gt.domhi(2))then
        do j = domhi(2)+1, hi(2)
          do i = lo(1), domlo(1)-1
            dqrad(i,j) = dqrad(domlo(1),domhi(2))
          enddo
        enddo
      endif

    endif            

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
      if(jlo.le.jhi)then
        do j = jlo,jhi
          do i = domhi(1)+1,hi(1)
            dqrad(i,j) = dqrad(domhi(1),j)
          enddo
        enddo
      endif
      if (lo(2).lt.domlo(2)) then
        do j = lo(2), domlo(2)-1
          do i = domhi(1)+1,hi(1)
            dqrad(i,j) = dqrad(domhi(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(2).gt.domhi(2))then
        do j = domhi(2)+1, hi(2)
          do i = domhi(1)+1,hi(1)
            dqrad(i,j) = dqrad(domhi(1),domhi(2))
          enddo
        enddo
      endif
    endif            

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
      if(ilo.le.ihi)then
        do j = lo(2), domlo(2)-1
          do i = ilo,ihi
            dqrad(i,j) = dqrad(i,domlo(2))
          enddo
        enddo
      endif
      if (lo(1).lt.domlo(1)) then
        do j = lo(2), domlo(2)-1
          do i = lo(1), domlo(1)-1
            dqrad(i,j) = dqrad(domlo(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(1).gt.domhi(1))then
        do j = lo(2), domlo(2)-1
          do i = domhi(1)+1, hi(1)
            dqrad(i,j) = dqrad(domhi(1),domlo(2))
          enddo
        enddo
      endif

    endif            

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
      if(ilo.le.ihi)then
        do j = domhi(2)+1, hi(2)
          do i = ilo,ihi
            dqrad(i,j) = dqrad(i,domhi(2))
          enddo
        enddo
      endif
      if (lo(1).lt.domlo(1)) then
        do j = domhi(2)+1, hi(2)
          do i = lo(1), domlo(1)-1
            dqrad(i,j) = dqrad(domlo(1),domhi(2))
          enddo
        enddo
      endif
      if(hi(1).gt.domhi(1))then
        do j = domhi(2)+1, hi(2)
          do i = domhi(1)+1, hi(1)
            dqrad(i,j) = dqrad(domhi(1),domhi(2))
          enddo
        enddo
      endif
    endif            

  end subroutine dqrad_fill

!--------------------------------------------------

  subroutine divu_fill (divu,DIMS(divu),domlo,domhi,delta, &
                        xlo,time,bc )bind(C, name="divu_fill")

    integer    DIMDEC(divu)
    integer    bc(SDIM,2)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     divu(DIMV(divu))

    integer    i, j
    integer    ilo, ihi, jlo, jhi

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(divu)
    hi(1) = ARG_H1(divu)
    lo(2) = ARG_L2(divu)
    hi(2) = ARG_H2(divu)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
      if(jlo.le.jhi)then
        do j = jlo, jhi
          do i = lo(1), domlo(1)-1
            divu(i,j) = divu(domlo(1),j)
          enddo
        enddo
      endif
      if (lo(2).lt.domlo(2)) then
        do j = lo(2), domlo(2)-1
          do i = lo(1), domlo(1)-1
            divu(i,j) = divu(domlo(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(2).gt.domhi(2))then
        do j = domhi(2)+1, hi(2)
          do i = lo(1), domlo(1)-1
            divu(i,j) = divu(domlo(1),domhi(2))
          enddo
        enddo
      endif

    endif            

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
      if(jlo.le.jhi)then
        do j = jlo,jhi
          do i = domhi(1)+1,hi(1)
            divu(i,j) = divu(domhi(1),j)
          enddo
        enddo
      endif
      if (lo(2).lt.domlo(2)) then
        do j = lo(2), domlo(2)-1
          do i = domhi(1)+1,hi(1)
            divu(i,j) = divu(domhi(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(2).gt.domhi(2))then
        do j = domhi(2)+1, hi(2)
          do i = domhi(1)+1,hi(1)
            divu(i,j) = divu(domhi(1),domhi(2))
          enddo
        enddo
      endif
    endif            

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
      if(ilo.le.ihi)then
        do j = lo(2), domlo(2)-1
          do i = ilo,ihi
            divu(i,j) = divu(i,domlo(2))
          enddo
        enddo
      endif
      if (lo(1).lt.domlo(1)) then
        do j = lo(2), domlo(2)-1
          do i = lo(1), domlo(1)-1
            divu(i,j) = divu(domlo(1),domlo(2))
          enddo
        enddo
      endif
      if(hi(1).gt.domhi(1))then
        do j = lo(2), domlo(2)-1
          do i = domhi(1)+1, hi(1)
            divu(i,j) = divu(domhi(1),domlo(2))
          enddo
        enddo
      endif

    endif            

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
      if(ilo.le.ihi)then
        do j = domhi(2)+1, hi(2)
          do i = ilo,ihi
            divu(i,j) = divu(i,domhi(2))
          enddo
        enddo
      endif
      if (lo(1).lt.domlo(1)) then
        do j = domhi(2)+1, hi(2)
          do i = lo(1), domlo(1)-1
            divu(i,j) = divu(domlo(1),domhi(2))
          enddo
        enddo
      endif
      if(hi(1).gt.domhi(1))then
        do j = domhi(2)+1, hi(2)
          do i = domhi(1)+1, hi(1)
            divu(i,j) = divu(domhi(1),domhi(2))
          enddo
        enddo
      endif
    endif            

  end subroutine divu_fill

!--------------------------------------------

  subroutine dsdt_fill (dsdt,DIMS(dsdt),domlo,domhi,delta, &
                              xlo,time,bc )bind(C, name="dsdt_fill")

    integer    DIMDEC(dsdt)
    integer    bc(SDIM,2)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     dsdt(DIMV(dsdt))

    integer    i, j
    integer    ilo, ihi, jlo, jhi

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(dsdt)
    hi(1) = ARG_H1(dsdt)
    lo(2) = ARG_L2(dsdt)
    hi(2) = ARG_H2(dsdt)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
      do i = lo(1), domlo(1)-1
        do j = lo(2), hi(2)
          dsdt(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
      do i = domhi(1)+1, hi(1)
        do j = lo(2), hi(2)
          dsdt(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
!                               inflow for burner in a can (bic, biac)

      do j = lo(2), domlo(2)-1
        do i = lo(1), hi(1)
          dsdt(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
      do j = domhi(2)+1, hi(2)
        do i = lo(1), hi(1)
          dsdt(i,j) = zero
        enddo
      enddo
    endif            

  end subroutine dsdt_fill

!--------------------------------------

  subroutine ydot_fill (ydot,DIMS(ydot),domlo,domhi,delta, &
                               xlo,time,bc)bind(C, name="ydot_fill")

    integer    DIMDEC(ydot), bc(SDIM,2)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     ydot(DIMV(ydot))

    integer    i, j
    integer    jlo, jhi, ilo, ihi
    REAL_T     x

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(ydot)
    hi(1) = ARG_H1(ydot)
    lo(2) = ARG_L2(ydot)
    hi(2) = ARG_H2(ydot)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
!                               inflow for pipe
      do j = lo(2), domlo(2)-1
        do i = lo(1), hi(1)
          ydot(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
      do j = domhi(2)+1, hi(2)
        do i = lo(1), hi(1)
          ydot(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
      do i = lo(1), domlo(1)-1
        do j = lo(2), hi(2)
          ydot(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
      do i = domhi(1)+1, hi(1)
        do j = lo(2),hi(2)
          x = (float(i-lo(1))+.5)*delta(1)+xlo(1)
          ydot(i,j) = zero
        enddo
      enddo
    endif            

  end subroutine ydot_fill

!-------------------------------------------

  subroutine rhoYdot_fill (rhoydot,DIMS(rhoydot),domlo,domhi,delta, &
                                  xlo,time,bc)bind(C, name="rhoYdot_fill")

    integer    DIMDEC(rhoydot), bc(SDIM,2)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     rhoydot(DIMV(rhoydot))

    integer    i, j
    integer    jlo, jhi, ilo, ihi
    REAL_T     x

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(rhoydot)
    hi(1) = ARG_H1(rhoydot)
    lo(2) = ARG_L2(rhoydot)
    hi(2) = ARG_H2(rhoydot)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (rhoydot,DIMS(rhoydot),domlo,domhi,delta,xlo,bc)

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
!                               inflow for pipe
      do j = lo(2), domlo(2)-1
        do i = lo(1), hi(1)
          rhoydot(i,j) = 0.0
        enddo
      enddo
    endif            

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
      do j = domhi(2)+1, hi(2)
        do i = lo(1), hi(1)
          rhoydot(i,j) = 0.0
        enddo
      enddo
    endif            

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
      do i = lo(1), domlo(1)-1
        do j = lo(2), hi(2)
          rhoydot(i,j) = zero
        enddo
      enddo
    endif            

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
      do i = domhi(1)+1, hi(1)
        do j = lo(2),hi(2)
          x = (float(i-lo(1))+.5)*delta(1)+xlo(1)
          rhoydot(i,j) = zero
        enddo
      enddo
    endif            

  end subroutine rhoYdot_fill

!------------------------------------------------

  subroutine fab_minmax (lo, hi, & 
                        fab, DIMS(fab), &
                        fmin, fmax, nc) &
                        bind(C, name="fab_minmax")

    integer lo(SDIM), hi(SDIM), nc
    integer DIMDEC(fab)
    REAL_T  fab(DIMV(fab),nc)
    REAL_T  fmin, fmax

    integer i,j,n

    do n = 1,nc
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          fab(i,j,n) = MAX( fmin, MIN( fmax, fab(i,j,n) ) )
        end do
      end do
    end do

  end subroutine fab_minmax

!--------------------------------------------------

  subroutine repair_flux (lo, hi, dlo, dhi, &
                          flux, DIMS(flux), &
                          RhoY, DIMS(RhoY), dir, Ybc)bind(C, name="repair_flux")

    use network,        only : nspec

    implicit none

    integer lo(SDIM), hi(SDIM), dlo(SDIM), dhi(SDIM), dir, Ybc(SDIM,2)
    integer DIMDEC(flux)
    integer DIMDEC(RhoY)
    REAL_T flux(DIMV(flux),Nspec)
    REAL_T RhoY(DIMV(RhoY),Nspec)
      
    integer i, j, n
    REAL_T sumFlux, RhoYe(Nspec), sumRhoYe


    if (dir.eq.0) then
!write(*,*) "DEBUG 1",lo,hi,dlo,dhi
!write(*,*) "DEBUG 2",DIMS(flux)
!     First, assume away from physical boundaries, then use boundary-aware version below if applicable
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          sumFlux = 0.d0
          sumRhoYe = 0.d0
          do n=1,Nspec
            sumFlux = sumFlux + flux(i,j,n)
            RhoYe(n) = 0.5d0*(RhoY(i-1,j,n) + RhoY(i,j,n))
            sumRhoYe = sumRhoYe + RhoYe(n)
          end do               
          do n=1,Nspec
            flux(i,j,n) = flux(i,j,n) - sumFlux*RhoYe(n)/sumRhoYe
          end do
        end do
      end do
!     xlo
      if (Ybc(1,1).eq.EXT_DIR.and.lo(1).le.dlo(1)) then
        do i = lo(1),dlo(1)
          do j = lo(2),hi(2)
            sumFlux = 0.d0
            sumRhoYe = 0.d0
            do n=1,Nspec
              sumFlux = sumFlux + flux(i,j,n)
              sumRhoYe = sumRhoYe + RhoY(i-1,j,n)
            enddo
            do n=1,Nspec
              flux(i,j,n) = flux(i,j,n) - sumFlux*RhoY(i-1,j,n)/sumRhoYe
            enddo
          enddo
        enddo
      endif
!     xhi
      if (Ybc(1,2).eq.EXT_DIR.and.hi(1).ge.dhi(1)) then
        do i = dhi(1),hi(1)
          do j = lo(2),hi(2)
            sumFlux = 0.d0
            sumRhoYe = 0.d0
            do n=1,Nspec
              sumFlux = sumFlux + flux(i,j,n)
              sumRhoYe = sumRhoYe + RhoY(i,j,n)
            enddo
            do n=1,Nspec
              flux(i,j,n) = flux(i,j,n) - sumFlux*RhoY(i,j,n)/sumRhoYe
            enddo
          enddo
        enddo
       endif

    else if (dir.eq.1) then

!     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          sumFlux = 0.d0
          sumRhoYe = 0.d0
          do n=1,Nspec
            sumFlux = sumFlux + flux(i,j,n)
            RhoYe(n) = 0.5d0*(RhoY(i,j-1,n) + RhoY(i,j,n))
            sumRhoYe = sumRhoYe + RhoYe(n)
          enddo
          do n=1,Nspec
            flux(i,j,n) = flux(i,j,n) - sumFlux*RhoYe(n)/sumRhoYe
          end do
        end do
      end do
!     ylo
      if (Ybc(2,1).eq.EXT_DIR.and.lo(2).le.dlo(2)) then
        do j = lo(2),dlo(2)
          do i = lo(1),hi(1)
            sumFlux = 0.d0
            sumRhoYe = 0.d0
            do n=1,Nspec
              sumFlux = sumFlux + flux(i,j,n)
              sumRhoYe = sumRhoYe + RhoY(i,j-1,n)
            enddo
            do n=1,Nspec
              flux(i,j,n) = flux(i,j,n) - sumFlux*RhoY(i,j-1,n)/sumRhoYe
            enddo
          enddo
        enddo
      endif
!     yhi
      if (Ybc(2,2).eq.EXT_DIR.and.hi(2).ge.dhi(2)) then
        do j = dhi(2),hi(2)
          do i = lo(1),hi(1)
            sumFlux = 0.d0
            sumRhoYe = 0.d0
            do n=1,Nspec
              sumFlux = sumFlux + flux(i,j,n)
              sumRhoYe = sumRhoYe + RhoY(i,j,n)
            enddo
            do n=1,Nspec
              flux(i,j,n) = flux(i,j,n) - sumFlux*RhoY(i,j,n)/sumRhoYe
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
                              stateo, DIMS(stateo), &
                              staten, DIMS(staten), &
                              vol,    DIMS(vol), &
                              nc, dt) &
                              bind(C, name="incrwext_flx_div")

    implicit none
    integer lo(SDIM), hi(SDIM), nc
    integer DIMDEC(xflux)
    integer DIMDEC(yflux)
    integer DIMDEC(stateo)
    integer DIMDEC(staten)
    integer DIMDEC(vol)
    REAL_T xflux(DIMV(xflux),nc)
    REAL_T yflux(DIMV(yflux),nc)
    REAL_T stateo(DIMV(stateo))
    REAL_T staten(DIMV(staten))
    REAL_T vol(DIMV(vol))
    REAL_T dt

    integer i, j, n
    REAL_T dF
      
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        dF = zero
        do n=1,nc
          dF = dF + ( xflux(i+1,j,n) - xflux(i,j,n) ) &
                   +    ( yflux(i,j+1,n) - yflux(i,j,n) ) 
        end do
        staten(i,j) = stateo(i,j) + dF*dt/vol(i,j)
      end do
    end do

  end subroutine incrwext_flx_div

!----------------------------------------------

  subroutine flux_div (lo, hi, &
                       update, DIMS(update), &
                       xflux,  DIMS(xflux), &
                       yflux,  DIMS(yflux), &
                       vol,    DIMS(vol), &
                       nc, scal)&
                       bind(C, name="flux_div")

    implicit none
    integer lo(SDIM), hi(SDIM), nc
    integer DIMDEC(update)
    integer DIMDEC(xflux)
    integer DIMDEC(yflux)
    integer DIMDEC(vol)
    REAL_T update(DIMV(update),nc)
    REAL_T xflux(DIMV(xflux),nc)
    REAL_T yflux(DIMV(yflux),nc)
    REAL_T vol(DIMV(vol))
    REAL_T scal

    integer i, j, n

    do n=1,nc
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          update(i,j,n)=scal*  &
                   ( (xflux(i+1,j,n)-xflux(i,j,n)) &
                   + (yflux(i,j+1,n)-yflux(i,j,n)) )/vol(i,j)
        end do
      end do
    end do

  end subroutine flux_div

!-----------------------------------------

  subroutine compute_ugradp (p, DIMS(p), ugradp, DIMS(ugp), &
                            umac,  DIMS(umac), &
                            vmac,  DIMS(vmac), &
                            lo, hi, dx) &
                            bind(C, name="compute_ugradp")

    implicit none
    integer lo(SDIM), hi(SDIM)
    integer DIMDEC(p)
    integer DIMDEC(ugp)
    integer DIMDEC(umac)
    integer DIMDEC(vmac)
    REAL_T  umac(DIMV(umac))
    REAL_T  vmac(DIMV(vmac))
    REAL_T      p(DIMV(p))
    REAL_T ugradp(DIMV(ugp))
    REAL_T dx(SDIM)

    integer i, j
    REAL_T uadv, vadv
    REAL_T p_x_lo, p_x_hi
    REAL_T p_y_lo, p_y_hi
    
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        uadv = half*(umac(i,j) + umac(i+1,j))
        vadv = half*(vmac(i,j) + vmac(i,j+1))
        p_x_hi = merge(p(i  ,j),p(i+1,j),umac(i+1,j)>=zero)
        p_x_lo = merge(p(i-1,j),p(i  ,j),umac(i  ,j)>=zero)
        p_y_hi = merge(p(i,j  ),p(i,j+1),vmac(i,j+1)>=zero)
        p_y_lo = merge(p(i,j-1),p(i,j  ),vmac(i,j  )>=zero)
        ugradp(i,j) = uadv * (p_x_hi - p_x_lo) / dx(1) + &
                         vadv * (p_y_hi - p_y_lo) / dx(2)      
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
    integer n,i,j, loC(SDIM),hiC(SDIM),ii,jj,iii,jjj,ncells
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

    do j=loC(2),hiC(2)
      do i=loC(1),hiC(1)
            
        bad_T = .false.
        do jj=0,ratio(2)-1
          jjj = ratio(2)*j + jj
          do ii=0,ratio(1)-1
            iii = ratio(1)*i + ii
            if (state(iii,jjj,Tcomp).lt.min_T) then
              bad_T = .true.
            endif
          enddo
        enddo
            
        if (bad_T .eqv. .true.) then

          tmp(Rcomp) = 0.d0
          do n=first_spec,last_spec
            tmp(n) = 0.d0
          enddo
          tmp(RhoH) = 0.d0

          do jj=0,ratio(2)-1
            jjj = ratio(2)*j + jj
            do ii=0,ratio(1)-1
              iii = ratio(1)*i + ii

              tmp(Rcomp) = tmp(Rcomp) + state(iii,jjj,Rcomp)
              do n=first_spec,last_spec
                tmp(n) = tmp(n) + state(iii,jjj,n)
              enddo
              tmp(RhoH) = tmp(RhoH) + state(iii,jjj,RhoH)
                     
            enddo
          enddo

          conservative_T_floor = conservative_T_floor + ncells

          tmp(Rcomp) = tmp(Rcomp) * ncellsInv
          do n=first_spec,last_spec
            tmp(n) = tmp(n) * ncellsInv
          enddo
          tmp(RhoH) = tmp(RhoH)* ncellsInv
               
          do jj=0,ratio(2)-1
            jjj = ratio(2)*j + jj
            do ii=0,ratio(1)-1
              iii = ratio(1)*i + ii
                     
              state(iii,jjj,Rcomp) = tmp(Rcomp)
              do n=first_spec,last_spec
                state(iii,jjj,n) = tmp(n)
              enddo
              state(iii,jjj,RhoH) = tmp(RhoH)
                     
            enddo
          enddo
               
        endif

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

  subroutine part_cnt_err(tag,tagl1,tagl2,tagh1,tagh2, &
                          set,clear, &
                          var,varl1,varl2,varh1,varh2, &
                          lo,hi,nd,domlo,domhi, &
                          delta,xlo,problo,time,level) &
                          bind(C, name="part_cnt_err")

    implicit none

    integer          :: set, clear, nd, level
    integer          :: tagl1,tagl2,tagh1,tagh2
    integer          :: varl1,varl2,varh1,varh2
    integer          :: lo(2), hi(2), domlo(2), domhi(2)
    integer          :: tag(tagl1:tagh1,tagl2:tagh2)
    double precision :: var(varl1:varh1,varl2:varh2,nd)
    double precision :: delta(2), xlo(2), problo(2), time
    integer          :: i,j

    do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (var(i,j,1) .gt. 0.0d0) tag(i,j) = set
         !if (var(i,j,1) .gt. 0.0d0) print *,'TAGGING ',i,j
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

    integer i,j
    REAL_T mag,gTx,gTy

!     Fill normal on nodes (assumes 1 grow cell properly filled)

    do j = lo(2),hi(2)+1
      do i = lo(1),hi(1)+1

        gTx = (T(i,j) + T(i,j-1) - T(i-1,j) - T(i-1,j-1))/(two*delta(1))
        gTy = (T(i,j) - T(i,j-1) + T(i-1,j) - T(i-1,j-1))/(two*delta(2))

        mag = SQRT(gTx*gTx+gTy*gTy)
            
!     Note: This sets the normal to point into the FUEL
        if (mag.eq.zero) then
          wrk(i,j,1) = zero
          wrk(i,j,2) = zero
        else
          wrk(i,j,1) = -gTx/mag
          wrk(i,j,2) = -gTy/mag
        endif
      end do
    end do
!     Get curvature at centers from nodal normal, get normal at centers
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        curv(i,j) = half*((wrk(i+1,j+1,1)-wrk(i,j+1,1))/delta(1) &
                +            (wrk(i+1,j  ,1)-wrk(i,j  ,1))/delta(1) &
                +            (wrk(i+1,j+1,2)-wrk(i+1,j,2))/delta(2) &
                +            (wrk(i  ,j+1,2)-wrk(i  ,j,2))/delta(2))
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

    integer i,j,ii,jj

    do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        Tout(i,j) = zero
        do jj=0,1
          do ii=0,1
            Tout(i,j) = Tout(i,j)  &
                      + Tin(i+ii,j+jj)   + Tin(i+ii-1,j+jj) &
                      + Tin(i+ii,j+jj-1) + Tin(i+ii-1,j+jj-1)
          end do
        end do
        Tout(i,j) = Tout(i,j) / 16.d0
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
    integer i,j

    fac = mult / dx

    if (inc .eq. 0) then

!     compute grad wbar fluxes

      if (dir.eq.0) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            Wgr =   fac*(Wbar(i,j) - Wbar(i-1,j))
            flux(i,j) = rDe(i,j) * Wgr * area(i,j)
          enddo
        enddo
      else if (dir.eq.1) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            Wgr =   fac*(Wbar(i,j) - Wbar(i,j-1))
            flux(i,j) = rDe(i,j) * Wgr * area(i,j)
          enddo
        enddo
      else
        call bl_pd_abort('Bad dir in grad_wbar')
      endif

    else

!     increment grad wbar fluxes by a factor of inc (can be negative)

      if (dir.eq.0) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            Wgr =   fac*(Wbar(i,j) - Wbar(i-1,j))
            flux(i,j) = flux(i,j) + inc * rDe(i,j) * Wgr * area(i,j)
          enddo
        enddo
      else if (dir.eq.1) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            Wgr =   fac*(Wbar(i,j) - Wbar(i,j-1))
            flux(i,j) = flux(i,j) + inc * rDe(i,j) * Wgr * area(i,j)
          enddo
        enddo
      else
            call bl_pd_abort('Bad dir in grad_wbar')
      endif
    end if
      
  end subroutine grad_wbar

!-----------------------------------

  subroutine recomp_update(lo, hi, &
                           update, DIMS(update), &
                           xflux,  DIMS(xflux), &
                           yflux,  DIMS(yflux), &
                           vol,    DIMS(vol), &
                           nc) &
                           bind(C, name="recomp_update")

    implicit none
    integer lo(SDIM), hi(SDIM), nc
    integer DIMDEC(update)
    integer DIMDEC(xflux)
    integer DIMDEC(yflux)
    integer DIMDEC(vol)
    REAL_T update(DIMV(update),nc)
    REAL_T xflux(DIMV(xflux),nc)
    REAL_T yflux(DIMV(yflux),nc)
    REAL_T vol(DIMV(vol))

    integer i, j, n

    do n=1,nc      
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          update(i,j,n)=-((xflux(i+1,j,n)-xflux(i,j,n)) &
                   +          (yflux(i,j+1,n)-yflux(i,j,n)))/vol(i,j)
        end do
      end do
    end do

  end subroutine recomp_update

!-----------------------------------

  subroutine valgt_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="valgt_error")

    implicit none

#include "probdata.H"
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    integer   tag(DIMV(tag))
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.value)
       end do
    end do

  end subroutine valgt_error

  subroutine vallt_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="vallt_error")

    implicit none

#include "probdata.H"
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    integer   tag(DIMV(tag))
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          tag(i,j) = merge(set,tag(i,j),adv(i,j,1).lt.value)
       end do
    end do

  end subroutine vallt_error
  

  subroutine magvort_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="magvort_error")

    implicit none

#include "probdata.H"

    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    integer   tag(DIMV(tag))
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)

         tag(i,j) = merge(set,tag(i,j), &
                    ABS(adv(i,j,1)).ge.value*2.d0**level)

       end do
    end do

  end subroutine magvort_error



  subroutine diffgt_error (tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="diffgt_error")

    implicit none

#include "probdata.H"
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    integer   tag(DIMV(tag))
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value
    REAL_T    axp, axm, ayp, aym, aerr

    integer   i, j


    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          axp = ABS(adv(i+1,j,1) - adv(i,j,1))
          axm = ABS(adv(i-1,j,1) - adv(i,j,1))
          ayp = ABS(adv(i,j+1,1) - adv(i,j,1))
          aym = ABS(adv(i,j-1,1) - adv(i,j,1))
          aerr = MAX(axp,MAX(axm,MAX(ayp,aym)))

          if (aerr.gt.value) then
             tag(i,j) = set
          endif
       end do
    end do
  end subroutine diffgt_error

end module PeleLM_2d
