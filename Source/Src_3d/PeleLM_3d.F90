
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

module PeleLM_3d

  use fuego_chemistry
  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public ::  calc_divu_fortran, calc_gamma_pinv, floor_spec, enth_diff_terms, &
             vel_visc, spec_temp_visc, beta_wbar, &
             est_divu_dt, check_divu_dt, dqrad_fill, divu_fill, &
             dsdt_fill, ydot_fill, rhoYdot_fill, fab_minmax, repair_flux, &
             incrwext_flx_div, flux_div, compute_ugradp, conservative_T_floor, &
             part_cnt_err, mcurve, smooth, grad_wbar, recomp_update, &
             valgt_error, vallt_error, magvort_error, diffgt_error, init_data_new_mech, &
             pphys_PfromRTY, pphys_mass_to_mole, pphys_massr_to_conc, pphys_HfromT, &
             pphys_HMIXfromTY, pphys_RHOfromPTY, FORT_AVERAGE_EDGE_STATES 

contains
             
  subroutine calc_divu_fortran(lo, hi, &
                              divu, DIMS(divu), rYdot, DIMS(rYdot), &
                              vtY,  DIMS(vtY),    vtT, DIMS(vtT), &
                              rhoY, DIMS(rhoY),     T, DIMS(T)) &
                              bind(C, name="calc_divu_fortran")

      use network,        only : nspecies

      implicit none

      integer lo(dim),hi(dim)
      integer DIMDEC(divu)
      integer DIMDEC(rYdot)
      integer DIMDEC(vtY)
      integer DIMDEC(vtT)
      integer DIMDEC(rhoY)
      integer DIMDEC(T)
      REAL_T  divu(DIMV(divu))
      REAL_T  rYdot(DIMV(rYdot),1:nspecies)
      REAL_T  vtY(DIMV(vtY),1:nspecies)
      REAL_T  vtT(DIMV(vtT))
      REAL_T  rhoY(DIMV(rhoY),1:nspecies)
      REAL_T  T(DIMV(T))
      
      integer i, j, k, n
      REAL_T Y(nspecies), H(nspecies), cpmix, rho, rhoInv, tmp, mmw, invmwt(nspecies)

      call CKWT(invmwt)
      do n=1,nspecies
         invmwt(n) = one / invmwt(n)

      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rho = 0.d0
               do n=1,nspecies
                  rho = rho + rhoY(i,j,k,n)
               enddo
               rhoInv = 1.d0 / rho
               do n=1,nspecies

                  Y(n) = rhoInv*rhoY(i,j,k,n)
               enddo
               CALL CKCPBS(T(i,j,k),Y,cpmix)
               CALL CKHMS(T(i,j,k),H)
               CALL CKMMWY(Y,mmw)

               cpmix = cpmix*1.d-4
               do n=1,nspecies
                  H(n) = H(n)*1.d-4
               enddo

               divu(i,j,k) = (divu(i,j,k) + vtT(i,j,k))/(rho*cpmix*T(i,j,k))
               do n=1,nspecies
                  divu(i,j,k) = divu(i,j,k) &
                       + (vtY(i,j,k,n) + rYdot(i,j,k,n)) &
                       *(invmwt(n)*mmw*rhoInv - H(n)/(rho*cpmix*T(i,j,k)))

               enddo
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
                                                        
      use network,        only : nspecies

      implicit none
      
      integer lo(dim),hi(dim)
      integer DIMDEC(theta)
      integer DIMDEC(rhoY)
      integer DIMDEC(T)
      REAL_T  theta(DIMV(theta))
      REAL_T  rhoY(DIMV(rhoY),1:nspecies)
      REAL_T  T(DIMV(T))
      REAL_T  Pamb_in
      
      integer i, j, k, n
      REAL_T Y(nspecies), cpmix, cvmix, rhoInv

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               rhoInv = 0.d0
               do n=1,nspecies
                  rhoInv = rhoInv + rhoY(i,j,k,n)
               enddo
               rhoInv = 1.d0 / rhoInv
               do n=1,nspecies
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

    use network,        only : nspecies

    implicit none
      
      integer lo(dim),hi(dim)
      integer DIMDEC(spec)
      REAL_T  spec(DIMV(spec),1:nspecies)
      
      integer i, j, k, n

      do n=1,nspecies
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
                              T, DIMS(T), &
                              rhoDx, DIMS(rhoDx), Fx, DIMS(Fx), Ax, DIMS(Ax), &
                              rhoDy, DIMS(rhoDy), Fy, DIMS(Fy), Ay, DIMS(Ay), &
                              rhoDz, DIMS(rhoDz), Fz, DIMS(Fz), Az, DIMS(Az), &
                              Tbc ) &
                              bind(C, name="enth_diff_terms")

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim), dlo(dim), dhi(dim), Tbc(dim,2)
      REAL_T  dx(dim)
      integer DIMDEC(T)
      REAL_T  T(DIMV(T))


      integer DIMDEC(rhoDx)
      REAL_T  rhoDx(DIMV(rhoDx))
      integer DIMDEC(Fx)
      REAL_T  Fx(DIMV(Fx),nspecies+3)
      integer DIMDEC(Ax)
      REAL_T  Ax(DIMV(Ax))

      integer DIMDEC(rhoDy)
      REAL_T  rhoDy(DIMV(rhoDy))
      integer DIMDEC(Fy)
      REAL_T  Fy(DIMV(Fy),nspecies+3)
      integer DIMDEC(Ay)
      REAL_T  Ay(DIMV(Ay))

      integer DIMDEC(rhoDz)
      REAL_T  rhoDz(DIMV(rhoDz))
      integer DIMDEC(Fz)
      REAL_T  Fz(DIMV(Fz),nspecies+3)
      integer DIMDEC(Az)
      REAL_T  Az(DIMV(Az))

      REAL_T, allocatable :: H(:,:,:,:)

      integer i, j, k, d, n
      integer lob(dim), hib(dim)
      REAL_T dxInv, dyInv, dzInv

!     Compute species enthalpies on box grown by one
      do d=1,dim
         lob(d) = lo(d)-1
         hib(d) = hi(d)+1
      enddo

!     Make space for Hi, use T box, since this better be big enough as well.
!     Note that any cells on a physical boundary with Dirichlet conditions will 
!     actually be centered on the edge, so the stencils below must reflect this

      allocate( H(DIMV(T),1:nspecies) )

      call pphys_HfromT(lob, hib, H, DIMS(T), T, DIMS(T))

!     On entry, Fx(1:nspecies) = spec flux, Fx(nspecies+1) = rhoh flux (both untouched)
!     On exit:
!     Fx(nspecies+2) = sum[ (species flux).(species enthalpy) ]
!     Fx(nspecies+3) = extensive heat conduction


!     Compute lambda.Grad(T)
      dxInv = 1.d0 / dx(1)
      dyInv = 1.d0 / dx(2)
      dzInv = 1.d0 / dx(3)

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
            Fx(i,j,k,nspecies+3) = - rhoDx(i,j,k)*(T(i,j,k) - T(i-1,j,k)) * dxInv * Ax(i,j,k)
          enddo
        enddo
      enddo

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
            Fy(i,j,k,nspecies+3) = - rhoDy(i,j,k)*(T(i,j,k) - T(i,j-1,k)) * dyInv * Ay(i,j,k)
          enddo
        enddo
      enddo

      do k=lo(3),hi(3)+1
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            Fz(i,j,k,nspecies+3) = - rhoDz(i,j,k)*(T(i,j,k) - T(i,j,k-1)) * dzInv * Az(i,j,k)
          enddo
        enddo
      enddo

!     xlo
      if (lo(1).eq.dlo(1) .and. Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         do k=lo(3),hi(3)
           do j=lo(2),hi(2)
             Fx(i,j,k,nspecies+3) = 2*Fx(i,j,k,nspecies+3)
           enddo
         enddo
      endif
!     xhi
      if (hi(1).eq.dhi(1) .and. Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         do k=lo(3),hi(3)
           do j=lo(2),hi(2)
             Fx(i,j,k,nspecies+3) = 2*Fx(i,j,k,nspecies+3)
           enddo
         enddo
      endif      
!     ylo
      if (lo(2).eq.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
         j=lo(2)
         do k=lo(3),hi(3)
           do i=lo(1),hi(1)
             Fy(i,j,k,nspecies+3) = 2*Fy(i,j,k,nspecies+3)
           enddo
         enddo
      endif
!     yhi
      if (hi(2).eq.dhi(2) .and. Tbc(2,2).eq.EXT_DIR) then
         j=hi(2)+1
         do k=lo(3),hi(3)
           do i=lo(1),hi(1)
             Fy(i,j,k,nspecies+3) = 2*Fy(i,j,k,nspecies+3)
           enddo
         enddo
      endif      
!     zlo
      if (lo(3).eq.dlo(3) .and. Tbc(3,1).eq.EXT_DIR) then
         k=lo(3)
         do j=lo(2),hi(2)
           do i=lo(1),hi(1)
             Fz(i,j,k,nspecies+3) = 2*Fz(i,j,k,nspecies+3)
           enddo
         enddo
      endif
!     zhi
      if (hi(3).eq.dhi(3) .and. Tbc(3,2).eq.EXT_DIR) then
         k=hi(3)+1
         do j=lo(2),hi(2)
           do i=lo(1),hi(1)
             Fz(i,j,k,nspecies+3) = 2*Fz(i,j,k,nspecies+3)
           enddo
         enddo
      endif

!     Compute hi*Fi

      Fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,nspecies+2) = 0.d0
      Fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,nspecies+2) = 0.d0
      Fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,nspecies+2) = 0.d0

      do n=1,nspecies
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
              Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*(H(i,j,k,n)+H(i-1,j,k,n))*0.5d0
            enddo
          enddo
        enddo
      enddo

      do n=1,nspecies
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
              Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n)*(H(i,j,k,n)+H(i,j-1,k,n))*0.5d0
            enddo
          enddo
        enddo
      enddo

      do n=1,nspecies
        do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
            do i=lo(1),hi(1)
              Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*(H(i,j,k,n)+H(i,j,k-1,n))*0.5d0
            enddo
          enddo
        enddo
      enddo

!     xlo
      if (lo(1).eq.dlo(1) .and. Tbc(1,1).eq.EXT_DIR) then
         i = dlo(1)
         Fx(i,lo(2):hi(2),lo(3):hi(3),nspecies+2) = 0.d0
         do n=1,nspecies
           do k=lo(3),hi(3)
             do j=lo(2),hi(2)
               Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*H(i-1,j,k,n)
             enddo
           enddo
         enddo
      endif      
!     xhi
      if (hi(1).eq.dhi(1) .and. Tbc(1,2).eq.EXT_DIR) then
         i = dhi(1)+1
         Fx(i,lo(2):hi(2),lo(3):hi(3),nspecies+2) = 0.d0
         do n=1,nspecies
           do k=lo(3),hi(3)
             do j=lo(2),hi(2)
               Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*H(i,j,k,n)
             enddo
           enddo
         enddo
      endif
!     ylo
      if (lo(2).eq.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
         j = dlo(2)
         Fy(lo(1):hi(1),j,lo(3):hi(3),nspecies+2) = 0.d0
         do n=1,nspecies
           do k=lo(3),hi(3)
             do i=lo(1),hi(1)
               Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n)*H(i,j-1,k,n)
             enddo
           enddo
         enddo
      endif
!     yhi
      if (hi(2).eq.dhi(2) .and. Tbc(2,2).eq.EXT_DIR) then
         j = dhi(2)+1
         Fy(lo(1):hi(1),j,lo(3):hi(3),nspecies+2) = 0.d0
         do n=1,nspecies
           do k=lo(3),hi(3)
             do i=lo(1),hi(1)
               Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n)*H(i,j,k,n)
             enddo
           enddo
         enddo
      endif
!     zlo
      if (lo(3).eq.dlo(3) .and. Tbc(3,1).eq.EXT_DIR) then
         k = dlo(3)
         Fz(lo(1):hi(1),lo(2):hi(2),k,nspecies+2) = 0.d0
         do n=1,nspecies
           do j=lo(2),hi(2)
             do i=lo(1),hi(1)
               Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*H(i,j,k-1,n)
             enddo
           enddo
         enddo
      endif
!     zhi
      if (hi(3).eq.dhi(3) .and. Tbc(3,2).eq.EXT_DIR) then
         k = dhi(3)+1
         Fz(lo(1):hi(1),lo(2):hi(2),k,nspecies+2) = 0.d0
         do n=1,nspecies
           do j=lo(2),hi(2)
             do i=lo(1),hi(1)
               Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*H(i,j,k,n)
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
      
      use network,        only : nspecies
      use PeleLM_F,       only : pphys_calc_src_sdc

      implicit none

      integer lo(dim)
      integer hi(dim)
      integer DIMDEC(RhoY)
      integer DIMDEC(RhoH)
      integer DIMDEC(T)
      integer DIMDEC(RhoYdot)
      REAL_T RhoY(DIMV(RhoY),nspecies)
      REAL_T RhoH(DIMV(RhoH))
      REAL_T T(DIMV(T))
      REAL_T RhoYdot(DIMV(RhoYdot),nspecies)

      REAL_T Zt(nspecies+1),Zdott(nspecies+1)
      REAL_T Temperature, TIME
      integer i,j,k,n

      TIME = 0.

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
 
               Zt(nspecies+1) = RhoH(i,j,k)
               do n=1,nspecies
                  Zt(n) = RhoY(i,j,k,n)
               end do
               Temperature = T(i,j,k)

               call pphys_calc_src_sdc(nspecies,TIME,Temperature,Zt,Zdott)

               do n=1,nspecies
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

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim)
      integer DIMDEC(P)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T P(DIMV(P))
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspecies), RHOt, SCAL, SCAL1
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
!           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
      SCAL = 1.d-1
      SCAL1 = SCAL**3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n=1,nspecies
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

      use network,        only : nspecies

      implicit none

      integer lo(dim)
      integer hi(dim)
      integer DIMDEC(Y)
      integer DIMDEC(X)
      REAL_T Y(DIMV(Y),*)
      REAL_T X(DIMV(X),*)

      REAL_T Xt(nspecies), Yt(nspecies)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,nspecies
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,Xt)
               do n = 1,nspecies
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
                                   
      use network,        only : nspecies

      implicit none

      integer lo(dim)
      integer hi(dim)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      integer DIMDEC(RHO)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      REAL_T RHO(DIMV(RHO))

      REAL_T Yt(nspecies), Ct(nspecies), rhoScl
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,nspecies
                  Yt(n) = Y(i,j,k,n)
               end do
               rhoScl = RHO(i,j,k)*1.e-3
               CALL CKYTCR(rhoScl,T(i,j,k),Yt,Ct)
               do n = 1,nspecies
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

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim)
      integer DIMDEC(HMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T HMIX(DIMV(HMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspecies), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n=1,nspecies
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

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(nspecies), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      CALL CKRP(RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                do n=1,nspecies
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

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim)
      integer DIMDEC(H)
      integer DIMDEC(T)
      REAL_T H(DIMV(H),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, Ht(nspecies)
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),Ht)
               do n=1,nspecies
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do

  end subroutine pphys_HfromT

!-------------------------------------

  subroutine pphys_MWMIXfromY(lo, hi, MWMIX, DIMS(MWMIX), Y, DIMS(Y))&
                        bind(C, name="pphys_MWMIXfromY")

      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim)
      integer DIMDEC(MWMIX)
      integer DIMDEC(Y)
      REAL_T MWMIX(DIMV(MWMIX))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(nspecies)

!     Returns mean molecular weight in kg/kmole

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,nspecies
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

      use network,        only : nspecies
                         
      implicit none

      integer         , intent(in   ) ::     lo(3),      hi(3)
      integer         , intent(in   ) ::   T_lo(3),    T_hi(3)
      integer         , intent(in   ) ::   Y_lo(3),    Y_hi(3)
      integer         , intent(in   ) ::cmix_lo(3),cmix_hi(3)
      REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
      REAL_T, intent(in   ) :: Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspecies)
      REAL_T, intent(out  ) :: CPMIX(cmix_lo(1):cmix_hi(1),cmix_lo(2):cmix_hi(2),cmix_lo(3):cmix_hi(3))
      
      integer i, j, k, n
      REAL_T Yt(nspecies), SCAL
      
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            do n=1,nspecies
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
                           
      use network,        only : nspecies
      use PeleLM_F,       only : pphys_TfromHYpt
      
      implicit none

      integer lo(dim), hi(dim)
      integer NiterMAX
      integer DIMDEC(T)
      integer DIMDEC(HMIX)
      integer DIMDEC(Y)
      REAL_T T(DIMV(T))
      REAL_T HMIX(DIMV(HMIX))
      REAL_T Y(DIMV(Y),*)
      REAL_T errMAX
      REAL_T res(0:NiterMAX-1)

      REAL_T Yt(nspecies), lres(0:NiterMAX-1)
      integer i, j, k, n, Niter,MAXiters

      MAXiters = 0
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,nspecies
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
  
! ::: -----------------------------------------------------------
      
  subroutine init_data_new_mech (level,time,lo,hi,nscal, &
                                 vel,scal,DIMS(state),press,DIMS(press), &
                                 delta,xlo,xhi)&
                                 bind(C, name="init_data_new_mech")
          
      use network,   only: nspecies
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac
      
      implicit none
      integer  level, nscal
      integer  lo(dim), hi(dim)
      integer  DIMDEC(state)
      integer  DIMDEC(press)
      REAL_T   xlo(dim), xhi(dim)
      REAL_T   time, delta(dim)
      REAL_T   vel(DIMV(state),dim)
      REAL_T   scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

      integer i, j, k, n
      REAL_T Patm
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               scal(i,j,k,Trac) = zero
            end do
         end do
      end do
 
      Patm = pamb / pphys_getP1atm_MKS()
 
      call pphys_RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state), &
          Patm)
      call pphys_HMIXfromTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state))
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n = 0,nspecies-1
                  scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
               enddo
               scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
            enddo
         enddo
      enddo
 
  end subroutine init_data_new_mech

!--------------------------------------------------------------

  subroutine vel_visc (lo,hi, &
                       T, T_lo, T_hi, &
                       RY, RY_lo, RY_hi, &
                       mu, mu_lo, mu_hi) &
                       bind(C, name="vel_visc")

    use network,          only : nspecies
    use transport_module, only: get_visco_coeffs
    use mod_Fvar_def, only : LeEQ1

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::   T_lo(3),  T_hi(3)
    integer         , intent(in   ) ::  RY_lo(3), RY_hi(3)
    integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
    REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    REAL_T, intent(in   ) :: RY(RY_lo(1):RY_hi(1),RY_lo(2):RY_hi(2),RY_lo(3):RY_hi(3),nspecies)
    REAL_T, intent(out  ) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))

    REAL_T Yt(lo(1):hi(1),nspecies), RHO_MKS, RHO_MKS_inv
    integer i, j, k, n

    integer bl(3), bh(3)

    bl = 1
    bh = 1
    bl(1) = lo(1)
    bh(1) = hi(1)

    if (.not. LeEQ1) then
       do k=lo(3), hi(3)
          do j=lo(2), hi(2)
             
             do i=lo(1), hi(1)
                RHO_MKS = 0.d0
                do n=1,nspecies
                   RHO_MKS = RHO_MKS + RY(i,j,k,n)
                end do
                RHO_MKS_inv = 1.d0 / RHO_MKS
                do n=1,nspecies
                   Yt(i,n) = RY(i,j,k,n) * RHO_MKS_inv
                end do
             enddo

             call get_visco_coeffs(bl, bh, &
                                   Yt, bl, bh, &
                                   T(lo(1),j,k), bl, bh, &
                                   mu(lo(1),j,k), bl, bh)

             mu(lo(1):hi(1),j,k) = mu(lo(1):hi(1),j,k) * 1.0d-1

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

  end subroutine vel_visc

!-------------------------------------------------

  subroutine beta_wbar (lo, hi, &
                        RD, RD_lo, RD_hi, &
                        RDW, RDW_lo, RDW_hi, &
                        Y, Y_lo, Y_hi) &
                        bind(C, name="beta_wbar")

    use network, only : nspecies
    use fuego_chemistry, only : CKMMWY

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  RD_lo(3), RD_hi(3)
    integer         , intent(inout) :: RDW_lo(3),RDW_hi(3)
    integer         , intent(in   ) ::   Y_lo(3),  Y_hi(3)
    REAL_T, intent(in   ) :: RD(RD_lo(1):RD_hi(1),RD_lo(2):RD_hi(2),RD_lo(3):RD_hi(3),*)
    REAL_T, intent(out  ) :: RDW(RDW_lo(1):RDW_hi(1),RDW_lo(2):RDW_hi(2),RDW_lo(3):RDW_hi(3),*)
    REAL_T, intent(in   ) :: Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspecies)

    integer i, j, k, n
    REAL_T Yt(nspecies), RHO, Wavg

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             RHO = 0.d0
             do n=1,nspecies
                RHO = RHO + Y(i,j,k,n)
             enddo

             do n=1,nspecies
                Yt(n) = Y(i,j,k,n) / RHO
             enddo

             CALL CKMMWY(Yt,Wavg)

             do n=1,nspecies
                RDW(i,j,k,n) = RD(i,j,k,n) * Yt(n) / Wavg
             enddo

          enddo
       enddo
    enddo

  end subroutine beta_wbar

!-------------------------------------------------

  subroutine spec_temp_visc(lo,hi, &
                            T, T_lo, T_hi, &
                            RY, RY_lo, RY_hi, &
                            rhoD, rhoD_lo, rhoD_hi, &
                            ncompd, P1ATM_MKS, do_temp, do_VelVisc, &
                            Pamb_in) &
                            bind(C, name="spec_temp_visc")

    use network, only : nspecies
    use transport_module, only : get_transport_coeffs
    use mod_Fvar_def, only : Pr, Sc, LeEQ1, thickFac
    
    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::   T_lo(3),  T_hi(3)
    integer         , intent(in   ) ::  RY_lo(3), RY_hi(3)
    integer         , intent(in   ) ::rhoD_lo(3),  rhoD_hi(3)
    REAL_T, intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    REAL_T, intent(in   ) :: RY(RY_lo(1):RY_hi(1),RY_lo(2):RY_hi(2),RY_lo(3):RY_hi(3),nspecies)
    REAL_T, intent(out  ) :: rhoD(rhoD_lo(1):rhoD_hi(1),rhoD_lo(2):rhoD_hi(2),rhoD_lo(3):rhoD_hi(3),ncompd)
    integer ncompd, do_temp, do_VelVisc
    REAL_T  Pamb_in, P1ATM_MKS
 
    integer i, j, k, n, nc, ncs
    REAL_T  Patm, Yl(nspecies)
    REAL_T  Yt(lo(1):hi(1),nspecies), invmwt(nspecies), Wavg(lo(1):hi(1))
    REAL_T  Tfac, Yfac, cpmix(1), RHO_MKS, RHO_MKS_inv, RHO_CGS(lo(1):hi(1))
    REAL_T  D(lo(1):hi(1),nspecies), mu(lo(1):hi(1)), lambda(lo(1):hi(1)), xi(lo(1):hi(1))

    integer bl(3), bh(3)

    bl = 1
    bh = 1
    bl(1) = lo(1)
    bh(1) = hi(1)

    ! nspecies+1 is Temp stuff, if requested
    ! nspecies+2 is Velocity stuff, if requested
    nc = nspecies
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

    ! Warning, FORT_VELVISC is called separately from this routine, so if there's
    ! any hacking to be done on viscosity, be sure to catch it there as well.
    Tfac = thickFac / Pr
    Yfac = thickFac / Sc

    call CKWT(invmwt)
    do n=1,nspecies
       invmwt(n) = one / invmwt(n)
    end do

    if (.not. LeEQ1) then
       Patm = Pamb_in / P1ATM_MKS

       do k=lo(3),hi(3) 
          do j=lo(2), hi(2)
             
             do i=lo(1), hi(1)                
                RHO_MKS = 0.d0
                do n=1,nspecies
                   RHO_MKS = RHO_MKS + RY(i,j,k,n)
                end do
                RHO_MKS_inv = 1.d0 / RHO_MKS
                do n=1,nspecies
                   Yt(i,n) = RY(i,j,k,n) * RHO_MKS_inv
                   Yl(n) = Yt(i,n)
                end do
                call CKMMWY(Yl,Wavg(i))
                RHO_CGS(i) = RHO_MKS * 1.d-3
             enddo

             ! Get transport coefficients over vector in i
             call get_transport_coeffs(bl, bh, &
                                       Yt,       bl, bh, &
                                       T(lo(1),j,k), bl, bh, &
                                       RHO_CGS,      bl, bh, &
                                       D,            bl, bh, &
                                       mu,           bl, bh, &
                                       xi,           bl, bh, &
                                       lambda,       bl, bh)
             
             do i=lo(1), hi(1)                
                do n=1,nspecies
                   rhoD(i,j,k,n) = Wavg(i) * invmwt(n) * D(i,n)  * 1.0d-1 
                end do
                if (do_temp .ne. 0) then 
                   rhoD(i,j,k,nspecies+1) = lambda(i) * 1.d-5
                end if
                if (thickFac.ne.1.d0) then
                   do n=1,nspecies+1
                      rhoD(i,j,k,n) = rhoD(i,j,k,n)*thickFac
                   enddo
                endif
                if (do_VelVisc .ne. 0) then 
                   rhoD(i,j,k,nspecies+2) = mu(i) * 1.0d-1
                end if
             end do
             
          end do
       end do
      
    else
    
       call vel_visc(lo, hi,&
                     T,    T_lo,    T_hi, &
                     RY,   RY_lo,   RY_hi, &
                     rhoD, rhoD_lo, rhoD_hi)

       ! Set rhoD[1] = mu * Yfac
       do k=lo(3),hi(3) 
          do j=lo(2), hi(2)
             do i=lo(1), hi(1)
                do n=nspecies+1,nc
                   rhoD(i,j,k,n) = rhoD(i,j,k,1)
                end do
                rhoD(i,j,k,1) = rhoD(i,j,k,1) * Yfac
             end do
          end do
       end do
       ! Set rhoD[2:N] = rhoD[1]
       do n=2,nspecies
          do k=lo(3),hi(3) 
             do j=lo(2), hi(2)
                do i=lo(1), hi(1)
                   rhoD(i,j,k,n) = rhoD(i,j,k,1)
                end do
             end do
          end do
       end do
       
       if (do_temp .ne. 0) then
          ! Set lambda = mu * cpmix, ptwise
          do k=lo(3),hi(3) 
             do j=lo(2), hi(2)               
                do i=lo(1), hi(1)
                   RHO_MKS = 0.d0
                   do n=1,nspecies
                      RHO_MKS = RHO_MKS + RY(i,j,k,n)
                   end do
                   RHO_MKS_inv = 1.d0 / RHO_MKS
                   do n=1,nspecies
                      Yl(n) = RY(i,j,k,n) * RHO_MKS_inv
                   end do
                   CALL pphys_CPMIXfromTY(bl, bl, &
                                          cpmix,    bl, bl, &
                                          T(i,j,k), bl, bl, &
                                          Yl,       bl, bl)
                   rhoD(i,j,k,nspecies+1) = rhoD(i,j,k,nspecies+1)*cpmix(1)*Tfac
                end do
             end do
          end do
       end if
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
      integer lo(dim), hi(dim)
      REAL_T  delta(dim)
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
      integer lo(dim), hi(dim)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(u)
      REAL_T  delta(dim)
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
      integer    bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     dqrad(DIMV(dqrad))

      call filcc (dqrad,DIMS(dqrad),domlo,domhi,delta,xlo,bc)
      call fillEdges(dqrad,DIMS(dqrad),domlo,domhi,bc)
      
  end subroutine dqrad_fill

!--------------------------------------------------

  subroutine divu_fill  (divu,DIMS(divu),domlo,domhi,delta, &
                         xlo,time,bc )bind(C, name="divu_fill")

      integer    DIMDEC(divu)
      integer    bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     divu(DIMV(divu))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      call fillEdges(divu,DIMS(divu),domlo,domhi,bc)

  end subroutine divu_fill

!--------------------------------------------

  subroutine dsdt_fill (dsdt,DIMS(dsdt),domlo,domhi,delta, &
                        xlo,time,bc )bind(C, name="dsdt_fill")

      integer    DIMDEC(dsdt)
      integer    bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     dsdt(DIMV(dsdt))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(dsdt,DIMS(dsdt),domlo,domhi,bc)

  end subroutine dsdt_fill

!--------------------------------------

  subroutine ydot_fill (ydot,DIMS(ydot),domlo,domhi,delta, &
                        xlo,time,bc)bind(C, name="ydot_fill")

      integer    DIMDEC(ydot), bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     ydot(DIMV(ydot))

      call filcc (ydot,DIMS(ydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(ydot,DIMS(ydot),domlo,domhi,bc)

  end subroutine ydot_fill

!-------------------------------------------

  subroutine rhoYdot_fill (rhoydot,DIMS(rhoydot),domlo,domhi,delta, &
                           xlo,time,bc)bind(C, name="rhoYdot_fill")

      integer    DIMDEC(rhoydot), bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     rhoydot(DIMV(rhoydot))

      call filcc (rhoydot,DIMS(rhoydot),domlo,domhi,delta,xlo,bc)

      call fillWithZeros(rhoydot,DIMS(rhoydot),domlo,domhi,bc)

  end subroutine rhoYdot_fill

!------------------------------------------------

  subroutine fab_minmax(lo, hi, &
                        fab, DIMS(fab), &
                        fmin, fmax, nc) &
                        bind(C, name="fab_minmax")
                        
      integer lo(dim), hi(dim), nc
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
                          
      use network,        only : nspecies

      implicit none

      integer lo(dim), hi(dim), dlo(dim), dhi(dim), dir, Ybc(dim,2)
      integer DIMDEC(flux)
      integer DIMDEC(RhoY)
      REAL_T flux(DIMV(flux),nspecies)
      REAL_T RhoY(DIMV(RhoY),nspecies)
      
      integer i, j, k, n
      REAL_T sumFlux, RhoYe(nspecies), sumRhoYe

      if (dir.eq.0) then

!     First, assume away from physical boundaries, then use boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i-1,j,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  end do
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i-1,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
                  do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j-1,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j-1,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
                  do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j,k-1,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k-1,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
                     do n=1,nspecies
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,nspecies
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
      integer lo(dim), hi(dim), nc
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
      integer lo(dim), hi(dim), nc
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
      integer    bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     dat(DIMV(dat))

      integer    lo(dim), hi(dim)
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
      integer    bc(dim,2)
      integer    domlo(dim), domhi(dim)
      REAL_T     dat(DIMV(dat))

      integer    lo(dim), hi(dim)
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
      integer lo(dim), hi(dim)
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
      REAL_T dx(dim)

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

      integer loF(dim), hiF(dim)
      integer DIMDEC(state)
      REAL_T  state(DIMV(state),0:*)      
      integer Tcomp, Rcomp, first_spec, last_spec, RhoH, ratio(dim), nt
      REAL_T  min_T, tmp(0:nt-1)
      integer n,i,j,k, loC(dim),hiC(dim),ii,jj,kk,iii,jjj,kkk,ncells
      Real ncellsInv
      logical bad_T

!     Returns the number of fine cells fixed up
      conservative_T_floor = 0

      ncells = 1
      do n=1,dim
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

      integer lo(dim), hi(dim)
      integer DIMDEC(T)
      integer DIMDEC(curv)
      integer DIMDEC(wrk)
      REAL_T    T(DIMV(T))
      REAL_T curv(DIMV(curv))
      REAL_T wrk(DIMV(wrk),dim)
      REAL_T delta(dim)

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
      integer lo(dim), hi(dim)
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
      
      integer lo(dim), hi(dim)
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
      
      integer lo(dim), hi(dim), nc
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

!-----------------------------------

  subroutine valgt_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="valgt_error")

    implicit none
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(dim), domhi(dim)
    integer   lo(dim), hi(dim)
    integer   tag(DIMV(tag))
    REAL_T    dx(dim), xlo(dim), problo(dim), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.value)
          end do
       end do
    end do

  end subroutine valgt_error

  subroutine vallt_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="vallt_error")

    implicit none
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(dim), domhi(dim)
    integer   lo(dim), hi(dim)
    integer   tag(DIMV(tag))
    REAL_T    dx(dim), xlo(dim), problo(dim), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).lt.value)
          end do
       end do
    end do

  end subroutine vallt_error

    subroutine magvort_error(tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="magvort_error")

    implicit none

    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(dim), domhi(dim)
    integer   lo(dim), hi(dim)
    integer   tag(DIMV(tag))
    REAL_T    dx(dim), xlo(dim), problo(dim), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value

    integer   i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

              tag(i,j,k) = merge(set,tag(i,j,k), &
                    ABS(adv(i,j,k,1)).ge.value*2.d0**level)

          end do
       end do
    end do

  end subroutine magvort_error  

  subroutine diffgt_error (tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),&
       lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level,value) bind(C, name="diffgt_error")

    implicit none
      
    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(dim), domhi(dim)
    integer   lo(dim), hi(dim)
    integer   tag(DIMV(tag))
    REAL_T    dx(dim), xlo(dim), problo(dim), time
    REAL_T    adv(DIMV(adv),1)
    REAL_T    value
    REAL_T    axp, axm, ayp, aym, azp, azm, aerr

    integer   i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             axp = ABS(adv(i+1,j,k,1) - adv(i,j,k,1))
             axm = ABS(adv(i-1,j,k,1) - adv(i,j,k,1))
             ayp = ABS(adv(i,j+1,k,1) - adv(i,j,k,1))
             aym = ABS(adv(i,j-1,k,1) - adv(i,j,k,1))
             azp = ABS(adv(i,j,k+1,1) - adv(i,j,k,1))
             azm = ABS(adv(i,j,k-1,1) - adv(i,j,k,1))
             aerr = MAX(azp,MAX(azm,MAX(axp,MAX(axm,MAX(ayp,aym)))))

             if (aerr.gt.value) then
                tag(i,j,k) = set
             endif
          end do
       end do
    end do
  end subroutine diffgt_error

  subroutine box_error (tag,DIMS(tag),set,clear,&
       boxlo,boxhi,lo,hi,&
       domlo,domhi,dx,xlo,&
       problo,time,level) bind(C, name="box_error")

    implicit none
      
    integer   DIMDEC(tag)
    integer   set, clear, level
    integer   domlo(dim), domhi(dim)
    integer   lo(dim), hi(dim)
    integer   tag(DIMV(tag))
    REAL_T    boxlo(dim),boxhi(dim),dx(dim), xlo(dim), problo(dim), time

    integer   i, j, k
    REAL_T    x, y, z


    do k = lo(3), hi(3)
       z = (float(k)+.5)*dx(3)+problo(3)
       if (z.ge.boxlo(3) .and. z.le.boxhi(3)) then
          do j = lo(2), hi(2)
             y = (float(j)+.5)*dx(2)+problo(2)
             if (y.ge.boxlo(2) .and. y.le.boxhi(2)) then
                do i = lo(1), hi(1)
                   x = (float(i)+.5)*dx(1)+problo(1)
                   if (x.ge.boxlo(1) .and. x.le.boxhi(1)) then
                      tag(i,j,k) = set
                   endif
                end do
             endif
          end do
       endif
    end do
  end subroutine box_error

!c     
!c     
!c     ::: -----------------------------------------------------------
!c     
!c     This routine averages the mac face velocities for makeforce at half time
!c
      subroutine FORT_AVERAGE_EDGE_STATES(vel,umacx,umacy,umacz, &
                                         DIMS(vel),DIMS(umacx),DIMS(umacy),DIMS(umacz),&
                                         getForceVerbose)&
                                         bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer    DIMDEC(vel)
      integer    DIMDEC(umacx)
      integer    DIMDEC(umacy)
      integer    DIMDEC(umacz)
      integer    getForceVerbose
      REAL_T     vel  (DIMV(vel),dim)
      REAL_T     umacx(DIMV(umacx))
      REAL_T     umacy(DIMV(umacy))
      REAL_T     umacz(DIMV(umacz))

      integer i,j,k,n
      integer ilo,jlo,klo
      integer ihi,jhi,khi

      integer isioproc

      REAL_T  velmin(3)
      REAL_T  velmax(3)

      do n = 1, 3
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      ilo = vel_l1
      jlo = vel_l2
      klo = vel_l3
      ihi = vel_h1
      jhi = vel_h2
      khi = vel_h3

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
               do n=1, 3
                  velmin(n)=min(velmin(n),vel(i,j,k,n))
                  velmax(n)=max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, 3
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

      end subroutine FORT_AVERAGE_EDGE_STATES

end module PeleLM_3d
