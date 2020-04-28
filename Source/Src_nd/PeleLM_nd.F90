#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif 

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module PeleLM_nd

  use fuego_chemistry
  use amrex_fort_module,  only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_filcc_module, only : amrex_filccn
  use network,            only : nspecies

  implicit none

  private

  public :: floor_spec, calc_divu_fortran, calc_gamma_pinv, &
            pphys_PfromRTY, pphys_mass_to_mole, pphys_massr_to_conc, pphys_HfromT, &
            pphys_HMIXfromTY, pphys_RHOfromPTY, pphys_CPMIXfromTY, init_data_new_mech, &
            spec_temp_visc, vel_visc, beta_wbar, est_divu_dt, check_divu_dt,&
            dqrad_fill, divu_fill, dsdt_fill, ydot_fill, rhoYdot_fill, &
            fab_minmax, repair_flux, incrwext_flx_div, flux_div, compute_ugradp, conservative_T_floor, &
            part_cnt_err, mcurve, smooth, grad_wbar, recomp_update, &
            valgt_error, vallt_error, magvort_error, diffgt_error, &
            FORT_AVERAGE_EDGE_STATES

contains
 
!=========================================================
!  Floor species mass fraction at 0.0
!=========================================================

   subroutine floor_spec(lo, hi, &
                         spec, s_lo, s_hi) &
                         bind(C, name="floor_spec")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: s_lo(3), s_hi(3)
      REAL_T, intent(inout), dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nspecies) :: spec

! Local
      integer :: i, j, k, n

      do n = 1, nspecies
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  spec(i,j,k,n) = max(0.d0,spec(i,j,k,n))
               end do
            enddo
         enddo
      enddo

   end subroutine floor_spec

!=========================================================
!  Compute DivU by gathering its components
!=========================================================

   subroutine calc_divu_fortran( lo, hi, &
                                 divu, d_lo, d_hi, &
                                 rYdot, rd_lo, rd_hi, &
                                 vtY, v_lo, v_hi, &
                                 vtT, vT_lo, vT_hi, &
                                 rhoY, rY_lo, rY_hi, &
                                 T, t_lo, t_hi) &
                                 bind(C, name="calc_divu_fortran")

      implicit none

! In/Out
      integer, intent(in) :: lo(3),hi(3)
      integer, intent(in) :: d_lo(3), d_hi(3)
      integer, intent(in) :: rd_lo(3), rd_hi(3)
      integer, intent(in) :: v_lo(3), v_hi(3)
      integer, intent(in) :: vT_lo(3), vT_hi(3)
      integer, intent(in) :: rY_lo(3), rY_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: divu
      REAL_T, dimension(rd_lo(1):rd_hi(1),rd_lo(2):rd_hi(2),rd_lo(3):rd_hi(3),nspecies) :: rYdot
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nspecies) :: vtY
      REAL_T, dimension(vT_lo(1):vT_hi(1),vT_lo(2):vT_hi(2),vT_lo(3):vT_hi(3)) :: vtT
      REAL_T, dimension(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies) :: rhoY
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T

! Local
      REAL_T, dimension(1:nspecies) :: Y, H, invmtw
      REAL_T :: cpmix, rho, rhoInv, tmp, mmw
      integer :: i, j, k, n

      call CKWT(invmtw)
      do n = 1,nspecies
         invmtw(n) = one / invmtw(n)
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

               tmp = 0.d0
               
               !write(*,*) "DEBUG IN calc_divu_fortran ",i,j,k,divu(i,j,k),vtT(i,j,k),rho,cpmix,T(i,j,k)
               
               divu(i,j,k) = (divu(i,j,k) + vtT(i,j,k))/(rho*cpmix*T(i,j,k))
               do n=1,nspecies
               !write(*,*) "DEBUG n ",n,divu(i,j,k),vtY(i,j,k,n),rYdot(i,j,k,n),invmtw(n),H(n)
                  divu(i,j,k) = divu(i,j,k) &
                              + (vtY(i,j,k,n) + rYdot(i,j,k,n)) &
                              *(invmtw(n)*mmw*rhoInv - H(n)/(rho*cpmix*T(i,j,k)))
               enddo
            enddo
         enddo
      enddo

   end subroutine calc_divu_fortran

!=========================================================
!  Compute theta = 1.0 / (\gamma * P )
!=========================================================

   subroutine calc_gamma_pinv(lo, hi, &
                              theta, th_lo, th_hi, &
                              rhoY, r_lo, r_hi, &
                              T, t_lo, t_hi, &
                              Pamb_in)&
                              bind(C, name="calc_gamma_pinv")

      implicit none

! In/Out
      integer, intent(in) :: lo(3),hi(3)
      integer, intent(in) :: th_lo(3), th_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      REAL_T, dimension(th_lo(1):th_hi(1),th_lo(2):th_hi(2),th_lo(3):th_hi(3)) :: theta
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspecies) :: rhoY
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T  :: Pamb_in

! Local
      REAL_T, dimension(1:nspecies) :: Y
      REAL_T :: cpmix, cvmix, rhoInv
      integer :: i, j, k, n

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

!=========================================================
!  Compute reaction rate rhoY source terms
!=========================================================

   subroutine pphys_RRATERHOY(lo, hi, &
                              RhoY, rY_lo, rY_hi, &
                              RhoH, rh_lo, rh_hi, &
                              T, t_lo, t_hi, &
                              mask, m_lo, m_hi, &
                              RhoYdot, rd_lo, rd_hi)&
                              bind(C, name="pphys_RRATERHOY")

      use PeleLM_F,       only : pphys_calc_src_sdc

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: rY_lo(3), rY_hi(3)
      integer, intent(in) :: rh_lo(3), rh_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: m_lo(3), m_hi(3)
      integer, intent(in) :: rd_lo(3), rd_hi(3)
      REAL_T, dimension(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies) :: RhoY
      REAL_T, dimension(rh_lo(1):rh_hi(1),rh_lo(2):rh_hi(2),rh_lo(3):rh_hi(3)) :: RhoH
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      INTEGER, dimension(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3)) :: mask
      REAL_T, dimension(rd_lo(1):rd_hi(1),rd_lo(2):rd_hi(2),rd_lo(3):rd_hi(3),nspecies) :: RhoYdot

! Local
      REAL_T  :: Zt(nspecies+1), Zdott(nspecies+1)
      REAL_T  :: Temperature, TIME

      integer :: i, j, k, n

      TIME = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( mask(i,j,k) == -1 ) then
                  RhoYdot(i,j,k,:) = 0.0d0
                  CYCLE
               end if

               Zt(nspecies+1) = RhoH(i,j,k)
               do n = 1,nspecies
                  Zt(n) = RhoY(i,j,k,n)
               end do
               Temperature = T(i,j,k)

               call pphys_calc_src_sdc(nspecies,TIME,Temperature,Zt,Zdott)

               do n = 1,nspecies
                  RhoYdot(i,j,k,n) = Zdott(n)
               end do
            end do
         end do
      end do

   end subroutine pphys_RRATERHOY

!=========================================================
!  Compute P from rho, rhoY and T
!=========================================================

   subroutine pphys_PfromRTY(lo, hi, &
                             P, p_lo, p_hi, &
                             Rho, r_lo, r_hi, &
                             T, t_lo, t_hi, &
                             Y, y_lo, y_hi)&
                             bind(C, name="pphys_PfromRTY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: p_lo(3), p_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: P
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: Rho
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T

! Local
      REAL_T :: Yt(nspecies), RHOt, SCAL, SCAL1
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
!           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
      SCAL = 1.d-1
      SCAL1 = SCAL**3

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n = 1, nspecies
                   Yt(n) = Y(i,j,k,n)
                end do

                RHOt = RHO(i,j,k) * SCAL1
                CALL CKPY(RHOt,T(i,j,k),Yt,P(i,j,k))

                P(i,j,k) = P(i,j,k) * SCAL

            end do
         end do
      end do

   end subroutine pphys_PfromRTY

!=========================================================
!  Compute Xm from Ym
!=========================================================

   subroutine pphys_mass_to_mole(lo, hi, &
                                 Y, y_lo, y_hi, &
                                 X, x_lo, x_hi) &
                                 bind(C, name="pphys_mass_to_mole")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      integer, intent(in) :: x_lo(3), x_hi(3)
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y
      REAL_T, dimension(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nspecies) :: X

! Local
      REAL_T  :: Xt(nspecies), Yt(nspecies)
      integer :: i, j, k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1, nspecies
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,Xt)
               do n = 1, nspecies
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      
   end subroutine pphys_mass_to_mole

!=========================================================
!  Compute Cm from Ym
!=========================================================

   subroutine pphys_massr_to_conc(lo, hi, &
                                  Y, y_lo, y_hi, &
                                  T, t_lo, t_hi, &
                                  Rho, r_lo, r_hi, &
                                  C, c_lo, c_hi)&
                                  bind(C, name="pphys_massr_to_conc")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: c_lo(3), c_hi(3)
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: Rho
      REAL_T, dimension(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nspecies) :: C

! Local
      REAL_T  :: Yt(nspecies), Ct(nspecies), rhoScl
      integer :: i, j, k, n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1, nspecies
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

!=========================================================
!  Compute mixture enthalpy from T and Y
!=========================================================

   subroutine pphys_HMIXfromTY(lo, hi, &
                               Hmix, h_lo, h_hi, & 
                               T, t_lo, t_hi, &
                               Y, y_lo, y_hi)&
                               bind(C, name="pphys_HMIXfromTY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: h_lo(3), h_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3)) :: Hmix
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y

! Local
      REAL_T  :: Yt(nspecies), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n = 1, nspecies
                   Yt(n) = Y(i,j,k,n)
                end do

                CALL CKHBMS(T(i,j,k),Yt,HMIX(i,j,k))

                HMIX(i,j,k) = HMIX(i,j,k) * SCAL

            end do
         end do
      end do

   end subroutine pphys_HMIXfromTY

!=========================================================
!  Compute rho from P, T and Y
!=========================================================

   subroutine pphys_RHOfromPTY(lo, hi, &
                               Rho, r_lo, r_hi, &
                               T, t_lo, t_hi, &
                               Y, y_lo, y_hi, Patm) &
                               bind(C, name="pphys_RHOfromPTY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: Rho
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y
      REAL_T  :: Patm

! Local
      REAL_T  :: RU, RUC, P1ATM, Ptmp, Yt(nspecies), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      CALL CKRP(RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                do n = 1, nspecies
                   Yt(n) = Y(i,j,k,n)
                end do
                CALL CKRHOY(Ptmp,T(i,j,k),Yt,RHO(i,j,k))
                RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do

   end subroutine pphys_RHOfromPTY

!=========================================================
!  Compute species H from T
!=========================================================

   subroutine pphys_HfromT(lo, hi, &
                           H, h_lo, h_hi, &
                           T, t_lo, t_hi)&
                           bind(C, name="pphys_HfromT")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: h_lo(3), h_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      REAL_T, dimension(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3),nspecies) :: H
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T

! Local
      REAL_T  :: SCAL, Ht(nspecies)
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),Ht)
               do n = 1, nspecies
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do

   end subroutine pphys_HfromT

!=========================================================
!  Compute mixture mean molecular weight from Y
!=========================================================

   subroutine pphys_MWMIXfromY(lo, hi, &
                               MWmix, m_lo, m_hi, &
                               Y, y_lo, y_hi)&
                               bind(C, name="pphys_MWMIXfromY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: m_lo(3), m_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3)) :: MWmix
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies) :: Y

! Local
      REAL_T  :: Yt(nspecies)
      integer :: i, j, k, n

!     Returns mean molecular weight in kg/kmole

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n = 1, nspecies
                  Yt(n) = Y(i,j,k,n)
               end do

               CALL CKMMWY(Yt,MWMIX(i,j,k))

            end do
         end do
      end do

   end subroutine pphys_MWMIXfromY

!=========================================================
!  Compute mixture mean heat capacity
!=========================================================

  subroutine pphys_CPMIXfromTY(lo, hi, &
                               CPmix, c_lo, c_hi, &
                               T, t_lo, t_hi, &
                               Y, y_lo, y_hi )&
                               bind(C,name="pphys_CPMIXfromTY")

      implicit none

! In/Out
      integer, intent(in) :: lo(3),   hi(3)
      integer, intent(in) :: c_lo(3), c_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3)), intent(out) :: CPmix
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)), intent(in) :: T
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies), intent(in) :: Y

! Local
      REAL_T  :: Yt(nspecies), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
      SCAL = 1.0d-4

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            do n = 1, nspecies
               Yt(n) = Y(i,j,k,n)
            end do
            CALL CKCPBS(T(i,j,k),Yt,CPMIX(i,j,k))
            CPMIX(i,j,k) = CPMIX(i,j,k) * SCAL
          end do
        end do
      end do

   end subroutine pphys_CPMIXfromTY

!=========================================================
!  Iterate on T until it matches Hmix and Ym
!=========================================================

   integer function pphys_TfromHY(lo, hi, &
                                  T, t_lo, t_hi, &
                                  Hmix, h_lo, h_hi, &
                                  Y, y_lo, y_hi, &
                                  errMax, NiterMAX, res) &
                                  bind(C, name="pphys_TfromHY")

      use PeleLM_F, only : pphys_TfromHYpt

      implicit none

! In/Out
      integer, intent(in) :: lo(3),   hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: h_lo(3), h_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      integer, intent(in) :: NiterMAX
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)), intent(out) :: T
      REAL_T, dimension(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3)), intent(in) :: Hmix
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nspecies), intent(in) :: Y
      REAL_T, intent(in)  :: errMAX
      REAL_T, dimension(0:NiterMAX-1), intent(out)  :: res

! Local
      REAL_T :: Yt(nspecies)
      integer :: i, j, k, n, Niter, MAXiters

      MAXiters = 0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               do n = 1, nspecies
                  Yt(n) = Y(i,j,k,n)
               end do

               CALL pphys_TfromHYpt( T(i,j,k), Hmix(i,j,k), Yt, errMax, NiterMAX, res, Niter)

               if (Niter < 0) then
                  write(*,*)"Error code ", Niter
                  call amrex_abort(" Something went wrong in pphys_TfromHYpt ")
               end if

               if (Niter > MAXiters) then
                  MAXiters = Niter
               end if

            end do
         end do
      end do

!     Set max iters taken during this solve, and exit
      pphys_TfromHY = MAXiters

      return

   end function pphys_TfromHY

!=========================================================
!  init data to transition between 2 chem mechanism.
!=========================================================

   subroutine init_data_new_mech (level, time, lo, hi, nscal, &
                                  vel, scal, s_lo, s_hi, &
                                  press, p_lo, p_hi, &
                                  delta, xlo, xhi)&
                                  bind(C, name="init_data_new_mech")

      use PeleLM_F,  only: pphys_getP1atm_MKS
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb

      implicit none

! In/Out
      integer, intent(in) :: level, nscal
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: s_lo(3), s_hi(3)
      integer, intent(in) :: p_lo(3), p_hi(3)
      REAL_T, intent(in)  :: xlo(3), xhi(3)
      REAL_T, intent(in)  :: time, delta(3)
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),dim), intent(out) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal), intent(out) :: scal
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)), intent(out) :: press

! Local
      REAL_T  :: Patm
      integer :: i, j, k, n

      Patm = pamb / pphys_getP1atm_MKS()

      call pphys_RHOfromPTY(lo,hi, &
                            scal(:,:,:,Density),   s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi, &
                            Patm)
      call pphys_HMIXfromTY(lo,hi, &
                            scal(:,:,:,RhoH),      s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi)
 
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

!=========================================================
!  Compute transport properties
!=========================================================

   subroutine spec_temp_visc(lo,hi, &
                             T, t_lo, t_hi, &
                             RhoY, rY_lo, rY_hi, &
                             RhoD, rd_lo, rd_hi, &
                             ncompd, P1ATM_MKS, do_temp, do_VelVisc, &
                             Pamb_in) &
                             bind(C, name="spec_temp_visc")

    use transport_module, only : get_transport_coeffs
    use mod_Fvar_def, only : Pr, Sc, LeEQ1, thickFac

    implicit none

! In/Out
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: rY_lo(3), rY_hi(3)
    integer, intent(in) :: rd_lo(3), rd_hi(3)
    integer, intent(in) :: ncompd, do_temp, do_VelVisc
    REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)), intent(in) :: T
    REAL_T, dimension(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies), intent(in) :: RhoY
    REAL_T, dimension(rd_lo(1):rd_hi(1),rd_lo(2):rd_hi(2),rd_lo(3):rd_hi(3),ncompd), intent(out) :: RhoD
    REAL_T, intent(in) :: Pamb_in, P1ATM_MKS

! Local
    integer :: i, j, k, n, nc, ncs
    REAL_T  :: Patm, Yl(nspecies)
    REAL_T  :: Yt(lo(1):hi(1),nspecies), invmwt(nspecies), Wavg(lo(1):hi(1))
    REAL_T  :: Tfac, Yfac, cpmix(1), RHO_MKS, RHO_MKS_inv, RHO_CGS(lo(1):hi(1))
    REAL_T  :: D(lo(1):hi(1),nspecies), mu(lo(1):hi(1)), lambda(lo(1):hi(1)), xi(lo(1):hi(1))

    integer :: bl(3), bh(3)

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
                   RHO_MKS = RHO_MKS + RhoY(i,j,k,n)
                end do
                RHO_MKS_inv = 1.d0 / RHO_MKS
                do n = 1, nspecies
                   Yt(i,n) = RhoY(i,j,k,n) * RHO_MKS_inv
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
                     T,    t_lo,  t_hi, &
                     RhoY, rY_lo, rY_hi, &
                     rhoD, rd_lo, rd_hi)

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
       
       if (do_temp /= 0) then
          ! Set lambda = mu * cpmix, ptwise
          do k=lo(3),hi(3) 
             do j=lo(2), hi(2)
                do i=lo(1), hi(1)
                   RHO_MKS = 0.d0
                   do n=1,nspecies
                      RHO_MKS = RHO_MKS + RhoY(i,j,k,n)
                   end do
                   RHO_MKS_inv = 1.d0 / RHO_MKS
                   do n=1,nspecies
                      Yl(n) = RhoY(i,j,k,n) * RHO_MKS_inv
                   end do
                   CALL pphys_CPMIXfromTY(bl, bh, &
                                          cpmix,    bl, bh, &
                                          T(i,j,k), bl, bh, &
                                          Yl,       bl, bh)
                   rhoD(i,j,k,nspecies+1) = rhoD(i,j,k,nspecies+1)*cpmix(1)*Tfac
                end do
             end do
          end do
       end if
    end if

   end subroutine spec_temp_visc

!=========================================================
!  Compute enthalpy diffusion term
!=========================================================

!   subroutine enth_diff_terms (lo, hi, &
!                               dlo, dhi, dx, &
!                               T, t_lo, t_hi, &
!
!                               rhoDx, rdx_lo, rdx_hi, Fx, fx_lo, fx_hi, Ax, ax_lo, ax_hi, &
!                               rhoDy, rdy_lo, rdy_hi, Fy, fy_lo, fy_hi, Ay, ay_lo, ay_hi, &
!#if ( AMREX_SPACEDIM == 3)
!                               rhoDz, rdz_lo, rdz_hi, Fz, fz_lo, fz_hi, Az, az_lo, az_hi, &
!#endif
!                               Tbc) &
!                               bind(C, name="enth_diff_terms")
!
!      use amrex_mempool_module
!
!      implicit none
!
!! In/Out
!      integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
!      REAL_T, intent(in)  :: dx(3)
!      integer, intent(in) :: t_lo(3), t_hi(3)
!      integer, intent(in) :: rdx_lo(3), rdx_hi(3)
!      integer, intent(in) :: fx_lo(3), fx_hi(3)
!      integer, intent(in) :: ax_lo(3), ax_hi(3)
!      integer, intent(in) :: rdy_lo(3), rdy_hi(3)
!      integer, intent(in) :: fy_lo(3), fy_hi(3)
!      integer, intent(in) :: ay_lo(3), ay_hi(3)
!#if ( AMREX_SPACEDIM == 3)
!      integer, intent(in) :: rdz_lo(3), rdz_hi(3)
!      integer, intent(in) :: fz_lo(3), fz_hi(3)
!      integer, intent(in) :: az_lo(3), az_hi(3)
!#endif
!      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)), intent(in) :: T
!      REAL_T, dimension(rdx_lo(1):rdx_hi(1),rdx_lo(2):rdx_hi(2),rdx_lo(3):rdx_hi(3)), intent(in) :: rhoDx
!      REAL_T, dimension(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nspecies+3), intent(out) :: Fx
!      REAL_T, dimension(ax_lo(1):ax_hi(1),ax_lo(2):ax_hi(2),ax_lo(3):ax_hi(3)), intent(in) :: Ax
!      REAL_T, dimension(rdy_lo(1):rdy_hi(1),rdy_lo(2):rdy_hi(2),rdy_lo(3):rdy_hi(3)), intent(in) :: rhoDy
!      REAL_T, dimension(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),nspecies+3), intent(out) :: Fy
!      REAL_T, dimension(ay_lo(1):ay_hi(1),ay_lo(2):ay_hi(2),ay_lo(3):ay_hi(3)), intent(in) :: Ay
!#if ( AMREX_SPACEDIM == 3)
!      REAL_T, dimension(rdz_lo(1):rdz_hi(1),rdz_lo(2):rdz_hi(2),rdz_lo(3):rdz_hi(3)), intent(in) :: rhoDz
!      REAL_T, dimension(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),nspecies+3), intent(out) :: Fz
!      REAL_T, dimension(az_lo(1):az_hi(1),az_lo(2):az_hi(2),az_lo(3):az_hi(3)), intent(in) :: Az
!#endif
!      integer, intent(in) :: Tbc(3,2)
!
!! Local
!      REAL_T, pointer, dimension(:,:,:,:) :: H
!      REAL_T :: dxInv, dyInv, dzInv
!      integer :: lob(3), hib(3)
!
!      integer :: i, j, k, d, n
!
!!     Compute species enthalpies on box grown by one
!!     init grown box in all dims then update only active dims
!      lob(:) = lo(:)
!      hib(:) = hi(:)
!      do d = 1, dim
!         lob(d) = lo(d)-1
!         hib(d) = hi(d)+1
!      enddo
!
!!     Make space for Hi, use T box, since this better be big enough as well.
!!     Note that any cells on a physical boundary with Dirichlet conditions will 
!!     actually be centered on the edge, so the stencils below must reflect this
!
!      call amrex_allocate(H, t_lo(1), t_hi(1), t_lo(2), t_hi(2), t_lo(3), t_hi(3), 1, nspecies) 
!
!      call pphys_HfromT( lob, hib, H(:,:,:,1), t_lo, t_hi, T(:,:,:), t_lo, t_hi )
!
!
!!   On entry, Fx(1:nspecies) = spec flux, Fx(nspecies+1) = rhoh flux (both untouched)
!!   On exit:
!!   Fx(nspecies+2) = sum[ (species flux).(species enthalpy) ]
!!   Fx(nspecies+3) = extensive heat conduction
!
!!     Compute lambda.Grad(T)
!      dxInv = 1.d0 / dx(1)
!      dyInv = 1.d0 / dx(2)
!#if ( AMREX_SPACEDIM == 3)
!      dzInv = 1.d0 / dx(3)
!#endif
!
!      do k=lo(3),hi(3)
!         do j=lo(2),hi(2)
!            do i=lo(1),hi(1)+1
!               Fx(i,j,k,nspecies+3) = - rhoDx(i,j,k)*(T(i,j,k) - T(i-1,j,k))* dxInv * Ax(i,j,k)
!            enddo
!         enddo
!      enddo
!
!      do k=lo(3),hi(3)
!         do j=lo(2),hi(2)+1
!            do i=lo(1),hi(1)
!               Fy(i,j,k,nspecies+3) = - rhoDy(i,j,k)*(T(i,j,k) - T(i,j-1,k)) * dyInv * Ay(i,j,k)
!            enddo
!         enddo
!      enddo
!
!#if ( AMREX_SPACEDIM == 3)
!      do k=lo(3),hi(3)+1
!         do j=lo(2),hi(2)
!            do i=lo(1),hi(1)
!               Fz(i,j,k,nspecies+3) = - rhoDz(i,j,k)*(T(i,j,k) - T(i,j,k-1)) * dzInv * Az(i,j,k)
!            enddo
!         enddo
!      enddo
!#endif
!
!!     xlo
!      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
!         i = dlo(1)
!         do k=lo(3),hi(3)
!            do j=lo(2),hi(2)
!               Fx(i,j,k,nspecies+3) = 2*Fx(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!!     xhi
!      if (hi(1)+1.eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
!         i = dhi(1)+1
!         do k=lo(3),hi(3)
!            do j=lo(2),hi(2)
!               Fx(i,j,k,nspecies+3) = 2*Fx(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!!     ylo
!      if (lo(2).eq.dlo(2) .and. Tbc(2,1).eq.EXT_DIR) then
!         j=lo(2)
!         do k=lo(3),hi(3)
!            do i=lo(1),hi(1)
!               Fy(i,j,k,nspecies+3) = 2*Fy(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!!     yhi
!      if (hi(2)+1.eq.dhi(2)+1 .and. Tbc(2,2).eq.EXT_DIR) then
!         j=hi(2)+1
!         do k=lo(3),hi(3)
!            do i=lo(1),hi(1)
!               Fy(i,j,k,nspecies+3) = 2*Fy(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!#if ( AMREX_SPACEDIM == 3)
!!     zlo
!      if (lo(3).eq.dlo(3) .and. Tbc(3,1).eq.EXT_DIR) then
!         k=lo(3)
!         do j=lo(2),hi(2)
!            do i=lo(1),hi(1)
!               Fz(i,j,k,nspecies+3) = 2*Fz(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!!     zhi
!      if (hi(3)+1.eq.dhi(3)+1 .and. Tbc(3,2).eq.EXT_DIR) then
!         k=hi(3)+1
!         do j=lo(2),hi(2)
!            do i=lo(1),hi(1)
!               Fz(i,j,k,nspecies+3) = 2*Fz(i,j,k,nspecies+3)
!            enddo
!         enddo
!      endif
!#endif
!
!!     Compute enthalpy flux as hi*Fi
!
!      Fx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),nspecies+2) = 0.d0
!      Fy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),nspecies+2) = 0.d0
!#if ( AMREX_SPACEDIM == 3)
!      Fz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,nspecies+2) = 0.d0
!#endif
!
!      do n = 1, nspecies
!         do k = lo(3), hi(3)
!            do j = lo(2), hi(2)
!               do i = lo(1), hi(1)+1
!                  Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*(H(i,j,k,n)+H(i-1,j,k,n))*0.5d0
!               enddo
!            enddo
!         enddo
!      enddo
!
!      do n = 1, nspecies
!         do k=lo(3),hi(3)
!            do j=lo(2),hi(2)+1
!               do i=lo(1),hi(1)
!                  Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n) *(H(i,j,k,n)+H(i,j-1,k,n))*0.5d0
!               enddo
!            enddo
!         enddo
!      enddo
!
!#if ( AMREX_SPACEDIM == 3)
!      do n = 1, nspecies
!         do k=lo(3),hi(3)+1
!            do j=lo(2),hi(2)
!               do i=lo(1),hi(1)
!                  Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*(H(i,j,k,n)+H(i,j,k-1,n))*0.5d0
!               enddo
!            enddo
!         enddo
!      enddo
!#endif
!
!!     xlo
!      if (lo(1).eq.dlo(1)  .and.  Tbc(1,1).eq.EXT_DIR) then
!         i = dlo(1)
!         Fx(i,lo(2):hi(2),lo(3):hi(3),nspecies+2) = 0.d0
!         do n=1,nspecies
!            do k=lo(3),hi(3)
!               do j=lo(2),hi(2)
!                  Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*H(i-1,j,k,n)
!               enddo
!            enddo
!         enddo
!      endif
!!     xhi
!      if (hi(1)+1.eq.dhi(1)+1  .and.  Tbc(1,2).eq.EXT_DIR) then
!         i = dhi(1)+1
!         Fx(i,lo(2):hi(2),lo(3):hi(3),nspecies+2) = 0.d0
!         do n=1,nspecies
!            do k=lo(3),hi(3)
!               do j=lo(2),hi(2)
!                  Fx(i,j,k,nspecies+2) = Fx(i,j,k,nspecies+2) + Fx(i,j,k,n)*H(i,j,k,n)
!               enddo
!            enddo
!         enddo
!      endif
!!     ylo
!      if (lo(2).eq.dlo(2)  .and.  Tbc(2,1).eq.EXT_DIR) then
!         j = dlo(2)
!         Fy(lo(1):hi(1),j,lo(3):hi(3),nspecies+2) = 0.d0
!         do n=1,nspecies
!            do k=lo(3),hi(3)
!               do i=lo(1),hi(1)
!                  Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n)*H(i,j-1,k,n)
!               enddo
!            enddo
!         enddo
!      endif
!!     yhi
!      if (hi(2)+1.eq.dhi(2)+1  .and.  Tbc(2,2).eq.EXT_DIR) then
!         j = dhi(2)+1
!         Fy(lo(1):hi(1),j,lo(3):hi(3),nspecies+2) = 0.d0
!         do n=1,nspecies
!            do k=lo(3),hi(3)
!               do i=lo(1),hi(1)
!                  Fy(i,j,k,nspecies+2) = Fy(i,j,k,nspecies+2) + Fy(i,j,k,n)*H(i,j,k,n)
!               enddo
!            enddo
!         enddo
!      endif
!#if ( AMREX_SPACEDIM == 3)
!!     zlo
!      if (lo(3).eq.dlo(3)  .and.  Tbc(3,1).eq.EXT_DIR) then
!         k = dlo(3)
!         Fz(lo(1):hi(1),lo(2):hi(2),k,nspecies+2) = 0.d0
!         do n=1,nspecies
!            do j=lo(2),hi(2)
!               do i=lo(1),hi(1)
!                  Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*H(i,j,k-1,n)
!               enddo
!            enddo
!         enddo
!      endif
!!     zhi
!      if (hi(3)+1.eq.dhi(3)+1  .and.  Tbc(3,2).eq.EXT_DIR) then
!         k = dhi(3)+1
!         Fz(lo(1):hi(1),lo(2):hi(2),k,nspecies+2) = 0.d0
!         do n=1,nspecies
!            do j=lo(2),hi(2)
!               do i=lo(1),hi(1)
!                  Fz(i,j,k,nspecies+2) = Fz(i,j,k,nspecies+2) + Fz(i,j,k,n)*H(i,j,k,n)
!               enddo
!            enddo
!         enddo
!      endif
!#endif
!
!      call amrex_deallocate(H)
!
!   end subroutine enth_diff_terms

!=========================================================
!  Fluid viscosity
!=========================================================

   subroutine vel_visc (lo, hi, &
                        T, t_lo, t_hi, &
                        RhoY, rY_lo, rY_hi, &
                        mu, mu_lo, mu_hi) &
                        bind(C, name="vel_visc")

      use transport_module, only: get_visco_coeffs
      use mod_Fvar_def, only : LeEQ1

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: rY_lo(3), rY_hi(3)
      integer, intent(in) :: mu_lo(3), mu_hi(3)
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)), intent(in) :: T
      REAL_T, dimension(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies), intent(in) :: RhoY
      REAL_T, dimension(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3)), intent(out) :: mu

! Local
      REAL_T  :: Yt(lo(1):hi(1),nspecies), RHO_MKS, RHO_MKS_inv
      integer :: i, j, k, n

      integer :: bl(3), bh(3)

      bl = 1
      bh = 1
      bl(1) = lo(1)
      bh(1) = hi(1)

      if (.not. LeEQ1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i=lo(1), hi(1)
                  RHO_MKS = 0.d0
                  do n = 1, nspecies
                     RHO_MKS = RHO_MKS + RhoY(i,j,k,n)
                  end do
                  RHO_MKS_inv = 1.d0 / RHO_MKS
                  do n = 1, nspecies
                     Yt(i,n) = RhoY(i,j,k,n) * RHO_MKS_inv
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
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  mu(i,j,k) = 1.85e-5*(T(i,j,k)/298.0)**.7
               end do
            end do
         end do
      end if

   end subroutine vel_visc

!=========================================================
!  Beta for Wbar diffusion
!=========================================================

   subroutine beta_wbar (lo, hi, &
                         RD, RD_lo, RD_hi, &
                         RDW, RDW_lo, RDW_hi, &
                         Y, Y_lo, Y_hi) &
                         bind(C, name="beta_wbar")

      use fuego_chemistry, only : CKMMWY

      implicit none

      integer, intent(in) ::     lo(3),    hi(3)
      integer, intent(in) ::  RD_lo(3), RD_hi(3)
      integer, intent(in) :: RDW_lo(3),RDW_hi(3)
      integer, intent(in) ::   Y_lo(3),  Y_hi(3)
      REAL_T, dimension(RD_lo(1):RD_hi(1),RD_lo(2):RD_hi(2),RD_lo(3):RD_hi(3),*), intent(in) :: RD
      REAL_T, dimension(RDW_lo(1):RDW_hi(1),RDW_lo(2):RDW_hi(2),RDW_lo(3):RDW_hi(3),*), intent(out) :: RDW
      REAL_T, dimension(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspecies), intent(in) :: Y

      integer :: i, j, k, n
      REAL_T :: Yt(nspecies), RHO, Wavg

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

!=========================================================
!  Estimate dt based on divU
!=========================================================

   subroutine est_divu_dt(flag, dtfactor, delta,&
                          divu, dv_lo,dv_hi,&
                          dsdt, ds_lo,ds_hi,&
                          Rho, r_lo, r_hi,&
                          U, u_lo, u_hi,&
                          volume, v_lo, v_hi,&
#ifdef AMREX_USE_EB
                          volfrac, vf_lo, vf_hi,&
#endif
                          areax, ax_lo, ax_hi,&
                          areay, ay_lo, ay_hi,&
#if ( AMREX_SPACEDIM == 3 )
                          areaz, az_lo, az_hi,&
#endif
                          lo, hi, dt, rhomin) &
                          bind(C, name="est_divu_dt")

      implicit none

      integer :: flag
      integer, intent(in) :: lo(3), hi(3)
      REAL_T , intent(in) :: delta(3)
      integer, intent(in) :: dv_lo(3), dv_hi(3)
      integer, intent(in) :: ds_lo(3), ds_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: u_lo(3), u_hi(3)
      integer, intent(in) :: v_lo(3), v_hi(3)
#ifdef AMREX_USE_EB
      integer, intent(in) :: vf_lo(3), vf_hi(3)
#endif
      integer, intent(in) :: ax_lo(3), ax_hi(3)
      integer, intent(in) :: ay_lo(3), ay_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer, intent(in) :: az_lo(3), az_hi(3)
#endif
      REAL_T, dimension(dv_lo(1):dv_hi(1),dv_lo(2):dv_hi(2),dv_lo(3):dv_hi(3)), intent(in) :: divu
      REAL_T, dimension(ds_lo(1):ds_hi(1),ds_lo(2):ds_hi(2),ds_lo(3):ds_hi(3)), intent(in) :: dsdt
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)), intent(in) :: Rho
      REAL_T, dimension(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3), AMREX_SPACEDIM), intent(in) :: U
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)), intent(in) :: volume
#ifdef AMREX_USE_EB
      REAL_T, dimension(vf_lo(1):vf_hi(1),vf_lo(2):vf_hi(2),vf_lo(3):vf_hi(3)), intent(in) :: volfrac
#endif
      REAL_T, dimension(ax_lo(1):ax_hi(1),ax_lo(2):ax_hi(2),ax_lo(3):ax_hi(3)), intent(in) :: areax
      REAL_T, dimension(ay_lo(1):ay_hi(1),ay_lo(2):ay_hi(2),ay_lo(3):ay_hi(3)), intent(in) :: areay
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(az_lo(1):az_hi(1),az_lo(2):az_hi(2),az_lo(3):az_hi(3)), intent(in) :: areaz
#endif
      REAL_T, intent(in) :: rhomin, dtfactor
      REAL_T, intent(out) :: dt

      integer :: i, j, k
      REAL_T  :: dtcell, dtcell2, denom, rhominij, rhoij
      REAL_T  :: fluxxlo, fluxxhi, fluxylo, fluxyhi, fluxzlo, fluxzhi
      REAL_T  :: a,b,c

      dt = 1.0D20

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dtcell = dt
               if ( flag == 1 ) then
                  if ( divu(i,j,k) > zero ) then
                     if ( rho(i,j,k) > rhomin) then
                        dtcell = dtfactor * ( one - rhomin / Rho(i,j,k) ) / divu(i,j,k)
                     else
                        dtcell = dtfactor * 0.5D0 / divu(i,j,k)
                     endif
                     if ( dsdt(i,j,k) > 1.0D-20 ) then
                        if ( abs(rho(i,j,k)) > rhomin ) then
                           rhominij = rhomin
                        else
                           rhominij = 0.9D0 * abs(rho(i,j,k))
                        endif
                        rhoij = abs(rho(i,j,k))
!
!     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!     We use the latter because it is more robust
!
                        a = rhoij * dsdt(i,j,k) * half
                        b = rhoij * divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two * c / (-b-sqrt(b**2-four*a*c))

                        dtcell2 = dtfactor * dtcell2
                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
                  if ( dtcell<=zero )then
                     write(6,*)'aha'
                  endif
               else if ( flag == 2 ) then
                  denom = u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1)&
                        + u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2)&
#if ( AMREX_SPACEDIM == 3 )
                        + u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)&
#endif
                        + rho(i,j,k) * divu(i,j,k)
                  if ( denom > zero ) then
                     if ( rho(i,j,k) > rhomin ) then
                        dtcell = dtfactor * ( rho(i,j,k) - rhomin ) / denom
                     else
                        dtcell = dtfactor * abs( rho(i,j,k) ) / denom
                     endif
                  endif
               else if (flag == 3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k)) &
                                   *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k)) &
                                   *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k)) &
                                   *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k)) &
                                   *(u(i,j,k,2)+u(i,j+1,k,2))
#if ( AMREX_SPACEDIM == 3 )
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k-1))  &
                                   *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k+1)) &
                                   *(u(i,j,k,3)+u(i,j,k+1,3))
#endif
                  denom = (  (areax(i+1,j,k) * fluxxhi - areax(i,j,k) * fluxxlo) &
                           + (areay(i,j+1,k) * fluxyhi - areay(i,j,k) * fluxylo) &
#if ( AMREX_SPACEDIM == 3 )
                           + (areaz(i,j,k+1) * fluxzhi - areaz(i,j,k) * fluxzlo) &
#endif
#ifdef AMREX_USE_EB
                          ) / (volume(i,j,k) * volfrac(i,j,k))
#else
                          ) / volume(i,j,k)
#endif
                  
                  if ( denom > zero ) then
                     if ( rho(i,j,k) > rhomin ) then
                        dtcell = dtfactor * ( rho(i,j,k) - rhomin ) / denom
                     else
                        dtcell = dtfactor * abs( rho(i,j,k) ) / denom
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

!=========================================================
!  Check that current time step is smaller than divU dt
!=========================================================

   subroutine check_divu_dt(flag, dtfactor, delta,&
                            divu, dv_lo,dv_hi,&
                            dsdt, ds_lo,ds_hi,&
                            Rho, r_lo, r_hi,&
                            U, u_lo, u_hi,&
                            volume, v_lo, v_hi,&
                            areax, ax_lo, ax_hi,&
                            areay, ay_lo, ay_hi,&
#if ( AMREX_SPACEDIM == 3 )
                            areaz, az_lo, az_hi,&
#endif
                            lo, hi, dt, rhomin) &
                            bind(C, name="check_divu_dt")

      implicit none

      integer :: flag
      integer, intent(in) :: lo(3), hi(3)
      REAL_T , intent(in) :: delta(3)
      integer, intent(in) :: dv_lo(3), dv_hi(3)
      integer, intent(in) :: ds_lo(3), ds_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: u_lo(3), u_hi(3)
      integer, intent(in) :: v_lo(3), v_hi(3)
      integer, intent(in) :: ax_lo(3), ax_hi(3)
      integer, intent(in) :: ay_lo(3), ay_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer, intent(in) :: az_lo(3), az_hi(3)
#endif
      REAL_T, dimension(dv_lo(1):dv_hi(1),dv_lo(2):dv_hi(2),dv_lo(3):dv_hi(3)), intent(in) :: divu
      REAL_T, dimension(ds_lo(1):ds_hi(1),ds_lo(2):ds_hi(2),ds_lo(3):ds_hi(3)), intent(in) :: dsdt
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)), intent(in) :: Rho
      REAL_T, dimension(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3), AMREX_SPACEDIM), intent(in) :: U
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)), intent(in) :: volume
      REAL_T, dimension(ax_lo(1):ax_hi(1),ax_lo(2):ax_hi(2),ax_lo(3):ax_hi(3)), intent(in) :: areax
      REAL_T, dimension(ay_lo(1):ay_hi(1),ay_lo(2):ay_hi(2),ay_lo(3):ay_hi(3)), intent(in) :: areay
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(az_lo(1):az_hi(1),az_lo(2):az_hi(2),az_lo(3):az_hi(3)), intent(in) :: areaz
#endif
      REAL_T, intent(in) :: rhomin, dtfactor
      REAL_T, intent(out) :: dt

      integer :: i, j, k
      REAL_T  :: dtcell, dtcell2, denom, rhominij, rhoij
      REAL_T  :: fluxxlo, fluxxhi, fluxylo, fluxyhi, fluxzlo, fluxzhi
      REAL_T  :: a,b,c

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dtcell = bigreal
               if ( flag == 1 ) then
                  if ( divu(i,j,k) > zero ) then
                     if ( rho(i,j,k) > rhomin) then
                        dtcell = ( one - rhomin / Rho(i,j,k) ) / divu(i,j,k)
                     else
                        dtcell = one / divu(i,j,k)
                     endif
                     if ( dsdt(i,j,k) > 1.0D-20 ) then
                        if ( abs(rho(i,j,k)) > rhomin ) then
                           rhominij = rhomin
                        else
                           rhominij = 0.9D0 * abs(rho(i,j,k))
                        endif
                        rhoij = abs(rho(i,j,k))
!
!     ... note: (-b+sqrt(b^2-4ac))/2a = 2c/(-b-sqrt(b^2-4ac))
!     We use the latter because it is more robust
!
                        a = rhoij * dsdt(i,j,k) * half
                        b = rhoij * divu(i,j,k)
                        c = rhominij - rhoij
                        dtcell2 = two * c / (-b-sqrt(b**2-four*a*c))

                        dtcell = min(dtcell,dtcell2)
                     endif
                  endif
               else if ( flag == 2 ) then
                  denom = u(i,j,k,1)*(rho(i+1,j,k)-rho(i-1,j,k))/delta(1)&
                        + u(i,j,k,2)*(rho(i,j+1,k)-rho(i,j-1,k))/delta(2)&
#if ( AMREX_SPACEDIM == 3 )
                        + u(i,j,k,3)*(rho(i,j,k+1)-rho(i,j,k-1))/delta(3)&
#endif
                        + rho(i,j,k) * divu(i,j,k)
                  if ( denom > zero ) then
                     if ( rho(i,j,k) > rhomin ) then
                        dtcell = ( rho(i,j,k) - rhomin ) / denom
                     else
                        dtcell = abs( rho(i,j,k) ) / denom
                     endif
                  endif
               else if (flag == 3) then
                  fluxxlo = fourth*(rho(i,j,k)+rho(i-1,j,k)) &
                                   *(u(i,j,k,1)+u(i-1,j,k,1))
                  fluxxhi = fourth*(rho(i,j,k)+rho(i+1,j,k)) &
                                   *(u(i,j,k,1)+u(i+1,j,k,1))
                  fluxylo = fourth*(rho(i,j,k)+rho(i,j-1,k)) &
                                   *(u(i,j,k,2)+u(i,j-1,k,2))
                  fluxyhi = fourth*(rho(i,j,k)+rho(i,j+1,k)) &
                                   *(u(i,j,k,2)+u(i,j+1,k,2))
#if ( AMREX_SPACEDIM == 3 )
                  fluxzhi = fourth*(rho(i,j,k)+rho(i,j,k-1))  &
                                   *(u(i,j,k,3)+u(i,j,k-1,3))
                  fluxzlo = fourth*(rho(i,j,k)+rho(i,j,k+1)) &
                                   *(u(i,j,k,3)+u(i,j,k+1,3))
#endif
                  denom = (  (areax(i+1,j,k) * fluxxhi - areax(i,j,k) * fluxxlo) &
                           + (areay(i,j+1,k) * fluxyhi - areay(i,j,k) * fluxylo) &
#if ( AMREX_SPACEDIM == 3 )
                           + (areaz(i,j,k+1) * fluxzhi - areaz(i,j,k) * fluxzlo) &
#endif
                          ) / volume(i,j,k)
                  
                  if ( denom > zero ) then
                     if ( rho(i,j,k) > rhomin ) then
                        dtcell = ( rho(i,j,k) - rhomin ) / denom
                     else
                        dtcell = abs( rho(i,j,k) ) / denom
                     endif
                  endif
               endif
               if (dt>dtcell) then
                  write(6,*)'ERROR: FORT_CHECK_DIVU_DT : i,j,k,dt>dtcell = ', &
                      i,j,k,dt,dtcell
               else if (dt>dtcell*dtfactor) then
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

   subroutine dqrad_fill ( dqrad, d_lo, d_hi,&
                           domlo , domhi, delta,&
                           xlo, time, bc ) bind(C, name="dqrad_fill")

      implicit none

      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: dqrad

      call amrex_filccn ( d_lo, d_hi, dqrad, d_lo, d_hi, 1, domlo, domhi, delta, xlo, bc)
      call fillEdges( dqrad, d_lo, d_hi, domlo, domhi, bc)

   end subroutine dqrad_fill

!--------------------------------------------------

   subroutine divu_fill ( divu, d_lo, d_hi,&
                          domlo , domhi, delta,&
                          xlo, time, bc ) bind(C, name="divu_fill")

      implicit none

      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: divu

      call amrex_filccn ( d_lo, d_hi, divu, d_lo, d_hi, 1, domlo, domhi, delta, xlo, bc)
      call fillEdges( divu, d_lo, d_hi, domlo, domhi, bc)

   end subroutine divu_fill

!--------------------------------------------------

   subroutine dsdt_fill ( dsdt, d_lo, d_hi,&
                          domlo , domhi, delta,&
                          xlo, time, bc ) bind(C, name="dsdt_fill")

      implicit none

      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: dsdt

      call amrex_filccn ( d_lo, d_hi, dsdt, d_lo, d_hi, 1, domlo, domhi, delta, xlo, bc)
      call fillEdges( dsdt, d_lo, d_hi, domlo, domhi, bc)

   end subroutine dsdt_fill

!--------------------------------------

   subroutine ydot_fill ( ydot, yd_lo, yd_hi,&
                          domlo, domhi, delta,&
                          xlo, time, bc ) bind(C, name="ydot_fill")

      implicit none

      integer :: yd_lo(3), yd_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(yd_lo(1):yd_hi(1),yd_lo(2):yd_hi(2),yd_lo(3):yd_hi(3)) :: ydot

      call amrex_filccn ( yd_lo, yd_hi, ydot, yd_lo, yd_hi, 1, domlo, domhi, delta, xlo, bc)
      call fillEdges( ydot, yd_lo, yd_hi, domlo, domhi, bc)

   end subroutine ydot_fill

!-------------------------------------------

   subroutine rhoYdot_fill ( rydot, ryd_lo, ryd_hi,&
                             domlo, domhi, delta,&
                             xlo, time, bc ) bind(C, name="rhoYdot_fill")

      implicit none

      integer :: ryd_lo(3), ryd_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(ryd_lo(1):ryd_hi(1),ryd_lo(2):ryd_hi(2),ryd_lo(3):ryd_hi(3)) :: rydot

      call amrex_filccn ( ryd_lo, ryd_hi, rydot, ryd_lo, ryd_hi, 1, domlo, domhi, delta, xlo, bc)
      call fillEdges( rydot, ryd_lo, ryd_hi, domlo, domhi, bc)

   end subroutine rhoYdot_fill

!=========================================================
!  Clip fab min/max
!=========================================================

   subroutine fab_minmax(lo, hi, &
                         fab, f_lo, f_hi,&
                         fmin, fmax, nc) &
                         bind(C, name="fab_minmax")

      implicit none

      integer :: lo(3), hi(3)
      integer :: nc
      integer :: f_lo(3), f_hi(3)
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc) :: fab
      REAL_T  :: fmin, fmax

      integer :: i, j, k, n

      do n = 1,nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  fab(i,j,k,n) = MAX( fmin, MIN( fmax, fab(i,j,k,n) ) )
               end do
            end do
         end do
      end do

   end subroutine fab_minmax

!=========================================================
!  Correct flux to ensure that sum of diffusicve flux is zero
!=========================================================

   subroutine repair_flux_eb (lo, hi, dlo, dhi, &
                           flux, f_lo, f_hi,&
                           RhoY, r_lo, r_hi,&
                           xstate, xstatelo, xstatehi, &
                           afrac_x, axlo, axhi, &
                           ystate, ystatelo, ystatehi, &
                           afrac_y, aylo, ayhi, &
#if ( AMREX_SPACEDIM == 3 )
                           zstate, zstatelo, zstatehi, &
                           afrac_z, azlo, azhi, &
#endif
                           dir, Ybc)&
                           bind(C, name="repair_flux_eb")

      implicit none

      integer :: lo(3), hi(3)
      integer :: dlo(3), dhi(3)
      integer :: dir, Ybc(dim,2)
      integer :: f_lo(3), f_hi(3)
      integer :: r_lo(3), r_hi(3)
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nspecies) :: flux
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspecies) :: RhoY

      integer,  intent(in   ) :: xstatelo(3), xstatehi(3)
      integer,  intent(in   ) :: ystatelo(3), ystatehi(3)
      integer,  intent(in   ) :: axlo(3), axhi(3)
      integer,  intent(in   ) :: aylo(3), ayhi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer,  intent(in   ) :: zstatelo(3), zstatehi(3)
      integer,  intent(in   ) :: azlo(3), azhi(3)
#endif

      REAL_T,  intent(in) ::  xstate(xstatelo(1):xstatehi(1),xstatelo(2):xstatehi(2),xstatelo(3):xstatehi(3),nspecies)
      REAL_T,  intent(in) ::  ystate(ystatelo(1):ystatehi(1),ystatelo(2):ystatehi(2),ystatelo(3):ystatehi(3),nspecies)
      REAL_T,  intent(in) ::  afrac_x(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
      REAL_T,  intent(in) ::  afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
#if ( AMREX_SPACEDIM == 3 ) 
      REAL_T,  intent(in) ::  zstate(zstatelo(1):zstatehi(1),zstatelo(2):zstatehi(2),zstatelo(3):zstatehi(3),nspecies)
      REAL_T,  intent(in) ::  afrac_z(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
#endif
      
      integer :: i, j, k, n
      REAL_T :: sumFlux, RhoYe(nspecies), sumRhoYe

!write(*,*) ' DEBUG LO HI',lo,hi      
!write(*,*) 'DEBUG AFRAC_X', lbound(afrac_x), ubound(afrac_x)
!write(*,*) 'DEBUG AFRAC_Y', lbound(afrac_y), ubound(afrac_y)
      if (dir.eq.0) then

!     First, assume away from physical boundaries, then use boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
        
                 if ( afrac_x(i,j,k) > zero ) then
               
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = xstate(i,j,k,n)
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  end do
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,nspecies
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                  end do
                  else
                  do n=1,nspecies
                    flux(i,j,k,n) = 0.0d0
                  end do
                  end if
                  
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
               
                 if ( afrac_y(i,j,k) > zero ) then
                   sumFlux = 0.d0
                   sumRhoYe = 0.d0
                   do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = ystate(i,j,k,n) 
                     sumRhoYe = sumRhoYe + RhoYe(n)
                   enddo
                   sumRhoYe = 1.0D0/sumRhoYe
                   do n=1,nspecies
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                   end do         
                 else
                   do n=1,nspecies
                     flux(i,j,k,n) = 0.0d0
                   end do
                 end if
                  
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

#if ( AMREX_SPACEDIM == 3 )
      else if (dir.eq.2) then

!     First, assume away from physical boundaries, then replace with boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
   
                 if ( afrac_z(i,j,k) > zero ) then
                   sumFlux = 0.d0
                   sumRhoYe = 0.d0
                   do n=1,nspecies
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = zstate(i,j,k,n) 
                     sumRhoYe = sumRhoYe + RhoYe(n)
                   enddo
                   sumRhoYe = 1.0D0/sumRhoYe
                   do n=1,nspecies
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                   end do
                  
                 else
                   do n=1,nspecies
                     flux(i,j,k,n) = 0.0d0
                   end do
                 end if
              
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
#endif

      endif

   end subroutine repair_flux_eb

!=========================================================
!  Correct flux to ensure that sum of diffusicve flux is zero
!=========================================================

   subroutine repair_flux (lo, hi, dlo, dhi, &
                           flux, f_lo, f_hi,&
                           RhoY, r_lo, r_hi, dir, Ybc)&
                           bind(C, name="repair_flux")

      implicit none

      integer :: lo(3), hi(3)
      integer :: dlo(3), dhi(3)
      integer :: dir, Ybc(dim,2)
      integer :: f_lo(3), f_hi(3)
      integer :: r_lo(3), r_hi(3)
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nspecies) :: flux
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspecies) :: RhoY
      
      integer :: i, j, k, n
      REAL_T :: sumFlux, RhoYe(nspecies), sumRhoYe

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

#if ( AMREX_SPACEDIM == 3 )
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
#endif

      endif

   end subroutine repair_flux


!=========================================================
! Increment old state with flux divergence
!=========================================================

   subroutine incrwext_flx_div(lo, hi, &
                               xflux,  xf_lo, xf_hi,&
                               yflux,  yf_lo, yf_hi,&
#if ( AMREX_SPACEDIM == 3 )
                               zflux,  zf_lo, zf_hi,&
#endif
                               stateo, so_lo, so_hi,&
                               staten, sn_lo, sn_hi,&
                               vol,    v_lo, v_hi,&
                               nc, dt) &
                               bind(C, name="incrwext_flx_div")

      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: nc
      integer, intent(in) :: xf_lo(3), xf_hi(3)
      integer, intent(in) :: yf_lo(3), yf_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer, intent(in) :: zf_lo(3), zf_hi(3)
#endif
      integer, intent(in) :: so_lo(3), so_hi(3)
      integer, intent(in) :: sn_lo(3), sn_hi(3)
      integer, intent(in) :: v_lo(3), v_hi(3)
      REAL_T, dimension(xf_lo(1):xf_hi(1),xf_lo(2):xf_hi(2),xf_lo(3):xf_hi(3),nc) ::  xflux
      REAL_T, dimension(yf_lo(1):yf_hi(1),yf_lo(2):yf_hi(2),yf_lo(3):yf_hi(3),nc) ::  yflux
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(zf_lo(1):zf_hi(1),zf_lo(2):zf_hi(2),zf_lo(3):zf_hi(3),nc) ::  zflux
#endif
      REAL_T, dimension(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3)) ::  stateo
      REAL_T, dimension(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3)) ::  staten
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) ::  vol
      REAL_T :: dt

      integer :: i, j, k, n
      REAL_T  :: dF

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dF = zero
               do n = 1, nc
                  dF = dF + ( ( xflux(i+1,j,k,n) - xflux(i,j,k,n) ) &
                          +   ( yflux(i,j+1,k,n) - yflux(i,j,k,n) ) &
#if ( AMREX_SPACEDIM == 3 )
                          +   ( zflux(i,j,k+1,n) - zflux(i,j,k,n) ) &
#endif
                            )
               end do
               staten(i,j,k) = stateo(i,j,k) + dF * dt / vol(i,j,k)
            end do
         end do
      end do

   end subroutine incrwext_flx_div

!=========================================================
! Compute flux divergence
!=========================================================

   subroutine flux_div (lo, hi, &
                        update, u_lo, u_hi,&
                        mask, m_lo, m_hi, &
                        xflux,  xf_lo, xf_hi,&
                        yflux,  yf_lo, yf_hi,&
#if ( AMREX_SPACEDIM == 3 )
                        zflux,  zf_lo, zf_hi,&
#endif
                        vol,    v_lo, v_hi, &
#ifdef AMREX_USE_EB
                        volfrac, vf_lo, vf_hi, &
#endif
                        nc, scal) &
                        bind(C, name="flux_div")

      use amrex_mempool_module

      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: nc
      integer, intent(in) :: u_lo(3), u_hi(3)
      integer, intent(in) :: m_lo(3), m_hi(3)
      integer, intent(in) :: xf_lo(3), xf_hi(3)
      integer, intent(in) :: yf_lo(3), yf_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer, intent(in) :: zf_lo(3), zf_hi(3)
#endif
      integer, intent(in) :: v_lo(3), v_hi(3)
      REAL_T, dimension(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc) :: update
      integer, dimension(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3)) :: mask
      REAL_T, dimension(xf_lo(1):xf_hi(1),xf_lo(2):xf_hi(2),xf_lo(3):xf_hi(3),nc) :: xflux
      REAL_T, dimension(yf_lo(1):yf_hi(1),yf_lo(2):yf_hi(2),yf_lo(3):yf_hi(3),nc) :: yflux
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(zf_lo(1):zf_hi(1),zf_lo(2):zf_hi(2),zf_lo(3):zf_hi(3),nc) :: zflux
#endif
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vol
#ifdef AMREX_USE_EB
      integer, intent(in) :: vf_lo(3), vf_hi(3)
      REAL_T, dimension(vf_lo(1):vf_hi(1),vf_lo(2):vf_hi(2),vf_lo(3):vf_hi(3)), intent(in) :: volfrac
#endif

      REAL_T :: scal

      REAL_T, pointer, dimension(:,:,:) :: ivol
      
      REAL_T :: infinity

      integer :: i, j, k, n

      infinity = HUGE(0.0)
      
      call amrex_allocate(ivol,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( mask(i,j,k) /= -1 ) then
#ifdef AMREX_USE_EB
                  ivol(i,j,k) = scal / ((vol(i,j,k) * volfrac(i,j,k)))
#else
                  ivol(i,j,k) = scal / vol(i,j,k)
#endif
               end if
            end do
         end do
      end do

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( mask(i,j,k) /= -1 ) then
                     update(i,j,k,n) = &
                        ( (xflux(i+1,j,k,n)-xflux(i,j,k,n)) &
                        + (yflux(i,j+1,k,n)-yflux(i,j,k,n)) & 
#if ( AMREX_SPACEDIM == 3 )
                        + (zflux(i,j,k+1,n)-zflux(i,j,k,n)) &
#endif
                        ) * ivol(i,j,k)
                  end if
               end do
            end do
         end do
      end do

      call amrex_deallocate(ivol)

   end subroutine flux_div

!=========================================================
! This routine fills ghost cells with the value from the nearest
! interior cell.
!=========================================================

   subroutine fillEdges( dat, d_lo, d_hi,&
                         domlo, domhi, bc )&
                         bind(C, name="fillEdges")

      implicit none

      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: dat

      integer :: ilo, ihi, jlo, jhi, klo, khi
      logical :: do_xlo,do_xhi,do_ylo,do_yhi,do_zlo,do_zhi
      integer :: i, j, k

      ! First pass to put all at input d_lo/d_hi
      ilo = d_lo(1)
      ihi = d_hi(1)
      jlo = d_lo(2)
      jhi = d_hi(2)
      klo = d_lo(3)
      khi = d_hi(3)

      ! Now only take min/max for the dimensions
      ilo = max(d_lo(1),domlo(1))
      ihi = min(d_hi(1),domhi(1))
      jlo = max(d_lo(2),domlo(2))
      jhi = min(d_hi(2),domhi(2))
#if ( AMREX_SPACEDIM == 3 )
      klo = max(d_lo(3),domlo(3))
      khi = min(d_hi(3),domhi(3))
#endif

      do_xlo = .false.
      do_xhi = .false.
      do_ylo = .false.
      do_yhi = .false.
      do_zlo = .false.
      do_zhi = .false.

      do_xlo = d_lo(1).lt.domlo(1)
      do_xhi = d_hi(1).gt.domhi(1)
      do_ylo = d_lo(2).lt.domlo(2)
      do_yhi = d_hi(2).gt.domhi(2)
#if ( AMREX_SPACEDIM == 3 )
      do_zlo = d_lo(3).lt.domlo(3)
      do_zhi = d_hi(3).gt.domhi(3)
#endif

      if (bc(1,1).eq.EXT_DIR.and. do_xlo) then
         do k = klo,khi
            do j = jlo, jhi
               do i = d_lo(1), domlo(1)-1
                  dat(i,j,k) = dat(domlo(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = d_lo(2), domlo(2)-1
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3),domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, d_hi(2)
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
#if ( AMREX_SPACEDIM == 3 )
         if (do_zlo) then
            do k = d_lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, d_hi(3)
               do j = jlo, jhi
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
#endif
      endif ! if on xho


      if (bc(1,2).eq.EXT_DIR.and. do_xhi) then
         do k = klo, khi
            do j = jlo, jhi
               do i = domhi(1)+1, d_hi(1)
                  dat(i,j,k) = dat(domhi(1),j,k)
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = klo, khi
               do j = d_lo(2), domlo(2)-1
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3),domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
         if (do_yhi) then
            do k = klo, khi
               do j = domhi(2)+1, d_hi(2)
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
#if ( AMREX_SPACEDIM == 3 )
         if (do_zlo) then
            do k = d_lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, d_hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
#endif
      endif  ! if on xhi

      if (bc(2,1).eq.EXT_DIR.and.do_ylo) then
         do k = klo, khi
            do j = d_lo(2), domlo(2)-1
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,domlo(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = d_lo(2), domlo(2)-1
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domlo(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = d_lo(2), domlo(2)-1
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),domlo(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
#if ( AMREX_SPACEDIM == 3 )
         if (do_zlo) then
            do k = d_lo(3), domlo(3)-1
               do j = d_lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, d_hi(3)
               do j = d_lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
#endif
      endif ! if on ylo

      if (bc(2,2).eq.EXT_DIR.and.do_yhi) then
         do k = klo, khi
            do j = domhi(2)+1, d_hi(2)
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,domhi(2),k)
               enddo
            enddo
         enddo
         if (do_xlo) then
            do k = klo, khi
               do j = domhi(2)+1, d_hi(2)
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),domhi(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
         if(do_xhi)then
            do k = klo, khi
               do j = domhi(2)+1, d_hi(2)
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),domhi(2),k)
                  enddo
               enddo
            enddo
#if ( AMREX_SPACEDIM == 3 )
            if (do_zlo) then
               do k = d_lo(3),domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_zhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
#endif
         endif
#if ( AMREX_SPACEDIM == 3 )
         if (do_zlo) then
            do k = d_lo(3), domlo(3)-1
               do j = domhi(2)+1, d_hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_zhi)then
            do k = domhi(3)+1, d_hi(3)
               do j = domhi(2)+1, d_hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
         endif
#endif
      endif ! if on yhi

#if ( AMREX_SPACEDIM == 3 )
      if (bc(3,1).eq.EXT_DIR.and. do_zlo) then
         do k = d_lo(3), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,j,domlo(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = d_lo(3), domlo(3)-1
               do j = d_lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = d_lo(3), domlo(3)-1
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = d_lo(3), domlo(3)-1
               do j = domhi(2)+1, d_hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domlo(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = d_lo(3), domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1),domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = d_lo(3), domlo(3)-1
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domlo(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = d_lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = d_lo(3), domlo(3)-1
               do j = jlo, jhi
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domlo(3))
                  enddo
               enddo
            enddo
         endif
      endif            

      if (bc(3,2).eq.EXT_DIR.and. do_zhi) then
         do k = domhi(3)+1, d_hi(3)
            do j = jlo, jhi
               do i = ilo, ihi
                  dat(i,j,k) = dat(i,j,domhi(3))
               enddo
            enddo
         enddo
         if (do_ylo) then
            do k = domhi(3)+1, d_hi(3)
               do j = d_lo(2), domlo(2)-1
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domlo(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = d_lo(2), domlo(2)-1
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domlo(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_yhi) then
            do k = domhi(3)+1, d_hi(3)
               do j = domhi(2)+1, d_hi(2)
                  do i = ilo, ihi
                     dat(i,j,k) = dat(i,domhi(2),domhi(3))
                  enddo
               enddo
            enddo
            if (do_xlo) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = d_lo(1), domlo(1)-1
                        dat(i,j,k) = dat(domlo(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
            if (do_xhi) then
               do k = domhi(3)+1, d_hi(3)
                  do j = domhi(2)+1, d_hi(2)
                     do i = domhi(1)+1, d_hi(1)
                        dat(i,j,k) = dat(domhi(1),domhi(2),domhi(3))
                     enddo
                  enddo
               enddo
            endif
         endif
         if (do_xlo) then
            do k = domhi(3)+1, d_hi(3)
               do j = jlo, jhi
                  do i = d_lo(1), domlo(1)-1
                     dat(i,j,k) = dat(domlo(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
         if(do_xhi)then
            do k = domhi(3)+1, d_hi(3)
               do j = jlo, jhi
                  do i = domhi(1)+1, d_hi(1)
                     dat(i,j,k) = dat(domhi(1),j,domhi(3))
                  enddo
               enddo
            enddo
         endif
      endif
#endif

   end subroutine fillEdges

!=========================================================
! This routine fills ghost cells with zeros
!=========================================================

   subroutine fillWithZeros ( dat, d_lo, d_hi,&
                              domlo, domhi, bc )&
                              bind(C, name="fillWithZeros")

      implicit none

      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: dat

      integer :: i, j, k

      if (bc(1,1).eq.EXT_DIR.and.d_lo(1).lt.domlo(1)) then
         do k = d_lo(3), d_hi(3)
            do j = d_lo(2), d_hi(2)
               do i = d_lo(1), domlo(1)-1
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif

      if (bc(1,2).eq.EXT_DIR.and.d_hi(1).gt.domhi(1)) then
         do k = d_lo(3), d_hi(3)
            do j = d_lo(2), d_hi(2)
               do i = domhi(1)+1, d_hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif

      if (bc(2,1).eq.EXT_DIR.and.d_lo(2).lt.domlo(2)) then
         do k = d_lo(3), d_hi(3)
            do j = d_lo(2), domlo(2)-1
               do i = d_lo(1), d_hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif

      if (bc(2,2).eq.EXT_DIR.and.d_hi(2).gt.domhi(2)) then
         do k = d_lo(3), d_hi(3)
            do j = domhi(2)+1, d_hi(2)
               do i = d_lo(1), d_hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )
      if (bc(3,1).eq.EXT_DIR.and.d_lo(3).lt.domlo(3)) then
         do k = d_lo(3), domlo(3)-1
            do j = d_lo(2), d_hi(2)
               do i = d_lo(1), d_hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif

      if (bc(3,2).eq.EXT_DIR.and.d_hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, d_hi(3)
            do j = d_lo(2), d_hi(2)
               do i = d_lo(1), d_hi(1)
                  dat(i,j,k) = zero
               enddo
            enddo
         enddo
      endif
#endif

   end subroutine fillWithZeros

!=========================================================
! Compute u times gradP 
!=========================================================

   subroutine compute_ugradp(p, p_lo, p_hi,&
                             ugradp, ugp_lo, ugp_hi,&
                             umac, u_lo, u_hi,&
                             vmac, v_lo, v_hi,&
#if ( AMREX_SPACEDIM == 3 )
                             wmac, w_lo, w_hi,&
#endif
                             lo, hi, dx)&
                             bind(C, name="compute_ugradp")

      implicit none

      integer :: lo(3), hi(3)
      integer :: p_lo(3), p_hi(3)
      integer :: ugp_lo(3), ugp_hi(3)
      integer :: u_lo(3), u_hi(3)
      integer :: v_lo(3), v_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer :: w_lo(3), w_hi(3)
#endif
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p
      REAL_T, dimension(ugp_lo(1):ugp_hi(1),ugp_lo(2):ugp_hi(2),ugp_lo(3):ugp_hi(3)) :: ugradp
      REAL_T, dimension(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3)) :: umac
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vmac
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3)) :: wmac
#endif
      REAL_T :: dx(3)

      integer :: i, j, k
      REAL_T :: uadv, vadv, wadv
      REAL_T :: p_x_lo, p_x_hi
      REAL_T :: p_y_lo, p_y_hi
      REAL_T :: p_z_lo, p_z_hi

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uadv = half*(umac(i,j,k) + umac(i+1,j,k))
               vadv = half*(vmac(i,j,k) + vmac(i,j+1,k))
#if ( AMREX_SPACEDIM == 3 )
               wadv = half*(wmac(i,j,k) + wmac(i,j,k+1))
#endif
               p_x_hi = merge(p(i  ,j,k),p(i+1,j,k),umac(i+1,j,k)>=zero)
               p_x_lo = merge(p(i-1,j,k),p(i  ,j,k),umac(i  ,j,k)>=zero)
               p_y_hi = merge(p(i,j  ,k),p(i,j+1,k),vmac(i,j+1,k)>=zero)
               p_y_lo = merge(p(i,j-1,k),p(i,j  ,k),vmac(i,j  ,k)>=zero)
#if ( AMREX_SPACEDIM == 3 )
               p_z_hi = merge(p(i,j,k  ),p(i,j,k+1),wmac(i,j,k+1)>=zero)
               p_z_lo = merge(p(i,j,k-1),p(i,j,k  ),wmac(i,j,k  )>=zero)
#endif
               ugradp(i,j,k) = (   uadv * (p_x_hi - p_x_lo) / dx(1) &
                                 + vadv * (p_y_hi - p_y_lo) / dx(2) &
#if ( AMREX_SPACEDIM == 3 )
                                 + wadv * (p_z_hi - p_z_lo) / dx(3) &
#endif
                               )
            end do
         end do
      end do

   end subroutine compute_ugradp

!=========================================================
! Floor T in a conservative manner
!=========================================================

   integer function conservative_T_floor ( lo, hi,&
                                           state, s_lo, s_hi,&
                                           min_T, Tcomp, Rcomp, first_spec, last_spec,&
                                           RhoH, ratio, tmp, nt ) &
                                           bind(C, name="conservative_T_floor")

      implicit none

      integer :: lo(3), hi(3)
      integer :: s_lo(3), s_hi(3)
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:*) :: state
      integer :: Tcomp, Rcomp, first_spec, last_spec, RhoH, ratio(3), nt
      REAL_T  :: min_T

      REAL_T  :: tmp(0:nt-1)
      integer :: loC(3), hiC(3)
      integer :: ii, jj, kk, iii, jjj, kkk, ncells
      REAL_T  :: ncellsInv
      logical :: bad_T

      integer :: n, i, j, k

!     Returns the number of fine cells fixed up
      conservative_T_floor = 0

      ncells = 1
      do n = 1, dim
         loC(n) = lo(n)/ratio(n)
         hiC(n) = (hi(n)+1)/ratio(n) - 1
         ncells = ncells*ratio(n)
      enddo
      ncellsInv = 1.d0 / ncells

      do k = loC(3), hiC(3)
         do j = loC(2), hiC(2)
            do i = loC(1), hiC(1)

               bad_T = .false.
               do kk = 0, ratio(3)-1
                  kkk = ratio(3)*k + kk
                  do jj = 0, ratio(2)-1
                     jjj = ratio(2)*j + jj
                     do ii = 0, ratio(1)-1
                        iii = ratio(1)*i + ii
                        if (state(iii,jjj,kkk,Tcomp).lt.min_T) then
                           bad_T = .true.
                        endif
                     enddo
                  enddo
               enddo

               if (bad_T .eqv. .true.) then

                  tmp(Rcomp) = 0.d0
                  do n = first_spec, last_spec
                     tmp(n) = 0.d0
                  enddo
                  tmp(RhoH) = 0.d0


                  do kk = 0,ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj = 0, ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii = 0, ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           tmp(Rcomp) = tmp(Rcomp) + state(iii,jjj,kkk,Rcomp)
                           do n = first_spec, last_spec
                              tmp(n) = tmp(n) + state(iii,jjj,kkk,n)
                           enddo
                           tmp(RhoH) = tmp(RhoH) + state(iii,jjj,kkk,RhoH)
                           
                        enddo
                     enddo
                  enddo

                  conservative_T_floor = conservative_T_floor + ncells
                  tmp(Rcomp) = tmp(Rcomp) * ncellsInv
                  do n = first_spec, last_spec
                     tmp(n) = tmp(n) * ncellsInv
                  enddo
                  tmp(RhoH) = tmp(RhoH)* ncellsInv
                  
                  do kk = 0, ratio(3)-1
                     kkk = ratio(3)*k + kk
                     do jj = 0, ratio(2)-1
                        jjj = ratio(2)*j + jj
                        do ii = 0, ratio(1)-1
                           iii = ratio(1)*i + ii
                           
                           state(iii,jjj,kkk,Rcomp) = tmp(Rcomp)
                           do n = first_spec, last_spec
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

   subroutine part_cnt_err(tag, t_lo, t_hi,&
                           set, clear,&
                           var, v_lo, v_hi,&
                           lo, hi, nd, domlo, domhi,&
                           delta, xlo, problo, time, level) &
                           bind(C, name="part_cnt_err")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nd) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( var(i,j,k,1) .gt. zero) tag(i,j,k) = set
               !if (var(i,j,k,1) .gt. zero) print *,'TAGGING ',i,j,k
            end do
         end do
      end do

   end subroutine part_cnt_err

!=========================================================
! Compute curvature from temperature
!=========================================================

   subroutine mcurve ( lo, hi, &
                       T, t_lo, t_hi,&
                       curv, c_lo, c_hi,&
                       wrk, w_lo, w_hi, delta )&
                       bind(C, name="mcurve")

      implicit none

      integer lo(3), hi(3)
      integer :: t_lo(3), t_hi(3)
      integer :: c_lo(3), c_hi(3)
      integer :: w_lo(3), w_hi(3)
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3)) :: curv
      REAL_T, dimension(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),dim) :: wrk
      REAL_T :: delta(3)

      REAL_T  :: mag
      REAL_T  :: Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Txz,tdxI(3),dxSqI(3)
      integer :: i, j, k

      ! Fill normal on nodes (assumes 1 grow cell properly filled)
      do i = 1, dim
         tdxI(i) =  one / (two*delta(i))
         dxSqI(i) = one / (delta(i)*delta(i))
      enddo
 
      ! Init those to zero. If AMREX_SPACEDIM /= 3, all z component remain zero.
      Tx  = 0.0d0
      Ty  = 0.0d0
      Tz  = 0.0d0
      Txx = 0.0d0
      Tyy = 0.0d0
      Tzz = 0.0d0
      Txy = 0.0d0
      Txz = 0.0d0
      Tyz = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               Tx = tdxI(1)*(T(i+1,j,k) - T(i-1,j,k))
               Ty = tdxI(2)*(T(i,j+1,k) - T(i,j-1,k))
#if ( AMREX_SPACEDIM == 3 )
               Tz = tdxI(3)*(T(i,j,k+1) - T(i,j,k-1))
#endif

               Txx = dxSqI(1)*(T(i+1,j,k) - two*T(i,j,k) + T(i-1,j,k))
               Tyy = dxSqI(2)*(T(i,j+1,k) - two*T(i,j,k) + T(i,j-1,k))
#if ( AMREX_SPACEDIM == 3 )
               Tzz = dxSqI(3)*(T(i,j,k+1) - two*T(i,j,k) + T(i,j,k-1))
#endif

               Txy = tdxI(1)*tdxI(2)*(T(i+1,j+1,k) - T(i-1,j+1,k) - T(i+1,j-1,k) + T(i-1,j-1,k))
#if ( AMREX_SPACEDIM == 3 )
               Txz = tdxI(1)*tdxI(3)*(T(i+1,j,k+1) - T(i-1,j,k+1) - T(i+1,j,k-1) + T(i-1,j,k-1))
               Tyz = tdxI(2)*tdxI(3)*(T(i,j+1,k+1) - T(i,j-1,k+1) - T(i,j+1,k-1) + T(i,j-1,k-1))
#endif

               mag = max(1.0d-12, SQRT(Tx*Tx + Ty*Ty + Tz*Tz))

               curv(i,j,k) = - half * (  Txx + Tyy + Tzz &
                                       - (  Tx*(Tx*Txx + Ty*Txy + Tz*Txz) &
                                          + Ty*(Tx*Txy + Ty*Tyy + Tz*Tyz) &
                                          + Tz*(Tx*Txz + Ty*Tyz + Tz*Tzz) ) / mag**2.0d0 ) / mag
            end do
         end do
      end do

   end subroutine mcurve

!=========================================================
! Smooth variable by averaging surrounding cells
!=========================================================

   subroutine smooth ( lo, hi,&
                       Tin, ti_lo, ti_hi,&
                       Tout, to_lo, to_hi )&
                       bind(C, name="smooth")

      implicit none

      integer :: lo(3), hi(3)
      integer :: ti_lo(3), ti_hi(3)
      integer :: to_lo(3), to_hi(3)
      REAL_T, dimension(ti_lo(1):ti_hi(1),ti_lo(2):ti_hi(2),ti_lo(3):ti_hi(3)) :: Tin
      REAL_T, dimension(to_lo(1):to_hi(1),to_lo(2):to_hi(2),to_lo(3):to_hi(3)) :: Tout
      REAL_T :: factor
      integer :: i, j, k, ii, jj, kk

#if ( AMREX_SPACEDIM == 2 )
      factor = 1.0d0 / 16.0d0
#elif ( AMREX_SPACEDIM == 3 )
      factor = 1.0d0 / 64.0d0
#endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               Tout(i,j,k) = zero
#if ( AMREX_SPACEDIM == 2 )
                do jj = 0, 1
                   do ii = 0, 1
                      Tout(i,j,k) = Tout(i,j,k) &
                          + Tin(i+ii,j+jj,  k)   + Tin(i+ii-1,j+jj,  k  ) &
                          + Tin(i+ii,j+jj-1,k  ) + Tin(i+ii-1,j+jj-1,k  )
                   end do
                end do
#elif ( AMREX_SPACEDIM == 3 )
               do kk = 0, 1
                  do jj = 0, 1
                     do ii = 0, 1
                        Tout(i,j,k) = Tout(i,j,k) &
                            + Tin(i+ii,j+jj,  k+kk-1) + Tin(i+ii-1,j+jj,  k+kk-1) &
                            + Tin(i+ii,j+jj-1,k+kk-1) + Tin(i+ii-1,j+jj-1,k+kk-1) &
                            + Tin(i+ii,j+jj,  k+kk)   + Tin(i+ii-1,j+jj,  k+kk  ) &
                            + Tin(i+ii,j+jj-1,k+kk  ) + Tin(i+ii-1,j+jj-1,k+kk  )
                     end do
                  end do
               end do
#endif
               Tout(i,j,k) = Tout(i,j,k) * factor
            end do
         end do
      end do

   end subroutine smooth

!=========================================================
!  Compute gradient of Wbar
!=========================================================

   subroutine grad_wbar ( lo, hi, &
                          Wbar, w_lo, w_hi,&
                          rDe, r_lo, r_hi,&
                          flux, f_lo, f_hi,&
                          area, a_lo, a_hi,&
                          dx, dir, mult, inc) &
                          bind(C, name="grad_wbar")

      implicit none

      integer :: lo(3), hi(3)
      integer :: w_lo(3), w_hi(3)
      integer :: r_lo(3), r_hi(3)
      integer :: f_lo(3), f_hi(3)
      integer :: a_lo(3), a_hi(3)
      integer :: dir
      REAL_T, dimension(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3)) :: Wbar
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: rDe
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3)) :: flux
      REAL_T, dimension(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3)) :: area
      REAL_T  :: dx, mult, inc

      REAL_T  :: Wgr, fac
      integer :: i, j, k

      fac = mult / dx

      if (inc == 0) then

!     compute grad wbar fluxes

         if (dir==0) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir==1) then

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

#if ( AMREX_SPACEDIM == 3 )
         else if (dir==2) then

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
#endif

         else
            call amrex_abort('Bad dir in grad_wbar')
         endif

      else

!     increment grad wbar fluxes by a factor of inc (can be negative)

         if (dir==0) then

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i-1,j,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

         else if (dir==1) then

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j-1,k))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo

#if ( AMREX_SPACEDIM == 3 )
         else if (dir==2) then

            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     Wgr =   fac*(Wbar(i,j,k) - Wbar(i,j,k-1))
                     flux(i,j,k) = flux(i,j,k) + inc * rDe(i,j,k) * Wgr * area(i,j,k)
                  enddo
               enddo
            enddo
#endif

         else
            call amrex_abort('Bad dir in grad_wbar')
         endif

      end if

   end subroutine grad_wbar

!=========================================================
!  Recompute cell update from fluxes
!=========================================================

   subroutine recomp_update(lo, hi, &
                            update, u_lo, u_hi,&
                            xflux,  xf_lo, xf_hi,&
                            yflux,  yf_lo, yf_hi,&
#if ( AMREX_SPACEDIM == 3 )
                            zflux,  zf_lo, zf_hi,&
#endif
                            vol,    v_lo, v_hi, &
                            nc) &
                            bind(C, name="recomp_update")

      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: nc
      integer, intent(in) :: u_lo(3), u_hi(3)
      integer, intent(in) :: xf_lo(3), xf_hi(3)
      integer, intent(in) :: yf_lo(3), yf_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer, intent(in) :: zf_lo(3), zf_hi(3)
#endif
      integer, intent(in) :: v_lo(3), v_hi(3)
      REAL_T, dimension(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc) :: update
      REAL_T, dimension(xf_lo(1):xf_hi(1),xf_lo(2):xf_hi(2),xf_lo(3):xf_hi(3),nc) :: xflux
      REAL_T, dimension(yf_lo(1):yf_hi(1),yf_lo(2):yf_hi(2),yf_lo(3):yf_hi(3),nc) :: yflux
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(zf_lo(1):zf_hi(1),zf_lo(2):zf_hi(2),zf_lo(3):zf_hi(3),nc) :: zflux
#endif
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vol

      integer :: i, j, k, n

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  update(i,j,k,n) = - ( (xflux(i+1,j,k,n)-xflux(i,j,k,n)) &
                                      + (yflux(i,j+1,k,n)-yflux(i,j,k,n)) &
#if ( AMREX_SPACEDIM == 3 )
                                      + (zflux(i,j,k+1,n)-zflux(i,j,k,n)) &
#endif
                                      ) / vol(i,j,k)
               end do
            end do
         end do
      end do

   end subroutine recomp_update

!=========================================================
!  Error tagging function : greater than 
!=========================================================

   subroutine valgt_error ( tag, t_lo, t_hi,&
                            set, clear,&
                            var, v_lo, v_hi,&
                            lo, hi, nd, domlo, domhi,&
                            delta, xlo, problo, time, level,value) &
                            bind(C, name="valgt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), var(i,j,k,1)>value)
            end do
         end do
      end do

   end subroutine valgt_error

!=========================================================
!  Error tagging function : lower than 
!=========================================================

   subroutine vallt_error ( tag, t_lo, t_hi,&
                            set, clear,&
                            var, v_lo, v_hi,&
                            lo, hi, nd, domlo, domhi,&
                            delta, xlo, problo, time, level,value) &
                            bind(C, name="vallt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), var(i,j,k,1)<value)
            end do
         end do
      end do

   end subroutine vallt_error

!=========================================================
!  Error tagging function : vorticity mag 
!  Value depends on level
!=========================================================

   subroutine magvort_error ( tag, t_lo, t_hi,&
                              set, clear,&
                              var, v_lo, v_hi,&
                              lo, hi, nd, domlo, domhi,&
                              delta, xlo, problo, time, level,value) &
                              bind(C, name="magvort_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), ABS(var(i,j,k,1)) >= value*2.d0**level)
            end do
         end do
      end do

   end subroutine magvort_error

!=========================================================
!  Error tagging function : diff adj. cells greater than 
!=========================================================

   subroutine diffgt_error ( tag, t_lo, t_hi,&
                             set, clear,&
                             var, v_lo, v_hi,&
                             lo, hi, nd, domlo, domhi,&
                             delta, xlo, problo, time, level,value) &
                             bind(C, name="diffgt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      REAL_T  :: axp, axm, ayp, aym, azp, azm, aerr
      integer :: i, j, k

      axp = 0.0d0
      axm = 0.0d0
      ayp = 0.0d0
      aym = 0.0d0
      azp = 0.0d0
      azm = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               axp = ABS(var(i+1,j,k,1) - var(i,j,k,1))
               axm = ABS(var(i-1,j,k,1) - var(i,j,k,1))
               ayp = ABS(var(i,j+1,k,1) - var(i,j,k,1))
               aym = ABS(var(i,j-1,k,1) - var(i,j,k,1))
#if ( AMREX_SPACEDIM == 3 )
               azp = ABS(var(i,j,k+1,1) - var(i,j,k,1))
               azm = ABS(var(i,j,k-1,1) - var(i,j,k,1))
#endif
               aerr = MAX(azp,MAX(azm,MAX(axp,MAX(axm,MAX(ayp,aym)))))

               if ( aerr>=value ) then
                  tag(i,j,k) = set
               endif
            end do
         end do
      end do

   end subroutine diffgt_error

!=========================================================
!  Error tagging function : tag box 
!=========================================================

   subroutine box_error ( tag, t_lo, t_hi,&
                          set, clear,&
                          boxlo, boxhi, lo, hi,&
                          domlo, domhi,&
                          delta, xlo, problo, time, level) &
                          bind(C, name="box_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      REAL_T :: boxlo(3), boxhi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T :: delta(3), xlo(3), problo(3), time

      integer :: i, j, k
      REAL_T  :: x, y, z

      do k = lo(3), hi(3)
#if ( AMREX_SPACEDIM == 3 )
         z = (float(k)+.5)*delta(3)+problo(3)
         if (z.ge.boxlo(3) .and. z.le.boxhi(3)) then
#endif
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+problo(2)
               if (y.ge.boxlo(2) .and. y.le.boxhi(2)) then
                  do i = lo(1), hi(1)
                     x = (float(i)+.5)*delta(1)+problo(1)
                     if (x.ge.boxlo(1) .and. x.le.boxhi(1)) then
                        tag(i,j,k) = set
                     endif
                  end do
               endif
            end do
#if ( AMREX_SPACEDIM == 3 )
         endif
#endif
      end do

   end subroutine box_error

!=========================================================
!  This routine averages the mac face velocities for makeforce at half time
!=========================================================

   subroutine FORT_AVERAGE_EDGE_STATES( vel, v_lo, v_hi,&
                                        umacx, ux_lo, ux_hi,&
                                        umacy, uy_lo, uy_hi,&
#if ( AMREX_SPACEDIM == 3 )
                                        umacz, uz_lo, uz_hi,&
#endif
                                        getForceVerbose)&
                                        bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer :: v_lo(3), v_hi(3)
      integer :: ux_lo(3), ux_hi(3)
      integer :: uy_lo(3), uy_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer :: uz_lo(3), uz_hi(3)
#endif
      integer :: getForceVerbose
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), dim) :: vel
      REAL_T, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      REAL_T, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(uz_lo(1):uz_hi(1),uz_lo(2):uz_hi(2),uz_lo(3):uz_hi(3)) :: umacz
#endif

      REAL_T  :: velmin(3)
      REAL_T  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, dim
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      do k = v_lo(3), v_hi(3)
         do j = v_lo(2), v_hi(2)
            do i = v_lo(1), v_hi(1)
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
#if ( AMREX_SPACEDIM == 3 )
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
#endif
               do n = 1, dim
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, dim
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

end module PeleLM_nd
