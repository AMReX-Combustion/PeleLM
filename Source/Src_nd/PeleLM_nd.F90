
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

module PeleLM_nd

  use fuego_chemistry
  use amrex_fort_module, only : dim=>amrex_spacedim
  use network,           only : nspecies

  implicit none

  private

  public :: floor_spec, calc_divu_fortran, calc_gamma_pinv

contains

!=========================================================
!  Floor species mass fraction at 0.0
!=========================================================

   subroutine floor_spec( lo, hi, spec, s_lo, s_hi) &
                        bind(C, name="floor_spec")

      implicit none

! In/Out
      integer :: lo(3), hi(3)
      integer :: s_lo(3), s_hi(3)
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
      integer :: lo(3),hi(3)
      integer :: d_lo(3), d_hi(3)
      integer :: rd_lo(3), rd_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: vT_lo(3), vT_hi(3)
      integer :: rY_lo(3), rY_hi(3)
      integer :: t_lo(3), t_hi(3)
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: divu
      REAL_T, dimension(rd_lo(1):rd_hi(1),rd_lo(2):rd_hi(2),rd_lo(3):rd_hi(3),nspecies) :: rYdot
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nspecies) :: vtY
      REAL_T, dimension(vT_lo(1):vT_hi(1),vT_lo(2):vT_hi(2),vT_lo(3):vT_hi(3)) :: vtT
      REAL_T, dimension(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies) :: rhoY
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T

! Local
      REAL_T, dimension(1:nspecies) :: Y, H, invmtw
      REAL_T :: cpmix, rhoInv, tmp, mmw
      integer :: i, j, k, n

      call CKWT(invmtw)
      do n = 1,nspecies
         invmtw(n) = one / invmtw(n)
      end do

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
               CALL CKHMS(T(i,j,k),H)
               CALL CKMMWY(Y,mmw)

               cpmix = cpmix*1.d-4
               do n=1,nspecies
                  H(n) = H(n)*1.d-4
               enddo

               tmp = 0.d0
               divu(i,j,k) = vtT(i,j,k)
               do n=1,nspecies
                  tmp = tmp + (rYdot(i,j,k,n)+vtY(i,j,k,n))*invmtw(n)
                  divu(i,j,k) = divu(i,j,k) - rYdot(i,j,k,n)*H(n)
               enddo
               divu(i,j,k) = ( divu(i,j,k)/(cpmix*T(i,j,k)) + tmp*mmw ) * rhoInv
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
      integer :: lo(3),hi(3)
      integer :: th_lo(3), th_hi(3)
      integer :: r_lo(3), r_hi(3)
      integer :: t_lo(3), t_hi(3)
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

!-----------------------------------  
  
end module PeleLM_nd
