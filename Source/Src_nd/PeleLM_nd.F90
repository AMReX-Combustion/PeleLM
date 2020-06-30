#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif 

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "PPHYS_CONSTANTS.H"
#include "mechanism.h"

module PeleLM_nd

  use amrex_fort_module,  only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_filcc_module, only : amrex_filccn

  implicit none

  private

  public :: pphys_HMIXfromTY, pphys_RHOfromPTY, pphys_TfromHY, &
            init_data_new_mech, &
            dqrad_fill, divu_fill, dsdt_fill, ydot_fill, rhoYdot_fill, &
            part_cnt_err, &
            valgt_error, vallt_error, magvort_error, diffgt_error, &
            FORT_AVERAGE_EDGE_STATES

contains
 
!=========================================================
!  Compute mixture enthalpy from T and Y
!=========================================================

   subroutine pphys_HMIXfromTY(lo, hi, &
                               Hmix, h_lo, h_hi, & 
                               T, t_lo, t_hi, &
                               Y, y_lo, y_hi)&
                               bind(C, name="pphys_HMIXfromTY")

      use fuego_chemistry, only : CKHBMS

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: h_lo(3), h_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3)) :: Hmix
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),NUM_SPECIES) :: Y

! Local
      REAL_T  :: Yt(NUM_SPECIES), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
      SCAL = 1.0d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                do n = 1, NUM_SPECIES
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

      use fuego_chemistry, only : CKRHOY

      implicit none

! In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: y_lo(3), y_hi(3)
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: Rho
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: T
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),NUM_SPECIES) :: Y
      REAL_T  :: Patm

! Local
      REAL_T  :: Ptmp, Yt(NUM_SPECIES), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      Ptmp = Patm * PP_PA_CGS
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                do n = 1, NUM_SPECIES
                   Yt(n) = Y(i,j,k,n)
                end do
                CALL CKRHOY(Ptmp,T(i,j,k),Yt,RHO(i,j,k))
                RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do

   end subroutine pphys_RHOfromPTY

!=========================================================
!  Iterate on T until it matches Hmix and Ym
!=========================================================

   function pphys_TfromHY(lo, hi, &
                          T, t_lo, t_hi, &
                          Hmix, h_lo, h_hi, &
                          Y, y_lo, y_hi, &
                          errMax, NiterMAX, res) &
                          bind(C, name="pphys_TfromHY") result(MAXiters)

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
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),NUM_SPECIES), intent(in) :: Y
      REAL_T, intent(in)  :: errMAX
      REAL_T, dimension(0:NiterMAX-1), intent(out)  :: res

! Local
      REAL_T :: Yt(NUM_SPECIES)
      integer :: i, j, k, n, Niter, MAXiters

      MAXiters = 0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               do n = 1, NUM_SPECIES
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
   end function pphys_TfromHY

!=========================================================
!  init data to transition between 2 chem mechanism.
!=========================================================

   subroutine init_data_new_mech (level, time, lo, hi, nscal, &
                                  vel, scal, s_lo, s_hi, &
                                  press, p_lo, p_hi, &
                                  delta, xlo, xhi)&
                                  bind(C, name="init_data_new_mech")

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

      Patm = pamb / PP_PA_MKS

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
               do n = 0,NUM_SPECIES-1
                  scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
               enddo
               scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
            enddo
         enddo
      enddo
 
   end subroutine init_data_new_mech

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

                     if (x.ge.boxlo(1) .and. x.le.boxhi(1) &
                          .and. y.ge.boxlo(2) .and. y.le.boxhi(2) &
#if ( AMREX_SPACEDIM == 3 )
                          .and. z.ge.boxlo(3) .and. z.le.boxhi(3) &
#endif
                          ) then
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
