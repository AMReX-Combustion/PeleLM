#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif 

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "mechanism.h"

module PeleLM_nd

  use fuego_chemistry
  use amrex_fort_module,  only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_filcc_module, only : amrex_filccn

  implicit none

  private

  public :: pphys_HMIXfromTY, pphys_RHOfromPTY, pphys_CPMIXfromTY, pphys_TfromHY, &
            init_data_new_mech, &
            dqrad_fill, divu_fill, dsdt_fill, ydot_fill, rhoYdot_fill, &
            repair_flux, conservative_T_floor, &
            part_cnt_err, mcurve, smooth, &
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
      REAL_T  :: RU, RUC, P1ATM, Ptmp, Yt(NUM_SPECIES), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
      SCAL = one * 1000
      CALL CKRP(RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
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
      REAL_T, dimension(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),NUM_SPECIES), intent(in) :: Y

! Local
      REAL_T  :: Yt(NUM_SPECIES), SCAL
      integer :: i, j, k, n

!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
      SCAL = 1.0d-4

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            do n = 1, NUM_SPECIES
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
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NUM_SPECIES) :: flux
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NUM_SPECIES) :: RhoY

      integer,  intent(in   ) :: xstatelo(3), xstatehi(3)
      integer,  intent(in   ) :: ystatelo(3), ystatehi(3)
      integer,  intent(in   ) :: axlo(3), axhi(3)
      integer,  intent(in   ) :: aylo(3), ayhi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer,  intent(in   ) :: zstatelo(3), zstatehi(3)
      integer,  intent(in   ) :: azlo(3), azhi(3)
#endif

      REAL_T,  intent(in) ::  xstate(xstatelo(1):xstatehi(1),xstatelo(2):xstatehi(2),xstatelo(3):xstatehi(3),NUM_SPECIES)
      REAL_T,  intent(in) ::  ystate(ystatelo(1):ystatehi(1),ystatelo(2):ystatehi(2),ystatelo(3):ystatehi(3),NUM_SPECIES)
      REAL_T,  intent(in) ::  afrac_x(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
      REAL_T,  intent(in) ::  afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
#if ( AMREX_SPACEDIM == 3 ) 
      REAL_T,  intent(in) ::  zstate(zstatelo(1):zstatehi(1),zstatelo(2):zstatehi(2),zstatelo(3):zstatehi(3),NUM_SPECIES)
      REAL_T,  intent(in) ::  afrac_z(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
#endif
      
      integer :: i, j, k, n
      REAL_T :: sumFlux, RhoYe(NUM_SPECIES), sumRhoYe

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
                  do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = xstate(i,j,k,n)
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  end do
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,NUM_SPECIES
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                  end do
                  else
                  do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i-1,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                   do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = ystate(i,j,k,n) 
                     sumRhoYe = sumRhoYe + RhoYe(n)
                   enddo
                   sumRhoYe = 1.0D0/sumRhoYe
                   do n=1,NUM_SPECIES
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                   end do         
                 else
                   do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j-1,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                   do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = zstate(i,j,k,n) 
                     sumRhoYe = sumRhoYe + RhoYe(n)
                   enddo
                   sumRhoYe = 1.0D0/sumRhoYe
                   do n=1,NUM_SPECIES
                     flux(i,j,k,n) = flux(i,j,k,n) - sumFlux*RhoYe(n)*sumRhoYe
                   end do
                  
                 else
                   do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k-1,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NUM_SPECIES) :: flux
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NUM_SPECIES) :: RhoY
      
      integer :: i, j, k, n
      REAL_T :: sumFlux, RhoYe(NUM_SPECIES), sumRhoYe

      if (dir.eq.0) then

!     First, assume away from physical boundaries, then use boundary-aware version below if applicable

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)               
                  sumFlux = 0.d0
                  sumRhoYe = 0.d0
                  do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i-1,j,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  end do
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i-1,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                  do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j-1,k,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j-1,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                  do n=1,NUM_SPECIES
                     sumFlux = sumFlux + flux(i,j,k,n)
                     RhoYe(n) = 0.5d0*(RhoY(i,j,k-1,n) + RhoY(i,j,k,n))
                     sumRhoYe = sumRhoYe + RhoYe(n)
                  enddo
                  sumRhoYe = 1.0D0/sumRhoYe
                  do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k-1,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
                     do n=1,NUM_SPECIES
                        sumFlux = sumFlux + flux(i,j,k,n)
                        sumRhoYe = sumRhoYe + RhoY(i,j,k,n)
                     enddo
                     sumRhoYe = 1.0D0/sumRhoYe
                     do n=1,NUM_SPECIES
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
