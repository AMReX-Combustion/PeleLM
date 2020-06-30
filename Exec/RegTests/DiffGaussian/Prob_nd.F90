
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>
#include "mechanism.h"
#include <PPHYS_CONSTANTS.H>

module prob_nd_module

  use fuego_chemistry
  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort

  implicit none

  private
  
  public :: amrex_probinit, init_data

contains

! ::: -----------------------------------------------------------
! ::: This routine is called at problem initialization time
! ::: and when restarting from a checkpoint file.
! ::: The purpose is (1) to specify the initial time value
! ::: (not all problems start at time=0.0) and (2) to read
! ::: problem specific data from a namelist or other input
! ::: files and possibly store them or derived information
! ::: in FORTRAN common blocks for later use.
! ::: 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: init      => TRUE if called at start of problem run
! :::              FALSE if called from restart
! ::: strttime <=  start problem with this time variable
! ::: 
! ::: -----------------------------------------------------------

   subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)



      use mod_Fvar_def, only : pamb

      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xgauss, ygauss, rgauss, ampgauss, gauss_type, start_time, &
                                 DiffCoeff   

      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i

      namelist /fortin/ meanFlowMag, meanFlowDir, T_mean, P_mean, &
                        xgauss, ygauss, rgauss, ampgauss, gauss_type, start_time, DiffCoeff
      namelist /heattransin/ pamb


!
!      Build `probin' filename -- the name of file containing fortin namelist.
!
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
!         call bl_abort('probinit called with init ne 1')
      end if

      if (namlen .gt. maxlen) then
         call bl_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         namlen = 6
         probin(1:namlen) = 'probin'
      else
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
!     Set defaults
      pamb = PP_PA_MKS

      meanFlowDir = 1.0
      meanFlowMag = 1.0d0
      T_mean = 298.0d0
      P_mean = pamb
      xgauss = 0.5
      ygauss = 0.5
      rgauss = 0.1
      ampgauss = 0.1
      start_time = 0.0d0
      DiffCoeff = 0.0d0
      gauss_type = "Spec"

      read(untin,fortin)
      
      read(untin,heattransin)
 
      close(unit=untin)

!     Do some checks on COVO params     
      IF (       meanFlowDir /= 1 .AND. meanFlowDir /= -1 &
           .AND. meanFlowDir /= 2 .AND. meanFlowDir /= -2 &
           .AND. meanFlowDir /= 3 .AND. meanFlowDir /= -3  ) THEN
          WRITE(*,*) " meanFlowDir should be either: "
          WRITE(*,*) " +/-1 for x direction" 
          WRITE(*,*) " +/-2 for y direction" 
          WRITE(*,*) " +/-3 for diagonal direction" 
          WRITE(*,*) " Note: the mean flow direction(s) must be periodic "
          CALL bl_abort('Correct meanFlowDir value !')
      END IF
      
      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if

  end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  The velocity field you
! ::: provide does not have to be divergence free and the pressure
! ::: field need not be set.  A subsequent projection iteration
! ::: will define aa divergence free velocity field along with a
! ::: consistant pressure.
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nscal     => number of scalar quantities.  You should know
! :::              this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::              ghost region).
! ::: -----------------------------------------------------------

   subroutine init_data(level, time, lo, hi, nscal, &
                        vel, scal, s_lo, s_hi, press, p_lo, p_hi, &
                        delta, xlo, xhi) &
                        bind(C, name="init_data")

      use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, domnlo


      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xgauss, ygauss, rgauss, ampgauss, gauss_type, start_time, DiffCoeff

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
      REAL_T :: dy, d_sq, r_sq
      REAL_T :: x, y, z, Yl(NUM_SPECIES), Patm
      REAL_T :: dx, y_lo, y_hi, x_lo, x_hi, dx_lo, dx_hi, dy_lo, dy_hi
      REAL_T :: pi_atan, inv_area, denom, sqrtdenom, scaleamp
      integer :: i, j, k, n

      pi_atan = 4.0d0*atan2(1.0,1.0)

      do k = lo(3), hi(3)
         z = (float(k)+.5d0)*delta(3)+domnlo(3)
         do j = lo(2), hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            y_lo = (float(j))*delta(2)+domnlo(2)
            y_hi = (float(j)+1.0d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5d0)*delta(1)+domnlo(1)
               x_lo = (float(i))*delta(1)+domnlo(1)
               x_hi = (float(i)+1.0d0)*delta(1)+domnlo(1)

               scal(i,j,k,Temp) = T_mean
               Yl(1) = 0.233d0
               Yl(2) = 0.767d0

               dx = x - xgauss
               dx_lo = x_lo - xgauss
               dx_hi = x_hi - xgauss
               dy = y - ygauss
               dy_lo = y_lo - ygauss
               dy_hi = y_hi - ygauss
               d_sq = dx*dx + dy*dy
               r_sq = rgauss*rgauss
               denom = r_sq + 4.0d0 * DiffCoeff * start_time
               sqrtdenom = sqrt(r_sq + 4.0d0 * DiffCoeff * start_time)
               scaleamp = 1.0d0 / (sqrt(4.0d0 * DiffCoeff * start_time / r_sq + 1.0d0))
               inv_area = 1.0d0 / ( delta(1) * delta(2) )

               IF ( gauss_type == "Spec" ) THEN
                  ! Cell center value
                  Yl(1) = Yl(1) + Yl(1) * ampgauss * exp(-d_sq/r_sq)
                  ! Cell average value
                  Yl(1) = Yl(1) + Yl(1) * ampgauss * scaleamp * 0.25d0 * pi_atan * denom &
                                                   * (erf(dx_lo/sqrtdenom) - erf(dx_hi/sqrtdenom)) &
                                                   * (erf(dy_lo/sqrtdenom) - erf(dy_hi/sqrtdenom)) &
                                                   * inv_area
                  Yl(2) = 1.0d0 - Yl(1)
               ELSE IF ( gauss_type == "Temp" ) THEN
                  ! Cell center value
                  scal(i,j,k,Temp) = T_mean + T_mean * ampgauss * exp(-d_sq/r_sq)
                  ! Cell average value
                  scal(i,j,k,Temp) = T_mean + T_mean * ampgauss * scaleamp * 0.25d0 * pi_atan * denom &
                                                     * (erf(dx_lo/sqrtdenom) - erf(dx_hi/sqrtdenom)) &
                                                     * (erf(dy_lo/sqrtdenom) - erf(dy_hi/sqrtdenom)) &
                                                     * inv_area
               ELSE
                  call bl_abort('gauss_type should be Spec or Temp')
               END IF

               do n = 1,NUM_SPECIES
                  scal(i,j,k,FirstSpec+n-1) = Yl(n)
               end do

               SELECT CASE ( meanFlowDir )
                  CASE (1)
                     vel(i,j,k,1) = meanFlowMag
                     vel(i,j,k,2) = 0.0d0
                  CASE (-1)
                     vel(i,j,k,1) = -meanFlowMag
                     vel(i,j,k,2) = 0.0d0
                  CASE (2)
                     vel(i,j,k,1) = 0.0d0
                     vel(i,j,k,2) = meanFlowMag
                  CASE (-2)
                     vel(i,j,k,1) = 0.0d0
                     vel(i,j,k,2) = -meanFlowMag
                  CASE (3)
                     vel(i,j,k,1) = meanFlowMag
                     vel(i,j,k,2) = meanFlowMag
                  CASE (-3)
                     vel(i,j,k,1) = -meanFlowMag
                     vel(i,j,k,2) = -meanFlowMag
               END SELECT

            end do

         end do
      end do

      Patm = P_mean / PP_PA_MKS

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

   end subroutine init_data

end module prob_nd_module
