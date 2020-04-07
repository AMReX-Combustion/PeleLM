
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

module prob_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  use fuego_chemistry

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


      use PeleLM_F,  only: pphys_getP1atm_MKS

      use mod_Fvar_def, only : pamb

      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xgauss, ygauss, rgauss, ampgauss, gauss_type

      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i

      namelist /fortin/ meanFlowMag, meanFlowDir, T_mean, P_mean, &
                        xgauss, ygauss, rgauss, ampgauss, gauss_type
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
      pamb = pphys_getP1atm_MKS()

      meanFlowDir = 1.0
      meanFlowMag = 1.0d0
      T_mean = 298.0d0
      P_mean = pamb
      xgauss = 0.5
      ygauss = 0.5
      rgauss = 0.1
      ampgauss = 0.1
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

      use network,   only: nspecies
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, domnlo


      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xgauss, ygauss, rgauss, ampgauss, gauss_type

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
      REAL_T :: x, y, z, Yl(nspecies), Patm
      REAL_T :: dx
      integer :: i, j, k, n

      do k = lo(3), hi(3)
         z = (float(k)+.5d0)*delta(3)+domnlo(3)
         do j = lo(2), hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5d0)*delta(1)+domnlo(1)

               scal(i,j,k,Temp) = T_mean
               Yl(1) = 0.233d0
               Yl(2) = 0.767d0

               dx = x - xgauss
               dy = y - ygauss
               d_sq = dx*dx + dy*dy
               r_sq = rgauss*rgauss

               IF ( gauss_type == "Spec" ) THEN
                  Yl(1) = Yl(1) + Yl(1) * ampgauss * exp(-d_sq/r_sq)
                  Yl(2) = 1.0d0 - Yl(1)
               ELSE IF ( gauss_type == "Temp" ) THEN
                  scal(i,j,k,Temp) = T_mean + T_mean * ampgauss * exp(-d_sq/r_sq)
               ELSE
                  call bl_abort('gauss_type should be Spec or Temp')
               END IF

               do n = 1,nspecies
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

      Patm = P_mean / pphys_getP1atm_MKS()

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

   end subroutine init_data

end module prob_nd_module
