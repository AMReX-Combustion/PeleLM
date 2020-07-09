#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>
#include <PPHYS_CONSTANTS.H>
#include "mechanism.h"

#include <Prob_F.H>
#include <PeleLM_F.H>


module prob_nD_module

   use amrex_fort_module, only : dim=>amrex_spacedim
   use amrex_error_module, only : amrex_abort
   use fuego_chemistry

  implicit none

  private
  
  public :: amrex_probinit, setupbc, init_data

contains

! ::: -----------------------------------------------------------
! ::: This routine is called at problem initialization time
! ::: and when restarting from a checkpoint file.
! ::: The purpose is (1) to specify the initial time value
! ::: (not all problems start at time=0.0) and (2) to read
! ::: problem specific data from a namelist or other input
! ::: files and possibly store them or derived information
! ::: in FORTRAN modules for later use.
! ::: 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: init      => TRUE if called at start of problem run
! :::              FALSE if called from restart
! ::: 
! ::: -----------------------------------------------------------

  subroutine amrex_probinit ( init, name, namlen, problo, probhi) bind(c)
  
      use mod_Fvar_def, only : pamb 
      use mod_Fvar_def, only : domnhi, domnlo
      use mod_Fvar_def, only : ac_hist_file, cfix, changemax_control, &
                               coft_old, controlvelmax, corr, dv_control, &
                               h_control, navg_pnts, scale_control, sest, &
                               tau_control, tbase_control, V_in, v_in_old, zbase_control, &
                               pseudo_gravity
      use probdata_module, only : T_in, splitx
      use probdata_module, only : midtanh, widthtanh, H2_enrich
      use user_defined_fcts_nd_module, only: set_Zst
      
      
      implicit none

! In/Out      
      integer, intent(in) :: init, namlen
      integer, intent(in) :: name(namlen)
      REAL_T, intent(in)  :: problo(dim), probhi(dim)

! Local     
      integer, parameter :: maxlen = 256
      integer :: realnamlen,untin, isioproc
      character(len=maxlen) :: probin
      integer :: i

! Namelists
      namelist /fortin/ V_in, T_in, midtanh, widthtanh, H2_enrich 

      namelist /heattransin/ pamb

      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
          zbase_control, pseudo_gravity, corr,controlVelMax,navg_pnts


      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
!         call bl_abort('probinit called with init ne 1')
      end if

!------------------------------------------      
!     Convert probin file name from C++ -> fortran
!------------------------------------------      
      if (namlen .gt. maxlen) then
         call amrex_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         realnamlen = 6
         probin(1:realnamlen) = 'probin'
      else
         realnamlen = namlen
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif

!------------------------------------------      
!     Set defaults
!------------------------------------------      
      pamb = PP_PA_MKS

!     Inflow/initial conditions parameters      
      T_in = 300.0d0                         ! Inflow temperature
      H2_enrich = 0.0d0                      ! Level of H2 enrichment on the fuel side
      midtanh = 0.6*(domnhi(1)+domnlo(1))    ! Mixing layer position: default around middomain
      widthtanh = 0.05*(domnhi(1)-domnlo(1)) ! Mixing layer thickness: default is 1/20th of domain width in x 
      splitx = 0.5d0*(domnhi(1)+ domnlo(1))  ! Half width position (in x direction)

!     Initialize inflow velocity control variables
      zbase_control = 0.d0
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      h_control = -one
      dV_control = zero
      tbase_control = zero
      h_control = -one
      pseudo_gravity = 0
      navg_pnts = 10

!------------------------------------------      
!     Read probin file 
!------------------------------------------      
      untin = 9
      open(untin,file=probin(1:realnamlen),form='formatted',status='old')
      read(untin,fortin)
      read(untin,heattransin)
      read(untin,control)
      close(unit=untin)

!------------------------------------------      
!     Initialize module data 
!------------------------------------------      
!     Initialize control variables that depend on fortin variables
      V_in_old = V_in

!     Set Zst
      call set_Zst()

!     Set up boundary functions
      call setupbc()

!     Pass inflow velocity control variables
      scale_control = 1.0
      if (h_control .gt. zero) then
         cfix = h_control
      endif

!     Dump read data on screen for debug purposes
      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
         write(6,control)
      end if

  end subroutine amrex_probinit
  
!------------------------------------
  
  subroutine setupbc()bind(C, name="setupbc")

    use mod_Fvar_def, only : V_in
    use probdata_module, only : T_bc, u_bc, v_bc
    use probdata_module, only : bcinit, T_in
  
    implicit none

    T_bc = T_in
    u_bc = zero
    v_bc = V_in
              
    bcinit = .true.

  end subroutine setupbc

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
! :::		   this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

  subroutine init_data(level, time, lo, hi, nscal, &
                       vel, scal, s_lo, s_hi, press, p_lo, p_hi, &
                       delta, xlo, xhi) &
                       bind(C, name="init_data")
                              

      use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb
      use mod_Fvar_def, only : bathID, oxidID, domnhi, domnlo, V_in
      use probdata_module, only : midtanh, widthtanh, splitx, T_in
      use user_defined_fcts_nd_module, only : set_Y_from_ksi

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

      integer :: nPMF

      integer :: i, j, k, n
      REAL_T  :: x, y, z, Yl(0:NUM_SPECIES-1,2), Xl(NUM_SPECIES), Patm
      REAL_T  :: pmf_vals(NUM_SPECIES+3), y1, y2
      REAL_T  :: pert,Lx, tanhval, rad, y_lo
      REAL_T  :: Gaussian_T, Gauss_T_width, Gauss_maxT
      REAL_T  :: Gaussian_Spec, Gaussian_Spec_width
      REAL_T  :: sumspec

!--------------------------------------------------------    
! Data initialization for triple flame case: 
!     - first set-up a mixing layer with a tanh in mixture fraction
!       between 0 and 1. It uses the fuel declared in input file. A
!       constant temperature is set first.
!     - Then a widening in the y-dir Gaussian profile of temperature
!       is super-imposed on the base flow. The composition in the hot
!       is replaced by air to avoid preheating premixed gas.
!         * y_lo : lower tip of the hot region
!         * Gauss_maxT : maximum temperature
!         * Gauss_T_width : width of the temperature profile
!         * Gaussian_Spec_width : width of the hot air profile
!
!     - Setting Gaussian_Spec_width > Gauss_T_width enable to avoid
!       heating fuel/air mixture at the border of the hot region.
!--------------------------------------------------------

      y_lo = 0.01d0
      Gauss_maxT = 1900.0d0
      Gauss_T_width = 0.0012d0
      Gaussian_Spec_width = 0.0012d0

!     write(6,*)" made it to initdata"
      if (bathID.lt.1 .or. bathID.gt.NUM_SPECIES) then
         call bl_pd_abort()
      endif

!     Set up hot air composition
      Yl(:,1) = 0.0d0
      Yl(bathID-1,1) = 0.767
      Yl(oxidID-1,1) = 0.233

      Patm = pamb / PP_PA_MKS

      do k = lo(3), hi(3)
         z = (float(k)+.5d0)*delta(3)+domnlo(3)
         do j = lo(2), hi(2)
           y = (float(j)+.5d0)*delta(2)+domnlo(2)
           do i = lo(1), hi(1)
             x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  
             ! Start with mixture fraction tanh between 0 and 1
             tanhval = 0.5d0*(1.0d0+TANH((x-midtanh)/widthtanh))
             call set_Y_from_Ksi(tanhval,Yl(0:NUM_SPECIES-1,2))

             ! Set background (inflow) T
             scal(i,j,k,Temp) = T_in

             ! Setup hot air region in the mixing layer
             rad = SQRT((x-splitx)**2.0+(MIN(y,y_lo)-y_lo)**2.0)
             Gaussian_T = EXP(-(rad)**2.0/(2.0*(Gauss_T_width + Gauss_T_width * (y-y_lo)/0.03)**2.0))
             scal(i,j,k,Temp) = T_in + (Gauss_maxT-T_in) * Gaussian_T

             Gaussian_spec = EXP(-(rad)**2.0/(2.0*(Gaussian_Spec_width + Gaussian_Spec_width * (y-y_lo)/0.03)**2.0))
             do n = 0, NUM_SPECIES-1
                scal(i,j,k,FirstSpec+n) = (1.0d0 - Gaussian_Spec) * Yl(n,2)  + Gaussian_spec * Yl(n,1)
             end do

             ! Maybe not necessary, but to be sure, normalized mass fraction ...
             sumspec = SUM(scal(i,j,k,FirstSpec:FirstSpec+NUM_SPECIES-1))
             do n = 0, NUM_SPECIES-1
                scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n) / sumspec
             end do

             vel(i,j,k,1) = 0.d0
             vel(i,j,k,2) = V_in 
#if ( AMREX_SPACEDIM == 3 ) 
             vel(i,j,k,3) = 0.d0
#endif

           end do
         end do
      end do

!     Compute rho and h from prescribed composition, temperature and pression
      call pphys_RHOfromPTY(lo,hi, &
                            scal(:,:,:,Density),   s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi, &
                            Patm)
      call pphys_HMIXfromTY(lo,hi, &
                            scal(:,:,:,RhoH),      s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi)

!     Yl -> rhoYl, h -> rhoYl
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

end module prob_nD_module
