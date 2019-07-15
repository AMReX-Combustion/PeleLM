
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

module prob_2D_module

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
      use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
      use mod_Fvar_def, only : dim
      use probdata_module, only: T_mean, P_mean, xblob, yblob, radblob
      
      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i
 
      namelist /fortin/ T_mean, P_mean, xblob, yblob, radblob
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber


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
      dpdt_factor = 0.3d0
      closed_chamber = 0

      T_mean = 298.0d0
      P_mean = pamb
      xblob = 0.
      yblob = 0.
      radblob = 0.4

      read(untin,fortin)
      
      read(untin,heattransin)
 
      close(unit=untin)

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
! :::		   this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

  subroutine init_data(level,time,lo,hi,nscal, &
     	 	                   vel,scal,DIMS(state),press,DIMS(press), &
                           delta,xlo,xhi) &
                           bind(C, name="init_data")
                              
      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, Trac, dim
      use mod_Fvar_def, only : domnlo, maxspec
      use probdata_module, only: T_mean, P_mean, xblob, yblob, radblob
      
      implicit none
      integer    level, nscal
      integer    lo(dim), hi(dim)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     time, delta(dim)
      REAL_T     vel(DIMV(state),dim)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))


      integer i, j, n
      REAL_T x, y, Yl(maxspec), Patm
      REAL_T :: dist

      do j = lo(2), hi(2)
         y = (float(j)+.5d0)*delta(2) +domnlo(2)
         do i = lo(1), hi(1)
            x = (float(i)+.5d0)*delta(1) +domnlo(1)
            
            vel(i,j,1) = tanh(30.*(.5-abs(y)))
            vel(i,j,2) = .05*sin(Pi*x)

            Yl(1) = 0.233
            Yl(2) = 0.767
            
            do n = 1,Nspec
               scal(i,j,FirstSpec+n-1) = Yl(n)
            end do

            dist = sqrt((x-xblob)**2 + (y-yblob)**2)
            scal(i,j,Trac) = merge(one,zero,dist.lt.radblob)
            scal(i,j,Temp) = T_mean

         end do
      end do

      Patm = P_mean / pphys_getP1atm_MKS()

      call pphys_RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
          Patm)

      call pphys_HMIXfromTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state)) 

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo
      
  end subroutine init_data
      




end module prob_2D_module
