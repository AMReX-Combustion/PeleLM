
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

  public :: amrex_probinit, setupbc, init_data


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
      use mod_Fvar_def, only : V_in
      use probdata_module, only : T_in, V_co, phi_in, T_co, &
                                   splitx, xfrontw, &
                                   blobr, bloby, blobx, Tfrontw, blobT  


      implicit none

      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i

      namelist /fortin/ V_in, T_in, V_co, phi_in, T_co, &
                        splitx, xfrontw, &
                        blobr, bloby, blobx, Tfrontw, blobT
      namelist /heattransin/ pamb, dpdt_factor


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

      read(untin,fortin)

      read(untin,heattransin)

      close(unit=untin)

!     Set up boundary functions
      call setupbc()


      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if

  end subroutine amrex_probinit

!------------------------------------

  subroutine setupbc()bind(C, name="setupbc")

    use network,   only: nspec
    use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, maxspec, maxspnml, V_in
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, Phi_in, T_in, V_co, T_co
    use user_defined_fcts_2d_module, only : getZone

    implicit none

    REAL_T Patm
    REAL_T Xt(maxspec), Yt(maxspec)

    integer b(2)
    integer zone, n, fuelZone, airZone
    integer num_zones_defined
    integer iN2, iO2,iCH3OCH3
    character*(maxspnml) name
    data  b / 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()
    num_zones_defined = 0


!     A diffusion flame
         fuelZone = getZone(0.d0, domnlo(2))
         airZone  = getZone(domnlo(1), domnlo(2))
         num_zones_defined = 2

!     Fuel
!     find O2 and DME
         iO2 = -1
         do n=1,Nspec
            call pphys_get_spec_name2(name,n)
            if (name .eq. 'O2') iO2 = n
         end do
         if (iO2.eq.-1) &
             write(6,*) '.....warning: no O2 in chemistry species list'

         iN2 = -1
         do n=1,Nspec
            call pphys_get_spec_name2(name,n)
            if (name .eq. 'N2') iN2 = n
         end do
         if (iN2.eq.-1) &
             write(6,*) '.....warning: no N2 in chemistry species list'

         iCH3OCH3 = -1
         do n=1,Nspec
            call pphys_get_spec_name2(name,n)
            if (name .eq. 'CH3OCH3') iCH3OCH3 = n
         end do
         if (iCH3OCH3.eq.-1) &
             write(6,*) '.....warning: no CH3OCH3 in chemistry species list'

!      set everything to zero
         do n = 1,Nspec
            Xt(n) = 0.d0
         end do 

!     fuel shall be half O2, half H2 by volume
         Xt(iCH3OCH3) = phi_in
         Xt(iN2) = 1.d0-Xt(iCH3OCH3)

!     convert mole frac to mass frac
         CALL CKXTY (Xt, Yt)

!     set bc's
         do n=1,Nspec
            Y_bc(n-1,fuelZone) = Yt(n)
         end do
         T_bc(fuelZone) = T_in
         u_bc(fuelZone) = 0.d0
         v_bc(fuelZone) = V_in

!     Air 
         do n=1,Nspec
            Xt(n) = zero
         enddo
         Xt(iN2) = 0.79d0
         Xt(iO2) = 0.21

         CALL CKXTY (Xt, Yt)         
         do n=1,Nspec
            Y_bc(n-1,airZone) = Yt(n)
         end do
         
         T_bc(airZone) = T_co
         u_bc(airZone) = 0.d0
         v_bc(airZone) = V_co

      do zone=1,num_zones_defined
!     Set density and hmix consistent with data

         call pphys_RHOfromPTY(b, b, &
                             rho_bc(zone), DIMARG(b), DIMARG(b), &
                             T_bc(zone),   DIMARG(b), DIMARG(b), &
                             Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)
         call pphys_HMIXfromTY(b, b, &
                             h_bc(zone),   DIMARG(b), DIMARG(b), &
                             T_bc(zone),   DIMARG(b), DIMARG(b), &
                             Y_bc(0,zone), DIMARG(b), DIMARG(b))
      enddo
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

  subroutine init_data(level,time,lo,hi,nscal, &
                       vel,scal,DIMS(state),press,DIMS(press), &
                       delta,xlo,xhi) &
                       bind(C, name="init_data")


      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : domnlo, maxspec, maxspnml
      use probdata_module, only : blobr, blobT, blobx, bloby, &
                                  splitx, tfrontw, xfrontw, &
                                  Y_bc, u_bc, v_bc, T_bc
      use user_defined_fcts_2d_module, only : getZone

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

      integer i, j, n, airZone, fuelZone
      REAL_T x, y, r, Patm, sigma, eta


         fuelZone = getZone(0.d0, domnlo(2))
         airZone  = getZone(domnlo(1), domnlo(2))

         sigma = 2.5d0*xfrontw*splitx

         do j = lo(2), hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)

!              if (blobr.lt.0) then

                  eta = 0.5d0 * (tanh((x + splitx)/sigma) &
                               - tanh((x - splitx)/sigma))
                  
                  do n=1,Nspec
                   scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*(1.d0-eta) &
                        + eta*Y_bc(n-1,fuelZone)
                  enddo
                  scal(i,j,Temp) = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)
                  vel(i,j,1) = u_bc(airZone)*(1.d0-eta) + eta*u_bc(fuelZone)
                  vel(i,j,2) = v_bc(airZone)*(1.d0-eta) + eta*v_bc(fuelZone)
                  scal(i,j,Trac) = 0.d0
!                 
!                 eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
!                 do n=1,Nspec
!                    scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta
!    &                    + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
!                 enddo

!              else
               if (blobr.gt.0) then

                  eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
                  do n=1,Nspec
                     scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta &
                         + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
                  enddo
               
                  
!     Superimpose blob of hot air
                  r = SQRT((x-blobx)**2 + (y-bloby)**2)
                  eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/Tfrontw))
                  do n=1,Nspec
                     scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta &
                         + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
                  enddo
                  scal(i,j,Temp) = blobT*eta + (1.d0-eta)*scal(i,j,Temp)
!                 vel(i,j,1) = 0.d0
!                 vel(i,j,2) = v_bc(airZone)
                  
                  vel(i,j,1) = u_bc(airZone)*eta + (1.d0-eta)*u_bc(fuelZone)
                  vel(i,j,2) = v_bc(airZone)*eta + (1.d0-eta)*v_bc(fuelZone)
!                 scal(i,j,Trac) = 0.d0
                  
               endif

            enddo
         enddo

       Patm = pamb / pphys_getP1atm_MKS()
!      write(6,*)"Patm",Patm

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



