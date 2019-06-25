
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

module prob_3D_module

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
                                   splitx, xfrontw, fuel_N2_vol_percent, &
                                   blobr, bloby, blobx, Tfrontw, blobT  


      implicit none

      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i

      namelist /fortin/ V_in, T_in, V_co, phi_in, T_co, &
                        splitx, xfrontw, fuel_N2_vol_percent, &
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
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, domnhi, maxspec, maxspnml, V_in
    use mod_Fvar_def, only : fuelID, oxidID
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, w_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, Phi_in, T_in, V_co, T_co, fuel_N2_vol_percent
    use user_defined_fcts_3d_module, only : getZone

    implicit none

    REAL_T Patm
    REAL_T Xt(maxspec), Yt(maxspec)

    integer b(3)
    integer zone, n, fuelZone, airZone
    integer num_zones_defined
    integer iN2, iO2,iCH3OCH3
    character*(maxspnml) name
    data  b / 1, 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()
    num_zones_defined = 0


!     A diffusion flame
    fuelZone = getZone(domnlo(1), 0.5*(domnlo(2)+domnhi(2)), domnlo(3))
    airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))
    num_zones_defined = 2

!     Fuel
    do n = 1,nspec
      Xt(n) = 0.d0
    end do 
    Xt(iN2) = fuel_N2_vol_percent*1.d-2
    Xt(fuelID) = 1.d0 - Xt(iN2)            

    CALL CKXTY (Xt, Yt)

    do n=1,Nspec
      Y_bc(n-1,fuelZone) = Yt(n)
    end do
    T_bc(fuelZone) = T_in
    u_bc(fuelZone) = 0.d0
    v_bc(fuelZone) = 0.d0
    w_bc(fuelZone) = V_in

!     Air
    do n=1,Nspec
      Xt(n) = zero
    enddo
    Xt(oxidID) = 0.21d0
    Xt(iN2)    = 1.d0 - Xt(oxidID)

    CALL CKXTY (Xt, Yt)         
    do n=1,Nspec
      Y_bc(n-1,airZone) = Yt(n)
    end do
         
    T_bc(airZone) = T_co
    u_bc(airZone) = 0.d0
    v_bc(airZone) = 0.d0
    w_bc(airZone) = V_co

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
      use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : domnlo, domnhi, maxspec, maxspnml
      use probdata_module, only : blobr, blobT, blobx, bloby, &
                                  splitx, tfrontw, xfrontw, &
                                  Y_bc, u_bc, v_bc, T_bc, IDX_COFLOW, &
                                  T_co
      use user_defined_fcts_3d_module, only : getZone, bcfunction

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

      integer i, j, k, n, airZone, fuelZone
      REAL_T x, y, z, ztemp, r, Patm, sigma, eta
      REAL_T pert, u,v,w,rho,T,h, Yl(maxspec), Xl(maxspec)

      fuelZone = getZone(domnlo(1), domnlo(2), domnlo(3))
      airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))

      do k = lo(3), hi(3)
        z = (float(k)+.5)*delta(3)+domnlo(3)
        eta = 0.5d0*(1.d0 - TANH(2.d0*(z-blobr)/Tfrontw))
        do j = lo(2), hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)

            ztemp = domnlo(3)
            call bcfunction(x,y,ztemp,1,1,time,u,v,w,rho,Yl,T,h,delta,.true.)
            do n = 1,Nspec
              scal(i,j,k,FirstSpec+n-1) = Yl(n)* eta + (1.d0-eta) *Y_bc(n-1,IDX_COFLOW)
            end do

            pert = dexp(-2.d0*((z-1.0d0*blobr)**2 + (y-blobr)**2)/blobr**2)
            pert = pert+ dexp(-2.d0*((z-1.0d0*blobr)**2 + (y+blobr)**2)/blobr**2)
            pert = pert+ dexp(-2.d0*((z-2.0d0*blobr)**2 + (y)**2)/blobr**2)

            scal(i,j,k,Trac) = 0.d0
            scal(i,j,k,Temp) = (T* eta + (1.d0-eta) *T_co) * (1.d0+ 1.d0*pert)

            vel(i,j,k,1) = u
            vel(i,j,k,2) = v
            vel(i,j,k,3) = w

          enddo
        enddo
      enddo

      Patm = pamb / pphys_getP1atm_MKS()
!      write(6,*)"Patm",Patm

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
               do n = 0,Nspec-1
                  scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
               enddo
               scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
            enddo
         enddo
      enddo

  end subroutine init_data
 
end module prob_3D_module



