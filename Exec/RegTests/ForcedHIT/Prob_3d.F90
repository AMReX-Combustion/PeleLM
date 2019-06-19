
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
      use mod_Fvar_def, only : domnhi, domnlo, dim
      use probdata_module, only : T_in, FTX, FTY, FTZ, TAT, TAP, FPX, FPY, FPZ, &
                                  FAX, FAY, FAZ, FPXX, FPYX, FPXY, FPYY, &
                                  FPZY, FPYZ, FPZX, FPXZ, FPZZ, &
                                  blrandseed, div_free_force, force_scale, &
                                  forcing_epsilon, forcing_time_scale_max, forcing_time_scale_min, &
                                  mode_start, moderate_zero_modes, nmodes, spectrum_type, time_offset
      
      implicit none
      
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)
      REAL_T     Lx, Ly, Lz
     double precision :: twicePi, kxd, kyd, kzd, rn
    double precision :: thetaTmp, phiTmp
    double precision :: cosThetaTmp, cosPhiTmp
    double precision :: sinThetaTmp, sinPhiTmp
    double precision :: px, py, pz, mp2, Ekh
    integer :: kx, ky, kz, mode_count, reduced_mode_count
    integer :: xstep, ystep, zstep
    integer :: nxmodes, nymodes, nzmodes

    double precision :: Lmin
    double precision :: kappa, kappaMax, freqMin, freqMax, freqDiff
      integer i

      namelist /fortin/ T_in
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber
      namelist /forcing/ spectrum_type, mode_start, nmodes, &
                      force_scale, forcing_time_scale_min, forcing_time_scale_max, &
                      forcing_time_scale_min, forcing_time_scale_max, &
                      time_offset, div_free_force

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

      div_free_force = 0
      spectrum_type          = 0
      mode_start             = 1
      nmodes                 = 4
      force_scale            = 1.0d10
      forcing_time_scale_min = 0.0
      forcing_time_scale_max = 0.0


      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
!     Set defaults
      pamb = pphys_getP1atm_MKS()
      closed_chamber = 1

      T_in = 300.0d0

      read(untin,fortin)
      
      read(untin,heattransin)
 
      read(untin,forcing)

      close(unit=untin)

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
         write(6,forcing)
      end if


    twicePi = two*Pi

    if (blrandseed.gt.0) then
       call blutilinitrand(blrandseed)
       call blutilrand(rn)
       call blutilinitrand(blrandseed)
       if (isioproc .eq. 1) then 
          write (*,*) "blrandseed = ",blrandseed
          write (*,*) "first random number = ",rn
       endif
    else
       call blutilinitrand(111397)
    endif

    Lx = domnhi(1)-domnlo(1)
    Ly = domnhi(2)-domnlo(2)
    Lz = domnhi(3)-domnlo(3)

    if (isioproc .eq. 1) then
       write(*,*) "Lx = ",Lx
       write(*,*) "Ly = ",Ly
       write(*,*) "Lz = ",Lz
    endif

    Lmin = min(Lx,Ly,Lz)
    kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
    nxmodes = nmodes*int(0.5+Lx/Lmin)
    nymodes = nmodes*int(0.5+Ly/Lmin)
    nzmodes = nmodes*int(0.5+Lz/Lmin)
    if (isioproc .eq. 1) then
      write(*,*) "Lmin = ",Lmin
      write(*,*) "kappaMax = ",kappaMax
      write(*,*) "nxmodes = ",nxmodes
      write(*,*) "nymodes = ",nymodes
      write(*,*) "nzmodes = ",nzmodes
    endif

    allocate(FTX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FTY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FTZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(TAT(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(TAP(1:nxmodes, 1:nymodes, 1:nzmodes))

    allocate(FPX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAZ(1:nxmodes, 1:nymodes, 1:nzmodes))

    allocate(FPXX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPXY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPXZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZZ(1:nxmodes, 1:nymodes, 1:nzmodes))


    if (forcing_time_scale_min.eq.zero) then
      forcing_time_scale_min = half
    endif
    if (forcing_time_scale_max.eq.zero) then
       forcing_time_scale_max = one
    endif

    freqMin = one/forcing_time_scale_max
    freqMax = one/forcing_time_scale_min
    freqDiff= freqMax-freqMin

    if (isioproc .eq. 1) then
       write(*,*) "forcing_time_scale_min = ",forcing_time_scale_min
       write(*,*) "forcing_time_scale_max = ",forcing_time_scale_max
       write(*,*) "freqMin = ",freqMin
       write(*,*) "freqMax = ",freqMax
       write(*,*) "freqDiff = ",freqDiff
    endif

    mode_count = 0

    xstep = int(Lx/Lmin+0.5)
    ystep = int(Ly/Lmin+0.5)
    zstep = int(Lz/Lmin+0.5)
    if (isioproc .eq. 1) then
      write (*,*) "Mode step ",xstep, ystep, zstep
    endif

    do kz = mode_start*zstep, nzmodes, zstep
       kzd = dfloat(kz)
       do ky = mode_start*ystep, nymodes, ystep
          kyd = dfloat(ky)
          do kx = mode_start*xstep, nxmodes, xstep
             kxd = dfloat(kx)

             kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

             if (kappa.le.kappaMax) then
                call blutilrand(rn)
                FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!     Translation angles, theta=0..2Pi and phi=0..Pi
                call blutilrand(rn)
                TAT(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                TAP(kx,ky,kz) = rn*Pi
!     Phases
                call blutilrand(rn)
                FPX(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPY(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPZ(kx,ky,kz) = rn*twicePi
                if (div_free_force.eq.1) then
                  call blutilrand(rn)
                  FPXX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPXY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPXZ(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYZ(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZZ(kx,ky,kz) = rn*twicePi
                endif

!     Amplitudes (alpha)
                call blutilrand(rn)
                thetaTmp      = rn*twicePi
                cosThetaTmp   = cos(thetaTmp)
                sinThetaTmp   = sin(thetaTmp)
                call blutilrand(rn)
                phiTmp        = rn*Pi
                cosPhiTmp     = cos(phiTmp)
                sinPhiTmp     = sin(phiTmp)

                px = cosThetaTmp * sinPhiTmp
                py = sinThetaTmp * sinPhiTmp
                pz =               cosPhiTmp

                mp2           = px*px + py*py + pz*pz
                if (kappa .lt. 0.000001) then
                  if (isioproc .eq. 1) then
                     write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                  endif
                  FAX(kx,ky,kz) = zero
                  FAY(kx,ky,kz) = zero
                  FAZ(kx,ky,kz) = zero
                else
!     Count modes that contribute
                   mode_count = mode_count + 1
!     Set amplitudes
                   if (spectrum_type.eq.1) then
                       Ekh        = one / kappa
                   else if (spectrum_type.eq.2) then
                       Ekh        = one / (kappa*kappa)
                   else
                       Ekh        = one
                   endif

                   if (div_free_force.eq.1) then
                       Ekh = Ekh / kappa
                   endif

                   if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                   endif

                   if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2
                   else
                           FAX(kx,ky,kz) = px * Ekh / mp2
                           FAY(kx,ky,kz) = py * Ekh / mp2
                           FAZ(kx,ky,kz) = pz * Ekh / mp2
                   endif

                   if (isioproc .eq. 1) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - A"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                    endif

                endif
             endif


          enddo
       enddo
    enddo

!    Now let's break symmetry, have to assume high aspect ratio in z for now
    reduced_mode_count = 0
    do kz = 1, zstep - 1
      kzd = dfloat(kz)
      do ky = mode_start, nymodes
         kyd = dfloat(ky)
         do kx = mode_start, nxmodes
            kxd = dfloat(kx)

            kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
            if (kappa.le.kappaMax) then
                call blutilrand(rn)
                FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!     Translation angles, theta=0..2Pi and phi=0..Pi
                call blutilrand(rn)
                TAT(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                TAP(kx,ky,kz) = rn*Pi
!     Phases
                call blutilrand(rn)
                FPX(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPY(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPZ(kx,ky,kz) = rn*twicePi
                if (div_free_force.eq.1) then
                   call blutilrand(rn)
                   FPXX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPXY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPXZ(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYZ(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZZ(kx,ky,kz) = rn*twicePi
                endif
!     Amplitudes (alpha)

                call blutilrand(rn)
                thetaTmp      = rn*twicePi
                cosThetaTmp   = cos(thetaTmp)
                sinThetaTmp   = sin(thetaTmp)
                call blutilrand(rn)
                phiTmp        = rn*Pi
                cosPhiTmp     = cos(phiTmp)
                sinPhiTmp     = sin(phiTmp)

                px = cosThetaTmp * sinPhiTmp
                py = sinThetaTmp * sinPhiTmp
                pz =               cosPhiTmp

                mp2           = px*px + py*py + pz*pz
                if (kappa .lt. 0.000001) then
                   if (isioproc .eq. 1) then
                      write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                   endif
                   FAX(kx,ky,kz) = zero
                   FAY(kx,ky,kz) = zero
                   FAZ(kx,ky,kz) = zero
                else
!     Count modes that contribute
                  reduced_mode_count = reduced_mode_count + 1
!     Set amplitudes
                  if (spectrum_type.eq.1) then
                    Ekh        = one / kappa
                  else if (spectrum_type.eq.2) then
                    Ekh        = one / (kappa*kappa)
                  else
                    Ekh        = one
                  endif

                  if (div_free_force.eq.1) then
                    Ekh = Ekh / kappa
                  endif

                  if (moderate_zero_modes.eq.1) then
                    if (kx.eq.0) Ekh = Ekh / two
                    if (ky.eq.0) Ekh = Ekh / two
                    if (kz.eq.0) Ekh = Ekh / two
                  endif
                  if (force_scale.gt.zero) then
                    FAX(kx,ky,kz) = forcing_epsilon * force_scale * px * Ekh / mp2
                    FAY(kx,ky,kz) = forcing_epsilon * force_scale * py * Ekh / mp2
                    FAZ(kx,ky,kz) = forcing_epsilon * force_scale * pz * Ekh / mp2
                  else
                    FAX(kx,ky,kz) = forcing_epsilon * px * Ekh / mp2
                    FAY(kx,ky,kz) = forcing_epsilon * py * Ekh / mp2
                    FAZ(kx,ky,kz) = forcing_epsilon * pz * Ekh / mp2
                  endif

                  if (isioproc .eq. 1) then
                     write (*,*) "Mode"
                     write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                     write (*,*) "Amplitudes - B"
                     write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                     write (*,*) "Frequencies"
                     write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                  endif
                endif
            endif


          enddo
       enddo
    enddo



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
      use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : domnlo, maxspec, maxspnml
      use probdata_module, only : T_in

      
      implicit none
      integer    level,nscal
      integer    lo(dim), hi(dim)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     time, delta(dim)
      REAL_T     vel(DIMV(state),dim)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

      integer i, j, k, n
      REAL_T x, y, z, Yl(maxspec), Patm

         do k = lo(3), hi(3)
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  
                  scal(i,j,k,Temp) = T_in
                  Yl(1) = 0.233d0
                  Yl(2) = 0.767d0

                  do n = 1,Nspec
                    scal(i,j,k,FirstSpec+n-1) = Yl(n)
                  end do

                  scal(i,j,k,Trac) = 0.d0
                 
                  
                  vel(i,j,k,1) = 0.d0
                  vel(i,j,k,2) = 0.d0
                  vel(i,j,k,3) = 0.0d0
                  
               end do
            end do
         end do
         
      Patm = pamb / pphys_getP1atm_MKS()

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
