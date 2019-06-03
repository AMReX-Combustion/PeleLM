
#include <MyProb_F.H>

module prob_3D_module

  implicit none

  private
  
  public :: amrex_probinit,setupbc, getZone, bcfunction, init_data_new_mech, init_data, &
            zero_visc, den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, FORT_ZVELFILL, chem_fill, press_fill, &
            FORT_MAKEFORCE 

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
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name="amrex_probinit")


    use chem_driver, only: P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
 
    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(SDIM), probhi(SDIM)

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif
#include <conp.H>

#ifdef DO_LMC_FORCE
#include <forcedata.H>
#endif

    integer i
    REAL_T area

    namelist /fortin/ vorterr, temperr, adverr, tempgrad,&
         flametracval, probtype,&
         max_temp_lev, max_vort_lev, max_trac_lev,&
         traceSpecVal,phi_in,T_in,&
         turb_scale, V_in, V_co, H2_frac,T_switch,&
         standoff, pertmag, pert_scale, nchemdiag, splitx, xfrontw,&
         splity, yfrontw, blobx, bloby, blobz, blobr,&
         blobT, Tfrontw, stTh, fuel_N2_vol_percent
    namelist /heattransin/ pamb, dpdt_factor, closed_chamber
#if defined(BL_DO_FLCT)
    namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi,&
         strmwse_dir, nCompInflow, flct_file, convVelInit
#endif
    namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
         zbase_control, pseudo_gravity, ac_hist_file,&
         corr,controlVelMax,navg_pnts,istemp

#ifdef DO_LMC_FORCE
    REAL_T  twicePi, kxd, kyd, kzd
    REAL_T  thetaTmp, phiTmp
    REAL_T  cosThetaTmp, cosPhiTmp
    REAL_T  sinThetaTmp, sinPhiTmp
    REAL_T  px, py, pz, mp2, Ekh
    integer kx, ky, kz, mode_count, reduced_mode_count
    integer xstep, ystep, zstep

    REAL_T  Lx, Ly, Lz, Lmin, rn 
    REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff

    namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz,&
         forcing_type, spectrum_type, time_offset, forcing_twice_wavelength,&
         forcing_xlength, forcing_ylength, forcing_zlength,&
         forcing_time_scale_min, forcing_time_scale_max, &
         force_scale, forcing_epsilon, blrandseed,&
         use_rho_in_forcing, do_mode_division,&
         div_free_force, moderate_zero_modes
#endif

    !
    !      Build `probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)
    character probin*(maxlen)

#if defined(BL_DO_FLCT)
    integer ierr, nCompFile
#endif


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

    !     Load domain dimensions into common (put something computable there for SDIM<3)
    do i=1,3
       domnlo(i) = 0.d0
       domnhi(i) = 0.d0
    enddo

    do i=1,SDIM
       domnlo(i) = problo(i)
       domnhi(i) = probhi(i)
    enddo

    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')

    !     Set defaults
    vorterr = 1.e20
    temperr = zero
    adverr = 1.e20
    tempgrad  = 50.0d0
    flametracval = 0.0001d0
    probtype = BL_PROB_UNDEFINED
    max_temp_lev = 0
    max_vort_lev = 0
    max_trac_lev = 100
    traceSpecVal = 1.d-10
    pamb = P1ATMMKS()
    dpdt_factor = 0.3d0
    closed_chamber = 0

    splitx = 0.5d0 * (domnhi(1) + domnlo(1))
    xfrontw = 0.05d0 * (domnhi(1) - domnlo(1))
    splity = 0.5d0 * (domnhi(2) + domnlo(2))
    yfrontw = 0.05d0 * (domnhi(2) - domnlo(2))
    blobx = 0.5d0 * (domnhi(1) + domnlo(1))
    bloby = 0.5d0 * (domnhi(2) + domnlo(2))
    blobz = 0.d0
    blobT = T_in
    Tfrontw = xfrontw
    stTh = -1.d0
    fuel_N2_vol_percent = 0.d0

#if defined(BL_DO_FLCT)
    !     add_turb = .FALSE.
    forceInflow = .FALSE.
    numInflPlanesStore = -1
    numInflPlanesStore = 30
    !
    !     Don't need to default 'nCompInflow' as it is block data'd to /3/
    !
    forceLo = .TRUE.
    forceHi = .FALSE.
    strmwse_dir = FLCT_ZVEL
    flct_file = ''
    turb_scale = 1
#endif

    zbase_control = 0.d0

    !     Note: for setup with no coflow, set Ro=Rf+wallth
    standoff = zero
    pertmag = 0.d0
    pert_scale = 1.0d0

    !     Initialize control variables
    tau_control = one
    sest = zero
    corr = one
    changeMax_control = .05
    coft_old = -one
    cfix = zero
    ac_hist_file = 'AC_History'
    h_control = -one
    nchemdiag = 1
    dV_control = zero
    tbase_control = zero
    h_control = -one
    pseudo_gravity = 0
    nchemdiag = 1
    controlVelMax = 5.d0

    if (isioproc .eq. 1) then
       write(6,*)"reading fortin"
    endif

    read(untin,fortin)

    if (isioproc .eq. 1) then
       write(6,*)"done reading fortin"
    endif

    !     Initialize control variables that depend on fortin variables
    V_in_old = V_in

    if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev

    read(untin,heattransin)

#if defined(BL_DO_FLCT)
    read(untin,flctin)
#endif

    read(untin,control)
    close(unit=untin)

#if defined(BL_DO_FLCT)
    if (forceInflow .eqv. .FALSE.) then
       forceLo = .FALSE.
       forceHi = .FALSE.
    else
       if (flct_file.ne.'') then
#define FF_UNIT 20
          ierr = 0
          write(6,*) '...initializing turbulence, reading header info'
          open(FF_UNIT, file=trim(flct_file)//'/HDR',form='formatted',status='old',iostat=ierr)
          if (ierr .ne. 0) then
             call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
          end if
          call RD_SCL_FLCTHD(FF_UNIT,nCompFile,dimFile,probSizeFile,dxFile)
          close(FF_UNIT)
       endif
    endif
    convVel = V_in
    if (probtype.ge.BL_PROB_PREMIXED_CONTROLLED_INFLOW) then
       convVel = one
    endif
#endif
    !     Set up boundary functions
    call setupbc()

    area = 1.d0
    do i=1,SDIM-1
       area = area*(domnhi(i)-domnlo(i))
    enddo
    if(istemp.eq.0)then
       scale_control = Y_bc(fuelID-1,BL_FUELPIPE) * rho_bc(BL_FUELPIPE) * area
    else
       scale_control = 1.d0
    endif
    if (isioproc .eq. 1) then
       write(6,*)" control setup", area, Y_bc(fuelID-1,BL_FUELPIPE), rho_bc(BL_FUELPIPE)
       write(6,*)" fuelpipe",BL_FUELPIPE
    endif

    if (h_control .gt. zero) then
       cfix = scale_control * h_control
    endif

#ifdef DO_LMC_FORCE
    if ( (probtype.eq.BL_PROB_PREMIXED_FIXED_INFLOW)&
         .or. (probtype.eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then

       if (isioproc .eq. 1) then
          write (*,*) "Initialising random number generator..."
       endif

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

       if (hack_lz.eq.1) then 
          Lz = Lz/two
       endif

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

                      if (isioproc.eq.1) then
                         write (*,*) "Mode"
                         write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                         write (*,*) "Amplitudes - A"
                         write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                         write (*,*) "Frequencies"
                         write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                         write (*,*) "TAT"
                         write (*,*) TAT(kx,ky,kz), TAP(kx,ky,kz)
                         write (*,*) "Amplitudes - AA"
                         write (*,*) FPXX(kx,ky,kz), FPYX(kx,ky,kz), FPZX(kx,ky,kz)
                         write (*,*) FPXY(kx,ky,kz), FPYY(kx,ky,kz), FPZY(kx,ky,kz)
                         write (*,*) FPXZ(kx,ky,kz), FPYZ(kx,ky,kz), FPZZ(kx,ky,kz)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

       !     Now break symmetry, have to assume high aspect ratio in z for now
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

                      if (isioproc.eq.1) then
                         write (*,*) "Mode"
                         write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                         write (*,*) "Amplitudes - B"
                         write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                         write (*,*) "Frequencies"
                         write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                         write (*,*) "TAT"
                         write (*,*) TAT(kx,ky,kz), TAP(kx,ky,kz)
                         write (*,*) "Amplitudes - BB"
                         write (*,*) FPXX(kx,ky,kz), FPYX(kx,ky,kz), FPZX(kx,ky,kz)
                         write (*,*) FPXY(kx,ky,kz), FPYY(kx,ky,kz), FPZY(kx,ky,kz)
                         write (*,*) FPXZ(kx,ky,kz), FPYZ(kx,ky,kz), FPZZ(kx,ky,kz)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

       if (isioproc .eq. 1) then
          write(*,*) "mode_count = ",mode_count
          write(*,*) "reduced_mode_count = ",reduced_mode_count
          if (spectrum_type.eq.1) then
             write (*,*) "Spectrum type 1"
          else if (spectrum_type.eq.2) then
             write (*,*) "Spectrum type 2"
          else
             write (*,*) "Spectrum type OTHER"
          endif
       endif
    endif
#endif

    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
#if defined(BL_DO_FLCT)
       write(6,flctin)
#endif
       write(6,control)
    end if

  end subroutine amrex_probinit

  subroutine setupbc() bind(C, name="setupbc")

    use chem_driver, only: P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    use probspec_module, only: set_Y_from_Phi
    use mod_Fvar_def, only : pamb

    implicit none


#include <cdwrk.H>
#include <conp.H>
#include <bc.H>
#include <probdata.H>

    REAL_T Patm, pmf_vals(maxspec+3), a
    REAL_T Xt(maxspec), Yt(maxspec), loc
    integer zone, n, fuelZone, airZone, region
    integer b(SDIM)
    integer num_zones_defined, len
    data  b / 1, 1, 1 /

    Patm = pamb / P1ATMMKS()
    num_zones_defined = 0
    len = len_trim(probtype)
    if ((probtype(1:len).eq.BL_PROB_PREMIXED_FREE) .or.  &
        (probtype(1:len).eq.BL_PROB_CHAMBER))  then

       !     Take fuel mixture from prob data

       region = BL_INTERIOR
       zone = getZone(domnlo(1), domnlo(2), domnlo(3))
       num_zones_defined = 1

       if (phi_in.gt.zero) then

          call set_Y_from_Phi(phi_in,Yt)
          do n=1,Nspec
             Y_bc(n-1,zone) = Yt(n)
          end do
          T_bc(zone) = T_in
          u_bc(zone) = zero
          v_bc(zone) = zero
          w_bc(zone) = V_in

       else 

          !     Take fuel mixture from pmf file
          loc = (domnlo(2)-standoff)*100.d0
          call pmf(loc,loc,pmf_vals,n)
          if (n.ne.Nspec+3) then
             call bl_pd_abort('INITDATA: n(pmf) .ne. Nspec+3')
          endif

          do n = 1,Nspec
             Xt(n) = pmf_vals(3+n)
          end do

          CALL CKXTY (Xt, Yt)

          do n=1,Nspec
             Y_bc(n-1,zone) = Yt(n)
          end do
          T_bc(zone) = pmf_vals(1)
          u_bc(zone) = zero
          v_bc(zone) = zero
          if (V_in .lt. 0) then
             w_bc(zone) = pmf_vals(2)*1.d-2
          else
             w_bc(zone) = V_in
          endif

       endif

    else if (probtype(1:len).eq.BL_PROB_DIFFUSION) then

       !     A diffusion flame
       fuelZone = getZone(domnlo(1), domnlo(2), domnlo(3))
       airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))
       num_zones_defined = 2

       !     Fuel
       do n = 1,Nspec
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
       v_bc(fuelZone) = V_in
       w_bc(fuelZone) = 0.d0

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

       T_bc(airZone) = T_in
       u_bc(airZone) = 0.d0
       v_bc(airZone) = V_in
       w_bc(airZone) = 0.d0

    else

       call bl_pd_abort()

    endif

    do zone=1,num_zones_defined
       !     Set density and hmix consistent with data

       call RHOfromPTY(b, b, &
            rho_bc(zone), DIMARG(b), DIMARG(b),&
            T_bc(zone),   DIMARG(b), DIMARG(b),&
            Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)
       call HMIXfromTY(b, b, &
            h_bc(zone),   DIMARG(b), DIMARG(b),&
            T_bc(zone),   DIMARG(b), DIMARG(b),&
            Y_bc(0,zone), DIMARG(b), DIMARG(b))
    enddo
    bcinit = .true.
  end subroutine setupbc

  ! ::: -----------------------------------------------------------

  integer function getZone(x, y, z) bind(C, name="getZone")
    implicit none
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
    REAL_T x, y, z
    integer len

#define BL_FUELPIPE 1
    !#define BL_COFLOW   2
    !#define BL_STICK    3
    !#define BL_WALL     4
    !#define BL_AMBIENT  5
    !#define BL_VOLUME   6
    !#define BL_PIPEEND  7

    getZone = BL_VOLUME
    len     = len_trim(probtype)

         if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
             .or. (probtype(1:len).eq.BL_PROB_CHAMBER) &
             .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
    !if  (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) then 

       getZone = BL_FUELPIPE

    else

       call bl_pd_abort('Unrecognized probtype')

    endif
  end function getZone

  ! ::: -----------------------------------------------------------

  subroutine bcfunction(RegionID,x,y,z,time,u,v,w,rho,Yl,T,h,dx,getuv)&
       bind(C, name="bcfunction")

    implicit none

    integer RegionID
    REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM)
    logical getuv

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer n, zone, len
    REAL_T eta, xmid, etamax

    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    len = len_trim(probtype)

    if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
         .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) &
         .or. (probtype(1:len).eq.BL_PROB_CHAMBER) &
         .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW)) then

     !  if (RegionID .eq. BL_ZLO) then

          zone = getZone(x,y,z)
          rho = rho_bc(zone)
          do n = 0, Nspec-1
             Yl(n) = Y_bc(n,zone)
          end do
          T = T_bc(zone)
          h = h_bc(zone)

          if (getuv .eqv. .TRUE.) then

             u = zero
             v = zero
             if  ((probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW)&
                  .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) ) then               
                w =  V_in + (time-tbase_control)*dV_control
             else if ((probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
                .or.  (probtype(1:len).eq.BL_PROB_CHAMBER)) then
                w = w_bc(zone)
             endif
          endif

     ! else
     !    write(6,*) 'No bcfunction instruction for RegionID = ', RegionID
     !    call bl_pd_abort(' ')
     ! endif

    elseif (probtype(1:len).eq.BL_PROB_DIFFUSION) then

       if (RegionID .eq. BL_YLO) then

          zone = getZone(x,y,z)
          rho = rho_bc(zone)
          do n = 0, Nspec-1
             Yl(n) = Y_bc(n,zone)
          end do
          T = T_bc(zone)
          h = h_bc(zone)

          if (getuv .eqv. .TRUE.) then
             u = zero
             v = zero
             w = V_in

             if (stTh.ge.0) then
                if (ABS(x - splitx).lt.0.5*stTh) then
                   v = 0.d0
                else
                   if (x .lt. splitx) then
                      etamax = splitx - 0.5*stTh
                      eta = (x - domnlo(1)) / etamax
                      v = V_in * (1.d0 - eta**2)
                   else
                      etamax = domnhi(1) - splitx - 0.5*stTh
                      eta = (domnhi(1) - x) / etamax
                      v = V_co * (1.d0 - eta**2)
                   endif
                endif
             endif

          endif

       else
          write(6,*) 'No bcfunction instruction for RegionID = ', RegionID
          call bl_pd_abort(' ')
       endif
    else
       write(6,*) 'No boundary condition for probtype = ', probtype(1:len)
       write(6,*) 'Available: '
       write(6,*) '            ',BL_PROB_PREMIXED_FIXED_INFLOW
       write(6,*) '            ',BL_PROB_PREMIXED_CONTROLLED_INFLOW
       write(6,*) '            ',BL_PROB_DIFFUSION
       call bl_pd_abort(' ')
    endif
  end subroutine bcfunction

  ! ::: -----------------------------------------------------------

  subroutine init_data_new_mech(level,time,lo,hi,nscal,&
       vel,scal,DIMS(state),press,DIMS(press),&
       delta,xlo,xhi)

    use chem_driver, only: P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac

    implicit none

    integer  level, nscal
    integer  lo(SDIM), hi(SDIM)
    integer  DIMDEC(state)
    integer  DIMDEC(press)
    REAL_T   xlo(SDIM), xhi(SDIM)
    REAL_T   time, delta(SDIM)
    REAL_T   vel(DIMV(state),SDIM)
    REAL_T   scal(DIMV(state),nscal)
    REAL_T   press(DIMV(press))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k, n
    REAL_T Patm

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             scal(i,j,k,Trac) = zero
          end do
       end do
    end do

    Patm = pamb / P1ATMMKS()
    call RHOfromPTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state),&
         Patm)
    call HMIXfromTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
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

  end subroutine init_data_new_mech

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
  subroutine init_data(level,time,lo,hi,nscal,&
       vel,scal,DIMS(state),press,DIMS(press),&
       delta,xlo,xhi) bind(C, name="init_data")

    use chem_driver, only: P1ATMMKS, get_spec_name
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac

    implicit none
    integer    level, nscal
    integer    lo(SDIM), hi(SDIM)
    integer    DIMDEC(state)
    integer    DIMDEC(press)
    REAL_T     xlo(SDIM), xhi(SDIM)
    REAL_T     time, delta(SDIM)
    REAL_T     vel(DIMV(state),SDIM)
    REAL_T    scal(DIMV(state),nscal)
    REAL_T   press(DIMV(press))
    integer tmpi, nPMF

#include <cdwrk.H>
#include <conp.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k, n, airZone, fuelZone, zone
    integer iO2,iH2,iCH4,len
    character*(maxspnml) name
    REAL_T x, y, z, r, Yl(maxspec), Xl(maxspec), Patm
    REAL_T Xlin(maxspec),alpha,beta,gamma,delt,factor
    REAL_T pmf_vals(maxspec+3), z1, z2, dx, Ly
    REAL_T pert,Lx,eta,u,v,w,rho,T,h

    if (iN2.lt.1 .or. iN2.gt.Nspec) then
       call bl_pd_abort()
    endif

    len = len_trim(probtype)

    if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) &
         .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW)&
         .or. (probtype(1:len).eq.BL_PROB_CHAMBER)&
         .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) ) then

       if(probtype(1:len).eq.BL_PROB_PREMIXED_FREE)then


          do n=1,Nspec

             call get_spec_name(name,n)
             if (name .eq. 'O2' ) iO2 = n
             if (name .eq. 'H2' ) iH2 = n
             if (name .eq. 'CH4' ) iCH4 = n

          enddo

       endif

       do k = lo(3), hi(3)
          z = (float(k)+.5d0)*delta(3)+domnlo(3)
          do j = lo(2), hi(2)
             y = (float(j)+.5d0)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5d0)*delta(1)+domnlo(1)

                pert = 0.d0
                if (pertmag .gt. 0.d0) then
                   Lx = (domnhi(1) - domnlo(1)) / pert_scale
                   Ly = (domnhi(2) - domnlo(2)) / pert_scale
                   pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)             * sin(2*Pi*5*y/Ly)&
                        + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)   * sin(2*Pi*4*(y-.0053765)/Ly)&
                        + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) * sin(2*Pi*3*(y-.02137)/Ly)&
                        + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)     * sin(2*Pi*6*(y-.018)/Ly)&
                        + .982 * sin(2*Pi*5*(x-.014234)/Lx) )

                endif

                z1 = (z - blobz - standoff - 0.5d0*delta(2) + pert )*100.d0
                z2 = (z - blobz - standoff + 0.5d0*delta(2) + pert )*100.d0 

                call pmf(z1,z2,pmf_vals,nPMF)               
                if (nPMF.ne.Nspec+3) then
                   call bl_abort('INITDATA: n .ne. Nspec+3')
                endif

                scal(i,j,k,Temp) = pmf_vals(1)
                do n = 1,Nspec
                   Xl(n) = pmf_vals(3+n)
                end do

                if(probtype(1:len).eq.BL_PROB_PREMIXED_FREE)then

                   iO2 = -1; iH2 = -1; iCH4 = -1

                   do n = 1,Nspec
                      Xlin(n) = 0.d0
                   end do

                   alpha = H2_frac
                   beta = 1.d0 - H2_frac
                   gamma = (0.5d0*alpha + 2.d0*beta) / phi_in
                   delt = gamma*.79d0/.21d0
                   factor = alpha+beta +gamma+delt

                   if (iH2  > 0) Xlin(iH2) = alpha / factor
                   if (iCH4 > 0) Xlin(iCH4) = beta / factor
                   if (iO2  > 0) Xlin(iO2) = gamma / factor
                   if (iN2  > 0) Xlin(iN2) = delt / factor

                   !                 blend here if needed
                   if( scal(i,j,k,Temp) .lt. T_switch)then

                      scal(i,j,k,TEMP) = T_in
                      do n = 1,Nspec
                         Xl(n) = Xlin(n)
                      end do
                   endif
                endif

                CALL CKXTY (Xl, Yl)

                do n = 1,Nspec
                   scal(i,j,k,FirstSpec+n-1) = Yl(n)
                end do

                scal(i,j,k,Trac) = 0.d0

                vel(i,j,k,1) = 0.d0
                vel(i,j,k,2) = 0.d0
                vel(i,j,k,3) = pmf_vals(2)*1.d-2

             end do
          end do
       end do

    else if ((probtype(1:len).eq.BL_PROB_DIFFUSION)) then

       fuelZone = getZone(domnlo(1), domnlo(2), domnlo(3))
       airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))

       do k = lo(3), hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                if (blobr.lt.0) then

                   eta = 0.5d0*(1.d0 - TANH(2.d0*(x-splitx)/xfrontw))

                   do n=1,Nspec
                      scal(i,j,k,FirstSpec-1+n) = Y_bc(n-1,airZone)*(1.d0-eta) + eta*Y_bc(n-1,fuelZone)
                   enddo
                   scal(i,j,k,Temp) = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)

                   eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
                   do n=1,Nspec
                      scal(i,j,k,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta&
                           + (1.d0-eta)*scal(i,j,k,FirstSpec-1+n)
                   enddo
                   scal(i,j,k,Temp) = blobT*eta + (1.d0-eta)*scal(i,j,k,Temp)

                   vel(i,j,k,1) = u_bc(airZone)*eta + (1.d0-eta)*u_bc(fuelZone)
                   vel(i,j,k,2) = v_bc(airZone)*eta + (1.d0-eta)*v_bc(fuelZone)
                   vel(i,j,k,3) = w_bc(airZone)*eta + (1.d0-eta)*w_bc(fuelZone)
                   scal(i,j,k,Trac) = 0.d0

                else

                   eta = 0.5d0*(1.d0 - TANH(2.d0*(x-splitx)/xfrontw))

                   do n=1,Nspec
                      scal(i,j,k,FirstSpec-1+n) = Y_bc(n-1,airZone)*(1.d0-eta) + eta*Y_bc(n-1,fuelZone)
                   enddo
                   scal(i,j,k,Temp) = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)

                   !     Superimpose blob of hot air
                   r = SQRT((x-blobx)**2 + (y-bloby)**2)
                   eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/Tfrontw))
                   do n=1,Nspec
                      scal(i,j,k,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta&
                           + (1.d0-eta)*scal(i,j,k,FirstSpec-1+n)
                   enddo
                   scal(i,j,k,Temp) = blobT*eta + (1.d0-eta)*scal(i,j,k,Temp)

                   vel(i,j,k,1) = u_bc(airZone)*eta + (1.d0-eta)*u_bc(fuelZone)
                   vel(i,j,k,2) = v_bc(airZone)*eta + (1.d0-eta)*v_bc(fuelZone)
                   vel(i,j,k,3) = w_bc(airZone)*eta + (1.d0-eta)*w_bc(fuelZone)
                   scal(i,j,k,Trac) = 0.d0

                   if (stTh.gt. 0.d0) then
                      zone = getZone(x,y,z)
                      call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                      vel(i,j,k,1) = u
                      vel(i,j,k,2) = v
                      vel(i,j,k,3) = w
                   endif
                endif

             enddo
          enddo
       enddo
    endif

    Patm = pamb / P1ATMMKS()

    call RHOfromPTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state),&
         Patm)

    call HMIXfromTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
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

  ! ::: -----------------------------------------------------------
  ! ::: This routine will zero out diffusivity on portions of the
  ! ::: boundary that are inflow, allowing that a "wall" block
  ! ::: the complement aperture
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: diff      <=> diffusivity on edges
  ! ::: DIMS(diff) => index extent of diff array
  ! ::: lo,hi      => region of interest, edge-based
  ! ::: domlo,hi   => index extent of problem domain, edge-based
  ! ::: dx         => cell spacing
  ! ::: problo     => phys loc of lower left corner of prob domain
  ! ::: bc         => boundary condition flag (on orient)
  ! :::                   in BC_TYPES::physicalBndryTypes
  ! ::: idir       => which face, 0=x, 1=y
  ! ::: isrz       => 1 if problem is r-z
  ! ::: id         => index of state, 0=u
  ! ::: ncomp      => components to modify
  ! ::: 
  ! ::: -----------------------------------------------------------
  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi,&
       dx,problo,bc,idir,isrz,id,ncomp) bind(C, name="zero_visc")

    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec

    implicit none

    integer DIMDEC(diff)
    integer lo(SDIM), hi(SDIM)
    integer domlo(SDIM), domhi(SDIM)
    integer bc(2*SDIM)
    integer idir, isrz, id, ncomp
    REAL_T  diff(DIMV(diff),*)
    REAL_T  dx(SDIM)
    REAL_T  problo(SDIM)

#include <probdata.H>
#include <cdwrk.H>

    integer i, j, k, n, Tid, RHid, YSid, YEid, ys, ye
    integer len
    logical do_T, do_RH, do_Y
    REAL_T xl, xr, xh, yb, yt, yh, z

    len = len_trim(probtype)

    if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW)&
         .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
       Tid  = Temp      - id + SDIM
       RHid = RhoH      - id + SDIM
       YSid = FirstSpec - id + SDIM
       YEid = LastSpec  - id + SDIM

       do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
       do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
       ys = MAX(YSid,1)
       ye = MIN(YEid,ncomp)
       do_Y = (ye - ys + 1) .GE. 1
       !     
       !     Do species, Temp, rhoH
       !     
       if ((idir.EQ.2) .AND. (lo(3) .LE. domlo(3))&
            .AND. (do_T .OR. do_RH .OR. do_Y) ) then

          k = lo(3)
          z = float(k)*dx(3)+domnlo(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                xl = float(i)*dx(1)+domnlo(1) 
                xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
                xh = 0.5d0*(xl+xr)
                yb = float(j)*dx(2)+domnlo(2) 
                yt = (float(j)+1.d0)*dx(2)+domnlo(2) 
                yh = 0.5d0*(yb+yt)

                if ( (getZone(xl,yb,z).eq.BL_STICK) .OR.&
                     (getZone(xh,yb,z).eq.BL_STICK) .OR.&
                     (getZone(xr,yb,z).eq.BL_STICK) .OR.&
                     (getZone(xl,yh,z).eq.BL_STICK) .OR.&
                     (getZone(xh,yh,z).eq.BL_STICK) .OR.&
                     (getZone(xr,yh,z).eq.BL_STICK) .OR.&
                     (getZone(xl,yt,z).eq.BL_STICK) .OR.&
                     (getZone(xh,yt,z).eq.BL_STICK) .OR.&
                     (getZone(xr,yt,z).eq.BL_STICK) .OR.&
                     (getZone(xl,yb,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xh,yb,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xr,yb,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xl,yh,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xh,yh,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xr,yh,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xl,yt,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xh,yt,z).eq.BL_PIPEEND) .OR.&
                     (getZone(xr,yt,z).eq.BL_PIPEEND) ) then

                   !                    if (do_T)  diff(i,j,k,Tid ) = 0.d0
                   !                    if (do_RH) diff(i,j,k,RHid) = 0.d0
                   if (do_Y) then
                      do n=ys,ye
                         diff(i,j,k,n) = 0.d0
                      enddo
                   endif

                endif
             end do
          end do
       endif
    end if
  end subroutine zero_visc

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data and that all non-interior cells have
  ! :::         have been filled with a large real number.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: den      <=  density array
  ! ::: DIMS(den) => index extent of den array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: dx        => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of den array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine den_fill(den,DIMS(den),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="den_fill")
    implicit none

    integer DIMDEC(den), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  den(DIMV(den))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(den)
    lo(2) = ARG_L2(den)
    lo(3) = ARG_L3(den)
    hi(1) = ARG_H1(den)
    hi(2) = ARG_H2(den)
    hi(3) = ARG_H3(den)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

    call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif

  end subroutine den_fill

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data and that all non-interior cells have
  ! :::        have been filled with a large real number.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: adv      <=  advected quantity array
  ! ::: DIMS(adv) => index extent of adv array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: dx        => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of adv array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine adv_fill (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc) bind(C, name="adv_fill")

    implicit none

    integer    DIMDEC(adv)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     delta(SDIM), xlo(SDIM), time
    REAL_T     adv(DIMV(adv))
    integer    bc(SDIM,2)

    integer    i,j,k
    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(adv)
    lo(2) = ARG_L2(adv)
    lo(3) = ARG_L3(adv)
    hi(1) = ARG_H1(adv)
    hi(2) = ARG_H2(adv)
    hi(3) = ARG_H3(adv)

    call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          do j = lo(2),hi(2)
             do i = lo(1), hi(1)
                adv(i,j,k) = 0.0d0
             enddo
          enddo
       enddo
    endif

  end subroutine adv_fill


  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.
  ! :::
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: temp     <=  temperature array
  ! ::: lo,hi     => index extent of adv array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::              corner of temperature array
  ! ::: time      => problem evolution time
  ! ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine temp_fill (temp,DIMS(temp),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="temp_fill")

    implicit none

    integer DIMDEC(temp), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  temp(DIMV(temp))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(temp)
    lo(2) = ARG_L2(temp)
    lo(3) = ARG_L3(temp)
    hi(1) = ARG_H1(temp)
    hi(2) = ARG_H2(temp)
    hi(3) = ARG_H3(temp)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

    call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif

  end subroutine temp_fill

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.
  ! :::
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: rhoh      <=  rho*h array
  ! ::: lo,hi     => index extent of adv array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::              corner of temperature array
  ! ::: time      => problem evolution time
  ! ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine rhoh_fill (rhoh,DIMS(rhoh),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="rhoh_fill")

    implicit none

    integer DIMDEC(rhoh), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  rhoh(DIMV(rhoh))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(rhoh)
    lo(2) = ARG_L2(rhoh)
    lo(3) = ARG_L3(rhoh)
    hi(1) = ARG_H1(rhoh)
    hi(2) = ARG_H2(rhoh)
    hi(3) = ARG_H3(rhoh)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

    call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

  end subroutine rhoh_fill

  !
  ! Fill x  y velocity at once.
  !
  subroutine vel_fill (vel,DIMS(vel),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="vel_fill")

    implicit none
    integer DIMDEC(vel), bc(SDIM,2,SDIM)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  vel(DIMV(vel),SDIM)

    call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1),&
         DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

    call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2),&
         DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

    call FORT_ZVELFILL (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3),&
         DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,3))

  end subroutine vel_fill

  !
  ! Fill all chem species at once
  !
  subroutine all_chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="all_chem_fill")

    implicit none
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer DIMDEC(rhoY), bc(SDIM,2,Nspec)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  rhoY(DIMV(rhoY),Nspec)

    integer n

    do n=1,Nspec
       call chem_fill (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),ARG_L3(rhoY),n),&
            DIMS(rhoY),domlo,domhi,delta,xlo,time,bc(1,1,n),n-1)
    enddo
  end subroutine all_chem_fill

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: xvel     <=  x velocity array
  ! ::: lo,hi     => index extent of xvel array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of rho array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="FORT_XVELFILL")

    implicit none

    integer DIMDEC(xvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  xvel(DIMV(xvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif      

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
    REAL_T  cs,sg,r,scaler,scalet
#if defined(BL_DO_FLCT)
    integer DIMDEC(uflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T uflct(:,:,:), vflct(:,:,:)
    allocatable uflct, vflct
#endif
    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(xvel)
    lo(2) = ARG_L2(xvel)
    lo(3) = ARG_L3(xvel)
    hi(1) = ARG_H1(xvel)
    hi(2) = ARG_H2(xvel)
    hi(3) = ARG_H3(xvel)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3)) 

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       do i = 1, SDIM
          loFlctArray(i) = lo(i)
          hiFlctArray(i) = hi(i)
       enddo
       loFlctArray(strmwse_dir) = 1
       hiFlctArray(strmwse_dir) = 1
       call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
       allocate(uflct(DIMV(uflct)))
       call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, time,
       $                  bc, domnlo, domnhi)
       allocate(vflct(DIMV(uflct)))
       call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, time,
       $                  bc, domnlo, domnhi)
    endif
#endif

    call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 1) then
                   xvel(i,j,k) = u + uflct(1,j,k)*turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 1) then
                   xvel(i,j,k) = u + uflct(1,j,k)*turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 2) then
                   xvel(i,j,k) = u + uflct(i,1,k)*turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 2) then
                   xvel(i,j,k) = u + uflct(i,1,k)*turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 3) then
                   xvel(i,j,k) = u + turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 3) then
                   xvel(i,j,k) = u + uflct(i,j,1)*turb_scale
                else
                   xvel(i,j,k) = u
                endif
#else
                xvel(i,j,k) = u
#endif
             enddo
          enddo
       enddo
    endif

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       deallocate(uflct)
       deallocate(vflct)
    endif
#endif

  end subroutine FORT_XVELFILL

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: yvel     <=  y velocity array
  ! ::: lo,hi     => index extent of yvel array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of rho array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="FORT_YVELFILL")

    implicit none

    integer DIMDEC(yvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  yvel(DIMV(yvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
    REAL_T  cs,sg,r,scaler,scalet

#if defined(BL_DO_FLCT)
    integer DIMDEC(vflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T uflct(:,:,:), vflct(:,:,:)
    allocatable uflct, vflct
#endif

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(yvel)
    lo(2) = ARG_L2(yvel)
    lo(3) = ARG_L3(yvel)
    hi(1) = ARG_H1(yvel)
    hi(2) = ARG_H2(yvel)
    hi(3) = ARG_H3(yvel)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       do i = 1, SDIM
          loFlctArray(i) = lo(i)
          hiFlctArray(i) = hi(i)
       enddo
       loFlctArray(strmwse_dir) = 1
       hiFlctArray(strmwse_dir) = 1
       call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
       allocate(uflct(DIMV(vflct)))
       call INFL_FILL(FLCT_XVEL, DIMS(vflct), uflct, xlo, delta, time,
       $                  bc, domnlo, domnhi)
       allocate(vflct(DIMV(vflct)))
       call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, time,
       $                  bc, domnlo, domnhi)
    endif
#endif

    call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 1) then
                   yvel(i,j,k) = v + vflct(1,j,k)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 1) then
                   yvel(i,j,k) = v + vflct(1,j,k)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 2) then
                   yvel(i,j,k) = v + vflct(i,1,k)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 2) then
                   yvel(i,j,k) = v + vflct(i,1,k)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 3) then
                   yvel(i,j,k) = v + vflct(i,j,1)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 3) then
                   yvel(i,j,k) = v + vflct(i,j,1)*turb_scale
                else
                   yvel(i,j,k) = v
                endif
#else
                yvel(i,j,k) = v
#endif
             enddo
          enddo
       enddo
    endif

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       deallocate(uflct)
       deallocate(vflct)
    endif
#endif

  end subroutine FORT_YVELFILL

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: zvel     <=  z velocity array
  ! ::: lo,hi     => index extent of zvel array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of rho array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: -----------------------------------------------------------

  subroutine FORT_ZVELFILL (zvel,DIMS(zvel),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="FORT_ZVELFILL")

    implicit none

    integer DIMDEC(zvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  zvel(DIMV(zvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h
#if defined(BL_DO_FLCT)
    integer DIMDEC(wflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T wflct(:,:,:)
    allocatable wflct
#endif

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(zvel)
    lo(2) = ARG_L2(zvel)
    lo(3) = ARG_L3(zvel)
    hi(1) = ARG_H1(zvel)
    hi(2) = ARG_H2(zvel)
    hi(3) = ARG_H3(zvel)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       do i = 1, SDIM
          loFlctArray(i) = lo(i)
          hiFlctArray(i) = hi(i)
       enddo
       loFlctArray(strmwse_dir) = 1
       hiFlctArray(strmwse_dir) = 1
       call SET_ARGS(DIMS(wflct), loFlctArray, hiFlctArray)
       allocate(wflct(DIMV(wflct)))
       call INFL_FILL(FLCT_ZVEL, DIMS(wflct), wflct, xlo, delta, time,
       $                  bc, domnlo, domnhi)
    endif
#endif

    call filcc (zvel,DIMS(zvel),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 1) then
                   zvel(i,j,k) = w + wflct(1,j,k)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 1) then
                   zvel(i,j,k) = w + wflct(1,j,k)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 2) then
                   zvel(i,j,k) = w + wflct(i,1,k)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 2) then
                   zvel(i,j,k) = w + wflct(i,1,k)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceLo .and. strmwse_dir .eq. 3) then
                   zvel(i,j,k) = w + wflct(i,j,1)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
#if defined(BL_DO_FLCT)
                if (forceHi .and. strmwse_dir .eq. 3) then
                   zvel(i,j,k) = w + wflct(i,j,1)*turb_scale
                else
                   zvel(i,j,k) = w
                endif
#else
                zvel(i,j,k) = w
#endif
             enddo
          enddo
       enddo
    endif

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       deallocate(wflct)
    endif
#endif

  end subroutine FORT_ZVELFILL

  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.
  ! :::
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: rhoY      <= rho*Y (Y=mass fraction) array
  ! ::: lo,hi     => index extent of adv array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::              corner of temperature array
  ! ::: time      => problem evolution time
  ! ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
  ! ::: stateID   => id index of state being filled
  ! ::: -----------------------------------------------------------

  subroutine chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta,&
       xlo,time,bc,id) bind(C, name="chem_fill")

    implicit none

    integer DIMDEC(rhoY), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM), id
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  rhoY(DIMV(rhoY))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    !      print *, 'FORT_CHEMFILL: ', domlo,domhi,delta,xlo,time

    lo(1) = ARG_L1(rhoY)
    lo(2) = ARG_L2(rhoY)
    lo(3) = ARG_L3(rhoY)
    hi(1) = ARG_H1(rhoY)
    hi(2) = ARG_H2(rhoY)
    hi(3) = ARG_H3(rhoY)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))

    call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(BL_XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

    if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(BL_ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif

  end subroutine chem_fill


  ! ::: -----------------------------------------------------------
  ! ::: This routine is called during a filpatch operation when
  ! ::: the patch to be filled falls outside the interior
  ! ::: of the problem domain.  You are requested to supply the
  ! ::: data outside the problem interior in such a way that the
  ! ::: data is consistant with the types of the boundary conditions
  ! ::: you specified in the C++ code.  
  ! ::: 
  ! ::: NOTE:  you can assume all interior cells have been filled
  ! :::        with valid data.
  ! ::: 
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: p        <=  pressure array
  ! ::: DIMS(p)   => index extent of p array
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: dx        => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::	           corner of rho array
  ! ::: time      => problem evolution time
  ! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
  ! ::: -----------------------------------------------------------


  subroutine press_fill (p,DIMS(p),domlo,domhi,dx,xlo,time,bc) &
       bind(C, name="press_fill")

    implicit none

    integer    DIMDEC(p)
    integer    domlo(SDIM), domhi(SDIM)
    REAL_T     dx(SDIM), xlo(SDIM), time
    REAL_T     p(DIMV(p))
    integer    bc(SDIM,2)

    integer    i, j, k
    integer    ilo, ihi, jlo, jhi, klo, khi
    logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi
    logical    per_xlo, per_xhi, per_ylo, per_yhi, per_zlo, per_zhi

    fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
    per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
    fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
    per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
    fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
    per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
    fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
    per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)
    fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
    per_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .eq. INT_DIR)
    fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)
    per_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .eq. INT_DIR)

    ilo = max(ARG_L1(p),domlo(1))
    jlo = max(ARG_L2(p),domlo(2))
    klo = max(ARG_L3(p),domlo(3))
    ihi = min(ARG_H1(p),domhi(1))
    jhi = min(ARG_H2(p),domhi(2))
    khi = min(ARG_H3(p),domhi(3))

    !***************
    !  SETTING BL_XLO
    !***************

    if (fix_xlo) then
       do i = ARG_L1(p), domlo(1)-1
          do k = klo, khi
             do j = jlo,jhi
                p(i,j,k) = p(ilo,j,k)
             end do
          end do
       end do

       if (fix_ylo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = klo, khi
                   p(i,j,k) = p(ilo,jlo,k)
                end do
             end do
          end do

          if (fix_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jlo,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jlo,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jlo,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jlo,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_yhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = klo, khi
                   p(i,j,k) = p(ilo,jhi,k)
                end do
             end do
          end do
          if (fix_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jhi,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jhi,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jhi,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jhi,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo, jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,klo)
                end do
             end do
          end do
          if (per_ylo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,j,klo)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,j,klo)
                   end do
                end do
             end do
          end if

       end if

       if (fix_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo, jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,khi)
                end do
             end do
          end do
          if (per_ylo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,j,khi)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,j,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_ylo) then
          do i = ARG_L1(p), domlo(1)-1
             do k = klo,khi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do i = ARG_L1(p), domlo(1)-1
             do k = klo,khi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo,jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo,jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_XHI
    !*****************************************************************************

    if (fix_xhi) then
       do i = domhi(1)+1, ARG_H1(p)
          do k = klo, khi
             do j = jlo,jhi
                p(i,j,k) = p(ihi,j,k)
             end do
          end do
       end do

       if (fix_ylo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = klo, khi
                   p(i,j,k) = p(ihi,jlo,k)
                end do
             end do
          end do

          if (fix_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jlo,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jlo,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jlo,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jlo,k)
                   end do
                end do
             end do
          end if
       end if
       if (fix_yhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = klo, khi
                   p(i,j,k) = p(ihi,jhi,k)
                end do
             end do
          end do
          if (fix_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jhi,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jhi,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jhi,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jhi,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo, jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,klo)
                end do
             end do
          end do
          if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,j,klo)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,j,klo)
                   end do
                end do
             end do
          end if

       end if

       if (fix_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo, jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,khi)
                end do
             end do
          end do
          if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,j,khi)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,j,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_ylo) then
          do i = domhi(1)+1, ARG_H1(p)
             do k = klo,khi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do k = klo,khi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo,jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo,jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if


       if (per_ylo .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_YLO
    !*****************************************************************************

    if (fix_ylo) then
       do j = ARG_L2(p), domlo(2)-1
          do k = klo, khi
             do i = ilo, ihi
                p(i,j,k) = p(i,jlo,k)
             end do
          end do
       end do

       if (fix_zlo) then
          do j = ARG_L2(p), domlo(2)-1
             do k = ARG_L3(p), domlo(3)-1
                do i = ilo, ihi
                   p(i,j,k) = p(i,jlo,klo)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jlo,klo)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jlo,klo)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zhi) then
          do j = ARG_L2(p), domlo(2)-1
             do k = domhi(3)+1, ARG_H3(p)
                do i = ilo, ihi
                   p(i,j,k) = p(i,jlo,khi)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jlo,khi)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jlo,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_xlo) then
          do j = ARG_L2(p), domlo(2)-1
             do k = klo,khi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do j = ARG_L2(p), domlo(2)-1
             do k = klo,khi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do j = ARG_L2(p), domlo(2)-1
             do i = ilo,ihi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do j = ARG_L2(p), domlo(2)-1
             do i = ilo,ihi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if


       if (per_xlo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_YHI
    !*****************************************************************************

    if (fix_yhi) then
       do j = domhi(2)+1, ARG_H2(p)
          do k = klo, khi
             do i = ilo, ihi
                p(i,j,k) = p(i,jhi,k)
             end do
          end do
       end do

       if (fix_zlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = ARG_L3(p), domlo(3)-1
                do i = ilo, ihi
                   p(i,j,k) = p(i,jhi,klo)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jhi,klo)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jhi,klo)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = domhi(3)+1, ARG_H3(p)
                do i = ilo, ihi
                   p(i,j,k) = p(i,jhi,khi)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jhi,khi)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jhi,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_xlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = klo,khi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = klo,khi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do i = ilo,ihi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do i = ilo,ihi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_ZLO
    !*****************************************************************************

    if (fix_zlo) then
       do k = ARG_L3(p), domlo(3)-1
          do j = jlo, jhi
             do i = ilo, ihi
                p(i,j,k) = p(i,j,klo)
             end do
          end do
       end do

       if (per_xlo) then
          do k = ARG_L3(p), domlo(3)-1
             do j = jlo,jhi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do k = ARG_L3(p), domlo(3)-1
             do j = jlo,jhi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ilo,ihi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ilo,ihi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_ZHI
    !*****************************************************************************

    if (fix_zhi) then
       do k = domhi(3)+1, ARG_H3(p)
          do j = jlo, jhi
             do i = ilo, ihi
                p(i,j,k) = p(i,j,khi)
             end do
          end do
       end do

       if (per_xlo) then
          do k = domhi(3)+1, ARG_H3(p)
             do j = jlo,jhi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do j = jlo,jhi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ilo,ihi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ilo,ihi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if


       if (per_xlo .and. per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

    end if

  end subroutine press_fill

  !
  !
  ! ::: -----------------------------------------------------------
  !
  !     This routine add the forcing terms to the momentum equation
  !
  subroutine FORT_MAKEFORCE(time,force,rho,&
       DIMS(istate),DIMS(state),&
       dx,xlo,xhi,gravity,scomp,ncomp)

    implicit none

    integer    DIMDEC(state)
    integer    DIMDEC(istate)
    integer    scomp, ncomp
    REAL_T     time, dx(SDIM)
    REAL_T     xlo(SDIM), xhi(SDIM)
    REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
    REAL_T     rho    (DIMV(state))
    REAL_T     gravity

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>

    integer i, j, k, n
    integer ilo, jlo, klo
    integer ihi, jhi, khi
    integer a2, a3, a4, a5
    REAL_T  x, y, z
    REAL_T  hx, hy, hz
    REAL_T  sga, cga
    integer isioproc
    integer nXvel, nYvel, nZvel, nRho, nTrac

    call bl_pd_is_ioproc(isioproc)

    if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
       write(*,*) "pseudo_gravity::dV_control = ",dV_control
    endif

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    ilo = istate_l1
    jlo = istate_l2
    klo = istate_l3
    ihi = istate_h1
    jhi = istate_h2
    khi = istate_h3

    !     Assumes components are in the following order
    nXvel = 1
    nYvel = 2
    nZvel = 3
    nRho  = 4
    nTrac = 5

    if (scomp.eq.0) then
       if (abs(gravity).gt.0.0001) then
          do k = klo, khi
             do j = jlo, jhi
                do i = ilo, ihi
                   force(i,j,k,nXvel) = zero
                   force(i,j,k,nYvel) = gravity*rho(i,j,k)
                   force(i,j,k,nZvel) = zero
                enddo
             enddo
          enddo
          !     else to zero
       else
          do k = klo, khi
             do j = jlo, jhi
                do i = ilo, ihi
                   force(i,j,k,nXvel) = zero
                   force(i,j,k,nYvel) = zero
                   force(i,j,k,nZvel) = zero
                enddo
             enddo
          enddo
       endif
       !     Add the pseudo gravity afterwards...
       if (pseudo_gravity.eq.1) then
          do k = klo, khi
             do j = jlo, jhi
                do i = ilo, ihi
                   force(i,j,k,nYvel) = force(i,j,k,nYvel) + dV_control*rho(i,j,k)
                enddo
             enddo
          enddo
       endif
       !     End of velocity forcing
    endif

    if ((scomp+ncomp).gt.BL_SPACEDIM) then
       !     Scalar forcing
       do n = max(scomp+1,nRho), scomp+ncomp
          if (n.eq.nRho) then
             !     Density
             do k = klo, khi
                do j = jlo, jhi
                   do i = ilo, ihi
                      force(i,j,k,n) = zero
                   enddo
                enddo
             enddo
          else
             !     Other scalar
             do k = klo, khi
                do j = jlo, jhi
                   do i = ilo, ihi
                      force(i,j,k,n) = zero
                   enddo
                enddo
             enddo
          endif
       enddo
    endif

  end subroutine FORT_MAKEFORCE
end module prob_3D_module
