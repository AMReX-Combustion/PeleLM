#include <MyProb_F.H>

module prob_3D_module

  implicit none

  private

  public :: setupbc, amrex_probinit, getZone, bcfunction, init_data_new_mech, init_data, &
            zero_visc, flame_tracer_error, adv_error, &
            temp_error, mv_error, den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, FORT_ZVELFILL, chem_fill, press_fill

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

    use chem_driver, only: get_spec_name, P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(SDIM), probhi(SDIM)

#include <probdata.H>
#include <htdata.H>
#include <cdwrk.H>
#include <bc.H>
#include <visc.H>
#include <conp.H>
#include <INFL_FORCE_F.H>

#ifdef DO_LMC_FORCE
#include <forcedata.H>
#endif

    integer i, j, zone, iH2
    character*(maxspnml) spName

    namelist /fortin/ vorterr, temperr, adverr, tempgrad, &
         flametracval, probtype,&
         max_temp_lev, max_vort_lev, max_trac_lev,&
         traceSpecVal,phi_in,T_in,&
         turb_scale, V_in, &
         zstandoff, pertmag, nchemdiag,&
         blobx,bloby,blobz,blobr
    namelist /heattransin/ pamb, dpdt_factor, closed_chamber
    namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi,&
         strmwse_dir, nCompInflow, flct_file
    namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, zbase_control

#ifdef DO_LMC_FORCE
    REAL_T  twicePi, kxd, kyd, kzd
    REAL_T  thetaTmp, phiTmp
    REAL_T  cosThetaTmp, cosPhiTmp
    REAL_T  sinThetaTmp, sinPhiTmp
    REAL_T  px, py, pz, mp2, Ekh
    integer kx, ky, kz, mode_count, reduced_mode_count
    integer xstep, ystep, zstep

    REAL_T  Lx, Ly, Lz, Lmin, rn 
    REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff, pdk

    namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz,&
         forcing_type, spectrum_type, time_offset,&
         forcing_xlength, forcing_ylength, forcing_zlength,&
         forcing_time_scale_min, forcing_time_scale_max, &
         force_scale, forcing_epsilon, blrandseed,&
         use_rho_in_forcing, do_mode_division,&
         AXY, BXY, CXY, DXY, PXY, QXY, RXY,&
         AZX, BZX, CZX, DZX, PZX, QZX, RZX,&
         AYZ, BYZ, CYZ, DYZ, PYZ, QYZ, RYZ,&
         FTX, FTY, FTZ, TAT, TAP,&
         FPX, FPY, FPZ, FAX, FAY, FAZ
#endif

    !
    !      Build 'probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)

    character probin*(maxlen)
    integer n

    call bl_pd_is_ioproc(isioproc)

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

    vorterr = 1.e20
    temperr = zero
    adverr = 1.e20
    tempgrad  = 50.0d0
    flametracval = 0.0001d0
    probtype = 1
    max_temp_lev = 0
    max_vort_lev = 0
    max_trac_lev = 100
    traceSpecVal = 1.d-10
    pamb = 101325.d0
    dpdt_factor = 0.3d0
    closed_chamber = 0

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
    zstandoff = zero
    pertmag = 0.d0

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

    nchemdiag = 1

    read(untin,fortin)

    !     Initialize control variables that depend on fortin variables
    V_in_old = V_in

    if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev

    read(untin,heattransin)

    read(untin,flctin)

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
    if (probtype.ge.3) then
       convVel = one
    endif
#endif
    !     Load domain dimensions into common, and set up boundary functions
    domnlo(1) = problo(1)
    domnlo(2) = problo(2)
    domnlo(3) = problo(3)
    domnhi(1) = probhi(1)
    domnhi(2) = probhi(2)
    domnhi(3) = probhi(3)

    call setupbc()
    bcinit = .true.

    iH2 = -1
    do n = 1,Nspec
       call get_spec_name(spName,n)
       if (spName .eq. 'H2') iH2 = n-1
    end do
    zone = getZone(half*(domnlo(1)+domnhi(1)),half*(domnlo(2)+domnhi(2)),domnlo(3))
    scale_control = Y_bc(iH2,zone)*rho_bc(zone)*(domnhi(1)-domnlo(1))*(domnhi(2)-domnlo(2))
    if (h_control .gt. zero) then
       cfix = scale_control * h_control
    endif

#ifdef DO_LMC_FORCE
    if ( probtype.eq.3 ) then

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

       if (hack_lz) then 
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
                   if (mp2 .lt. 0.000001) then
                      write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
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
                         write (*,*) "kappa = ",kx,ky,kz,kappa
                         write (*,*) "Amplitudes"
                         write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                         write (*,*) "Frequencies"
                         write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

       !     Now lets break symmetry, have to assume high aspect ratio in z for now
       reduced_mode_count = 0
       do kz = mode_start, zstep - 1
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
                   if (mp2 .lt. 0.000001) then
                      write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
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
                         write (*,*) "kappa = ",kx,ky,kz,kappa
                         write (*,*) "Amplitudes"
                         write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                         write (*,*) "Frequencies"
                         write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
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
    !      
    !     Set random numbers for Ylm perts
    !

    !
    !     In parallel these must be consistent, so set the seed to have
    !     the same value across all processors
    !
    call blutilinitrand(134527)
    call blutilinitrand(111397)

    do i = lmodemin,lmodemax
       do j = 0,i
          call blutilrand(alphalm(i,j))
          call blutilrand(betalm(i,j))
          call blutilrand(gammalm(i,j))
          !           write(6,*)"numbers",i,j,alphalm(i,j),betalm(i,j),gammalm(i,j)
       end do
    end do

    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
    end if

  end subroutine amrex_probinit

  subroutine setupbc() bind(C, name="setupbc")

    use chem_driver, only: P1ATMMKS, get_spec_name
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    
    implicit none

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <htdata.H>
      
    REAL_T Patm
    REAL_T Yfu(maxspec)
    REAL_T Y_FU_H2, Y_FU_O2
    REAL_T X_FU_H2, X_FU_O2, X_FU_N2, Xsum
    REAL_T Xfu(maxspec)
    integer zone, n, b(SDIM),iN2, iH2, iO2

    character*(maxspnml) name
    data b / 1, 1, 1 /

    REAL_T X_H2(Nzones),X_O2(Nzones),T(Nzones)
#ifdef SWIRL
    !     phi = 0.8
    data X_H2 /  0.0774907749077d0,   zero, 0.0774907749077d0,   zero,    zero,   zero ,   zero /
    data X_O2  /   0.193726937269d0,  0.21d0, 0.193726937269d0,  0.21d0,  0.21d0, 0.21d0,  0.21d0 /
    data T     /        3.d2,          3.d2,         3.d2,        3.d2,    3.d2,   3.d2 ,   3.d2 /
#else
    !     phi = 0.7
    !      data X_H2 /  0.0684676292501d0,  zero, 0.0684676292501d0,  zero, zero /
    !      data X_O2  /  0.195621797857d0, 0.21d0, 0.195621797857d0, 0.21d0, zero /
    !     phi = 1.0
    data X_H2 /  0.095023d0,   zero,  0.095023d0,   zero,   zero /
    data X_O2  /  0.190045d0, 0.21d0,  0.190045d0,  0.21d0,  zero /
    data T     / 3.d2,         3.d2,         3.d2,      3.d2 , zero/
#endif
    Patm = pamb / P1ATMMKS()
    iN2 = -1
    do n = 1,Nspec
       call get_spec_name(name, n)
       if (name .eq. 'N2' ) iN2 = n
    end do

    !     For probtype == 3, set boundary conditions as constants
    iN2 = -1
    iH2 = -1
    iO2 = -1
    do n = 1,Nspec
       call get_spec_name(name,n)
       if (name .eq. 'N2' ) iN2  = n
       if (name .eq. 'H2') iH2 = n
       if (name .eq. 'O2' ) iO2  = n
    end do
 
    if (phi_in.ge.zero) then
 
       do n=1,Nspec
          Xfu(n) = zero
       enddo
 
       X_FU_H2 = 2.d0*phi_in
       X_FU_O2 = 1.d0
       X_FU_N2 = .79d0/.21d0
 
       Xsum = X_FU_H2 + X_FU_O2 + X_FU_N2
 
       Xfu(iH2) = X_FU_H2/ Xsum
       Xfu(iO2) = X_FU_O2/ Xsum
       Xfu(iN2) = X_FU_N2/ Xsum
 
       !         write(6,*)" in bc routine", Xfu(iH2),Xfu(iO2) ,Xfu(iN2)
 
 
       call CKXTY(Xfu,Yfu)

    else
       !     phi = 1.0
       Y_FU_H2 = 0.0551670662321d0
       Y_FU_O2  = 0.220068142419d0
       !     phi = 0.7
       !         Y_FU_H2 = 0.0392668168889d0
       !         Y_FU_O2  = 0.223771589041d0
 
       do n=1,Nspec
          Yfu(n) = zero
       enddo
       Yfu(iH2) = Y_FU_H2
       Yfu(iO2)  = Y_FU_O2
       Yfu(iN2)  = one - Y_FU_H2 - Y_FU_O2
 
    endif
 
    zone = 1
    do n=1,Nspec
       Y_bc(n-1,zone) = Yfu(n)
    end do
    T_bc(zone) = T_in
    u_bc(zone) = zero
    v_bc(zone) = V_in

    call RHOfromPTY(b,b,&
         rho_bc(zone),DIMARG(b),DIMARG(b),&
         T_bc(zone),DIMARG(b),DIMARG(b),&
         Y_bc(0,zone),DIMARG(b),DIMARG(b), Patm)

    call HMIXfromTY(b,b,&
         h_bc(zone),DIMARG(b),DIMARG(b),&
         T_bc(zone),DIMARG(b),DIMARG(b),&
         Y_bc(0,zone),DIMARG(b),DIMARG(b))
 
  end subroutine setupbc

  integer function getZone(x, y, z) bind(C, name="getZone")
    
    implicit none
    
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
    REAL_T x, y, z
    !
    !     ZONE DEF:
    !     BL_FUELPIPE, BL_COFLOW, BL_STICK, BL_WALL, BL_AMBIENT, BL_VOLUME
    !
    !     getZone = BL_VOLUME
    !     if ( ((direction.eq."z") .and. (z.le.domnlo(3))).or.
    !    &     ((direction.eq."y") .and. (y.le.domnlo(2))).or.
    !    &     ((direction.eq."x") .and. (x.le.domnlo(1)))) then

    getZone = 1

    !      endif

  end function getZone
      
  !
  !     Compute the average velocity in a dx(1)*dx(2) square
  !     x,y  position on entry plane
  !     Rf     fuel tube radius
  !     wallTh fuel tube thickness
  !     Ro     oxidizer tube radius
  !     stTh   thickness of stick crossing origin
  !     ypls   dimensionless laminar sublayer thickness
  !     ypt    dimensionless height of transition layer
  !     nu     viscosity of fluid
  !     F      flux of flux (cm**2/s)
  !     r0     radius where log-law is replaced by parabola
  !
  subroutine VAvg(V,x,y,dx,Rf,wallTh,Ro,stTh,dBL,stBL,&
       Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW,probtype,coCUTOFFp) bind(C, name="VAvg")

    implicit none
    REAL_T V(SDIM)
    REAL_T x,y,dx(SDIM),Rf,wallTh,Ro,stTh,dBL,stBL
    REAL_T Vf,Vco_l,Vco_r,tVco_l,tVco_r,time,swK,swW,coCUTOFFp
    integer probtype

    REAL_T incx,incy,x1,x2,r
    REAL_T Vincx,Vincy,Vincz
    integer i,j,M
    parameter (M=10)
    !
    !     size of subvolumes, there are M*M of them
    !
    incx = dx(1) / M
    incy = dx(2) / M
    V(1) = 0.d0
    V(2) = 0.d0
    V(3) = 0.d0
    do i=1,M
       do j=1,M
          x1 = x - 0.5d0*dx(1) + (i-0.5d0)*incx
          x2 = y - 0.5d0*dx(2) + (j-0.5d0)*incy
          r = sqrt(x1*x1 + x2*x2)
          Vincx = 0.d0
          Vincy = 0.d0
          Vincz = 0.d0
          if (r.le.Ro) then
#ifdef SWIRL
!               eta = r/Ro
!               Vincz = Vf*(0.63336d0 - 2.3529d0*eta + 25.595d0*eta**2 - 115.84d0*eta**3
!     $              + 255.47d0*eta**4 - 260.06d0*eta**5+ 99.786d0*eta**6)
#else
             Vincz = Vf*TANH(two*(one - r/Ro)/dBL)
#endif
          else if (r.ge.Ro+wallTh) then
             if (time .le. tVco_l) then
                Vincz = Vco_l
             else if (time .ge. tVco_r) then
                Vincz = Vco_r
             else
                Vincz = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
             endif

             !     For stagnation case, ramp co-flow to zero at coCUTOFF
             !              if (probtype.eq.3) then
             !                 coCUTOFF = coCUTOFFp - Ro - wallTh
             !                 Vincz = Vincz*MAX(zero,(coCUTOFF-r+Ro+wallTh)/coCUTOFF)
             !              endif
          endif

          if (stTh.gt.zero) then
#ifdef SWIRL
             call bl_abort('better not be here')
#endif
             if (ABS(x1).le.half*stTh.and.r.le.Ro+wallTh) then
                Vincz = 0.d0
             elseif(r.le.Ro+wallTh)then
                Vincz = Vincz * TANH(2.d0*(ABS(x1)-half*stTh)/stBL)
             end if
          end if

#ifdef SWIRL
          if (r.le.Ro) then
             !
             !     JBB fit
             !              Vswirl = swK*1.1d0*EXP(-20*(.95-r/Rf)**2)*TANH(40*(1.d0-r/Rf))
             !
             !     Original swirl parameterization
             !              eta = r/Ro
             !              if ((r.gt.Rf*(one-swW)).and.(r.le.Rf)) then
             !                 Vswirl = swK*(SIN(Pi*(one-r/Rf)/swW))**2
             !                 Vincx = Vswirl*(-x2/r)
             !                 Vincy = Vswirl*( x1/r)
             !              end if
             !
             !     Fit from averaged EBHC calc
             !              Vswirl = Vf*(-0.067641d0 + 2.0113d0*eta - 19.305d0*eta**2 + 72.404d0*eta**3
             !     &                     - 126.14d0*eta**4 + 111.66d0*eta**5 - 31.098d0*eta**6)
             !              Vradial =Vf*(-0.0022569d0 + 0.3562d0*eta - 0.84751d0*eta**2 + 1.8533d0*eta**3
             !     &                    - 5.1925d0*eta**4 + 7.7047d0*eta**5 - 3.8714d0*eta**6)
             !              Vincx = Vswirl*(-x2/r) + Vradial*x1/r
             !              Vincy = Vswirl*( x1/r) + Vradial*x2/r
          end if
#endif
              
          V(1) = V(1) + Vincx
          V(2) = V(2) + Vincy
          V(3) = V(3) + Vincz
       end do
    end do
    V(1) = V(1) / (M*M)
    V(2) = V(2) / (M*M)
    V(3) = V(3) / (M*M)
  end subroutine VAvg

! ::: -----------------------------------------------------------
      
  subroutine bcfunction(orient,x,y,z,time,u,v,w,rho,Yl,T,h,dx,getuvw) bind(C, name="bcfunction")

    implicit none

    integer orient
    REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM)
    logical getuvw

#include <htdata.H>
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer n, zone

    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    if (orient .lt. 6) then
       zone = getZone(x,y,z)
       u = u_bc(zone)
       v = v_bc(zone)
       w = 0.d0
       rho = rho_bc(zone)
       do n = 0, Nspec-1
          Yl(n) = Y_bc(n,zone)
       end do
       T = T_bc(zone)
       h = h_bc(zone)

       if (getuvw) then
          if (probtype.eq.3 .or. probtype.eq.1) then
             u = zero
             v = zero
             w =  V_in + (time-tbase_control)*dV_control
          endif
       endif
    else
       write(6,*) 'No boundary condition for orientation = ', orient
       call bl_abort(" ")
    endif

  end subroutine bcfunction

  REAL_T function shapet(r) bind(C,name="shapet")
    REAL_T r
    if(r.le.1.d0)then
       shapet = tanh(100.d0*(1-r))
    else
       shapet = 0.d0
    endif
  end function shapet

  REAL_T function shapen(r) bind(C,name="shapen")
    REAL_T r
    if(r.le.1.d0)then
       shapen = tanh(100.d0*(1-r))
    else
       shapen = 0.d0
    endif
  end function shapen

  REAL_T function turbSclT(r) bind(C,name="turbSclT")
    REAL_T r, eta
#include <probdata.H>
    turbSclT = 0.d0
    eta = r/Ro
    if(r.lt.Ro)then
#ifdef SWIRL
       turbSclT = V_in*(0.10278d0 + 0.34958d0 *eta - 0.86045d0 *eta**2&
            - 6.9187d0 *eta**3 + 30.903d0 *eta**4 - 31.34d0 *eta**5 + 8.588d0 *eta**6)

#else
       turbSclT = 1.d0
#endif
    endif
  end function turbSclT

  REAL_T function turbSclR(r) bind(C,name="turbSclR")
    REAL_T r, eta
#include <probdata.H>
    turbSclR = 0.d0
    eta = r/Ro
    if(r.lt.Ro)then
#ifdef SWIRL
       turbSclR=V_in*(0.086871d0 + 1.5445d0*eta - 17.316d0*eta**2 + 79.469d0*eta**3&
            - 164.43d0*eta**4 + 160.85d0*eta**5 - 60.233d0*eta**6)
       
#else
       turbSclR = (0.25d0 + 0.95d0*tanh((1.d0-r/Ro)/.0502362048))/1.2d0
#endif
    endif
  end function turbSclR

  REAL_T function turbSclX(x,y,z) bind(C,name="turnSclX")
#include <probdata.H>
    REAL_T x, y, z, shapet, eta
    eta = SQRT(x*x+y*y)/Ro
    turbSclX = shapet(eta)
    if (stBL.gt.zero) then
       turbSclX = turbSclX*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
    end if
  end function turbSclX

  REAL_T function turbSclY(x,y,z) bind(C,name="turbSclY")
    REAL_T x, y, z
    turbSclY = turbSclX(x,y,z)
  end function turbSclY

  REAL_T function turbSclZ(x,y,z) bind(C,name="turbSclZ")
#include <probdata.H>
    REAL_T x, y, z, eta
    turbSclZ = zero
    eta = SQRT(x*x+y*y)/Ro
#ifdef SWIRL
    turbSclZ = one
    return
!      if(eta.lt.one)then
!         turbSclZ = V_in*(-0.020318d0 + 2.18d0*eta - 18.157d0*eta**2 + 66.538d0*eta**3 - 
!     &         112.55d0*eta**4 + 94.931d0*eta**5 - 32.329d0*eta**6)
!      endif
#else
    if(eta.lt.one)then
       turbSclZ = 2.d0 -1.d0 * tanh((1-eta)/.0428222042d0)
    endif
    if (stBL.gt.zero) then
       turbSclZ = turbSclZ*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
    end if
#endif
  end function turbSclZ

  subroutine init_data_new_mech(level,time,lo,hi,nscal,&
       vel,scal,DIMS(state),press,DIMS(press),&
       delta,xlo,xhi) bind(C, name="init_data_new_mech")

    use chem_driver, only: P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    
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
 
#include <htdata.H>
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

    use chem_driver, only: P1ATMMKS
    use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
    use chem_driver, only: get_spec_name

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
    
#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>
#include <INFL_FORCE_F.H>

    integer i, j, k, n, iO2
    REAL_T x, y, z, Yl(maxspec), Xl(maxspec), Patm
    REAL_T pmf_vals(maxspec+3), r, theta
    REAL_T  xloc,yloc,zloc,xtemp
    character*(maxspnml) name
    REAL_T   phi,rfront
    integer l,m,ctr


#if defined(BL_DO_FLCT)
    integer DIMDEC(uflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T uflct(:,:,:)
    allocatable uflct
    REAL_T fileDataZLO, fileDataZHI, gridLo(SDIM)
    integer bc(SDIM,2)

    do i = 1, SDIM
       gridLo(i) = xlo(i) - half*probSizeFile(i)
    end do
    fileDataZLO = half*(domnhi(3)+domnlo(3)-probSizeFile(3))
    fileDataZHI = half*(domnhi(3)+domnlo(3)+probSizeFile(3))
    !         
    !     Make a platter of data in XY, will march over z
    !
    do i = 1, SDIM
       loFlctArray(i) = lo(i)
       hiFlctArray(i) = hi(i)
    enddo

    strmwse_dir = 3
    loFlctArray(strmwse_dir) = 1
    hiFlctArray(strmwse_dir) = 1
    call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
    allocate(uflct(DIMV(uflct)))

    do i = 1, SDIM
       bc(i,1) = EXT_DIR
       bc(i,2) = EXT_DIR
    end do
    convVelSAVE = convVel
    convVel = 15.d0
#endif 
           
    if (iN2.lt.1 .or. iN2.gt.Nspec) then
       call bl_pd_abort()
    endif
    do n = 1,Nspec
       call get_spec_name(name, n)
       if (name .eq. 'O2' ) iO2 = n
    end do

    do k= lo(3),hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n=1,SDIM
                vel(i,j,k,n) = zero
             end do
          end do
       end do
    enddo

    do k= lo(3),hi(3)
       z = (dfloat(k)+.5)*delta(3)+domnlo(3)            
       do j = lo(2), hi(2)
          y = (dfloat(j)+0.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (dfloat(i)+.5)*delta(1)+domnlo(1)
             r = SQRT((x-blobx)**2+(y-bloby)**2+(z-blobz)**2)

             !              if (r.LE.delta(1)) then
             !                 theta = zero
             !                 phi = zero
             !              else
             !                 theta = ATAN2(SQRT((x-blobx)**2+(y-bloby)**2)/r,(z-blobz)/r)-Pi
             !                 phi = ATAN2((y-bloby)/r,(x-blobx)/r)
             !              endif
               
             rfront = zero
             ctr = 0
             do l = lmodemin,lmodemax
                do m = 0,l
                   ctr = ctr + 1

                   ! set local coordinates
                   if (r.LE.delta(1)) then
                      theta = zero
                      phi = zero
                   else
                      xtemp = cos(2.d0*Pi*alphalm(l,m))*x + sin(2.d0*Pi*alphalm(l,m))*y
                      yloc = -sin(2.d0*Pi*alphalm(l,m))*x + cos(2.d0*Pi*alphalm(l,m))*y
                      
                      xloc = cos(2.d0*Pi*gammalm(l,m))*xtemp + sin(2.d0*Pi*gammalm(l,m))*z
                      zloc = -sin(2.d0*Pi*gammalm(l,m))*xtemp + cos(2.d0*Pi*gammalm(l,m))*z

                      theta = ATAN2(SQRT((xloc-blobx)**2+(yloc-bloby)**2)/r,(zloc-blobz)/r)-Pi
                      phi = ATAN2((yloc-bloby)/r,(xloc-blobx)/r)
                   endif
                   rfront = rfront +&
                        betalm(l,m)&
                        * Ylm(l,m,phi,COS(theta))

                end do
             end do
             if (ctr.GT.0) rfront = rfront/ctr

             rfront = (blobr + pertmag* rfront*blobr - r + zstandoff) *100.d0
               
             call pmf(rfront,rfront,pmf_vals,n)
               
             if (n.ne.Nspec+3) then
                call bl_abort('INITDATA: n .ne. Nspec+3')
             endif

             scal(i,j,k,Temp) = pmf_vals(1)
             do n = 1,Nspec
                Xl(n) = pmf_vals(3+n)
             end do
             
             call CKXTY(Xl,YL)
               
             do n = 1,Nspec
                scal(i,j,k,FirstSpec+n-1) = Yl(n)
             end do

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

#if defined(BL_DO_FLCT)          
    convVel = convVelSAVE
    deallocate(uflct)
#endif
  end subroutine init_data
 
  REAL_T function Ylm(l,m,phi,costheta) bind(C, name="Ylm")
    REAL_T phi, costheta, facN, facD, plgndr
    integer m, l, n
    facN = l - m
    if (facN .EQ. 0) facN=1
    do n = l-m-1,2,-1
       facN = facN*n
    end do
    facD = l+m
    if (facD.EQ.0) facD=1
    do n = l+m-1,2,-1
       facD = facD*n
    end do
    Ylm = SQRT((2*l+1)*facN/(4*Pi*facD))*plgndr(l,m,costheta)*COS(m*phi)
  end function Ylm
 
  REAL_T FUNCTION plgndr(l,m,x) bind(C, name="plgndr")
    INTEGER l,m
    REAL_T x
    INTEGER i,ll
    REAL_T fact,pll,pmm,pmmp1,somx2
    if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) then
       call bl_abort('bad arguments in plgndr')
    endif
    pmm=1.d0
    if(m.gt.0) then
       somx2=sqrt((1.d0-x)*(1.d0+x))
       fact=1.d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       enddo
    endif
    if(l.eq.m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          enddo
          plgndr=pll
       endif
    endif
  end FUNCTION plgndr

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
#include <htdata.H>

    integer i, j, k, n, Tid, RHid, YSid, YEid, ys, ye
    logical do_T, do_RH, do_Y
    REAL_T xl, xr, xh, yb, yt, yh, z

    if (probtype.lt.4) then
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

  subroutine flame_tracer_error (tag,DIMS(tag),set,clear,&
       ftrac,DIMS(ftrac),lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level) bind(C, name="flame_tracer_error")

    implicit none

    integer   DIMDEC(ftrac)
    integer   DIMDEC(tag)
    integer   lo(SDIM), hi(SDIM)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    integer   tag(DIMV(tag))
    REAL_T    ftrac(DIMV(ftrac), nvar)

    integer   i, j, k

#include <probdata.H>

    if (level.lt.max_trac_lev) then
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j,k) = merge(set,tag(i,j,k),&
                     ftrac(i,j,k,1).gt.flametracval)
             enddo
          enddo
       enddo
!$omp end parallel do
    endif
  end subroutine flame_tracer_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the 
! ::: density gradient
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: DIMS(tag) => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: adv       => scalar array
! ::: DIMS(adv) => index extent of scalar array
! ::: lo,hi     => index extent of grid
! ::: nvar      => number of components in rho array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: -----------------------------------------------------------
  subroutine adv_error (tag,DIMS(tag),set,clear,&
       adv,DIMS(adv),lo,hi,nvar,&
       domlo,domhi,delta,xlo,&
       problo,time,level) bind(C, name="adv_error")

    implicit none

    integer   DIMDEC(tag)
    integer   DIMDEC(adv)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    REAL_T    delta(SDIM), xlo(SDIM), problo(SDIM), time
    integer   tag(DIMV(tag))
    REAL_T    adv(DIMV(adv),nvar)

#include <probdata.H>
      
    if ((probtype .eq. 8 .or. probtype .eq. 28 .or. probtype .eq. 29) &
         .and. (time .eq. zero)) then
       call mv_error(tag,DIMS(tag),set,clear,&
            adv,DIMS(adv),lo,hi,nvar,&
            domlo,domhi,delta,xlo,&
            problo,time,level)
    endif
    
  end subroutine adv_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the
! ::: temperature gradient
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: DIMS(tag) => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => density array
! ::: DIMS(temp)=> index extent of temp array
! ::: lo,hi     => index extent of grid
! ::: nvar      => number of components in rho array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: -----------------------------------------------------------
  subroutine temp_error (tag,DIMS(tag),set,clear,&
       temperature,DIMS(temp),lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level) bind(C, name="temp_error")    

    implicit none

    integer   DIMDEC(tag)
    integer   DIMDEC(temp)
    integer   nvar, set, clear, level
    integer   domlo(SDIM), domhi(SDIM)
    integer   lo(SDIM), hi(SDIM)
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    integer   tag(DIMV(tag))
    REAL_T    temperature(DIMV(temp),nvar)

    integer   i, j, k, ng

#include <probdata.H>

    ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),ARG_H3(temp)-hi(3),&
         lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp),lo(3)-ARG_L3(temp))

    if (ng .lt. 1) then
       write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
       call bl_abort(" ")
    endif

    if (level .lt. max_temp_lev) then
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j,k) = merge(set,tag(i,j,k),temperature(i,j,k,1).lt.temperr)
             enddo
          enddo
       enddo
!$omp end parallel do
    endif
  end subroutine temp_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the 
! ::: magnitude of vorticity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: DIMS(tag) => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vort      => array of vorticity values
! ::: DIMS(vor) => index extent of vort array
! ::: nvar      => number of components in vort array (should be 1)
! ::: lo,hi     => index extent of grid
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: -----------------------------------------------------------
  subroutine mv_error (tag,DIMS(tag),set,clear,&
       vort,DIMS(vort),lo,hi,nvar,&
       domlo,domhi,dx,xlo,&
       problo,time,level) bind(C, name="mv_error")

    implicit none

    integer   DIMDEC(tag)
    integer   DIMDEC(vort)
    integer   nvar, set, clear, level
    integer   lo(SDIM), hi(SDIM)
    integer   domlo(SDIM), domhi(SDIM)
    REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
    integer   tag(DIMV(tag))
    REAL_T    vort(DIMV(vort),nvar)

    integer   i, j, k

#include <probdata.H>

    if (level .lt. max_vort_lev) then
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j,k) = merge(set,tag(i,j,k),&
                     ABS(vort(i,j,k,1)).ge.vorterr*2.d0**level)
             enddo
          enddo
       enddo
!$omp end parallel do
    end if
  end subroutine mv_error

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
! ::         have been filled with a large real number.
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

  subroutine den_fill (den,DIMS(den),domlo,domhi,delta,&
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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
! ::         have been filled with a large real number.
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

  subroutine temp_fill (temp,DIMS(temp),domlo,domhi,delta,xlo,time,bc) bind(C, name="temp_fill")

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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

  subroutine rhoh_fill (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,time,bc) bind(C, name="rhoh_fill")

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

  end subroutine rhoh_fill
  !
  ! Fill x, y & z velocity at once.
  !
  subroutine vel_fill (vel,DIMS(vel),domlo,domhi,delta,xlo,time,bc) bind(C, name="vel_fill")

    implicit none

    integer DIMDEC(vel), bc(SDIM,2,SDIM)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  vel(DIMV(vel),SDIM)

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <INFL_FORCE_F.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

#if defined(BL_DO_FLCT)
    integer DIMDEC(uflct)
    REAL_T t_flct, dt_flct
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T  uflct(:,:,:), vflct(:,:,:),  wflct(:,:,:)
    allocatable uflct, vflct, wflct
#endif
    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(vel)
    lo(2) = ARG_L2(vel)
    lo(3) = ARG_L3(vel)
    hi(1) = ARG_H1(vel)
    hi(2) = ARG_H2(vel)
    hi(3) = ARG_H3(vel)

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
       allocate(vflct(DIMV(uflct)))
       allocate(wflct(DIMV(uflct)))
       dt_flct = time - tbase_control
       t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
       call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct,&
            bc(1,1,1), domnlo, domnhi)
       call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, t_flct,&
            bc(1,1,2), domnlo, domnhi)
       call INFL_FILL(FLCT_ZVEL, DIMS(uflct), wflct, xlo, delta, t_flct,&
            bc(1,1,3), domnlo, domnhi)
    endif
#endif

    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1),&
         DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2),&
         DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3),&
         DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))

    if (lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)

                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(1,1,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 1) then
                      !                       vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                      vel(i,j,k,1) = u + uflct(1,j,k)*turb_scale
                   endif
#endif
                endif

                if (bc(1,1,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 1) then
                      !                       vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                      vel(i,j,k,2) = v + vflct(1,j,k)*turb_scale
                   endif
#endif
                endif

                if (bc(1,1,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 1) then
                      !                       vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                      vel(i,j,k,3) = w + wflct(1,j,k)*turb_scale
                   endif
#endif
                endif
             enddo
          enddo
       enddo
    endif

    if (hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)

                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(1,2,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 1) then
                      vel(i,j,k,1) = u + uflct(1,j,k)*turbSclX(x,y,z)*turb_scale
                   endif
#endif
                endif

                if (bc(1,2,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 1) then
                      vel(i,j,k,2) = v + vflct(1,j,k)*turbSclY(x,y,z)*turb_scale
                   endif
#endif
                endif

                if (bc(1,2,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 1) then
                      vel(i,j,k,3) = w + wflct(1,j,k)*turbSclZ(x,y,z)*turb_scale
                   endif
#endif
                endif
             enddo
          enddo
       enddo
    endif

    if (lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(2,1,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                      vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                   endif
#endif
                endif

                if (bc(2,1,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                      vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                   endif
#endif
                endif

                if (bc(2,1,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                      vel(i,j,k,3) = w + wflct(i,1,k)*turb_scale
                   endif

#endif
                endif
             enddo
          enddo
       enddo
    endif

    if (hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(2,2,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,1) = u + uflct(i,1,k)*turbSclX(x,y,z)*turb_scale
                      vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                   endif
#endif
                endif

                if (bc(2,2,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,2) = v + vflct(i,1,k)*turbSclY(x,y,z)*turb_scale
                      vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                   endif
#endif
                endif

                if (bc(2,2,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 2) then
                      !                       vel(i,j,k,3) = w + wflct(i,1,k)*turbSclZ(x,y,z)*turb_scale
                      vel(i,j,k,3) = w + wflct(i,1,k)*turb_scale
                   endif
#endif
                endif
             enddo
          enddo
       enddo
    endif

    if (lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(3,1,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 3) then
                      vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                   endif
#endif
                endif

                if (bc(3,1,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 3) then
                      vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                   endif
#endif
                endif

                if (bc(3,1,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceLo .and. strmwse_dir .eq. 3) then
                      vel(i,j,k,3) = w + wflct(i,j,1)*turb_scale
                   endif
#endif
                endif
             enddo
          enddo
       enddo
    endif

    if (hi(3).gt.domhi(3)) then
       do k = domhi(3)+1, hi(3)
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(3,2,1).eq.EXT_DIR) then
                   vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 3) then
                      !                       vel(i,j,k,1) = u + uflct(i,j,1)*turbSclX(x,y,z)*turb_scale
                      vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                   endif
#endif
                endif

                if (bc(3,2,2).eq.EXT_DIR) then
                   vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 3) then
                      !                       vel(i,j,k,2) = v + vflct(i,j,1)*turbSclY(x,y,z)*turb_scale
                      vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                   endif
#endif
                endif

                if (bc(3,2,3).eq.EXT_DIR) then
                   vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                   if (forceHi .and. strmwse_dir .eq. 3) then
                      !                       vel(i,j,k,3) = w + wflct(i,j,1)*turbSclZ(x,y,z)*turb_scale
                      vel(i,j,k,3) = w + wflct(i,j,1)*turb_scale
                   endif
#endif
                endif
             enddo
          enddo
       enddo
    endif

#if defined(BL_DO_FLCT)
    if (forceInflow) then
       deallocate(uflct)
       deallocate(vflct)
       deallocate(wflct)
    endif
#endif

  end subroutine vel_fill

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

  subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta,xlo,time,bc) bind(C,name="FORT_XVELFILL")

    implicit none

    integer DIMDEC(xvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  xvel(DIMV(xvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <INFL_FORCE_F.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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

  subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,xlo,time,bc) bind(C,name="FORT_YVELFILL")

    implicit none

    integer DIMDEC(yvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  yvel(DIMV(yvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <INFL_FORCE_F.H>

    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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

  subroutine FORT_ZVELFILL (zvel,DIMS(zvel),domlo,domhi,delta,xlo,time,bc) bind(C,name="FORT_ZVELFILL")

    implicit none

    integer DIMDEC(zvel), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  zvel(DIMV(zvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <INFL_FORCE_F.H>

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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

  subroutine all_chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,time,bc) bind(C, name="all_chem_fill")

    implicit none

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    integer DIMDEC(rhoY), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  rhoY(DIMV(rhoY),Nspec)

    integer i, j, k, n
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(maxspec), T, h

    integer lo(SDIM), hi(SDIM)

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

    do n = 1,Nspec
       call filcc (rhoY(lo(1),lo(2),lo(3),n),&
            DIMS(rhoY),domlo,domhi,delta,xlo,bc)
    end do

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do k = lo(3),hi(3)
             z = (float(k)+.5)*delta(3)+domnlo(3)
             do j = lo(2), hi(2)
                y = (float(j)+.5)*delta(2)+domnlo(2)
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,Nspec
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
             enddo
          enddo
       enddo
    endif

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

  subroutine chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,time,bc,id) bind(C, name="chem_fill")

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
                call bcfunction(XLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(XHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(YHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZLO,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                call bcfunction(ZHI,x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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

  subroutine press_fill (p,DIMS(p),domlo,domhi,dx,xlo,time,bc) bind(C, name="press_fill")

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
    !  SETTING XLO
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
    ! SETTING XHI
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
    ! SETTING YLO
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
    ! SETTING YHI
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
    ! SETTING ZLO
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
    ! SETTING ZHI
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

end module prob_3D_module
