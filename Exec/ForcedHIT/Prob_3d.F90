
#include <MyProb_F.H>

module prob_3D_module

  implicit none

  private
  
  public :: amrex_probinit,setupbc, getZone, bcfunction, init_data_new_mech, init_data, &
            zero_visc, FORT_DENERROR, flame_tracer_error, adv_error, &
            temp_error, mv_error, den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, FORT_ZVELFILL, chem_fill, press_fill, &
#ifdef MOREGENGETFORCE
            FORT_AVERAGE_EDGE_STATES, &
#endif
            FORT_MAKEFORCE, read_csv, locate

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
  
      
      use chem_driver, only: P1ATMMKS
      use amrex_constants_module, only: HALF
      use amrex_fort_module
      use probdata_module
      
      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(SDIM), probhi(SDIM)

#include <probdata.H>
#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif
#include <visc.H>
#include <conp.H>

#ifdef DO_LMC_FORCE
#include <forcedata.H>
#endif

      integer i,istemp
      REAL_T FORT_P1ATMMKS, area

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, &
                       flametracval, probtype, &
     		        max_temp_lev, max_vort_lev, max_trac_lev, &
                       traceSpecVal,phi_in,T_in,  &
                       turb_scale, V_in, V_co, &
                       standoff, pertmag, nchemdiag, splitx, xfrontw, &
                       splity, yfrontw, blobx, bloby, blobz, blobr, &
                       blobT, Tfrontw, stTh, fuel_N2_vol_percent, iname, inres, &
                       urms0, uin_norm, siname, sinres, Gain, sGain
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber
#if defined(BL_DO_FLCT)
      integer nCompFile
      namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi, &
          strmwse_dir, nCompInflow, flct_file, convVelInit
#endif
      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
          zbase_control, pseudo_gravity, istemp,corr,controlVelMax,navg_pnts
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

      namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz, &
                       forcing_type, spectrum_type, time_offset, forcing_twice_wavelength, &
                       forcing_xlength, forcing_ylength, forcing_zlength, &
                       forcing_time_scale_min, forcing_time_scale_max,  &
                       force_scale, forcing_epsilon, blrandseed, &
                       use_rho_in_forcing, do_mode_division, &
                       div_free_force, moderate_zero_modes
#endif

!
!      Build `probin' filename -- the name of file containing fortin namelist.
!
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)

#if defined(BL_DO_FLCT)
      integer ierr
#endif

      ! Local variables
      REAL_T, dimension(:), allocatable :: data
      REAL_T lint, ubar, vbar, wbar, nu, dudx2, dudy2, dudz2, dvdx2, dvdy2, dvdz2, dwdx2, dwdy2, dwdz2
      REAL_T linf, uinf, einf, tinf
      REAL_T Z2mean, Zmean
      integer(kind=8) :: nnx, nny, nnz
      integer len, j, k
      len = len_trim(probtype)


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
      istemp = 0
      navg_pnts = 10

      read(untin,fortin)
      
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
      scale_control = Y_bc(fuelID-1,BL_FUELPIPE) * rho_bc(BL_FUELPIPE) * area

      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif

#ifdef DO_LMC_FORCE
      if ( (probtype.eq.BL_PROB_FORCEDHIT) ) then

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

!     Now let's break symmetry, have to assume high aspect ratio in z for now
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

if ((probtype(1:len).eq.BL_PROB_FORCEDHIT) .or. (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then


   ! Read in the velocity data
   nnx = int8(inres)
   nny = int8(inres)
   nnz = int8(inres)

   allocate(data(0:nnx*nny*nnz*6-1))
   allocate(xinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(yinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(zinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(uinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(vinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(winput(0:nnx-1,0:nny-1,0:nnz-1))
   ! Only read in CSV files for now
   call read_csv(iname, nnx, nny, nnz, data, int8(6))

   uinput = urms0 / uin_norm * reshape(data(3::6), (/nnx, nny, nnz/))
   vinput = urms0 / uin_norm * reshape(data(4::6), (/nnx, nny, nnz/))
   winput = urms0 / uin_norm * reshape(data(5::6), (/nnx, nny, nnz/))
   xinput = reshape(data(0::6), (/nnx, nny, nnz/))
   yinput = reshape(data(1::6), (/nnx, nny, nnz/))
   zinput = reshape(data(2::6), (/nnx, nny, nnz/))

   ! Get the xarray table and the differences
   allocate(xarray(0:nnx-1))
   allocate(xdiff(0:nnx-1))
   xarray(0:nnx-1) = xinput(:,0,0)
   xdiff(:nnx-2) = xarray(1:) - xarray(:nnx-2)
   xdiff(nnx-1) = xarray(nnx-1) - xarray(nnx-2)

   ! Dimensions of the input box
   Linput = maxval(xinput(:,0,0)) + HALF*xdiff(nnx-1)

   ! Compute Aforce and tke0 for use in the forcing scheme
   ! Assume integral length scale is about 20% of computational domain
   lint = 0.2*(domnhi(1)-domnlo(1))
   ! Compute ubar, vbar, and wbar
   ubar = 0.d0
   vbar = 0.d0
   wbar = 0.d0
   do k = 0, nnz-1
      do j = 0, nny-1
         do i = 0, nnx-1
            ubar = ubar + uinput(i,j,k)
            vbar = vbar + vinput(i,j,k)
            wbar = wbar + winput(i,j,k)
         enddo
      enddo
   enddo
   ubar = ubar / (nnx*nny*nnz)
   vbar = vbar / (nnx*nny*nnz)
   wbar = wbar / (nnx*nny*nnz)
   write(*,*) "ubar = ", ubar, "vbar = ", vbar, "wbar = ", wbar
   ! Compute tke0
   tke0 = 0.d0
   do k = 0, nnz-1
      do j = 0, nny-1
         do i = 0, nnx-1
            tke0 = tke0 + 0.5*( (uinput(i,j,k)-ubar)**2 + (vinput(i,j,k)-vbar)**2 + (winput(i,j,k)-wbar)**2 )
         enddo
      enddo
   enddo
   tke0 = tke0 / (nnx*nny*nnz)

   ! Compute dissipation, assume periodicity
   nu = 1.5d-5
   ! Compute derivatives in x
   dudx2 = 0.0
   dvdx2 = 0.0
   dwdx2 = 0.0
   do k = 0, nnz-1
      do j = 0, nny-1
         do i = 1, nnx-2
            dudx2 = dudx2 + ((uinput(i+1,j,k) - uinput(i-1,j,k)) / (xinput(i+1,j,k) - xinput(i-1,j,k)))**2
            dvdx2 = dvdx2 + ((vinput(i+1,j,k) - vinput(i-1,j,k)) / (xinput(i+1,j,k) - xinput(i-1,j,k)))**2
            dwdx2 = dwdx2 + ((winput(i+1,j,k) - winput(i-1,j,k)) / (xinput(i+1,j,k) - xinput(i-1,j,k)))**2
         enddo
         dudx2 = dudx2 + ((uinput(1,j,k) - uinput(nnx-1,j,k)) / (xinput(1,j,k) - xinput(nnx-1,j,k)))**2
         dudx2 = dudx2 + ((uinput(0,j,k) - uinput(nnx-2,j,k)) / (xinput(0,j,k) - xinput(nnx-2,j,k)))**2
         dvdx2 = dvdx2 + ((vinput(1,j,k) - vinput(nnx-1,j,k)) / (xinput(1,j,k) - xinput(nnx-1,j,k)))**2
         dvdx2 = dvdx2 + ((vinput(0,j,k) - vinput(nnx-2,j,k)) / (xinput(0,j,k) - xinput(nnx-2,j,k)))**2
         dwdx2 = dwdx2 + ((winput(1,j,k) - winput(nnx-1,j,k)) / (xinput(1,j,k) - xinput(nnx-1,j,k)))**2
         dwdx2 = dwdx2 + ((winput(0,j,k) - winput(nnx-2,j,k)) / (xinput(0,j,k) - xinput(nnx-2,j,k)))**2
      enddo
   enddo
   dudx2 = dudx2 / (nnx*nny*nnz)
   dvdx2 = dvdx2 / (nnx*nny*nnz)
   dwdx2 = dwdx2 / (nnx*nny*nnz)

   ! Compute derivatives in y
   dudy2 = 0.0
   dvdy2 = 0.0
   dwdy2 = 0.0
   do k = 0, nnz-1
      do i = 0, nnx-1
         do j = 1, nny-2
            dudy2 = dudy2 + ((uinput(i,j+1,k) - uinput(i,j-1,k)) / (yinput(i,j+1,k) - yinput(i,j-1,k)))**2
            dvdy2 = dvdy2 + ((vinput(i,j+1,k) - vinput(i,j-1,k)) / (yinput(i,j+1,k) - yinput(i,j-1,k)))**2
            dwdy2 = dwdy2 + ((winput(i,j+1,k) - winput(i,j-1,k)) / (yinput(i,j+1,k) - yinput(i,j-1,k)))**2
         enddo
         dudy2 = dudy2 + ((uinput(i,1,k) - uinput(i,nny-1,k)) / (yinput(i,1,k) - yinput(i,nny-1,k)))**2
         dudy2 = dudy2 + ((uinput(i,0,k) - uinput(i,nny-2,k)) / (yinput(i,0,k) - yinput(i,nny-2,k)))**2
         dvdy2 = dvdy2 + ((vinput(i,1,k) - vinput(i,nny-1,k)) / (yinput(i,1,k) - yinput(i,nny-1,k)))**2
         dvdy2 = dvdy2 + ((vinput(i,0,k) - vinput(i,nny-2,k)) / (yinput(i,0,k) - yinput(i,nny-2,k)))**2
         dwdy2 = dwdy2 + ((winput(i,1,k) - winput(i,nny-1,k)) / (yinput(i,1,k) - yinput(i,nny-1,k)))**2
         dwdy2 = dwdy2 + ((winput(i,0,k) - winput(i,nny-2,k)) / (yinput(i,0,k) - yinput(i,nny-2,k)))**2
      enddo
   enddo
   dudy2 = dudy2 / (nnx*nny*nnz)
   dvdy2 = dvdy2 / (nnx*nny*nnz)
   dwdy2 = dwdy2 / (nnx*nny*nnz)

   ! Compute derivatives in z
   dudz2 = 0.0
   dvdz2 = 0.0
   dwdz2 = 0.0
   do j = 0, nny-1
      do i = 0, nnx-1
         do k = 1, nnz-2
            dudz2 = dudz2 + ((uinput(i,j,k+1) - uinput(i,j,k-1)) / (zinput(i,j,k+1) - zinput(i,j,k-1)))**2
            dvdz2 = dvdz2 + ((vinput(i,j,k+1) - vinput(i,j,k-1)) / (zinput(i,j,k+1) - zinput(i,j,k-1)))**2
            dwdz2 = dwdz2 + ((winput(i,j,k+1) - winput(i,j,k-1)) / (zinput(i,j,k+1) - zinput(i,j,k-1)))**2
         enddo
         dudz2 = dudz2 + ((uinput(i,j,1) - uinput(i,j,nnz-1)) / (zinput(i,j,1) - zinput(i,j,nnz-1)))**2
         dudz2 = dudz2 + ((uinput(i,j,0) - uinput(i,j,nnz-2)) / (zinput(i,j,0) - zinput(i,j,nnz-2)))**2
         dvdz2 = dvdz2 + ((vinput(i,j,1) - vinput(i,j,nnz-1)) / (zinput(i,j,1) - zinput(i,j,nnz-1)))**2
         dvdz2 = dvdz2 + ((vinput(i,j,0) - vinput(i,j,nnz-2)) / (zinput(i,j,0) - zinput(i,j,nnz-2)))**2
         dwdz2 = dwdz2 + ((winput(i,j,1) - winput(i,j,nnz-1)) / (zinput(i,j,1) - zinput(i,j,nnz-1)))**2
         dwdz2 = dwdz2 + ((winput(i,j,0) - winput(i,j,nnz-2)) / (zinput(i,j,0) - zinput(i,j,nnz-2)))**2
      enddo
   enddo
   dudz2 = dudz2 / (nnx*nny*nnz)
   dvdz2 = dvdz2 / (nnx*nny*nnz)
   dwdz2 = dwdz2 / (nnx*nny*nnz)

   ! nu = 1.5d-5
   ! dudx2 = 0.0
   ! dvdx2 = 0.0
   ! dwdx2 = 0.0
   ! do i = 1, nnx-2
   !    dudx2 = dudx2 + sum( ((uinput(i+1,:,:) - uinput(i-1,:,:)) / (xinput(i+1,:,:) - xinput(i-1,:,:)))**2 )
   !    dvdx2 = dvdx2 + sum( ((vinput(i+1,:,:) - vinput(i-1,:,:)) / (xinput(i+1,:,:) - xinput(i-1,:,:)))**2 )
   !    dwdx2 = dwdx2 + sum( ((winput(i+1,:,:) - winput(i-1,:,:)) / (xinput(i+1,:,:) - xinput(i-1,:,:)))**2 )
   ! enddo
   ! dudx2 = dudx2 + sum( ((uinput(1,:,:) - uinput(nnx-1,:,:)) / (xinput(1,:,:) - xinput(nnx-1,:,:)))**2 )
   ! dudx2 = dudx2 + sum( ((uinput(0,:,:) - uinput(nnx-2,:,:)) / (xinput(0,:,:) - xinput(nnx-2,:,:)))**2 )
   ! dvdx2 = dvdx2 + sum( ((vinput(1,:,:) - vinput(nnx-1,:,:)) / (xinput(1,:,:) - xinput(nnx-1,:,:)))**2 )
   ! dvdx2 = dvdx2 + sum( ((vinput(0,:,:) - vinput(nnx-2,:,:)) / (xinput(0,:,:) - xinput(nnx-2,:,:)))**2 )
   ! dwdx2 = dwdx2 + sum( ((winput(1,:,:) - winput(nnx-1,:,:)) / (xinput(1,:,:) - xinput(nnx-1,:,:)))**2 )
   ! dwdx2 = dwdx2 + sum( ((winput(0,:,:) - winput(nnx-2,:,:)) / (xinput(0,:,:) - xinput(nnx-2,:,:)))**2 )
   ! dudx2 = dudx2 / (nnx*nny*nnz)
   ! dvdx2 = dvdx2 / (nnx*nny*nnz)
   ! dwdx2 = dwdx2 / (nnx*nny*nnz)

   ! dudy2 = 0.0
   ! dvdy2 = 0.0
   ! dwdy2 = 0.0
   ! do j = 1, nny-2
   !    dudy2 = dudy2 + sum( ((uinput(:,j+1,:) - uinput(:,j-1,:)) / (yinput(:,j+1,:) - yinput(:,j-1,:)))**2 )
   !    dvdy2 = dvdy2 + sum( ((vinput(:,j+1,:) - vinput(:,j-1,:)) / (yinput(:,j+1,:) - yinput(:,j-1,:)))**2 )
   !    dwdy2 = dwdy2 + sum( ((winput(:,j+1,:) - winput(:,j-1,:)) / (yinput(:,j+1,:) - yinput(:,j-1,:)))**2 )
   ! enddo
   ! dudy2 = dudy2 + sum( ((uinput(:,1,:) - uinput(:,nny-1,:)) / (yinput(:,1,:) - yinput(:,nny-1,:)))**2 )
   ! dudy2 = dudy2 + sum( ((uinput(:,0,:) - uinput(:,nny-2,:)) / (yinput(:,0,:) - yinput(:,nny-2,:)))**2 )
   ! dvdy2 = dvdy2 + sum( ((vinput(:,1,:) - vinput(:,nny-1,:)) / (yinput(:,1,:) - yinput(:,nny-1,:)))**2 )
   ! dvdy2 = dvdy2 + sum( ((vinput(:,0,:) - vinput(:,nny-2,:)) / (yinput(:,0,:) - yinput(:,nny-2,:)))**2 )
   ! dwdy2 = dwdy2 + sum( ((winput(:,1,:) - winput(:,nny-1,:)) / (yinput(:,1,:) - yinput(:,nny-1,:)))**2 )
   ! dwdy2 = dwdy2 + sum( ((winput(:,0,:) - winput(:,nny-2,:)) / (yinput(:,0,:) - yinput(:,nny-2,:)))**2 )
   ! dudy2 = dudy2 / (nnx*nny*nnz)
   ! dvdy2 = dvdy2 / (nnx*nny*nnz)
   ! dwdy2 = dwdy2 / (nnx*nny*nnz)

   ! dudz2 = 0.0
   ! dvdz2 = 0.0
   ! dwdz2 = 0.0
   ! do k = 1, nnz-2
   !    dudz2 = dudz2 + sum( ((uinput(:,:,k+1) - uinput(:,:,k-1)) / (zinput(:,:,k+1) - zinput(:,:,k-1)))**2 )
   !    dvdz2 = dvdz2 + sum( ((vinput(:,:,k+1) - vinput(:,:,k-1)) / (zinput(:,:,k+1) - zinput(:,:,k-1)))**2 )
   !    dwdz2 = dwdz2 + sum( ((winput(:,:,k+1) - winput(:,:,k-1)) / (zinput(:,:,k+1) - zinput(:,:,k-1)))**2 )
   ! enddo
   ! dudz2 = dudz2 + sum( ((uinput(:,:,1) - uinput(:,:,nnz-1)) / (zinput(:,:,1) - zinput(:,:,nnz-1)))**2 )
   ! dudz2 = dudz2 + sum( ((uinput(:,:,0) - uinput(:,:,nnz-2)) / (zinput(:,:,0) - zinput(:,:,nnz-2)))**2 )
   ! dvdz2 = dvdz2 + sum( ((vinput(:,:,1) - vinput(:,:,nnz-1)) / (zinput(:,:,1) - zinput(:,:,nnz-1)))**2 )
   ! dvdz2 = dvdz2 + sum( ((vinput(:,:,0) - vinput(:,:,nnz-2)) / (zinput(:,:,0) - zinput(:,:,nnz-2)))**2 )
   ! dwdz2 = dwdz2 + sum( ((winput(:,:,1) - winput(:,:,nnz-1)) / (zinput(:,:,1) - zinput(:,:,nnz-1)))**2 )
   ! dwdz2 = dwdz2 + sum( ((winput(:,:,0) - winput(:,:,nnz-2)) / (zinput(:,:,0) - zinput(:,:,nnz-2)))**2 )
   ! dudz2 = dudz2 / (nnx*nny*nnz)
   ! dvdz2 = dvdz2 / (nnx*nny*nnz)
   ! dwdz2 = dwdz2 / (nnx*nny*nnz)


   ! lint = (u'^2)^(3/2)/epsilon, epsilon = mean(nu*duidxj*duidxj)
   ! lint = (2.0/3.0*tke0)**(3.0/2.0) / ( nu*(dudx2 + dvdx2 + dwdx2 + dudy2 + dvdy2 + dwdy2 + dudz2 + dvdz2 + dwdz2) )
   ! einf = nu*(dudx2 + dvdx2 + dwdx2 + dudy2 + dvdy2 + dwdy2 + dudz2 + dvdz2 + dwdz2)
   uinf = (2.0*TKE0/3.0)**0.5
   ! linf = uinf**3.0 / einf
   ! tinf = linf / uinf
   tinf = lint / uinf
   ! write(*,*) "einf = ", einf
   write(*,*) "uinf = ", uinf
   write(*,*) "lint = ", lint
   write(*,*) "tinf = ", tinf

   ! Aforce = einf / (2.0*TKE0)
   ! Aforce = 1.0*sqrt(2.0/27.0*tke0/lint**2)     ! 1.5 times stuff works when lint = 20% of domain
   write(*,*) "TKE0 = ", tke0


   ! Deallocate
   deallocate(data)

if ((probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then
   ! Read in the scalar data
   nnx = int8(sinres)
   nny = int8(sinres)
   nnz = int8(sinres)

   ! only one scalar for now
   allocate(data(0:nnx*nny*nnz*4-1))
   allocate(sxinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(syinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(szinput(0:nnx-1,0:nny-1,0:nnz-1))
   allocate(sinput(0:nnx-1,0:nny-1,0:nnz-1))

   ! Only read in CSV files for now
   call read_csv(siname, nnx, nny, nnz, data, int8(4))

   sinput = reshape(data(3::4), (/nnx, nny, nnz/))
   sxinput = reshape(data(0::4), (/nnx, nny, nnz/))
   syinput = reshape(data(1::4), (/nnx, nny, nnz/))
   szinput = reshape(data(2::4), (/nnx, nny, nnz/))

   ! Get the xarray table and the differences
   allocate(sxarray(0:nnx-1))
   allocate(sxdiff(0:nnx-1))
   sxarray(0:nnx-1) = sxinput(:,0,0)
   sxdiff(:nnx-2) = sxarray(1:) - sxarray(:nnx-2)
   sxdiff(nnx-1) = sxarray(nnx-1) - sxarray(nnx-2)

   ! Dimensions of the input box
   sLinput = maxval(sxinput(:,0,0)) + HALF*sxdiff(nnx-1)

   ! Deallocate again
   deallocate(data)

   ! Compute varience of scalar field
   Z2mean = 0.d0
   Zmean = 0.d0
   do k = 0, nnz-1
      do j= 0, nny-1
         do i = 0, nnx-1
            if (sinput(i,j,k) .gt. 0) then
               Z2mean = Z2mean + sinput(i,j,k)**2
               Zmean = Zmean + sinput(i,j,k)
            endif
         enddo
      enddo
   enddo

   Z2mean = Z2mean / (nnx*nny*nnz)
   Zmean = Zmean / (nnx*nny*nnz)
   Svar0 = Z2mean - Zmean**2
   Zm0 = Zmean

   write(*,*) "Svar0 = ", Svar0
   write(*,*) "Zmean0 = ", Zm0

endif
endif



      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
#if defined(BL_DO_FLCT)
         write(6,flctin)
#endif
         write(6,control)
      end if

  end subroutine amrex_probinit
  
!------------------------------------

  subroutine setupbc()bind(C, name="setupbc")
  
      use chem_driver, only: P1ATMMKS
      use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
      use probspec_module, only: set_Y_from_Phi
  
      implicit none
      
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <htdata.H>
#include <conp.H>

      REAL_T Patm, pmf_vals(maxspec+3)
      REAL_T Xt(maxspec), Yt(maxspec), loc
      integer zone, n, len, b(SDIM), num_zones_defined
      integer fuelZone, airZone
      data  b / 1, 1, 1 /

      Patm = pamb / 101325.0d0
      num_zones_defined = 0
      len = len_trim(probtype)

      if ((probtype(1:len).eq.BL_PROB_FORCEDHIT) .or. &
          (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then

      fuelZone = BL_FUELPIPE
      airZone = BL_AMBIENT
      num_zones_defined = 2

!     Fuel
      do n = 1,Nspec
         Xt(n) = 0.d0
      end do
      Xt(iN2) = fuel_N2_vol_percent
      Xt(fuelID) = 1.d0 - Xt(iN2)
      
      CALL CKXTY (Xt, Yt)

      do n=1,Nspec
         Y_bc(n-1,fuelZone) = Yt(n)
      end do 
      T_bc(fuelZone) = T_in
      u_bc(fuelZone) = 0.d0
      v_bc(fuelZone) = 0.d0

!     Air
      do n = 1,Nspec
         Xt(n) = 0.d0
      end do
      Xt(oxidID) = 0.21d0
      Xt(iN2) = 1.d0 - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)

      do n=1,Nspec
         Y_bc(n-1,airZone) = Yt(n)
      end do 
      T_bc(airZone) = T_in
      u_bc(airZone) = 0.d0
      v_bc(airZone) = 0.d0
      w_bc(airZone) = 0.d0

      else 

         print *,'Unrecognized probtype'
         call bl_pd_abort()

      endif

      do zone=1,num_zones_defined
!     Set density and hmix consistent with data

         call RHOfromPTY(b, b, &
                             rho_bc(zone), DIMARG(b), DIMARG(b), &
                             T_bc(zone),   DIMARG(b), DIMARG(b), &
                             Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)
         call HMIXfromTY(b, b, &
                             h_bc(zone),   DIMARG(b), DIMARG(b), &
                             T_bc(zone),   DIMARG(b), DIMARG(b), &
                             Y_bc(0,zone), DIMARG(b), DIMARG(b))
      enddo
      bcinit = .true.

  end subroutine setupbc

! ::: -----------------------------------------------------------
      
  integer function getZone(x, y,z)bind(C, name="getZone")

      implicit none
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

      REAL_T x, y, z
      integer len

!#define BL_FUELPIPE 1
!#define BL_COFLOW   2
!#define BL_STICK    3
!#define BL_WALL     4
!#define BL_AMBIENT  5
!#define BL_VOLUME   6
!#define BL_PIPEEND  7

      getZone = BL_VOLUME
      len     = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_FORCEDHIT) &
          .or. (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR) ) then

         getZone = BL_FUELPIPE
         
      else

         call bl_pd_abort('Unrecognized probtype')

      endif
  end function getZone
      
!-----------------------

  subroutine bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,dx,getuvw) &
                        bind(C, name="bcfunction")

      implicit none

      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(SDIM)
      logical getuvw

#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

      integer n, zone, len

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      len = len_trim(probtype)
      
      if ((probtype(1:len).eq.BL_PROB_FORCEDHIT) .or. &
          (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then

         ! No boundary conditions needed yet

      else

         write(6,*) 'No boundary condition for probtype = ', probtype(1:len)
         write(6,*) 'Available: '
         write(6,*) '            ',BL_PROB_FORCEDHIT
         write(6,*) '            ',BL_PROB_FORCEDHITWSCALAR
         call bl_pd_abort(' ')
      endif
  end subroutine bcfunction

! ::: -----------------------------------------------------------
      
  subroutine init_data_new_mech (level,time,lo,hi,nscal, &
                                 vel,scal,DIMS(state),press,DIMS(press), &
                                 delta,xlo,xhi)&
                                 bind(C, name="init_data_new_mech")
          
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
 
#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>
 
      integer i, j, k, n
      REAL_T Patm
 
      Patm = pamb / 101325.0d0
 
      call RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state), &
          Patm)
      call HMIXfromTY(lo,hi, &
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
 
  end subroutine init_data_new_mech

!c ::: -----------------------------------------------------------
!c ::: This routine is called at problem setup time and is used
!c ::: to initialize data on each grid.  The velocity field you
!c ::: provide does not have to be divergence free and the pressure
!c ::: field need not be set.  A subsequent projection iteration
!c ::: will define aa divergence free velocity field along with a
!c ::: consistant pressure.
!c ::: 
!c ::: NOTE:  all arrays have one cell of ghost zones surrounding
!c :::        the grid interior.  Values in these cells need not
!c :::        be set here.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: level     => amr level of grid
!c ::: time      => time at which to init data             
!c ::: lo,hi     => index limits of grid interior (cell centered)
!c ::: nscal     => number of scalar quantities.  You should know
!c :::		   this already!
!c ::: vel      <=  Velocity array
!c ::: scal     <=  Scalar array
!c ::: press    <=  Pressure array
!c ::: delta     => cell size
!c ::: xlo,xhi   => physical locations of lower left and upper
!c :::              right hand corner of grid.  (does not include
!c :::		   ghost region).
!c ::: -----------------------------------------------------------

  subroutine init_data(level,time,lo,hi,nscal, &
     	 	               vel,scal,DIMS(state),press,DIMS(press), &
                       delta,xlo,xhi) &
                       bind(C, name="init_data")
                       
      use chem_driver, only: P1ATMMKS
      use chem_driver_3D, only: RHOfromPTY, HMIXfromTY
      use chem_driver, only: get_spec_name
      use probdata_module
      
      implicit none
      integer    level,nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     time, delta(SDIM)
      REAL_T     vel(DIMV(state),1:SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>
#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
#endif

      integer i, j, k, n, etZone, nPMF, len, airZone, fuelZone
      REAL_T x, y, z, Yl(maxspec), Xl(maxspec), Patm
      REAL_T pmf_vals(maxspec+3), y1, y2
      REAL_T pert,Lx,Ly,FORT_P1ATMMKS
      REAL_T xmod, ymod, zmod
      REAL_T r, s, t, uinterp, vinterp, winterp, sinterp
      REAL_T f0, f1, f2, f3, f4, f5, f6, f7
      integer im, imp1, in, inp1, ip, ipp1
      integer nnx, nny, nnz
      REAL_T ubar, vbar, wbar, lint, eta

      len = len_trim(probtype)
           
      if ((probtype(1:len).eq.BL_PROB_FORCEDHIT) .or. &
          (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then

         airZone = BL_AMBIENT
         fuelZone = BL_FUELPIPE
         ubar = 0.d0
         vbar = 0.d0
         wbar = 0.d0

         ! Interpolate HIT data previously laoded
         do k = lo(3), hi(3)
            z = (float(k)+0.5)*delta(3)+domnlo(3)
            zmod = mod(z,Linput)
            call locate(xarray, inres, zmod, ip)
            ipp1 = mod(ip+1,inres)
            t = (zmod - xarray(ip)) / xdiff(ip)

            do j = lo(2), hi(2)
               y = (float(j)+0.5)*delta(2)+domnlo(2)
               ymod = mod(y,Linput)
               call locate(xarray, inres, ymod, in)
               inp1 = mod(in+1,inres)
               s = (ymod - xarray(in)) / xdiff(in)

               do i = lo(1), hi(1)
                  x = (float(i)+0.5)*delta(1)+domnlo(1)
                  xmod = mod(x,Linput)
                  call locate(xarray, inres, xmod, im)
                  imp1 = mod(im+1,inres)
                  r = (xmod - xarray(im)) / xdiff(im)

               
                  ! Factors for trilinear interpolation
                  f0 = (1-r) * (1-s) * (1-t)
                  f1 = r * (1-s) * (1-t)
                  f2 = (1-r) * s * (1-t)
                  f3 = (1-r) * (1-s) * t
                  f4 = r * (1-s) * t
                  f5 = (1-r) * s * t
                  f6 = r * s * (1-t)
                  f7 = r * s * t
                  uinterp = uinput(im,in,ip)*f0 + uinput(imp1,in,ip)*f1 +   &
                            uinput(im,inp1,ip)*f2 + uinput(im,in,ipp1)*f3 + &
                            uinput(imp1,in,ipp1)*f4 + uinput(im,inp1,ipp1)*f5 + &
                            uinput(imp1,inp1,ip)*f6 + uinput(imp1,inp1,ipp1)*f7
                  vinterp = vinput(im,in,ip)*f0 + vinput(imp1,in,ip)*f1 +   &
                            vinput(im,inp1,ip)*f2 + vinput(im,in,ipp1)*f3 + &
                            vinput(imp1,in,ipp1)*f4 + vinput(im,inp1,ipp1)*f5 + &
                            vinput(imp1,inp1,ip)*f6 + vinput(imp1,inp1,ipp1)*f7
                  winterp = winput(im,in,ip)*f0 + winput(imp1,in,ip)*f1 +   &
                            winput(im,inp1,ip)*f2 + winput(im,in,ipp1)*f3 + &
                            winput(imp1,in,ipp1)*f4 + winput(im,inp1,ipp1)*f5 + &
                            winput(imp1,inp1,ip)*f6 + winput(imp1,inp1,ipp1)*f7


                  do n=1,Nspec
                     scal(i,j,k,FirstSpec-1+n) = Y_bc(n-1,airZone)
                  enddo
                  scal(i,j,k,Temp) = T_bc(airZone)
                  scal(i,j,k,Trac) = 0.d0

                  vel(i,j,k,1) = uinterp
                  vel(i,j,k,2) = vinterp
                  vel(i,j,k,3) = winterp
                  ubar = ubar + uinterp
                  vbar = vbar + vinterp
                  wbar = wbar + winterp
               enddo
            enddo
         enddo

      if (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR) then
         ! Interpolate scalar HIT data previously laoded
         do k = lo(3), hi(3)
            z = (float(k)+0.5)*delta(3)+domnlo(3)
            zmod = mod(z,sLinput)
            call locate(sxarray, sinres, zmod, ip)
            ipp1 = mod(ip+1,sinres)
            t = (zmod - sxarray(ip)) / sxdiff(ip)

            do j = lo(2), hi(2)
               y = (float(j)+0.5)*delta(2)+domnlo(2)
               ymod = mod(y,sLinput)
               call locate(sxarray, sinres, ymod, in)
               inp1 = mod(in+1,sinres)
               s = (ymod - sxarray(in)) / sxdiff(in)

               do i = lo(1), hi(1)
                  x = (float(i)+0.5)*delta(1)+domnlo(1)
                  xmod = mod(x,sLinput)
                  call locate(sxarray, sinres, xmod, im)
                  imp1 = mod(im+1,sinres)
                  r = (xmod - sxarray(im)) / sxdiff(im)

               
                  ! Factors for trilinear interpolation
                  f0 = (1-r) * (1-s) * (1-t)
                  f1 = r * (1-s) * (1-t)
                  f2 = (1-r) * s * (1-t)
                  f3 = (1-r) * (1-s) * t
                  f4 = r * (1-s) * t
                  f5 = (1-r) * s * t
                  f6 = r * s * (1-t)
                  f7 = r * s * t
                  sinterp = sinput(im,in,ip)*f0 + sinput(imp1,in,ip)*f1 +   &
                            sinput(im,inp1,ip)*f2 + sinput(im,in,ipp1)*f3 + &
                            sinput(imp1,in,ipp1)*f4 + sinput(im,inp1,ipp1)*f5 + &
                            sinput(imp1,inp1,ip)*f6 + sinput(imp1,inp1,ipp1)*f7

                  if (sinterp .gt. 0) then
                     if (sinterp .ge. 1.0) then
                        sinterp = 1.0
                     endif
                     do n=1,Nspec
                        scal(i,j,k,FirstSpec-1+n) = sinterp*Y_bc(n-1,fuelZone) + (1.0-sinterp)*Y_bc(n-1,airZone)
                     enddo
                  endif

               enddo
            enddo
         enddo
      endif


      else

         write(6,*) 'No initial condition for probtype = ', probtype(1:len)
         write(6,*) 'Available: '
         write(6,*) '            ',BL_PROB_FORCEDHIT
         write(6,*) '            ',BL_PROB_FORCEDHITWSCALAR
         call bl_pd_abort(' ')

      endif

      Patm = pamb / P1ATMMKS()

      call RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state), &
          Patm)

      call HMIXfromTY(lo,hi, &
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

! ::: ------------------------------------------------------------
! ::: Read a csv file
! :::
! ::: INPUTS/OUTPUTS
! :::
! ::: iname => filename
! ::: nx    => input resolution
! ::: ny    => input resolution
! ::: nz    => input resolution
! ::: data  <= output data
! ::: -----------------------------------------------------------
subroutine read_csv(iname,nx,ny,nz,data,ncomp)

    implicit none

    character(len=255), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny, nz, ncomp
    REAL_T, intent(out) :: data(0:nx*ny*nz*ncomp-1)

    integer :: i
    integer, parameter :: in_unit=1
    integer :: nlines = 0, ios = 0

    ! Get number of lines in file
    open(in_unit,file=trim(iname), access='sequential', form='formatted', status='old', action='read')
    read(in_unit,*) ! skip header
    nlines = 0
    do
       read(in_unit,*,iostat=ios)
       if (ios .ne. 0) exit
       nlines = nlines + 1
    enddo
    ios = 0

    ! Quick sanity check
    if (nlines .ne. nx*ny*nz) then
       write(*,'("Number of lines in the input file (=",I0,") does not ")')nlines
       write(*,'("  match the input resolution (n=",I0,") in the probin file")')nx
       stop 99
    endif

    ! Read the data from the file
    rewind(in_unit)
    read(in_unit,*) ! skip header
    do i = 0, nlines-1
       read(in_unit, *, iostat=ios)data(i*ncomp:(i+1)*ncomp-1)
       if (ios .ne. 0) then
          write(*,*)'Error in CSV input file read. Exiting with read error', ios
          stop 99
       endif
    enddo
    close(in_unit)

  end subroutine read_csv

! ::: -----------------------------------------------------------
  ! ::: Search for the closest index in an array to a given value
  ! ::: using the bisection technique.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: xtable(0:n-1) => array to search in (ascending order)
  ! ::: n             => number of elements in array
  ! ::: x             => x location
  ! ::: idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
  ! ::: -----------------------------------------------------------
  subroutine locate(xtable, n, x, idxlo)

    implicit none

    REAL_T, intent(in) :: xtable(0:n-1)
    integer, intent(in) :: n
    REAL_T, intent(in) :: x
    integer, intent(out) :: idxlo

    ! Local variables
    integer :: idxhi, idxmid
    logical :: notdone

    ! If x is out of bounds, return boundary index
    if (x >= xtable(n-1)) then
       idxlo=n-1
       return
    elseif (x <= xtable(0)) then
       idxlo=0
       return
    endif

    ! Make sure the search array is increasing
    if (xtable(0) > xtable(n-1)) then
       write(*,'("Error in locate: non ascending input search array.")')
       stop 99
    endif

    ! Do the bisection
    idxlo = 0
    idxhi = n-1
    notdone = .true.
    do while (notdone)
       if (idxhi-idxlo <= 1) then
          notdone = .false.
       else
          idxmid = (idxhi+idxlo)/2
          if (x >= xtable(idxmid)) then
             idxlo = idxmid
          else
             idxhi = idxmid
          endif
       endif
    enddo
    return
  end subroutine locate
  
!c ::: -----------------------------------------------------------
!c ::: This routine will zero out diffusivity on portions of the
!c ::: boundary that are inflow, allowing that a "wall" block
!c ::: the complement aperture
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: diff      <=> diffusivity on edges
!c ::: DIMS(diff) => index extent of diff array
!c ::: lo,hi      => region of interest, edge-based
!c ::: domlo,hi   => index extent of problem domain, edge-based
!c ::: dx         => cell spacing
!c ::: problo     => phys loc of lower left corner of prob domain
!c ::: bc         => boundary condition flag (on orient)
!c :::                   in BC_TYPES::physicalBndryTypes
!c ::: idir       => which face, 0=x, 1=y
!c ::: isrz       => 1 if problem is r-z
!c ::: id         => index of state, 0=u
!c ::: ncomp      => components to modify
!c ::: 
!c ::: -----------------------------------------------------------

subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                     dx,problo,bc,idir,isrz,id,ncomp)  &
                     bind(C, name="zero_visc")   
                                          
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
#include <htdata.H>
      integer i, j, k, n, Tid, RHid, YSid, YEid, ys, ye
      integer len
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, yb, yt, yh, z

      len = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
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
         if ((idir.EQ.2) .AND. (lo(3) .LE. domlo(3)) &
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
                  
                  if ( (getZone(xl,yb,z).eq.BL_STICK) .OR. &
                      (getZone(xh,yb,z).eq.BL_STICK) .OR. &
                      (getZone(xr,yb,z).eq.BL_STICK) .OR. &
                      (getZone(xl,yh,z).eq.BL_STICK) .OR. &
                      (getZone(xh,yh,z).eq.BL_STICK) .OR. &
                      (getZone(xr,yh,z).eq.BL_STICK) .OR. &
                      (getZone(xl,yt,z).eq.BL_STICK) .OR. &
                      (getZone(xh,yt,z).eq.BL_STICK) .OR. &
                      (getZone(xr,yt,z).eq.BL_STICK) .OR. &
                      (getZone(xl,yb,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xh,yb,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xr,yb,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xl,yh,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xh,yh,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xr,yh,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xl,yt,z).eq.BL_PIPEEND) .OR. &
                      (getZone(xh,yt,z).eq.BL_PIPEEND) .OR. &
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

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the 
!c ::: density gradient
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: rho       => density array
!c ::: DIMS(rho) => index extent of rho array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------

  subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
                            rho,DIMS(rho),lo,hi,nvar, &
                            domlo,domhi,dx,xlo, &
                		        problo,time,level) &
                             bind(C, name="FORT_DENERROR")
                               
      implicit none
      integer   DIMDEC(rho)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho), nvar)

#include <probdata.H>

      call bl_abort('DENERROR: should no be here')
      
  end subroutine FORT_DENERROR

! ::: -----------------------------------------------------------

  subroutine flame_tracer_error(tag,DIMS(tag),set,clear, &
                                ftrac,DIMS(ftrac),lo,hi,nvar, &
                                domlo,domhi,dx,xlo, &
                 	              problo,time,level) &
                                bind(C, name="flame_tracer_error")
                                
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
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k), &
                      ftrac(i,j,k,1).gt.flametracval)
               enddo
            enddo
         enddo
      endif

  end subroutine flame_tracer_error

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the 
!c ::: density gradient
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: adv       => scalar array
!c ::: DIMS(adv) => index extent of scalar array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------


  subroutine adv_error (tag,DIMS(tag),set,clear, &
                        adv,DIMS(adv),lo,hi,nvar, &
                        domlo,domhi,delta,xlo, &
          			        problo,time,level)&
                        bind(C, name="adv_error")
                  
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    delta(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag)), len
      REAL_T    adv(DIMV(adv),nvar)

#include <probdata.H>

      len = len_trim(probtype)
      
      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
          .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) &
          .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW_FORCED) &
          .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW_FORCED) ) then

         call mv_error(tag,DIMS(tag),set,clear, &
                       adv,DIMS(adv),lo,hi,nvar, &
                       domlo,domhi,delta,xlo, &
                       problo,time,level)

      endif
      
  end subroutine adv_error
  
!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the
!c ::: temperature gradient
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: temp      => density array
!c ::: DIMS(temp)=> index extent of temp array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::              corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------

  subroutine temp_error(tag,DIMS(tag),set,clear, &
                        temperature,DIMS(temp),lo,hi,nvar, &
                        domlo,domhi,dx,xlo, &
                        problo,time,level)&
                        bind(C, name="temp_error")
                               
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

      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),ARG_H3(temp)-hi(3), &
              lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp),lo(3)-ARG_L3(temp))

      if (ng .lt. 1) then
         write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
         call bl_abort(" ")
      endif
!
!     ::::: refine where there is temperature gradient
!
      if (level .lt. max_temp_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
!c                 ax = abs(temperature(i+1,j,k,1) - temperature(i,j,k,1))
!c                 ay = abs(temperature(i,j+1,k,1) - temperature(i,j,k,1))
!c                 az = abs(temperature(i,j,k+1,1) - temperature(i,j,k,1))
!c                 ax = MAX(ax,abs(temperature(i,j,k,1) - temperature(i-1,j,k,1)))
!c                 ay = MAX(ay,abs(temperature(i,j,k,1) - temperature(i,j-1,k,1)))
!c                 az = MAX(az,abs(temperature(i,j,k,1) - temperature(i,j,k-1,1)))
!c                 aerr = max(ax,ay,az)
!c                 tag(i,j,k) = merge(set,tag(i,j,k),aerr.ge.tempgrad)
                  tag(i,j,k) = merge(set,tag(i,j,k),temperature(i,j,k,1).lt.temperr)
               enddo
            enddo
         enddo
      endif

  end subroutine temp_error

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the 
!c ::: magnitude of vorticity
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: vort      => array of vorticity values
!c ::: DIMS(vor) => index extent of vort array
!c ::: nvar      => number of components in vort array (should be 1)
!c ::: lo,hi     => index extent of grid
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------

   subroutine mv_error (tag,DIMS(tag),set,clear, &
                        vort,DIMS(vort),lo,hi,nvar,  &
                        domlo,domhi,dx,xlo, &
     			              problo,time,level)&
                        bind(C, name="mv_error")
                 
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
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k), &
                              ABS(vort(i,j,k,1)).ge.vorterr*2.d0**level)
               enddo
            enddo
         enddo
      end if

  end subroutine mv_error 

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: den      <=  density array
!c ::: DIMS(den) => index extent of den array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of den array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

  subroutine den_fill (den,DIMS(den),domlo,domhi,delta, &
                       xlo,time,bc) &
                       bind(C, name="den_fill")
                       
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif

  end subroutine den_fill

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: adv      <=  advected quantity array
!c ::: DIMS(adv) => index extent of adv array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of adv array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

  subroutine adv_fill (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)&
                           bind(C, name="adv_fill")

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

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: temp     <=  temperature array
!c ::: lo,hi     => index extent of adv array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::              corner of temperature array
!c ::: time      => problem evolution time
!c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine temp_fill  (temp,DIMS(temp),domlo,domhi,delta, &
                             xlo,time,bc) &
                             bind(C, name="temp_fill")

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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif

  end subroutine temp_fill

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: rhoh      <=  rho*h array
!c ::: lo,hi     => index extent of adv array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::              corner of temperature array
!c ::: time      => problem evolution time
!c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

  subroutine rhoh_fill  (rhoh,DIMS(rhoh),domlo,domhi,delta, &
                         xlo,time,bc)&
                         bind(C, name="rhoh_fill")

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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif

  end subroutine rhoh_fill
!
! Fill x, y & z velocity at once.
!
  subroutine vel_fill  (vel,DIMS(vel),domlo,domhi,delta, &
                        xlo,time,bc)&
                        bind(C, name="vel_fill")

      implicit none

      integer DIMDEC(vel), bc(SDIM,2,SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  vel(DIMV(vel),SDIM)

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
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct, &
                       bc(1,1,1), domnlo, domnhi)
         call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, t_flct, &
                       bc(1,1,2), domnlo, domnhi)
         call INFL_FILL(FLCT_ZVEL, DIMS(uflct), wflct, xlo, delta, t_flct, &
                       bc(1,1,3), domnlo, domnhi)
      endif
#endif

      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1), &
                 DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2), &
                 DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
      call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3), &
                 DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))

      if (lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  
                  if ((bc(1,1,1).eq.EXT_DIR) &
                      .or. (bc(1,1,2).eq.EXT_DIR) &
                      .or. (bc(1,1,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

                  if (bc(1,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 1) then
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

                  if ((bc(1,2,1).eq.EXT_DIR) &
                      .or. (bc(1,2,2).eq.EXT_DIR) &
                      .or. (bc(1,2,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

                  if (bc(1,2,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,1) = u + uflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,2) = v + vflct(1,j,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(1,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 1) then
                        vel(i,j,k,3) = w + wflct(1,j,k)*turb_scale
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

                  if ((bc(2,1,1).eq.EXT_DIR) &
                      .or. (bc(2,1,2).eq.EXT_DIR) &
                      .or. (bc(2,1,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

                  if (bc(2,1,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,1,3).eq.EXT_DIR) then
                        vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceLo .and. strmwse_dir .eq. 2) then
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

                  if ((bc(2,2,1).eq.EXT_DIR) &
                      .or. (bc(2,2,2).eq.EXT_DIR) &
                      .or. (bc(2,2,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

                  if (bc(2,2,1).eq.EXT_DIR) then
                        vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,1) = u + uflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
                        vel(i,j,k,2) = v + vflct(i,1,k)*turb_scale
                     endif
#endif
                  endif

                  if (bc(2,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 2) then
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

                  if ((bc(3,1,1).eq.EXT_DIR) &
                      .or. (bc(3,1,2).eq.EXT_DIR) &
                      .or. (bc(3,1,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

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

                  if ((bc(3,2,1).eq.EXT_DIR) &
                      .or. (bc(3,2,2).eq.EXT_DIR) &
                      .or. (bc(3,2,3).eq.EXT_DIR)) then 
                     call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  endif

                  if (bc(3,2,1).eq.EXT_DIR) then
                     vel(i,j,k,1) = u
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,1) = u + uflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
                        vel(i,j,k,2) = v + vflct(i,j,1)*turb_scale
                     endif
#endif
                  endif

                  if (bc(3,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
#if defined(BL_DO_FLCT)
                     if (forceHi .and. strmwse_dir .eq. 3) then
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

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: xvel     <=  x velocity array
!c ::: lo,hi     => index extent of xvel array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta, &
                                xlo,time,bc)&
                                bind(C, name="FORT_XVELFILL")

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
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, time, &
                       bc, domnlo, domnhi)
         allocate(vflct(DIMV(uflct)))
         call INFL_FILL(FLCT_YVEL, DIMS(uflct), vflct, xlo, delta, time, &
                       bc, domnlo, domnhi)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: yvel     <=  y velocity array
!c ::: lo,hi     => index extent of yvel array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,&
                                xlo,time,bc)&
                                bind(C, name="FORT_YVELFILL")

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
         call INFL_FILL(FLCT_XVEL, DIMS(vflct), uflct, xlo, delta, time, &
                       bc, domnlo, domnhi)
         allocate(vflct(DIMV(vflct)))
         call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, time, &
                       bc, domnlo, domnhi)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: zvel     <=  z velocity array
!c ::: lo,hi     => index extent of zvel array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_ZVELFILL (zvel,DIMS(zvel),domlo,domhi,delta, &
                                xlo,time,bc) &
                                bind(C, name="FORT_ZVELFILL")

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
         call INFL_FILL(FLCT_ZVEL, DIMS(wflct), wflct, xlo, delta, time, &
                       bc, domnlo, domnhi)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
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
  
  subroutine all_chem_fill(rhoY,DIMS(rhoY),domlo,domhi,delta, &
                           xlo,time,bc) &
                           bind(C, name="all_chem_fill")

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

!      print *, 'FORT_ALLCHEMFILL: ', domlo,domhi,delta,xlo,time

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
         call filcc (rhoY(lo(1),lo(2),lo(3),n), &
                    DIMS(rhoY),domlo,domhi,delta,xlo,bc)
      end do

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  do n = 1,Nspec
                     rhoY(i,j,k,n) = rho*Yl(n)
                  end do
               enddo
            enddo
         enddo
      endif

  end subroutine all_chem_fill

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: rhoY      <= rho*Y (Y=mass fraction) array
!c ::: lo,hi     => index extent of adv array
!c ::: domlo,hi  => index extent of problem domain
!c ::: delta     => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::              corner of temperature array
!c ::: time      => problem evolution time
!c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: stateID   => id index of state being filled
!c ::: -----------------------------------------------------------

  subroutine chem_fill  (rhoY,DIMS(rhoY),domlo,domhi,delta, &
                         xlo,time,bc,id) &
                         bind(C, name="chem_fill")

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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
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
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.false.)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif

  end subroutine chem_fill

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.  
!c ::: 
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: p        <=  pressure array
!c ::: DIMS(p)   => index extent of p array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!c ::: -----------------------------------------------------------

  subroutine press_fill  (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)&
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
#ifdef GENGETFORCE
  subroutine FORT_MAKEFORCE(time,force,rho, &
                            DIMS(istate),DIMS(state), &
                            dx,xlo,xhi,gravity,scomp,ncomp) &
                            bind(C,name="FORT_MAKEFORCE")

      implicit none

      integer    DIMDEC(state)
      integer    DIMDEC(istate)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
      REAL_T     rho    (DIMV(state))
      REAL_T     gravity

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>

#ifdef DO_LMC_FORCE
#include <forcedata.H>
#endif

!
!     ::::: local variables
!
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1, f2, f3
      REAL_T  twicePi
      REAL_T  force_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc, len
      integer nXvel, nYvel, nZvel, nRho

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      len = len_trim(probtype)

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

      if (scomp.eq.0) then
!     Do velocity forcing
         if ( (probtype(1:len).eq.BL_PROB_FORCEDHIT) ) then
#ifdef DO_LMC_FORCE
!     Homogeneous Isotropic Turbulence
            twicePi=two*Pi
            
!     Adjust time offset if controlled:  FIXME: What is this really doing??
            if (time_offset.gt.(-half)) then
               force_time = time + time_offset
            else
               force_time = time
            endif

            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)

            if (hack_lz.eq.1) then 
               Lz = Lz/two
            endif
         
            Lmin = min(Lx,Ly,Lz)
            kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
            nxmodes = nmodes*int(0.5+Lx/Lmin)
            nymodes = nmodes*int(0.5+Ly/Lmin)
            nzmodes = nmodes*int(0.5+Lz/Lmin)

            xstep = int(Lx/Lmin+0.5)
            ystep = int(Ly/Lmin+0.5)
            zstep = int(Lz/Lmin+0.5)

            if (forcing_twice_wavelength.eq.1) then
               HLx = Lx/two
               HLy = Ly/two
               HLz = Lz/two
            else
               HLx = Lx
               HLy = Ly
               HLz = Lz
            endif

            do k = klo, khi
               z = xlo(3) + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     f1 = zero
                     f2 = zero
                     f3 = zero
                     do kz = mode_start*zstep, nzmodes, zstep
                        kzd = dfloat(kz)
                        do ky = mode_start*ystep, nymodes, ystep
                           kyd = dfloat(ky)
                           do kx = mode_start*xstep, nxmodes, xstep
                              kxd = dfloat(kx)
                              kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                              if (kappa.le.kappaMax) then
                                 xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))

                                 if (div_free_force.eq.1) then

                                    f1 = f1 + xT * &
                                        ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                        *   sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz))  &
                                        *   cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz))  &
                                        - FAY(kx,ky,kz)*twicePi*(kzd/HLz) &
                                        *   sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz))  &
                                        *   cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                                    f2 = f2 + xT * &
                                        ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                        *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                        *   cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz))   &
                                        - FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                        *   cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz))  &
                                        *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                                    f3 = f3 + xT * &
                                        ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                        *   cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz))  &
                                        - FAX(kx,ky,kz)*twicePi*(kyd/HLy) &
                                        *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                        *   cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                                 else

                                    f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                                    f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                                    f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 endif

                              endif
                           enddo
                        enddo
                     enddo
                     do kz = 1, zstep - 1
                        kzd = dfloat(kz)
                        do ky = mode_start, nymodes
                           kyd = dfloat(ky)
                           do kx = mode_start, nxmodes
                              kxd = dfloat(kx)
                              kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                              if (kappa.le.kappaMax) then
                                 xT = cos(FTX(kx,ky,kz)*force_time+TAT(kx,ky,kz))
                                 if (div_free_force.eq.1) then
                                    f1 = f1 + xT * &
                                        ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                        *   sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                        *   cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz))  &
                                        - FAY(kx,ky,kz)*twicePi*(kzd/HLz)  &
                                        *   sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                        *   cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                                    f2 = f2 + xT * &
                                        ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                        *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                        *   cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz))  &
                                        - FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                        *   cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                                    f3 = f3 + xT * &
                                        ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                        *   cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                        *   sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz))  &
                                        - FAX(kx,ky,kz)*twicePi*(kyd/HLy) &
                                        *   sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                        *   cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                        *   sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                                 else
                                    f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                                    f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz)) 

                                    f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                        *                     sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                        *                     cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                     if (use_rho_in_forcing.eq.1) then
                        force(i,j,k,nXvel) = f1*rho(i,j,k)
                        force(i,j,k,nYvel) = f2*rho(i,j,k)
                        force(i,j,k,nZvel) = f3*rho(i,j,k)
                     else
                        force(i,j,k,nXvel) = f1
                        force(i,j,k,nYvel) = f2
                        force(i,j,k,nZvel) = f3
                     endif
                  enddo
               enddo
            enddo
#endif
!     Default to gravity...
         elseif (abs(gravity).gt.0.0001) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*rho(i,j,k)
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
!     Let's add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nZvel) = force(i,j,k,nZvel)  + dV_control*rho(i,j,k)
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

#elif LINEARFORCING

  subroutine FORT_MAKEFORCE(time,force,vel,scal, &
                            DIMS(force),DIMS(vel),DIMS(scal), &
                            Ubar, TKEmean, Zmean, Svarmean, &
                            dx,xlo,xhi,gravity,scomp,ncomp,nscal,getForceVerbose) &
                            bind(C,name="FORT_MAKEFORCE")

      use probdata_module

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(vel)
      integer    DIMDEC(scal)
      integer    scomp, ncomp, nscal
      integer    getForceVerbose
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(force),scomp+1:scomp+ncomp)
      REAL_T     vel    (DIMV(vel),1:SDIM)
      REAL_T     scal   (DIMV(scal),1:nscal)
      REAL_T     Ubar   (1:SDIM+1)
      REAL_T     gravity
      REAL_T     TKEmean
      REAL_T     Svarmean(1:nscal)
      REAL_T     Zmean(1:nscal)


#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>
#include <conp.H>

!
!     ::::: local variables
!
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer len
      integer nXvel, nYvel, nZvel

      len = len_trim(probtype)

      ilo = force_l1
      jlo = force_l2
      klo = force_l3
      ihi = force_h1
      jhi = force_h2
      khi = force_h3


!     Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nZvel = 3
      nRho  = 4

      if (scomp.eq.0) then
!     Do velocity forcing
         if ((probtype(1:len).eq.BL_PROB_FORCEDHIT) .or. &
          (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR)) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     ! Forcing scheme inspired by Bassenne et al. (2016)
                     force(i,j,k,nXvel) = Gain*((tke0-TKEmean)/(2.0*TKEmean))*(vel(i,j,k,nXvel)-Ubar(nXvel))*scal(i,j,k,1)
                     force(i,j,k,nYvel) = Gain*((tke0-TKEmean)/(2.0*TKEmean))*(vel(i,j,k,nYvel)-Ubar(nYvel))*scal(i,j,k,1)
                     force(i,j,k,nZvel) = Gain*((tke0-TKEmean)/(2.0*TKEmean))*(vel(i,j,k,nZvel)-Ubar(nZvel))*scal(i,j,k,1)
                  enddo
               enddo
            enddo
         endif
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
                        if (probtype(1:len).eq.BL_PROB_FORCEDHITWSCALAR) then
                           if ((n-nRho).eq.fuelID) then
                              ! Force the x-component of velocity according to the magnitude of the fluctuations
                              ! Since this will never be negative, shift the function by one standard deviation
                              ! The gain term adjusts the magnitude of the forcing to maintain constant variance in fuel
                              force(i,j,k,n) = sGain*((Svar0-Svarmean(fuelID))/(2.0*Svarmean(fuelID)))*(abs(vel(i,j,k,nXvel)-Ubar(nXvel)) - Ubar(nZvel+1))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

   end subroutine FORT_MAKEFORCE


!c     ::: -----------------------------------------------------------
!c     
!c     This routine averages the mac face velocities for makeforce at half time
!c
   subroutine FORT_AVERAGE_EDGE_STATES(vel,umacx,umacy,umacz, &
                                      DIMS(vel),DIMS(umacx),DIMS(umacy),DIMS(umacz),&
                                      getForceVerbose)&
                                      bind(C, name="FORT_AVERAGE_EDGE_STATES")

   implicit none

   integer    DIMDEC(vel)
   integer    DIMDEC(umacx)
   integer    DIMDEC(umacy)
   integer    DIMDEC(umacz)
   integer    getForceVerbose
   REAL_T     vel  (DIMV(vel),SDIM)
   REAL_T     umacx(DIMV(umacx))
   REAL_T     umacy(DIMV(umacy))
   REAL_T     umacz(DIMV(umacz))

   integer i,j,k,n
   integer ilo,jlo,klo
   integer ihi,jhi,khi

   integer isioproc

   REAL_T  velmin(3)
   REAL_T  velmax(3)

   do n = 1, 3
      velmin(n) = 1.d234
      velmax(n) = -1.d234
   enddo

   ilo = vel_l1
   jlo = vel_l2
   klo = vel_l3
   ihi = vel_h1
   jhi = vel_h2
   khi = vel_h3

   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
            vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
            vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
            do n=1, 3
               velmin(n)=min(velmin(n),vel(i,j,k,n))
               velmax(n)=max(velmax(n),vel(i,j,k,n))
            enddo
         enddo
      enddo
   enddo
   
   if (getForceVerbose.gt.0) then
      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         do n = 1, 3
            write (6,*) "mac velmin (",n,") = ",velmin(n)
            write (6,*) "mac velmax (",n,") = ",velmax(n)
         enddo
      endif
   endif

   end subroutine FORT_AVERAGE_EDGE_STATES


#endif

end module prob_3D_module
