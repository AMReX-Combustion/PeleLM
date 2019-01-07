
#include <MyProb_F.H>

module prob_2D_module

  implicit none

  private
  
  public :: amrex_probinit, setupbc, getZone, bcfunction, init_data_new_mech, init_data, &
            zero_visc, FORT_DENERROR, flame_tracer_error, adv_error, &
            temp_error, mv_error, den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, chem_fill, press_fill, &
            FORT_MAKEFORCE 

contains

! ----- original -----
!      subroutine FORT_PROBINIT (init,name,namlen,problo,probhi)

      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
    
      use chem_driver, only: P1ATMMKS

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

      integer i, istemp
      REAL_T FORT_P1ATMMKS, area

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, &
                        flametracval, probtype, &
      		        max_temp_lev, max_vort_lev, max_trac_lev, &
                        phi_in, T_fu, T_ox, T_air, phi_fu, phi_ox, phi_air, phi_pe, phi_out, pipeTh, pipeBL, &
                        V_in, V_fu, V_ox, V_air, fuel_ox_split, ox_air_split, &
                        standoff, pertmag, nchemdiag, &
                        blobx, bloby, blobr, blobT, blobw, max_nozzle_lev, &
                        refine_nozzle_x, refine_nozzle_y, refine_nozzle_w
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber
#if defined(BL_DO_FLCT)
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
                        forcing_time_scale_min, forcing_time_scale_max, &
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
      max_nozzle_lev = 0
      refine_nozzle_x = 0.d0
      refine_nozzle_y = 0.d0
      refine_nozzle_w = 0.d0
!      pamb = FORT_P1ATMMKS()
      pamb = P1ATMMKS()
      dpdt_factor = 0.3d0
      closed_chamber = 0

      blobr = 0.d0
      blobx = 0.5d0 * (domnhi(1) + domnlo(1))
      bloby = 0.5d0 * (domnhi(2) + domnlo(2))
      blobw = 1.d0
      fuel_ox_split = .001d0
      ox_air_split =  .0125d0
      pipeTh = 0.d0
      pipeBL = 1.d0
      V_fu = -0.2d0
      V_ox = 0.2d0
      V_air = 0.2d0
      T_fu = 300.d0
      T_ox = 300.d0
      T_air = 300.d0
      phi_fu = 0
      phi_ox = 0
      phi_air = 0
      phi_pe = 10.0
      phi_out = 0.0
      phi_in = -1

#if defined(BL_DO_FLCT)
!     add_turb = .FALSE.
      forceInflow = .FALSE.
      numInflPlanesStore = -1
c
!     Don't need to default 'nCompInflow' as it is block data'd to /3/
c
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = BL_SPACEDIM
      flct_file = ''
      turb_scale = 1.d0
      nCompInFlow = BL_SPACEDIM
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

      if (isioproc .eq. 1) write(6,*)"reading fortin"
      read(untin,fortin)
      if (isioproc .eq. 1) write(6,*)"done reading fortin"
      
!     Initialize control variables that depend on fortin variables
      V_in_old = V_fu

      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
      read(untin,heattransin)
 
#if defined(BL_DO_FLCT)
      read(untin,flctin)
      convVel = one
#endif

      if (isioproc .eq. 1) write(6,*)"reading control"
      read(untin,control)
      close(unit=untin)
      if (isioproc .eq. 1) write(6,*)"done reading control"

#if defined(BL_DO_FLCT)
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file.ne.'') then
#define FF_UNIT 20
            ierr = 0
            if (isioproc .eq. 1) then
               write(6,*) '...initializing turbulence'
            end if
            open(FF_UNIT, file=trim(flct_file)//'/HDR',form='formatted',status='old',iostat=ierr)
            if (ierr .ne. 0) then
               call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
            end if
            call RD_SCL_FLCTHD(FF_UNIT,nCompFile,dimFile,probSizeFile,dxFile)
            close(FF_UNIT)
         endif
      endif
#endif
!     Set up boundary functions
      if (isioproc .eq. 1) write(6,*)" setup bc"
      call setupbc()
      if (isioproc .eq. 1) write(6,*)" done setup bc"
      
      area = 1.d0
      do i=1,SDIM-1
         area = area*(domnhi(i)-domnlo(i))
      enddo
      if(istemp.eq.0)then
      scale_control = Y_bc(fuelID-1,BL_FUELPIPE) * rho_bc(BL_FUELPIPE) * area
      else
      scale_control = 1.d0
      endif
      if (isioproc.eq.1) then
         write(6,*)" control setup", area, Y_bc(fuelID-1,BL_FUELPIPE), rho_bc(BL_FUELPIPE)
         write(6,*)" fuelpipe",BL_FUELPIPE
      endif

      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif

#ifdef DO_LMC_FORCE
      if ( (probtype.eq.BL_PROB_JET_DIFFUSION) ) then

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

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
#if defined(BL_DO_FLCT)
         write(6,flctin)
#endif
         write(6,control)
      end if

! ----- original -----
!      end

  end subroutine amrex_probinit

!------------------------------------

! ----- original
!      subroutine setupbc()

  subroutine setupbc()bind(C, name="setupbc")

    use chem_driver, only: P1ATMMKS
    use chem_driver_2D, only: RHOfromPTY, HMIXfromTY
    use probspec_module, only: set_Y_from_Phi
    use chem_driver, only: get_spec_name

      implicit none
#include <cdwrk.H>
#include <conp.H>
#include <bc.H>
#include <probdata.H>
#include <htdata.H>
      
!      REAL_T Patm, FORT_P1ATMMKS, pmf_vals(maxspec+3), a
      REAL_T Patm, pmf_vals(maxspec+3), a
      REAL_T Xt(maxspec), Yt(maxspec), loc
      REAL_T rho	! added
!      integer zone, n, getZone, fuelZone, airZone, oxZone, volZone, peZone, ofZone, region, len
      integer zone, n, fuelZone, airZone, oxZone, volZone, peZone, ofZone, region, len
      integer b(SDIM)
      integer num_zones_defined
      integer iO2,iH2,iCH4, iCO2
      character*(maxspnml) name
      data  b / 1, 1 /
      
!      Patm = pamb / FORT_P1ATMMKS()
    Patm = pamb / P1ATMMKS()
      num_zones_defined = 0
      len = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW)  &
           .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE)  &
           .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
         
         region = BL_INTERIOR
         num_zones_defined = 2

         fuelZone = getZone(0.5*(domnlo(1)+domnhi(1)), domnlo(2))

         if (phi_in.gt.zero) then

!      Take fuel mixture from prob data
            call set_Y_from_Phi(phi_in,Yt)
            do n=1,Nspec
               Y_bc(n-1,fuelZone) = Yt(n)
            end do
            T_bc(fuelZone) = T_fu

         else 

!      Take fuel mixture from pmf file
            loc = (domnlo(2)-standoff)*100.d0
            call pmf(loc,loc,pmf_vals,n)
            if (n.ne.Nspec+3) then
               call bl_pd_abort('setupbc: n(pmf) .ne. Nspec+3')
            endif

!      If starting from pmf file generated by 1D code, you get directly mass fractions (Y); otherwise you get mole fractions (X), then convert them into Y
           if (.false.) then
!            rho = 0
!            do n = 1,Nspec
!               rho = rho + pmf_vals(3+n)	! get density
!            end do
              do n = 1,Nspec
                 Yt(n) = pmf_vals(3+n)	! get mass fractions
              end do
           else            
              do n = 1,Nspec			! original
                 Xt(n) = pmf_vals(3+n)
              end do
              CALL CKXTY (Xt, IWRK, RWRK, Yt)
           endif 

            do n=1,Nspec
               Y_bc(n-1,fuelZone) = Yt(n)
            end do
            T_bc(fuelZone) = pmf_vals(1)
         endif

         phi_bc(fuelZone) = phi_fu

         u_bc(fuelZone) = zero
         if (V_fu .lt. 0) then
            v_bc(fuelZone) = pmf_vals(2)*1.d-2
         else
            v_bc(fuelZone) = V_fu
         endif

         ofZone = getZone(0.5*(domnlo(1)+domnhi(1)), domnhi(2))
         do n=1,Nspec
            Y_bc(n-1,ofZone) = Y_bc(n-1,fuelZone)
         end do
         T_bc(ofZone) = T_bc(fuelZone)
         u_bc(ofZone) = u_bc(fuelZone)
         v_bc(ofZone) = v_bc(fuelZone)
         phi_bc(ofZone) = phi_out

      else if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then

         do n=1,Nspec

            call get_spec_name(name,n)
            if (name .eq. 'N2' ) iN2 = n
            if (name .eq. 'O2' ) iO2 = n
!            if (name .eq. 'CO2' ) iCO2 = n
            if (name .eq. 'CH4' ) iCH4 = n

         enddo

!      A diffusion flame

         fuelZone = BL_FUELPIPE
         oxZone   = BL_OXIDIZER
         airZone  = BL_AIR
         volZone  = BL_VOLUME
         peZone   = BL_PIPEEND
         ofZone   = BL_OUTFLOW
         num_zones_defined = 6

!      Fuel
         do n = 1,Nspec
            Yt(n) = 0.d0
         end do 

         Yt(iCH4) = 1.d0
         Yt(iO2) = 0.d0
         Yt(iN2) = 1.d0 -  Yt(iCH4) -  Yt(iO2)

         do n=1,Nspec
            Y_bc(n-1,fuelZone) = Yt(n)
         end do
         T_bc(fuelZone) = T_fu
         u_bc(fuelZone) = 0.d0
         v_bc(fuelZone) = V_fu
	 phi_bc(fuelZone) = phi_fu

!      Oxidizer
         do n = 1,Nspec
            Yt(n) = 0.d0
         end do 

!         Yt(iO2) = 0.2395d0
!         Yt(iCO2) = 0.0005d0
!         Yt(iN2) = 1.d0 -  Yt(iO2) -  Yt(iCO2)
         Yt(iO2) = 0.233d0
         Yt(iN2) = 1.d0 -  Yt(iO2)

         do n=1,Nspec
            Y_bc(n-1,oxZone) = Yt(n)
         end do
         T_bc(oxZone) = T_ox
         u_bc(oxZone) = 0.d0
         v_bc(oxZone) = V_ox
	 phi_bc(oxZone) = phi_ox

!      Air
         do n=1,Nspec
            Xt(n) = zero
         enddo
         Xt(iO2) = 0.21d0
         Xt(iN2) = 1.d0 - Xt(iO2)

         CALL CKXTY (Xt, IWRK, RWRK, Yt)         
         do n=1,Nspec
            Y_bc(n-1,airZone) = Yt(n)
         end do

         T_bc(airZone) = T_air
         u_bc(airZone) = 0.d0
         v_bc(airZone) = V_air
	 phi_bc(airZone) = phi_air

!      Pipeend (as oxidizer,except for  phiV)
         do n=1,Nspec
            Y_bc(n-1,peZone) = Y_bc(n-1,oxZone)
         end do
         
         T_bc(peZone) = T_bc(oxZone)
         u_bc(peZone) = u_bc(oxZone)
         v_bc(peZone) = 0.d0
	 phi_bc(peZone) = phi_pe

!      Volume (as air)
         do n=1,Nspec
            Y_bc(n-1,volZone) = Y_bc(n-1,airZone)
         end do
         
         T_bc(volZone) = T_bc(airZone)
         u_bc(volZone) = u_bc(airZone)
         v_bc(volZone) = v_bc(airZone)
	 phi_bc(volZone) = phi_bc(airZone)

!	 print*, "Y_BC zone: ", Y_BC(:,volZone), sum(Y_BC(:,volZone))

!      Outflow (as air, except phi)
         do n=1,Nspec
            Y_bc(n-1,ofZone) = Y_bc(n-1,airZone)
         end do
         
         T_bc(ofZone) = T_bc(airZone)
         u_bc(ofZone) = u_bc(airZone)
         v_bc(ofZone) = v_bc(airZone)
	 phi_bc(ofZone) = phi_out

      else
         
         call bl_pd_abort('Unknown probtype')

      endif

      do zone=1,num_zones_defined
!      Set density and hmix consistent with data

! ----- original -----
!         call FORT_RHOfromPTY(b, b, &
!                              rho_bc(zone), DIMARG(b), DIMARG(b), &
!                              T_bc(zone),   DIMARG(b), DIMARG(b), &
!                              Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)

!         call FORT_HMIXfromTY(b, b, &
!                              h_bc(zone),   DIMARG(b), DIMARG(b), &
!                              T_bc(zone),   DIMARG(b), DIMARG(b), &
!                              Y_bc(0,zone), DIMARG(b), DIMARG(b))

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
      end

! ::: -----------------------------------------------------------

! ----- original -----      
!      integer function getZone(x, y)

  integer function getZone(x, y)bind(C, name="getZone")

      implicit none

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

      REAL_T x, y
      integer len

      getZone = BL_VOLUME
      len     = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) &
            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
         if (y.le.domnlo(2)) then
            getZone = BL_FUELPIPE
         else
            getZone = BL_OUTFLOW
         endif
      else if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then
         if (y.le.domnlo(2)) then

!            if (ABS(x).le.fuel_ox_split) then
!               getZone = BL_FUELPIPE
!            else if (ABS(x).le.fuel_ox_split+pipeTh) then
!               getZone = BL_PIPEEND

            if (ABS(x).le.fuel_ox_split) then
               getZone = BL_PIPEEND
            else if (ABS(x).le.fuel_ox_split+pipeTh) then
               getZone = BL_FUELPIPE

            else if (ABS(x).le.ox_air_split) then
               getZone = BL_OXIDIZER
            else
               getZone = BL_AIR
            endif
         else if (y.ge.domnhi(2)) then
            getZone = BL_OUTFLOW
         endif
      else
         call bl_pd_abort('Unrecognized probtype')
      endif

! ----- original -----
!      end

  end function getZone
      
! ::: -----------------------------------------------------------

! ----- original -----      
!      subroutine bcfunction(RegionID,x,y,time,u,v,rho,Yl,T,phiV,h,dx,getuv)

  subroutine bcfunction(RegionID,x,y,time,u,v,rho,Yl,T,phiV,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use chem_driver, only: P1ATMMKS
      use chem_driver_2D, only: RHOfromPTY, HMIXfromTY

      implicit none

      integer RegionID
      REAL_T x, y, time, u, v, rho, Yl(0:*), T, phiV, h, dx(SDIM), r
      logical getuv
      integer b(SDIM)
      data  b / 1, 1 /

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>
!      integer n, getZone, zone, len
!      REAL_T eta, eta1, xmid, etamax, Patm, FORT_P1ATMMKS
      integer n, zone, len
      REAL_T eta1, eta1a, eta2, eta3, eta4, xmid, etamax, Patm
      REAL_T h_fu(0:maxspec-1), h_ox(0:maxspec-1), hmix
         
      REAL_T Wf, Wa, Wm, mf, Yf

      REAL_T rho_V(1), T_V(1), h_V(1)			! it doesn't like usual rho, T, h

      REAL_T,  parameter :: HtoTerrMAX = BL_REAL_E(7.8,-12)
      integer, parameter :: HtoTiterMAX = 20
      REAL_T res(0:HtoTiterMAX-1)
      integer Niter
      integer iO2,iCH4, zoneL, zoneR, zoneT
      character*(maxspnml) name

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      len = len_trim(probtype)
      
      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) &
            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) &
           ) then

         zone = getZone(x,y)
         rho = rho_bc(zone)
         do n = 0, Nspec-1
            Yl(n) = Y_bc(n,zone)
         end do
         T = T_bc(zone)
         phiV = phi_bc(zone)
!       phiV = MIN(1.d0, (time - 0.023100422277313303d0)/1.d-4) * phi_bc(zone)
!       print *,"*********************************** ", time, MIN(1.d0, (time - 0.023100422277313303d0)/4.8d-4)
         h = h_bc(zone)
      
         if (getuv .eqv. .TRUE.) then

            u = zero
            if  ((probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) &
                  .or. (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) ) then               
               v =  V_fu + (time-tbase_control)*dV_control
            else if (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) then
               v = v_bc(zone)
            endif
         endif

      elseif (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then

         zone = getZone(x,y)
         if (zone .eq. BL_OUTFLOW) then
            rho = rho_bc(zone)
            do n = 0, Nspec-1
               Yl(n) = Y_bc(n,zone)
            end do
            T = T_bc(zone)
            phiV = phi_bc(zone)
            h = h_bc(zone)
            if (getuv .eqv. .TRUE.) then
               u = zero
               v = v_bc(zone)
            endif
         else
            zoneL = BL_FUELPIPE
            zoneR = BL_OXIDIZER
            zoneT = BL_VOLUME         
            eta1 = 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-fuel_ox_split)/blobw))
            eta1a = 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-(fuel_ox_split+pipeTh))/blobw))
            do n = 0, Nspec-1
!               Yl(n) = Y_bc(n,zoneL)*eta1 + Y_bc(n,zoneR)*(1.d0-eta1)
               Yl(n) =  Y_bc(n,zoneR)*((1.d0-eta1a)+eta1) + Y_bc(n,zoneL)*((1.d0-eta1)*eta1a)
            end do
            T = T_bc(zoneL)*eta1 + T_bc(zoneR)*(1.d0-eta1)
            phiV = phi_bc(zoneL)*eta1 + phi_bc(zoneR)*(1.d0-eta1)


            eta2 = 0.5d0*(1.d0 - TANH(2.d0*(y-bloby)/blobw))
            do n = 0, Nspec-1
               Yl(n) = Yl(n)*eta2 + Y_bc(n,zoneT)*(1.d0-eta2)
            end do
!            T = T*eta2 + T_bc(zoneT)*(1.d0-eta2)
            T = T_bc(zoneT)

!! ----- original -----
!!            Patm = pamb / FORT_P1ATMMKS()
!!!            call FORT_RHOfromPTY(b, b, &
!!                  rho, DIMARG(b), DIMARG(b), &
!!                  T,   DIMARG(b), DIMARG(b), &
!!                  Yl,  DIMARG(b), DIMARG(b), Patm)

!!            call FORT_HMIXfromTY(b, b, &
!!                  h,   DIMARG(b), DIMARG(b), &
!!                  T,   DIMARG(b), DIMARG(b), &
!!                  Yl,  DIMARG(b), DIMARG(b))

    	    Patm = pamb / P1ATMMKS()
	    rho_V(1) = rho
	    T_V(1) = T
	    h_V(1) = h
      	    call RHOfromPTY(b,b, &
          	 rho_V,  1, 1, 1, 1, &
          	 T_V,     1, 1, 1, 1, &
          	 Yl,1, 1, 1, 1, &
          	 Patm)

      	    call HMIXfromTY(b,b, &
          	 h_V,     1, 1, 1, 1, &
          	 T_V,     1, 1, 1, 1, &
          	 Yl, 1, 1, 1, 1) 

	    rho = rho_V(1)
	    T = T_V(1)
	    h = h_V(1)

!	    print*, Yl(0:Nspec-1)
!	    print*, rho, T, h

            if (getuv .eqv. .TRUE.) then
               u = zero
               zone = getZone(x,y)
               if (zone.eq.BL_FUELPIPE) then
                  zoneL = BL_PIPEEND
                  zoneR = BL_FUELPIPE
                  r = fuel_ox_split - ABS(x)
                  eta3 = 0.5d0*(1.d0 - TANH(2.d0*r/pipeBL))
                  v = v_bc(zoneL)*eta3 + v_bc(zoneR)*(1.d0-eta3)
                  if (eta3.lt.0.d0 .or. eta3.gt.1.d0) then
                     call bl_abort()
                  endif
               else if (zone.eq.BL_OXIDIZER) then
                  zoneL = BL_PIPEEND
                  zoneR = BL_OXIDIZER
                  r = ABS(x) - (fuel_ox_split+pipeTh)
                  eta4 = 0.5d0*(1.d0 - TANH(2.d0*r/pipeBL))
                  v = v_bc(zoneL)*eta4 + v_bc(zoneR)*(1.d0-eta4)
                  if (eta4.lt.0.d0 .or. eta4.gt.1.d0) then
                     call bl_abort()
                  endif
               else
                  v = v_bc(zone)
               endif
            endif
         endif
      else
         write(6,*) 'No boundary condition for probtype = ', probtype(1:len)
         write(6,*) 'Available: '
         write(6,*) '            ',BL_PROB_PREMIXED_FIXED_INFLOW
         write(6,*) '            ',BL_PROB_PREMIXED_CONTROLLED_INFLOW
         write(6,*) '            ',BL_PROB_JET_DIFFUSION
         call bl_pd_abort(' ')
      endif

! ----- original -----
!      end

  end subroutine bcfunction

! ::: -----------------------------------------------------------

! ----- original -----      
!      subroutine FORT_INITDATANEWMECH(level,time,lo,hi,nscal, &
!            vel,scal,DIMS(state),press,DIMS(press), &
!            delta,xlo,xhi)
!      implicit none
!      integer  level, nscal
!      integer  lo(SDIM), hi(SDIM)
!      integer  DIMDEC(state)
!      integer  DIMDEC(press)
!      REAL_T   xlo(SDIM), xhi(SDIM)
!      REAL_T   time, delta(SDIM)
!      REAL_T   vel(DIMV(state),SDIM)
!      REAL_T   scal(DIMV(state),nscal)
!      REAL_T   press(DIMV(press))
! 
!#include <cdwrk.H>
!#include <htdata.H>
!#include <bc.H>
!#include <probdata.H>
! 
!      integer i, j, n
!      REAL_T Patm, FORT_P1ATMMKS
! 
!      do j = lo(2), hi(2)
!         do i = lo(1), hi(1)
!            scal(i,j,Trac) = zero
!         end do
!      end do
! 
!      Patm = pamb / FORT_P1ATMMKS()
!      call FORT_RHOfromPTY(lo,hi, &
!            scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
!            Patm)
!      call FORT_HMIXfromTY(lo,hi, &
!            scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))
! 
!      do j = lo(2), hi(2)
!         do i = lo(1), hi(1)
!            do n = 0,Nspec-1
!               scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
!            enddo
!            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
!         enddo
!      enddo
! 
!      end

  subroutine init_data_new_mech (level,time,lo,hi,nscal, &
          vel,scal,DIMS(state),press,DIMS(press), &
          delta,xlo,xhi)&
          bind(C, name="init_data_new_mech")
          
      use chem_driver, only: P1ATMMKS
      use chem_driver_2D, only: RHOfromPTY, HMIXfromTY
      
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
 
      integer i, j, n
      REAL_T Patm
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            scal(i,j,Trac) = zero
         end do
      end do
 
      Patm = pamb / P1ATMMKS()
      call RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
          Patm)
      call HMIXfromTY(lo,hi, &
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

! ----- original -----
!      subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
!       	 	               vel,scal,DIMS(state),press,DIMS(press), &
!                                delta,xlo,xhi)

  subroutine init_data(level,time,lo,hi,nscal, &
     	 	                   vel,scal,DIMS(state),press,DIMS(press), &
                           delta,xlo,xhi) &
                           bind(C, name="init_data")
                              
      use chem_driver, only: P1ATMMKS
      use chem_driver_2D, only: RHOfromPTY, HMIXfromTY
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
      integer tmpi, nPMF

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

!      integer i, j, n, airZone, fuelZone, getZone, zone, len
      integer i, j, n, airZone, fuelZone, zone, len
      REAL_T x, y, r, Yl(maxspec), Xl(maxspec), Patm
      REAL_T pmf_vals(maxspec+3), y1, y2, dx
!      REAL_T pert,Lx,FORT_P1ATMMKS,eta,u,v,rho,T,h
      REAL_T pert,Lx,eta,u,v,rho,T,h
      REAL_T phiV

      if (iN2.lt.1 .or. iN2.gt.Nspec) then
         call bl_pd_abort('N2 id not sest prior to calling INITDATA')
      endif

      len = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then

         do j = lo(2), hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5d0)*delta(1)+domnlo(1)
               
               pert = 0.d0
               if (pertmag .gt. 0.d0) then
                  Lx = domnhi(1) - domnlo(1)
                  pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx) &
                       + 1.023 * sin(2*Pi*2*(x-.004598)/Lx) &
                           + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) &
                               + 1.017 * sin(2*Pi*5*(x-.0033)/Lx) &
                                    + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
               endif
                  
               y1 = (y - standoff - 0.5d0*delta(2) + pert)*100.d0
               y2 = (y - standoff + 0.5d0*delta(2) + pert)*100.d0
               
!               call pmf(y1,y2,pmf_vals,nPMF)    
               call pmf(0.5*(y1+y2),0.5*(y1+y2),pmf_vals,nPMF)                          
               if (nPMF.ne.Nspec+3) then
                  call bl_abort('INITDATA: n .ne. Nspec+3')
               endif 

!      If starting from pmf file generated by 1D code, you get directly mass fractions (Y); otherwise you get mole fractions (X), then convert them into Y
               if (.false.) then
!                rho = 0
!                do n = 1,Nspec
!                  rho = rho + pmf_vals(3+n)	! get density
!                end do
                  do n = 1,Nspec
                     Yl(n) = pmf_vals(3+n)	! get mass fractions
                  end do
               else            
                  do n = 1,Nspec			! original
                     Xl(n) = pmf_vals(3+n)
                  end do
                  CALL CKXTY (Xl, IWRK, RWRK, Yl)
               endif 
               
               scal(i,j,Temp) = pmf_vals(1)

               do n = 1,Nspec
                  scal(i,j,FirstSpec+n-1) = Yl(n)
               end do

               scal(i,j,Trac) = 0.d0

               vel(i,j,1) = 0.d0
               vel(i,j,2) = pmf_vals(2)*1.d-2

            end do
         end do

      else if ((probtype(1:len).eq.BL_PROB_JET_DIFFUSION)) then

         do j = lo(2), hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)

               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.true.)

               do n=1,Nspec
                  scal(i,j,FirstSpec-1+n) = Yl(n)
               enddo
               scal(i,j,Temp) = T

               if (blobr.gt.0.d0) then
                  r = SQRT((x-blobx)**2 + (y-bloby)**2)
                  eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/blobw))
                  scal(i,j,Temp) = blobT*eta + scal(i,j,Temp)*(1.d0-eta)
               endif

               vel(i,j,1) = u
               vel(i,j,2) = v
               scal(i,j,Trac) = 0.d0
               
!		if(i.eq.lo(1) .and. j.eq.lo(2)) then
!	          print*, i, j, scal(i,j,:)
!		  print*, getZone(x,y), BL_VOLUME
!		endif

            enddo
         enddo

      endif

! ----- original -----
!      Patm = pamb / FORT_P1ATMMKS()

!!      call FORT_RHOfromPTY(lo,hi, &
!           scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
!            Patm)

!      call FORT_HMIXfromTY(lo,hi, &
!            scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
!            scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))

      Patm = pamb / P1ATMMKS()

      call RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
          Patm)

      call HMIXfromTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state)) 

!      Update typical values
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               typVal_Y(n+1) = MAX(typVal_Y(n+1),scal(i,j,firstSpec+n))
            enddo
            typVal_Density = MAX(scal(i,j,Density),typVal_Density)
            typVal_Temp    = MAX(scal(i,j,Temp),   typVal_Temp)
            typVal_Trac    = MAX(scal(i,j,Trac),   typVal_Trac)
            typVal_RhoH = MAX(ABS(scal(i,j,RhoH)*scal(i,j,Density)),typVal_RhoH)
            do n = 1,BL_SPACEDIM
               typVal_Vel  = MAX(ABS(vel(i,j,n)),typVal_Vel)
            enddo
         enddo
      enddo
      do n = 1,Nspec
         typVal_Y(n) = MIN(MAX(typVal_Y(n),typVal_YMIN),typVal_YMAX)
      enddo

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo

! ----- original -----
!      end

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

! ----- original -----
!      subroutine FORT_ZEROVISC(diff,DIMS(diff),lo,hi,domlo,domhi, &
!                                dx,problo,bc,idir,isrz,id,ncomp)

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                           dx,problo,bc,idir,isrz,id,ncomp) &
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
      integer i, j, n, Tid, RHid, YSid, YEid, ys, ye
!      integer getZone, len
      integer len
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, yb, yt, yh, y

      len = len_trim(probtype)

      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then
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
!      Do species, Temp, rhoH
!      
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2)) &
                  .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
            y = float(j)*dx(2)+domnlo(2)
            j = lo(2)
            do i = lo(1), hi(1)
               
               xl = float(i)*dx(1)+domnlo(1) 
               xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
               xh = 0.5d0*(xl+xr)
                  
               if ( (getZone(xl,y).eq.BL_PIPEEND) .OR. &
                     (getZone(xh,y).eq.BL_PIPEEND) .OR. &
                     (getZone(xr,y).eq.BL_PIPEEND)  ) then

!               if (do_T)  diff(i,j,Tid ) = 0.d0
!               if (do_RH) diff(i,j,RHid) = 0.d0
                  if (do_Y) then
                     do n=ys,ye
                        diff(i,j,n) = 0.d0
                     enddo
                  endif
                     
               endif
            end do
         endif
      end if

! ----- original -----
!      end

  end subroutine zero_visc

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
      subroutine FORT_ZEROPHID(diff,DIMS(diff),lo,hi,domlo,domhi, &
                                dx,problo,bc,idir,isrz,id,ncomp)
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
      integer i, j, n
!      integer getZone, len
      integer len
      REAL_T x, y
      REAL_T etaFP, etaPO, eta

      len = len_trim(probtype)

      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))) then
            y = float(j)*dx(2)+domnlo(2)
            j = lo(2)
            do i = lo(1), hi(1)
               x = (i+0.5d0)*dx(1)+domnlo(1)
               etaFP = 1.d0 - 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-fuel_ox_split)/blobw))
               etaPO = 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-(fuel_ox_split+pipeTh))/blobw))
               eta = MIN(etaFP,etaPO)
               diff(i,j,1) = diff(i,j,1) * eta
            end do
         endif
      end if
      end

      subroutine FORT_MASKSOLIDWALL(efab,DIMS(efab),lo,hi,domlo,domhi, &
                                     dx,problo,bc,idir,isrz,ncomp,mask_is_zero)
      implicit none
      integer DIMDEC(efab)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(2*SDIM)
      integer idir, isrz, id, ncomp, mask_is_zero
      REAL_T  efab(DIMV(efab),ncomp)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM)

#include <probdata.H>
#include <cdwrk.H>
#include <htdata.H>
      integer i, j, n
!      integer getZone, len
      integer len
      REAL_T x, y
      REAL_T etaFP, etaPO, eta

      len = len_trim(probtype)

      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))) then
            y = float(j)*dx(2)+domnlo(2)
            j = lo(2)
            do i = lo(1), hi(1)
               x = (i+0.5d0)*dx(1)+domnlo(1)
               etaFP = 1.d0 - 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-fuel_ox_split)/blobw))
               etaPO = 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-(fuel_ox_split+pipeTh))/blobw))
               eta = MAX(0.d0, MIN(1.d0,etaFP,etaPO))

               if (mask_is_zero.eq.1) then
                  eta = 1.d0 - eta
               endif
               do n=1,ncomp
                  efab(i,j,n) = efab(i,j,n) * eta
               end do
            end do
         endif
      end if
      end

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
! ::: rho       => density array
! ::: DIMS(rho) => index extent of rho array
! ::: lo,hi     => index extent of grid
! ::: nvar      => number of components in rho array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: -----------------------------------------------------------

! ----- original
!      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
!                                 rho,DIMS(rho),lo,hi,nvar, &
!                                 domlo,domhi,dx,xlo, &
!       			        problo,time,level)

  subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
                               rho,DIMS(rho),lo,hi,nvar,  &
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

! ----- original      
!      end

  end subroutine FORT_DENERROR


! ::: -----------------------------------------------------------

! ----- original -----
!      subroutine FORT_FLAMETRACERROR (tag,DIMS(tag),set,clear, &
!                                       ftrac,DIMS(ftrac),lo,hi,nvar, &
!                                       domlo,domhi,dx,xlo, &
!       			              problo,time,level)

  subroutine flame_tracer_error (tag,DIMS(tag),set,clear, &
                                 ftrac,DIMS(ftrac),lo,hi,nvar,  &
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

      integer   i, j

#include <probdata.H>

      if (level.lt.max_trac_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
!               tag(i,j) = cvmgt(set,tag(i,j), &
               tag(i,j) = merge(set,tag(i,j), &
                     ftrac(i,j,1).gt.flametracval)
            enddo
         enddo
      endif

! ----- original -----
!      end

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

! ----- original -----
!      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear, &
!                                 adv,DIMS(adv),lo,hi,nvar, &
!                                 domlo,domhi,delta,xlo, &
!       			        problo,time,level)
!      implicit none
!      integer   DIMDEC(tag)
!      integer   DIMDEC(adv)
!      integer   nvar, set, clear, level
!      integer   domlo(SDIM), domhi(SDIM)
!      integer   lo(SDIM), hi(SDIM)
!      REAL_T    delta(SDIM), xlo(SDIM), problo(SDIM), time
!      integer   tag(DIMV(tag)), len
!      REAL_T    adv(DIMV(adv),nvar)

!#include <probdata.H>

!      len = len_trim(probtype)

!      if ( (probtype(1:len).eq.BL_PROB_PREMIXED_FIXED_INFLOW) &
!            .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
!         call FORT_MVERROR(tag,DIMS(tag),set,clear, &
!                            adv,DIMS(adv),lo,hi,nvar, &
!                            domlo,domhi,delta,xlo, &
!                            problo,time,level)
!      endif
!      
!      end

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
          .or. (probtype(1:len).eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then
         call mv_error(tag,DIMS(tag),set,clear, &
                          adv,DIMS(adv),lo,hi,nvar, &
                          domlo,domhi,delta,xlo, &
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

! ----- original -----
!      subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear, &
!                                 temperature,DIMS(temp),lo,hi,nvar, &
!                                 domlo,domhi,dx,xlo, &
!                                 problo,time,level)

  subroutine temp_error (tag,DIMS(tag),set,clear, &
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

      REAL_T    ax, ay, aerr, x, y
      integer   i, j, ng

#include <probdata.H>

      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2), &
                lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp))

      if (ng .lt. 1) then
         write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
         call bl_abort(" ")
      endif
! 
!   refine where there is temperature gradient
! 
      if (level .lt. max_temp_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
! $$$              ax = abs(temperature(i+1,j,1) - temperature(i,j,1))
! $$$              ay = abs(temperature(i,j+1,1) - temperature(i,j,1))
! $$$              ax = MAX(ax,abs(temperature(i,j,1) - temperature(i-1,j,1)))
! $$$              ay = MAX(ay,abs(temperature(i,j,1) - temperature(i,j-1,1)))
! $$$              aerr = max(ax,ay)
! $$$              tag(i,j) = cvmgt(set,tag(i,j),aerr.ge.tempgrad)
!               tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).gt.temperr)
               tag(i,j) = merge(set,tag(i,j),temperature(i,j,1).gt.temperr)
            enddo
         enddo
      endif

      if (level .lt. max_nozzle_lev) then
         do j = lo(2), hi(2)
            y = (float(j)+0.5d0)*dx(2)+domnlo(2) 
            do i = lo(1), hi(1)
               x = (float(i)+0.5d0)*dx(1)+domnlo(1) 
               if (ABS(x-refine_nozzle_x).le.refine_nozzle_w  &
                     .and. y.le.refine_nozzle_y) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

! ----- original -----
!      end

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

! ----- original -----
!      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear, &
!                                vort,DIMS(vort),lo,hi,nvar, &
!                                domlo,domhi,dx,xlo, &
!       			       problo,time,level)

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

      integer   i, j

#include <probdata.H>

      if (level .lt. max_vort_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
!               tag(i,j) = cvmgt(set,tag(i,j), &
               tag(i,j) = merge(set,tag(i,j), &
                     ABS(vort(i,j,1)).ge.vorterr*2.d0**level)
            enddo
         enddo
      end if

! ----- original -----
!      end

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

! ----- original -----
!      subroutine FORT_DENFILL (den,DIMS(den),domlo,domhi,delta, &
!                                xlo,time,bc)

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
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      lo(2) = ARG_L2(den)
      hi(1) = ARG_H1(den)
      hi(2) = ARG_H2(den)

      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif

! ----- original -----
!      end

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

! ----- original
!      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)

  subroutine adv_fill (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)&
                           bind(C, name="adv_fill")

      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i,j
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(adv)
      lo(2) = ARG_L2(adv)
      hi(1) = ARG_H1(adv)
      hi(2) = ARG_H2(adv)

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do i = lo(1), hi(1)
               adv(i,j) = 0.0d0
            enddo
         enddo
      endif

! ----- original -----
!      end

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

! ----- original -----
!      subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,delta, &
!                                xlo,time,bc)

      subroutine temp_fill (temp,DIMS(temp),domlo,domhi,delta, &
                              xlo,time,bc)&
                              bind(C, name="temp_fill")

      implicit none

      integer DIMDEC(temp), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  temp(DIMV(temp))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(temp)
      lo(2) = ARG_L2(temp)
      hi(1) = ARG_H1(temp)
      hi(2) = ARG_H2(temp)

      call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif

! ----- original -----      
!      end

  end subroutine temp_fill

      subroutine FORT_PHIEFIELDFILL (phiEField,DIMS(phiEField),domlo,domhi,delta, &
                                      xlo,time,bc)

      implicit none

      integer DIMDEC(phiEField), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  phiEField(DIMV(phiEField))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(phiEField)
      lo(2) = ARG_L2(phiEField)
      hi(1) = ARG_H1(phiEField)
      hi(2) = ARG_H2(phiEField)

      call filcc (phiEField,DIMS(phiEField),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               phiEField(i,j) = phiV
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               phiEField(i,j) = phiV
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               phiEField(i,j) =  phiV
            enddo
         enddo
      endif    

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               phiEField(i,j) =  phiV
            enddo
         enddo
      endif
      
      end

      subroutine FORT_EFIELDXFILL (Ex,DIMS(Ex),domlo,domhi,delta, &
                                    xlo,time,bc)

      implicit none

      integer DIMDEC(Ex), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  Ex(DIMV(Ex))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(Ex)
      lo(2) = ARG_L2(Ex)
      hi(1) = ARG_H1(Ex)
      hi(2) = ARG_H2(Ex)

      call filcc (Ex,DIMS(Ex),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ex(i,j) = 0.0d0 
	enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ex(i,j) = 0.0d0 
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ex(i,j) = 0.0d0 
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ex(i,j) = 0.0d0 
            enddo
         enddo
      endif
      
      end

      subroutine FORT_EFIELDYFILL (Ey,DIMS(Ey),domlo,domhi,delta, &
                                    xlo,time,bc)

      implicit none

      integer DIMDEC(Ey), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  Ey(DIMV(Ey))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(Ey)
      lo(2) = ARG_L2(Ey)
      hi(1) = ARG_H1(Ey)
      hi(2) = ARG_H2(Ey)

      call filcc (Ey,DIMS(Ey),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ey(i,j) = 0.0d0 
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ey(i,j) = 0.0d0
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ey(i,j) = 0.0d0 
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               Ey(i,j) = 0.0d0 
            enddo
         enddo
      endif
      
      end

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

! ----- original -----
!      subroutine FORT_RHOHFILL (rhoh,DIMS(rhoh),domlo,domhi,delta, &
!                                xlo,time,bc)

  subroutine rhoh_fill (rhoh,DIMS(rhoh),domlo,domhi,delta, &
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
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoh)
      lo(2) = ARG_L2(rhoh)
      hi(1) = ARG_H1(rhoh)
      hi(2) = ARG_H2(rhoh)

      call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,phiV,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif

! ----- original -----
!      end

  end subroutine rhoh_fill

! 
! Fill x & y velocity at once.
! 

! ----- original -----
!      subroutine FORT_VELFILL (vel,DIMS(vel),domlo,domhi,delta, &
!                                xlo,time,bc)

  subroutine vel_fill (vel,DIMS(vel),domlo,domhi,delta, &
                              xlo,time,bc)&
                              bind(C, name="vel_fill")

      implicit none
      integer DIMDEC(vel), bc(SDIM,2,SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  vel(DIMV(vel),SDIM)

      call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),1), &
        DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

      call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),2), &
        DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

! ----- original -----
!      end
  end subroutine vel_fill

! 
! Fill all chem species at once
! 

! ----- original -----
!      subroutine FORT_ALLCHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta, &
!                                    xlo,time,bc)

!      implicit none
!#include <cdwrk.H>
!#include <bc.H>
!#include <probdata.H>

!      integer DIMDEC(rhoY), bc(SDIM,2,Nspec)
!      integer domlo(SDIM), domhi(SDIM)
!      REAL_T  delta(SDIM), xlo(SDIM), time
!      REAL_T  rhoY(DIMV(rhoY),Nspec)

!      integer n
!      
!      do n=1,Nspec
!         call FORT_CHEMFILL (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),n), &
!               DIMS(rhoY),domlo,domhi,delta,xlo,time,bc(1,1,n),n-1)
!      enddo
!      end

  subroutine all_chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta, &
                            xlo,time,bc)&
                            bind(C, name="all_chem_fill")

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
         call chem_fill (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),n), &
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

! ----- original -----
!      subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta, &
!                                 xlo,time,bc)

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

      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(xvel)
      hi(1) = ARG_H1(xvel)
      lo(2) = ARG_L2(xvel)
      hi(2) = ARG_H2(xvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)
      
!      NOTE:
!      In order to set Dirichlet boundary conditions in a mulitspecies
!      problem, we have to know all the state values, in a sense.  For
!      example, the total density rho = sum_l(rho.Yl).  So to compute any
!      rho.Yl, we need all Yl's...also need to evaluate EOS since we
!      really are specifying T and Yl's.  so, all this is centralized
!      here.  Finally, a layer of flexibilty is added to for the usual case
!      that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif

! ----- original -----      
!      end

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

! ----- original -----
!      subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta, &
!                                 xlo,time,bc)

  subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta, &
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
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(yvel)
      hi(1) = ARG_H1(yvel)
      lo(2) = ARG_L2(yvel)
      hi(2) = ARG_H2(yvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

!      NOTE:
!      In order to set Dirichlet boundary conditions in a mulitspecies
!      problem, we have to know all the state values, in a sense.  For
!      example, the total density rho = sum_l(rho.Yl).  So to compute any
!      rho.Yl, we need all Yl's...also need to evaluate EOS since we
!      really are specifying T and Yl's.  so, all this is centralized
!      here.  Finally, a layer of flexibilty is added to for the usual case
!      that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif

! ----- original -----
!      end

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

! ----- original -----      
!      subroutine FORT_CHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta, &
!                                 xlo,time,bc,id )

 subroutine chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta, &
                            xlo,time,bc,id ) &
                            bind(C, name="chem_fill")
      implicit none
      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx, xhi(SDIM)
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h
      REAL_T  phiV
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoY)
      hi(1) = ARG_H1(rhoY)
      lo(2) = ARG_L2(rhoY)
      hi(2) = ARG_H2(rhoY)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)
      
!      NOTE:
!      In order to set Dirichlet boundary conditions in a mulitspecies
!      problem, we have to know all the state values, in a sense.  For
!      example, the total density rho = sum_l(rho.Yl).  So to compute any
!      rho.Yl, we need all Yl's...also need to evaluate EOS since we
!      really are specifying T and Yl's.  so, all this is centralized
!      here.  Finally, a layer of flexibilty is added to for the usual case
!      that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, phiV, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif

! ----- original -----      
!      end

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

! ----- original -----
!      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)

  subroutine press_fill (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)&
                            bind(C, name="press_fill")

      implicit none
      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi
      logical    per_xlo, per_xhi, per_ylo, per_yhi

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
! 
!   left side
! 
      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do j = jlo,jhi
               p(i,j) = p(ilo,j)
            end do
         end do
         if (fix_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,jlo)
               end do
            end do
         else if (per_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
         if (fix_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,jhi)
               end do
            end do
         else if (per_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
      end if
! 
!    right side
! 
      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do j = jlo,jhi
               p(i,j) = p(ihi,j)
            end do
	 end do
	 if (fix_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,jlo)
               end do
	    end do
	 else if (per_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
	 if (fix_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,jhi)
               end do
	    end do
	 else if (per_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
      end if
      
      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do i = ilo, ihi
               p(i,j) = p(i,jlo)
            end do
	 end do
	 if (per_xlo) then
          do j = ARG_L2(p), domlo(2)-1
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
	 if (per_xhi) then
           do j = ARG_L2(p), domlo(2)-1
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
      end if

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ilo, ihi
               p(i,j) = p(i,jhi)
            end do
	 end do
	 if (per_xlo) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jhi)
               end do
	    end do
         end if
	 if (per_xhi) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jhi)
               end do
	    end do
         end if
      end if

! ----- original -----
!      end

  end subroutine press_fill

      subroutine FORT_RADLOSS(lo,hi,rad,DIMS(rad), &
                               T,DIMS(T),Y,DIMS(Y),dx,Patm,time)
      implicit none

#include <cdwrk.H>
#include <probdata.H>

      integer DIMDEC(rad)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer lo(SDIM), hi(SDIM)
      REAL_T  rad(DIMV(rad))
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),1)
      REAL_T  dx(SDIM), Patm, time

      integer i, j
      
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rad(i,j) = zero
         end do
      end do
      end

! 
! 
! ::: -----------------------------------------------------------
! 
!     This routine add the forcing terms to the momentum equation
! 
      subroutine FORT_MAKEFORCE(time,force,rho, &
                                 DIMS(istate),DIMS(state), &
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

      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      integer a2, a3, a4, a5
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  sga, cga
      integer isioproc
      integer nXvel, nYvel, nRho, nTrac

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)

      ilo = istate_l1
      jlo = istate_l2
      ihi = istate_h1
      jhi = istate_h2

!      Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nRho  = 3
      nTrac = 4

      if (scomp.eq.0) then
         if (abs(gravity).gt.0.0001) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = gravity*rho(i,j)
               enddo
            enddo
!      else to zero
         else
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = zero
               enddo
            enddo
         endif
!      Add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nYvel) = force(i,j,nYvel) + dV_control*rho(i,j)
               enddo
            enddo
         endif
!      End of velocity forcing
      endif
      
      if ((scomp+ncomp).gt.BL_SPACEDIM) then
!      Scalar forcing
         do n = max(scomp+1,nRho), scomp+ncomp
            if (n.eq.nRho) then
!      Density
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
!      Tracer
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
!      Other scalar
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif

      end

end module prob_2D_module

