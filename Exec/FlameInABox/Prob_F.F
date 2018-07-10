
#include <MyProb_F.H>

      subroutine FORT_SET_PROB_SPEC(fuel, oxid, prod, numspec)
      implicit none
#include <cdwrk.H>
#include <probdata.H>
      integer fuel, oxid, prod, numspec
      fuelID = fuel + 1
      oxidID = oxid + 1
      prodID = prod + 1

      if (numspec .ne. Nspec) then
         call bl_pd_abort('number of species not consistent')
      endif
      end

c ::: -----------------------------------------------------------
c ::: This routine is called at problem initialization time
c ::: and when restarting from a checkpoint file.
c ::: The purpose is (1) to specify the initial time value
c ::: (not all problems start at time=0.0) and (2) to read
c ::: problem specific data from a namelist or other input
c ::: files and possibly store them or derived information
c ::: in FORTRAN common blocks for later use.
c ::: 
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: init      => TRUE if called at start of problem run
c :::              FALSE if called from restart
c ::: strttime <=  start problem with this time variable
c ::: 
c ::: -----------------------------------------------------------
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
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

      namelist /fortin/ vorterr, temperr, adverr, tempgrad, 
     &                  flametracval, probtype,
     &		        max_temp_lev, max_vort_lev, max_trac_lev,
     &                  traceSpecVal,phi_in,T_in,
     $                  turb_scale, V_in, V_co,
     $                  standoff, pertmag, nchemdiag, splitx, xfrontw,
     &                  splity, yfrontw, blobx, bloby, blobz, blobr,
     &                  blobT, Tfrontw, stTh, fuel_N2_vol_percent
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber
#if defined(BL_DO_FLCT)
      integer nCompFile
      namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi,
     &     strmwse_dir, nCompInflow, flct_file, convVelInit
#endif
      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, 
     &     zbase_control, pseudo_gravity, istemp,corr,controlVelMax,navg_pnts
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

      namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz,
     &                  forcing_type, spectrum_type, time_offset, forcing_twice_wavelength,
     &                  forcing_xlength, forcing_ylength, forcing_zlength,
     &                  forcing_time_scale_min, forcing_time_scale_max, 
     &                  force_scale, forcing_epsilon, blrandseed,
     &                  use_rho_in_forcing, do_mode_division,
     &                  div_free_force, moderate_zero_modes
#endif

c
c      Build `probin' filename -- the name of file containing fortin namelist.
c
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)

#if defined(BL_DO_FLCT)
      integer ierr
#endif


      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
c         call bl_abort('probinit called with init ne 1')
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

c     Load domain dimensions into common (put something computable there for SDIM<3)
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
      
c     Set defaults
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
      pamb = FORT_P1ATMMKS()
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
c     add_turb = .FALSE.
      forceInflow = .FALSE.
      numInflPlanesStore = -1
      numInflPlanesStore = 30
c
c     Don't need to default 'nCompInflow' as it is block data'd to /3/
c
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = ''
      turb_scale = 1
#endif

      zbase_control = 0.d0

c     Note: for setup with no coflow, set Ro=Rf+wallth
      standoff = zero
      pertmag = 0.d0

c     Initialize control variables
      tau_control = one
      controlVelMax = 2.d0
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
      
c     Initialize control variables that depend on fortin variables
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
c     Set up boundary functions
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
      if ( (probtype.eq.BL_PROB_PREMIXED_FIXED_INFLOW)
     &     .or. (probtype.eq.BL_PROB_PREMIXED_CONTROLLED_INFLOW) ) then

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
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
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
c     Amplitudes (alpha)
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
c     Count modes that contribute
                        mode_count = mode_count + 1
c     Set amplitudes
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

c     Now let's break symmetry, have to assume high aspect ratio in z for now
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
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
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
c     Amplitudes (alpha)
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
c     Count modes that contribute
                        reduced_mode_count = reduced_mode_count + 1
c     Set amplitudes
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

      end


      subroutine set_Y_from_Phi(phi,Yt)
      implicit none
#include <probdata.H>
#include <cdwrk.H>
#include <conp.H>
      REAL_T a, phi
      REAL_T Xt(maxspec), Yt(maxspec)
      integer n
      character*(maxspnml) name
      do n=1,Nspec
         Xt(n) = zero
      enddo
      
c     Set "a" for computing X from phi
c     hc + a.O2 -> b.CO2 + c.H2O
      
      call get_spec_name(name,fuelID)

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'H2') then
         a = .5d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'CH3OCH3') then
         a = 3.0d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      Xt(oxidID) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * Xt(oxidID) / a
      Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(oxidID)
      
      CALL CKXTY (Xt, IWRK, RWRK, Yt)
      end
