#include <MyProb_F.H>

module prob_2D_module
  
  implicit none

  private
  
  public :: setupbc, amrex_probinit, getZone, bcfunction, init_data_new_mech, init_data, &
            zero_visc, flame_tracer_error, adv_error, &
            temp_error, mv_error, den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, chem_fill, press_fill

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
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c, name="amrex_probinit")

    use chem_driver, only: P1ATMMKS
    use chem_driver_2D, only: RHOfromPTY, HMIXfromTY

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(SDIM), probhi(SDIM)

#include <probdata.H>
#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <visc.H>
#include <conp.H>
#if defined(BL_DO_FLCT)
#include <infl_frc.H>
#include <INFL_FORCE_F.H>
#endif

    integer zone,i

    namelist /fortin/ vorterr, temperr, adverr, tempgrad, flametracval, probtype,&
         max_temp_lev, max_vort_lev, splitx, splity, traceSpecVal,&
         xfrontw, yfrontw, refine_nozzle_x,&
         refine_nozzle_y, blobx, bloby,&
         blobr, xcen, v_strength,v_width,v_xcen,v_ycen,v_cl_x,&
         y_offset, reflect, &
         v_blob_r, v_blob_T, v_blob_airfrac, V_in, V_co,&
         stTh, Rf, R_hot, R_hotBL, twall, turb_scale, stBL,&
         max_trac_lev, max_nozzle_lev, pipeTh, pipeBL,&
         tV_in_l,tV_in_r,V_in_new,tV_co_l,tV_co_r,V_co_new,&
         T_stick, T_in,phi_in, thickFacTR, thickFacCH,nchemdiag,&
         pertmag
    namelist /heattransin/ pamb, dpdt_factor, closed_chamber
#if defined(BL_DO_FLCT)
    namelist /flctin/ forceInflow, numInflPlanesStore, forceLo, forceHi,&
         strmwse_dir, nCompInflow, flct_file
#endif
    namelist /control/ tau_control, sest, cfix, changeMax_control, h_control,&
         zbase_control, tbase_control
    !
    !      Build 'probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)

    character probin*(maxlen)

    call bl_pd_is_ioproc(isioproc)

    if (namlen .gt. maxlen) then
       write(6,*) 'probin file name too long'
       stop
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
    splitx = zero
    splity = zero
    traceSpecVal = 1.d-14
    xfrontw = fourth*half*(problo(1)+probhi(1))
    yfrontw = xfrontw
    refine_nozzle_y = zero
    refine_nozzle_x = zero
    max_nozzle_lev = 0
    blobx = -one
    bloby = -one
    blobr = -one
    xcen = problo(1)
    v_strength = zero
    v_width = zero
    v_xcen = zero
    v_ycen = zero
    v_cl_x = probhi(1)
    y_offset = zero
    reflect = zero
    pamb = 101325.d0
    closed_chamber = 0
    dpdt_factor = 0.3d0
    v_blob_r = zero
    v_blob_T = 2000.d0
    v_blob_airfrac = zero
    V_in = zero
    V_co = zero
    stTh = zero
    Rf = zero
    R_hot = half*probhi(1)
    R_hotBL = .1*probhi(1)
    twall = 298.d0
    stBL = R_hotBL
    pipeTh = zero
    tV_in_l = 1.e6
    tV_in_r = 2.e6
    V_in_new = -100.d0
    tV_co_l = 1.e6
    tV_co_r = 2.e6
    V_co_new = -100.d0
    T_stick(1) = -300.d0
    T_in = 300.d0
    phi_in = 1.d0
    zbase_control = 0.d0
    pertmag = 0.d0

#if defined(BL_DO_FLCT)
    !
    !     Don't need to default 'nCompInflow' as it is block data'd to /3/
    !
    forceInflow = .FALSE.
    numInflPlanesStore = -1
    forceLo = .TRUE.
    forceHi = .FALSE.
    strmwse_dir = FLCT_YVEL
    flct_file = ""
    turb_scale = 1
    nCompInFlow = 2
#endif

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
    tbase_control = zero
    
    read(untin,fortin)

    !     Initialize control variables that depend on fortin variables
    V_in_old = V_in
      
    !      if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev
      
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
       if (flct_file.ne."") then
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
    !      convVel = V_in
    convVel = one
#endif

    !     Load domain dimensions into common, and set up boundary functions
    domnlo(1) = problo(1)
    domnlo(2) = problo(2)
    domnhi(1) = probhi(1)
    domnhi(2) = probhi(2)
    
    call setupbc()
    bcinit = .true.

    zone = getZone(half*(domnlo(1)+domnhi(1)),domnlo(2))
    scale_control = Y_bc(fuelID-1,zone)*rho_bc(zone)*(domnhi(1)-domnlo(1))
    if (h_control .gt. zero) then
       cfix = scale_control * h_control
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

  !
  ! This is a version of INITDATA in which the velocities,
  ! the temperature & the mass fractions are passed in.
  !
  subroutine init_data_new_mech(level,time,lo,hi,nscal, &
       vel,scal,DIMS(state),press,DIMS(press), &
       delta,xlo,xhi) bind(C, name="init_data_new_mech")
    
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

    Patm = pamb / P1ATMMKS()
     
    call RHOfromPTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),Density), DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state),&
         Patm)
    call HMIXfromTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
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
  

  subroutine setupbc() bind(C,name='setupbc')

    use chem_driver, only: P1ATMMKS, get_spec_name
    use chem_driver_2D, only: RHOfromPTY, HMIXfromTY

    implicit none

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
#include <htdata.H>
      
    REAL_T Patm, y, pmf_vals(maxspec+3), sum, Y_pmf(maxspec), yEval
    REAL_T Yfu(maxspec), Yox(maxspec)
    REAL_T Y_FU_H2, Y_FU_O2, Y_OX_H2, Y_OX_O2
    REAL_T Xfu(maxspec), a
    integer zone, n, b(SDIM), iN2, iH2, iO2
    character*(maxspnml) name
    data b / 1, 1 /

    data Y_FU_H2 / 0.0392668168889d0 /
    data Y_FU_O2  / 0.223771589041d0 /
    data Y_OX_H2 / zero /
    data Y_OX_O2  / 0.232917518594d0 /
      
    Patm = pamb / P1ATMMKS()

    if (probtype .eq. 1) then
       !     For probtype == 1, set boundary conditions from premixed flame profile in y
       y = domnlo(2)
       zone = getZone(y,domnlo(2))
       iN2 = -1
       do n = 1,Nspec
          call get_spec_name(name,n)
          if (name .eq. 'N2' ) iN2 = n
       end do

       yEval = y*100. - y_offset*100.
       call pmf(yEval,yEval,pmf_vals,n)

       if (n.ne.Nspec+3) then
          print *,'n=',n
          print *,'Nspec=',Nspec
          call bl_abort('SETUPBC: n .ne. Nspec+3')
       endif

       CALL CKXTY (pmf_vals(4), Y_pmf)
         
       do n=0,Nspec-1
          Y_bc(n,zone) = Y_pmf(n+1)
       end do
       if (iN2.gt.0) then
          sum = zero
          do n = 1,Nspec
             if (n.ne.iN2) sum = sum+Y_bc(n-1,zone)
          end do
          Y_bc(iN2-1,zone) = one - sum
       end if

       T_bc(zone) = pmf_vals(1)
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

    else if (probtype .eq. 2) then
       iN2 = -1
       iH2 = -1
       iO2 = -1
       do n = 1,Nspec
          call get_spec_name(name,n)
          if (name .eq. 'N2' ) iN2  = n
          if (name .eq. 'H2') iH2 = n
          if (name .eq. 'O2' ) iO2  = n
       end do

       do n=1,Nspec
          Yfu(n) = zero
          Yox(n) = zero
       enddo
       Yfu(iH2) = Y_FU_H2
       Yfu(iO2)  = Y_FU_O2
       Yfu(iN2)  = one - Y_FU_H2 - Y_FU_O2

       Yox(iH2) = Y_OX_H2
       Yox(iO2)  = Y_OX_O2
       Yox(iN2)  = one - Y_OX_H2 - Y_OX_O2

       zone = M_FUEL
       do n=1,Nspec
          Y_bc(n-1,zone) = Yfu(n)
       end do
       T_bc(zone) = T_in
       u_bc(zone) = zero
       v_bc(zone) = V_in

       zone = M_CO
       do n=1,Nspec
          Y_bc(n-1,zone) = Yox(n)
       end do
       T_bc(zone) = T_in
       u_bc(zone) = zero
       v_bc(zone) = V_co

       zone = M_WALL
       do n=1,Nspec
          Y_bc(n-1,zone) = Yox(n)
       end do
       T_bc(zone) = twall
       u_bc(zone) = zero
       v_bc(zone) = V_co

       do zone=1,3
          call RHOfromPTY(b,b,&
               rho_bc(zone),DIMARG(b),DIMARG(b),&
               T_bc(zone),DIMARG(b),DIMARG(b),&
               Y_bc(0,zone),DIMARG(b),DIMARG(b), Patm)

          call HMIXfromTY(b,b,&
               h_bc(zone),DIMARG(b),DIMARG(b),&
               T_bc(zone),DIMARG(b),DIMARG(b),&
               Y_bc(0,zone),DIMARG(b),DIMARG(b))
       enddo

       !     Get something reasonable for stick
       if (T_stick(1).gt.zero) then

          call RHOfromPTY(b,b,&
               rho_stick,DIMARG(b),DIMARG(b),&
               T_stick,DIMARG(b),DIMARG(b),&
               Y_bc(0,M_FUEL),DIMARG(b),DIMARG(b), Patm)

          call HMIXfromTY(b,b,&
               h_stick,DIMARG(b),DIMARG(b),&
               T_stick,DIMARG(b),DIMARG(b),&
               Y_bc(0,M_FUEL),DIMARG(b),DIMARG(b))
       endif

    elseif (probtype .eq. 3) then
       !     For probtype == 3, set boundary conditions as constants

       iN2 = -1
       do n = 1,Nspec
          call get_spec_name(name,n)
          if (name .eq. 'N2' ) iN2  = n
       end do

       call get_spec_name(name,fuelID)
       if (name .eq. 'H2' ) then
          a = half
          print *,'this is hydrogen',fuelID,oxidID,iN2
       else if (name .eq. 'CH4' ) then
          a = two
          print *,'this is methane',fuelID,oxidID,iN2
       else if (name .eq. 'C3H8' ) then
          a = five
          print *,'this is propane',fuelID,oxidID,iN2
       else if (phi_in .ge. zero) then
          call bl_abort('Unknown fuel species')
       endif

       if (phi_in .ge. zero) then

          do n=1,Nspec
             Xfu(n) = zero
          enddo

          Xfu(oxidID) = one/(one + phi_in/a  + 0.79d0/0.21d0)
          Xfu(fuelID) = phi_in * Xfu(oxidID) / a
          Xfu(iN2) = (0.79d0/0.21d0) * Xfu(oxidID)

          call CKXTY(Xfu,Yfu)
            
       else

          !     This is for hydrogen
          !     phi = 1.0
          Y_FU_H2 = 0.0551670662321d0
          Y_FU_O2  = 0.220068142419d0
          !     phi = 0.7
          !     Y_FU_H2 = 0.0392668168889d0
          !     Y_FU_O2  = 0.223771589041d0
          do n=1,Nspec
             Yfu(n) = zero
          enddo
          Yfu(fuelID) = Y_FU_H2
          Yfu(oxidID) = Y_FU_O2
          Yfu(iN2)  = one - Y_FU_H2 - Y_FU_O2

       endif

       zone = M_FUEL
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

    else
       call bl_abort('No boundary init yet for this probtype')
    end if
      
  end subroutine setupbc

         
! ::: -----------------------------------------------------------
      
  integer function getZone(x, y) bind(C, name="getZone")

    implicit none

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

    REAL_T x, y, eps, absx
    data eps / 1.d-8 /

    if ((probtype.eq.1).or.(probtype.eq.3)) then
       getZone = 1
    else
       if (y.gt.domnlo(2)+eps) then
          getZone = M_WALL
       else
          absx = ABS(x)
          if (absx .lt. Rf + pipeTh) then
             getZone = M_FUEL
          else
             getZone = M_CO
          endif
       endif
    endif

  end function getZone
      

  subroutine VAvg(V,x,y,dx) bind(C, name="VAvg")

    implicit none

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>

    REAL_T V(SDIM)
    REAL_T x,y,dx(SDIM)

    REAL_T incx,x1
    REAL_T Vincx,Vincy
    integer i,M
    parameter (M=10)

    !     size of sub-lengths, there are M of them
    incx = dx(1) / M
    V(1) = zero
    V(2) = zero
    do i=1,M
       x1 = ABS(x - half*dx(1) + (i-half)*incx)
       Vincx = zero
       Vincy = zero

       if (x1.lt.Rf) then
          Vincy = v_bc(M_FUEL)*TANH(two*(one - x1/Rf)/pipeBL)
          if (stTh.gt.0) then
             if (x1.le.half*stTh) then
                Vincy = zero
             else
                Vincy = Vincy*TANH(two*(one - (Rf-x1)/((Rf-half*stTh)))/stBL)
             endif
          endif
       else if (x1.gt.Rf+pipeTh) then
          Vincy = v_bc(M_CO)*TANH(two*(one-MAX(zero,2*Rf+pipeTh-x1)/Rf)/pipeBL)
       endif

       V(1) = V(1) + Vincx
       V(2) = V(2) + Vincy
    enddo
    V(1) = V(1)/M
    V(2) = V(2)/M
  end subroutine VAvg

! ::: -----------------------------------------------------------
      
  subroutine bcfunction(orient, x, y, time, u, v, rho, Yl, T, h, dx)  bind(C, name="bcfunction")

    implicit none
    integer orient
    REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(SDIM)

#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

    REAL_T Va(SDIM)
    integer n, zone
    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    if (orient .lt. 4) then
       zone = getZone(x,y)
       u = u_bc(zone)
       v = v_bc(zone)
       rho = rho_bc(zone)
       do n = 0,Nspec-1
          Yl(n) = Y_bc(n,zone)
       end do
       T = T_bc(zone)
       h = h_bc(zone)

       if ((probtype.eq.1).or.(probtype.eq.3)) then

          v_bc(M_FUEL) = V_in + (time-tbase_control)*dV_control

       else if (probtype.eq.2) then

          !     Set fueling velocity as a function of time
          if (time .le. tV_in_l) then
             v_bc(M_FUEL) = V_in
          else if (time .ge. tV_in_r) then
             v_bc(M_FUEL) = V_in_new
          else
             v_bc(M_FUEL) = V_in+(time-tV_in_l)*(V_in_new-V_in)/(tV_in_r-tV_in_l)
          endif

          !     Set coflow velocity as a function of time
          if (time .le. tV_co_l) then
             v_bc(M_CO) = V_co
          else if (time .ge. tV_co_r) then
             v_bc(M_CO) = V_co_new
          else
             v_bc(M_CO) = V_co+(time-tV_co_l)*(V_co_new-V_co)/(tV_co_r-tV_co_l)
          endif

          call VAvg(Va,x,y,dx)
          u = Va(1)
          v = Va(2)

!     Set stick temperature
          if (zone.eq.M_FUEL  .and.  ABS(x).le..5*stTh .and. T_stick(1).gt.zero ) then
             T = T_stick(1)
             rho = rho_stick(1)
             h = h_stick(1)
          endif
       endif
    else
       print*, 'No boundary condition for orientation = ', orient
       call bl_abort(' ')
    end if
  end subroutine bcfunction
      
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

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

    integer i, j, n
    REAL_T  x, y
    character*(maxspnml) name
    REAL_T Patm
      
    REAL_T pmf_vals(maxspec+3), theta
    REAL_T  r,rfront,Xl(maxspec),Yl(maxspec), Xp(maxspec)
    REAL_T  a
    integer iAr, iO2, iH2, iH2O
    integer b(SDIM)
    data b /1,1/

    if ((Temp.gt.0).neqv.(RhoH.gt.0)) then
       call bl_abort('Need both Temp and RhoH, or neither')
    end if
      
    if ((Temp .LT. 0) .OR. (RhoH .LT. 0)) then
       call bl_abort('No IC''s for system without T, RhoH')
    endif

    iAr = -1
    iO2 = -1
    iH2 = -1
    if (iN2.lt.1 .or. iN2.gt.Nspec) then
       call bl_pd_abort()
    endif
    do n = 1,Nspec
       call get_spec_name(name,n)
       if (name .eq. 'AR' ) iAr = n
       if (name .eq. 'O2' ) iO2 = n
       if (name .eq. 'H2' ) iH2 = n
       if (name .eq. 'H2O' ) iH2O = n
    end do
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          do n=1,SDIM
             vel(i,j,n) = zero
          end do
       end do
    end do
 
    do j = lo(2), hi(2)
       y = (dfloat(j)+0.5)*delta(2)+domnlo(2)
       do i = lo(1), hi(1)
          x = (dfloat(i)+.5)*delta(1)+domnlo(1)
          r = SQRT((x-blobx)**2+(y-bloby)**2)
 
          if (r.LE.delta(1)) then
             theta = zero
          else
             theta = ATAN2((y-bloby)/r,(x-blobx)/r)
          endif

          rfront = sin(3.d0*theta) + 0.583d0*sin(5.d0*theta - 1.077)&
               + .634d0*sin(10.d0*theta - 2.0546d0)&
               + .364d0*sin(15.d0*theta - 0.0546d0)

          if(probtype .eq. 1) then
              
             if (reflect .eq. 0.d0) then
                rfront = (blobr + pertmag* rfront * blobr &
                     - r + y_offset) *100.d0
             else
                rfront = (y_offset - abs(y)) * 100.d0
             end if
 
             call pmf(rfront,rfront,pmf_vals,n)
 
             if (n.ne.Nspec+3) then
                call bl_abort('INITDATA: n .ne. Nspec+3')
             endif
               
             scal(i,j,Temp) = pmf_vals(1)
             do n = 1,Nspec
                Xl(n) = pmf_vals(3+n)
             end do

          else

             rfront = (blobr + pertmag* rfront * blobr - r ) / (.4*blobr) 

             scal(i,j,Temp) = T_in + (T_stick(1) - T_in ) *&
                  0.5d0*(1.d0+tanh(rfront))

             do n=1,Nspec
                Xl(n) = zero
                Xp(n) = zero
             enddo

             Xl(iO2) = one/(one + phi_in/half  + 0.79d0/0.21d0)
             Xl(iH2) = phi_in * Xl(iO2) / half
             Xl(iN2) = (0.79d0/0.21d0) * Xl(iO2)

             Xp(iO2) = Xl(iO2)-Xl(iH2)/2.d0
             Xp(iH2O) = Xl(iH2)
             Xp(iN2) = Xl(iN2)
  
             a = Xp(iO2) + Xp(iH2O) + Xp(iN2)
            
             Xp(iO2) = Xp(iO2)/ a
             Xp(iH2O) = Xp(iH2O)/ a
             Xp(iN2) = Xp(iN2)/ a

             a = zero
             do n=1,Nspec
                Xl(n) = Xl(n) + (Xp(n) - Xl(n))* 0.5d0*(1.d0+tanh(rfront))
                a = a + Xl(n)
             enddo

             if(abs(a-1.d0).gt.1.d-12)then
                write(6,*)i,j,Xp,Xl
                stop
             endif

          endif
 
          call CKXTY(Xl,Yl)
 
          do n = 1,Nspec
             scal(i,j,FirstSpec+n-1) = Yl(n)
          end do
          scal(i,j,Trac) = zero
 
       end do
    end do

    Patm = pamb / P1ATMMKS()
    
    call RHOfromPTY(b,b,&
         scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),firstSpec),DIMS(state),&
         Patm)

    call HMIXfromTY(b,b,&
         scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),firstSpec),DIMS(state))

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          do n = 0,Nspec-1
             scal(i,j,firstSpec+n) = scal(i,j,firstSpec+n)*scal(i,j,Density)
          enddo
          scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
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
! ::: lo,hi      => region of interest
! ::: domlo,hi   => index extent of problem domain
! ::: dx         => cell spacing
! ::: bc         => boundary condition flag (on orient)
! :::                   in BC_TYPES::physicalBndryTypes
! ::: problo     => phys loc of lower left corner of prob domain
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
#include <cdwrk.H>
#include <htdata.H>

    integer i, j, n, Tid, RHid, YSid, YEid, ys, ye
    logical do_T, do_RH, do_Y
    REAL_T xl, xr, xh, y, ap

    if (probtype.eq.2) then
       Tid  = Temp      - id + SDIM
       RHid = RhoH      - id + SDIM
       YSid = FirstSpec - id + SDIM
       YEid = LastSpec  - id + SDIM
         
       do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
       do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
       ys = MAX(YSid,1)
       ye = MIN(YEid,ncomp)
       do_Y = (ye - ys + 1) .GE. 1

       !     Do species, Temp, rhoH
       if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))&
            .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
          j = lo(2)
          y = float(j)*dx(2)+domnlo(2)
          do i = lo(1), hi(1)

             xl = float(i)*dx(1)+domnlo(1) 
             xr = (float(i)+1.d0)*dx(1)+domnlo(1) 
             xh = 0.5d0*(xl+xr)
               
             if ( (getZone(xl,y).eq.M_FUEL) .OR.&
                  (getZone(xh,y).eq.M_FUEL) .OR.&
                  (getZone(xr,y).eq.M_FUEL) ) then
                  
                ap = zero
#if 0
                if (xl.ge.-half*stTh) then
                   ap = ap + (MIN(Rf,xr)-MAX(xl,half*stTh))/dx(1)
                else if (xr.le.half*stTh) then
                   ap = ap + (MIN(xr,-half*stTh)-MAX(-Rf,xl))/dx(1)
                else
                   call bl_abort("cell spans zero")
                endif
#endif

#if 0
                if (do_T)  diff(i,j,Tid ) = diff(i,j,Tid )*ap
                if (do_RH) diff(i,j,RHid) = diff(i,j,RHid)*ap
#else
                !     Conducting stick: only zero conduction if not on stick
                if ( (T_stick(1) .le. zero)  .or. (&
                     (ABS(xl).gt..5*stTh) .and.&
                     (ABS(xh).gt..5*stTh) .and.&
                     (ABS(xr).gt..5*stTh) ) ) then

                   if (do_T)  diff(i,j,Tid ) = diff(i,j,Tid )*ap
                   if (do_RH) diff(i,j,RHid) = diff(i,j,RHid)*ap

                endif
#endif
                if (do_Y) then
                   do n=ys,ye
                      diff(i,j,n) = diff(i,j,n)*ap
                   enddo
                endif
             endif
          enddo
       endif
    endif

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

    integer   i, j
    REAL_T    x, y
    logical   in_refine_zone

#include <probdata.H>

    if (level .lt. max_trac_lev) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             tag(i,j) = merge(set,tag(i,j),ftrac(i,j,1).gt.flametracval)
          enddo
       enddo
    endif

    if (level .lt. max_nozzle_lev) then
       do j = lo(2), hi(2)
          y = (float(j)+.5)*dx(2)+problo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*dx(1)+problo(1)
             in_refine_zone = (y - problo(2)) .le. refine_nozzle_y
             in_refine_zone = in_refine_zone .and. &
                  (ABS(x-xcen) .le. refine_nozzle_x)
             tag(i,j) = merge(set,tag(i,j),in_refine_zone)
          end do
       end do
    end if
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
      
    if (time .eq. zero) then
       call mv_error(tag,DIMS(tag),set,clear,&
            adv,DIMS(adv),lo,hi,nvar,&
            domlo,domhi,delta,xlo,&
            problo,time,level)
    end if

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

    REAL_T    ax, ay, aerr
    integer   i, j, ng

#include <probdata.H>

    ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp))

    if (ng .lt. 1) then
       print*, "TEMPERR cannot compute gradient, ng = ",ng
       call bl_abort(' ')
    endif
    !
    !     refine where there is temperature gradient
    !
    if (level .lt. max_temp_lev) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
             ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
             aerr = max(ax,ay)
             tag(i,j) = merge(set,tag(i,j),aerr.ge.tempgrad)
             !              tag(i,j) = merge(set,tag(i,j),temperature(i,j,1).le.tempgrad)
          enddo
       enddo
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

    REAL_T    x, y
    integer   i, j
    REAL_T    v, r1sq, r2sq

#include <probdata.H>

    if (time.eq.zero .and. level.lt.max_vort_lev) then
       do j = lo(2), hi(2)
          y = float(j)*dx(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = float(i)*dx(1)+domnlo(1)
             r1sq = (x-v_xcen)**2              + (y-v_ycen)**2
             r2sq = (x-(two*v_cl_x-v_xcen))**2 + (y-v_ycen)**2
             v = v_strength*&
                  (dexp(-r1sq/v_width**2)-dexp(-r2sq/v_width**2)) 
             tag(i,j) = merge(set,tag(i,j),ABS(v).ge.vorterr)
          end do
       end do
    endif

    if (level .lt. max_vort_lev) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             tag(i,j) = merge(set,tag(i,j),ABS(vort(i,j,1)).ge.vorterr*2.d0**level)
          end do
       end do
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
      
    integer i, j
    integer ilo, ihi, jlo, jhi
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(den)
    hi(1) = ARG_H1(den)
    lo(2) = ARG_L2(den)
    hi(2) = ARG_H2(den)

    hx  = delta(1)
    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))
      
    call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)
      
    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             den(i,j) = rho
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             den(i,j) = rho
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             den(i,j) = rho
          enddo
       enddo
    endif
      
    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             den(i,j) = rho
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

    integer    i, j
    integer    ilo, ihi, jlo, jhi

#include <probdata.H>

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(adv)
    hi(1) = ARG_H1(adv)
    lo(2) = ARG_L2(adv)
    hi(2) = ARG_H2(adv)

    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))

    call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          do j = lo(2), hi(2)
             adv(i,j) = zero
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          do j = lo(2), hi(2)
             adv(i,j) = zero
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       !                               inflow for burner in a can (bic, biac)

       do j = lo(2), domlo(2)-1
          do i = lo(1), hi(1)
             adv(i,j) = zero
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          do i = lo(1), hi(1)
             adv(i,j) = zero
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
      
    integer i, j
    integer ilo, ihi, jlo, jhi
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(temp)
    hi(1) = ARG_H1(temp)
    lo(2) = ARG_L2(temp)
    hi(2) = ARG_H2(temp)

    hx  = delta(1)
    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))
      
    call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             temp(i,j) = T
          enddo
       enddo
    endif
      
    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             temp(i,j) = T
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             temp(i,j) = T
          enddo
       enddo
    endif
      
    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             temp(i,j) = T
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
       xlo,time,bc,id ) bind(C, name="rhoh_fill")
    
    implicit none

    integer DIMDEC(rhoh), bc(SDIM,2)
    integer domlo(SDIM), domhi(SDIM), id
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  rhoh(DIMV(rhoh))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
    integer i, j
    integer ilo, ihi, jlo, jhi
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

    integer lo(SDIM), hi(SDIM)

    lo(1) = ARG_L1(rhoh)
    hi(1) = ARG_H1(rhoh)
    lo(2) = ARG_L2(rhoh)
    hi(2) = ARG_H2(rhoh)

    hx  = delta(1)
    ilo = max(lo(1),domlo(1))
    ihi = min(hi(1),domhi(1))
    jlo = max(lo(2),domlo(2))
    jhi = min(hi(2),domhi(2))
      
    call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)
      
    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoh(i,j) = rho*h
          enddo
       enddo
    endif
      
    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoh(i,j) = rho*h
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoh(i,j) = rho*h
          enddo
       enddo
    endif
      
    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoh(i,j) = rho*h
          enddo
       enddo
    endif
      
  end subroutine rhoh_fill
!
! Fill x & y velocity at once.
!
  subroutine vel_fill (vel,DIMS(vel),domlo,domhi,delta,&
       xlo,time,bc) bind(C, name="vel_fill")

    implicit none

    integer DIMDEC(vel), bc(SDIM,2,SDIM)
    integer domlo(SDIM), domhi(SDIM)
    REAL_T  delta(SDIM), xlo(SDIM), time
    REAL_T  vel(DIMV(vel),SDIM)

    call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),1),&
         DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

    call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),2),&
         DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

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
       call chem_fill (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),n),&
            DIMS(rhoY),domlo,domhi,delta,xlo,time,bc(1,1,n),n-1)
    enddo
  end subroutine all_chem_fill


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

  REAL_T function turbSclX(x,y) bind(C,name="turbSclX")
#include <probdata.H>
    REAL_T x, y, shapet, eta
    if ((probtype.eq.1).or.(probtype.eq.3)) then
       turbSclX = one
    else
       eta = ABS(x)/Rf
       turbSclX = shapet(eta)
       if (stBL.gt.zero) then
          turbSclX = turbSclX*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
       end if
    end if
  end function turbSclX

  REAL_T function turbSclY(x,y) bind(C,name="turbSclY")
#include <probdata.H>
    REAL_T x, y, eta
    turbSclY = zero
    if ((probtype.eq.1).or.(probtype.eq.3)) then
       turbSclY = one
    else
       eta = ABS(x)/Rf
       if(eta.lt.one)then
          turbSclY = 2.d0 -1.d0 * tanh((1-eta)/.0428222042d0)
       endif
       if (stBL.gt.zero) then
          turbSclY = turbSclY*TANH(4.d0*MAX(0.d0,ABS(x)-0.5d0*stTh)/stBL)
       end if
    endif
  end function turbSclY


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
    REAL_T turbSclX, t_flct, dt_flct
    integer DIMDEC(uflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T uflct(:,:)
    allocatable uflct
#endif
      
    integer i, j
    integer ilo, ihi, jlo, jhi
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

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
       dt_flct = time - tbase_control
       t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
       call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct,&
            bc, domnlo, domnhi)
    endif
#endif

    call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)

    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             xvel(i,j) = u
#if defined(BL_DO_FLCT)
             if (forceLo .and. strmwse_dir .eq. 1) then
                xvel(i,j) = xvel(i,j) + uflct(1,j)*turb_scale
             endif
#endif
          enddo
       enddo
    endif
      
    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             xvel(i,j) = u
#if defined(BL_DO_FLCT)
             if (forceHi .and. strmwse_dir .eq. 1) then
                xvel(i,j) = xvel(i,j) + uflct(1,j)*turb_scale
             endif
#endif
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             xvel(i,j) = u
#if defined(BL_DO_FLCT)
             if (forceLo .and. strmwse_dir .eq. 2) then
                xvel(i,j) = xvel(i,j) + uflct(i,1)*turb_scale*turbSclX(x,y)
             endif
#endif
          enddo
       enddo
    endif
      
    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             xvel(i,j) = u
#if defined(BL_DO_FLCT)
             if (forceHi .and. strmwse_dir .eq. 2) then
                xvel(i,j) = xvel(i,j) + uflct(i,1)*turb_scale
             endif
#endif
          enddo
       enddo
    endif
      
#if defined(BL_DO_FLCT)
    if (forceInflow) then
       deallocate(uflct)
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
    REAL_T turbSclY, t_flct, dt_flct
    integer DIMDEC(vflct)
    integer loFlctArray(SDIM), hiFlctArray(SDIM)
    REAL_T vflct(:,:)
    allocatable vflct
#endif

    integer i, j
    integer ilo, ihi, jlo, jhi
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

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
      
#if defined(BL_DO_FLCT)
    if (forceInflow) then
       do i = 1, SDIM
          loFlctArray(i) = lo(i)
          hiFlctArray(i) = hi(i)
       enddo
       loFlctArray(strmwse_dir) = 1
       hiFlctArray(strmwse_dir) = 1
       call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
       allocate(vflct(DIMV(vflct)))
       dt_flct = time - tbase_control
       t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
       call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, t_flct,&
            bc, domnlo, domnhi)
    endif
#endif

    call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             yvel(i,j) = v
#if defined(BL_DO_FLCT)
             if (forceLo .and. strmwse_dir .eq. 1) then
                yvel(i,j) = yvel(i,j) + vflct(1,j)*turb_scale
             endif
#endif
          enddo
       enddo
    endif
      
    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             yvel(i,j) = v
#if defined(BL_DO_FLCT)
             if (forceHi .and. strmwse_dir .eq. 1) then
                yvel(i,j) = yvel(i,j) + vflct(1,j)*turb_scale
             endif
#endif
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             yvel(i,j) = v
#if defined(BL_DO_FLCT)
             if (forceLo .and. strmwse_dir .eq. 2) then
                yvel(i,j) = yvel(i,j) + vflct(i,1)*turb_scale*turbSclY(x,y)
             endif
#endif
          enddo
       enddo
    endif
      
    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             yvel(i,j) = v
#if defined(BL_DO_FLCT)
             if (forceHi .and. strmwse_dir .eq. 2) then
                yvel(i,j) = yvel(i,j) + vflct(i,1)*turb_scale
             endif
#endif
          enddo
       enddo
    endif
      
#if defined(BL_DO_FLCT)
    if (forceInflow) then
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
       xlo,time,bc,id ) bind(C, name="chem_fill")

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
    REAL_T  y, x, hx
    REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

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

    !     NOTE:
    !     In order to set Dirichlet boundary conditions in a mulitspecies
    !     problem, we have to know all the state values, in a sense.  For
    !     example, the total density rho = sum_l(rho.Yl).  So to compute any
    !     rho.Yl, we need all Yl's...also need to evaluate EOS since we
    !     really are specifying T and Yl's.  so, all this is centralized
    !     here.  Finally, a layer of flexibilty is added to for the usual case
    !     that the bc values may often be set up ahead of time.

    if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
       do i = lo(1), domlo(1)-1
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XLO, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoY(i,j) = rho*Yl(id)
          enddo
       enddo
    endif

    if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)
          do j = lo(2), hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             call bcfunction(XHI, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoY(i,j) = rho*Yl(id)
          enddo
       enddo
    endif

    if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
       do j = lo(2), domlo(2)-1
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YLO, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoY(i,j) = rho*Yl(id)
          enddo
       enddo
    endif

    if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
       do j = domhi(2)+1, hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)
             call bcfunction(YHI, x, y, time, u, v, rho, Yl, T, h, delta)
             rhoY(i,j) = rho*Yl(id)
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
    !     ::::: left side
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
    !     ::::: right side
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

  end subroutine press_fill

end module prob_2D_module
