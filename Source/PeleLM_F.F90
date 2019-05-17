
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#define SDIM 2

module PeleLM_F

  use, intrinsic :: iso_c_binding

  use fuego_chemistry

  implicit none

  integer :: use_eg

  private

  public :: set_scal_numb, get_typical_vals, set_typical_vals, &
            set_ht_visc_common, init_typcals_common, get_pamb, &
            get_closed_chamber, get_dpdt, set_common, active_control, &
            pphys_calc_src_sdc, pphys_getP1atm_MKS, &
            pphys_get_spec_name2, pphys_TfromHYpt

contains

! Init/Close PelePhysics network, transport, reaction. Similar to PeleC.  

  subroutine pphys_network_init() bind(C, name="pphys_network_init")                                                                                         

     use network, only: network_init
     
     call network_init()
       
  end subroutine pphys_network_init 

  subroutine pphys_network_close() bind(C, name="pphys_network_close")

     use network, only: network_close

     call network_close()

  end subroutine pphys_network_close

  subroutine pphys_transport_init(ieg) bind(C, name="pphys_transport_init")

     use transport_module, only: transport_init

     implicit none
     integer(c_int), intent(in   ) :: ieg

     call transport_init()
     use_eg = ieg

  end subroutine pphys_transport_init  

  subroutine pphys_transport_close() bind(C, name="pphys_transport_close")

      use transport_module, only: transport_close

      call transport_close()

  end subroutine pphys_transport_close

!  subroutine pphys_reactor_close() bind(C, name="pphys_reactor_close")
!
!      use reactor_module, only: reactor_close
!
!      call reactor_close()
!
!  end subroutine pphys_reactor_close

subroutine plm_extern_init(name,namlen) bind(C, name="plm_extern_init")

  ! initialize the external runtime parameters in
  ! extern_probin_module
  !use extern_probin_module, only: runtime_init 

  integer :: namlen
  integer :: name(namlen)

  call runtime_init(name,namlen)

end subroutine plm_extern_init

  subroutine pphys_get_num_spec(nspec_out) bind(C, name="pphys_get_num_spec")

      use network, only : nspec

      implicit none
      integer(c_int), intent(out) :: nspec_out

      nspec_out = nspec

  end subroutine pphys_get_num_spec  

  subroutine pphys_get_spec_name(spec_names_out,ispec,len) &
             bind(C, name="pphys_get_spec_name")

      use network, only : spec_names

      implicit none
      integer(c_int), intent(in   ) :: ispec
      integer(c_int), intent(inout) :: len
      integer(c_int), intent(inout) :: spec_names_out(len)

      ! Local variables
      integer   :: i

      len = len_trim(spec_names(ispec+1))

      do i = 1,len
         spec_names_out(i) = ichar(spec_names(ispec+1)(i:i))
      end do

  end subroutine pphys_get_spec_name  

  subroutine pphys_get_spec_name2(name, j)
  
    implicit none

#include "cdwrk.H"

    integer i, j
    integer coded(maxspnml), len
    character*(maxspnml) name

    print *, "In 'pphys_get_spec_name2' bef call to pphys_getckspecname",maxspnml, len

    len = pphys_getckspecname(j, coded)
    do i = 1, maxspnml
      name(i:i) = ' '
    end do
    do i = 1, len
      name(i:i) = char(coded(i))
    end do
    
  end subroutine pphys_get_spec_name2

  integer function pphys_getckspecname(i, coded)
  
    implicit none
      
#include "cdwrk.H"

    integer i
    integer coded(*)
    integer names(maxspec*maxspnml)
    integer j, str_len
    str_len = 0
    call cksyms(names, maxspnml)
    do j = 1, maxspnml
      coded(j) = names(maxspnml*(i-1)+j)
    end do
    do j = 1, maxspnml
      if (coded(j).eq.ICHAR(' ')) then
        str_len = j
        exit
      endif 
    end do
    pphys_getckspecname = str_len - 1
 
  end function pphys_getckspecname

  function pphys_getRuniversal() bind(C, name="pphys_getRuniversal") result(RUNIV)

    implicit none
    double precision Ruc, Pa, RUNIV

    call CKRP(RUNIV,Ruc,Pa)
!     1 erg/(mole.K) = 1.e-4 J/(kmole.K)
    RUNIV = RUNIV*1.d-4

  end function pphys_getRuniversal

  function pphys_getP1atm_MKS() bind(C, name="pphys_getP1atm_MKS") result(P1ATM)

    implicit none
    double precision Ru, Ruc, P1ATM

    call CKRP(Ru,Ruc,P1ATM)
!     1 N/(m.m) = 0.1 dyne/(cm.cm)
    P1ATM = P1ATM*1.d-1

  end function pphys_getP1atm_MKS

  function pphys_numReactions() bind(C, name="pphys_numReactions") result(NR)

    implicit none
    integer Nelt,Nspec,NR,Nfit
      
    call CKINDX(Nelt,Nspec,NR,Nfit)

  end function pphys_numReactions

  subroutine pphys_set_verbose_vode() bind(C, name="pphys_set_verbose_vode")

    use vode_module, only: verbose

    implicit none

    verbose = 5

  end subroutine pphys_set_verbose_vode

  subroutine pphys_calc_src_sdc(N, TIME, TEMP, Z, ZP) bind(C, name="pphys_calc_src_sdc")
!
!     Variables in Z are:  Z(1:K) = rhoY(K) [MKS]
!                          Z(K+1) = RhoH    [MKS]
    use network, only : nspec

    implicit none

    integer(c_int), intent(inout)   :: N
    double precision, intent(in   ) :: TIME
    double precision, intent(inout) :: TEMP
    double precision, intent(in   ) :: Z(nspec+1)
    double precision, intent(out  ) :: ZP(nspec+1)

    double precision WDOT_CGS(nspec)
    double precision Y(nspec), CONC_CGS(nspec), MWT(nspec)
    double precision RHO_MKS, RINV_MKS, THFAC, HMIX_MKS, HMIX_CGS
    integer :: lierr, K


    ! RHO MKS
    RHO_MKS  = sum(Z(1:nspec))
    RINV_MKS = 1.d0 / RHO_MKS
    ! MW CGS
    call CKWT(MWT);
      
    do K=1,nspec
      CONC_CGS(K) = Z(K)/MWT(K)*1.d-3
      Y(K) = Z(K) * RINV_MKS
      !print *," Y, WT ", Z(K), 1./MWT(K)
    enddo

    HMIX_MKS = (Z(nspec+1) + 0.0d0*TIME) * RINV_MKS
    HMIX_CGS = HMIX_MKS * 1.0d4
    call get_t_given_hY(HMIX_CGS, Y, TEMP, lierr);
    call CKWC(TEMP,CONC_CGS,WDOT_CGS)

    ZP(Nspec+1) = 0.0d0
    THFAC = 1.d3
    do k= 1, Nspec
      ZP(k) = WDOT_CGS(k) * MWT(k) * THFAC + 0.0d0
      !print *," RHO, C(CGS), H, T",RHO_MKS, CONC_CGS(k), HMIX_MKS, TEMP
      !print *," wdot(CGS), wdot", WDOT_CGS(k), ZP(k)
    end do
      
  end subroutine pphys_calc_src_sdc

!------------------------------------  

  subroutine set_scal_numb(DensityIn, TempIn, TracIn, RhoHIn, &
                           FirstSpecIn, LastSpecIn) &
                           bind(C, name="set_scal_numb")

    implicit none

#include <htdata.H>

    integer DensityIn, TempIn, TracIn, RhoHIn, FirstSpecIn, LastSpecIn

!
! ::: Remove SPACEDIM from the counter, since those spots contain the
! ::: velocity, and our INITDATA function below fills the scalar state
! ::: However, add one since the C++ is 0-based      
!     
    Density = DensityIn - BL_SPACEDIM + 1
    Temp = TempIn - BL_SPACEDIM + 1
    Trac = TracIn - BL_SPACEDIM + 1
    RhoH = RhoHIn - BL_SPACEDIM + 1
    FirstSpec = FirstSpecIn - BL_SPACEDIM + 1
    LastSpec = LastSpecIn - BL_SPACEDIM + 1

  end subroutine set_scal_numb

!------------------------------------------

  subroutine get_typical_vals(typ_vals,nVals)bind(C, name="get_typical_vals")

    use network, only : nspec

    implicit none

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>

    integer nVals,n,nVals1
    REAL_T typ_vals(nVals)
    nVals1 = nVals-BL_SPACEDIM
!     Note: typical values are defaulted to zero, and may be left that way

    if (Density.gt.nVals1 & 
          .or. Temp.gt.nVals1 &
          .or. RhoH.gt.nVals1 &
          .or. Trac.gt.nVals1 &
          .or. LastSpec.gt.nVals) then

      call bl_pd_abort('cannot write typical values')
    endif

    do n=1,BL_SPACEDIM
      typ_vals(n) = typVal_Vel
    enddo

    typ_vals(Density+BL_SPACEDIM) = typVal_Density
    typ_vals(Temp+BL_SPACEDIM)    = typVal_Temp
    typ_vals(RhoH+BL_SPACEDIM)    = typVal_RhoH
    typ_vals(Trac+BL_SPACEDIM)    = typVal_Trac
    do n=1,Nspec
      typ_vals(FirstSpec+n-1+BL_SPACEDIM) = typVal_Y(n)
    enddo

  end subroutine get_typical_vals

!-------------------------------------------------

  subroutine set_typical_vals(typ_vals,nVals)bind(C, name="set_typical_vals")

    use network, only : nspec

    implicit none

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>

    integer nVals,n,nVals1
    REAL_T typ_vals(nVals)
    nVals1 = nVals-BL_SPACEDIM
!     Note: typical values are defaulted to zero, and may be left that way

    if (Density.gt.nVals1 & 
          .or. Temp.gt.nVals1 &
          .or. RhoH.gt.nVals1 &
          .or. Trac.gt.nVals1 &
          .or. LastSpec.gt.nVals) then
      call bl_pd_abort('cannot write typical values')
    endif

    do n=1,BL_SPACEDIM
      typVal_Vel = typ_vals(n)
    enddo

    typVal_Density = typ_vals(Density+BL_SPACEDIM)
    typVal_Temp    = typ_vals(Temp+BL_SPACEDIM)
    typVal_RhoH    = typ_vals(RhoH+BL_SPACEDIM)
    typVal_Trac    = typ_vals(Trac+BL_SPACEDIM)

    do n=1,Nspec
      typVal_Y(n) = typ_vals(FirstSpec+n-1+BL_SPACEDIM)
    enddo

  end subroutine set_typical_vals

!------------------------------------------

  subroutine set_ht_visc_common(muIsVar,     muVal, &
                                lambdaIsVar, lambdaVal, &
                                rhoDIsVar,   rhoDVal, &
                                prandtl, schmidt, unityLe) &
                                bind(C, name="set_ht_visc_common") 

    implicit none
    integer muIsVar, lambdaIsVar, rhoDIsVar
    REAL_T muVal, lambdaVal, rhoDVal, prandtl, schmidt
    integer unityLe

#include <visc.H>

    if (muIsVar .EQ. 1) then
      use_constant_mu = .false.
      constant_mu_val = - one
    else
      use_constant_mu = .true.
      constant_mu_val = muVal
    end if

    if (lambdaIsVar .EQ. 1) then
      use_constant_lambda = .false.
      constant_lambda_val = - one
    else
      use_constant_lambda = .true.
      constant_lambda_val = lambdaVal
    end if

    if (rhoDIsVar .EQ. 1) then
      use_constant_rhoD = .false.
      constant_rhoD_val = - one
    else
      use_constant_rhoD = .true.
      constant_rhoD_val = rhoDVal
    end if

    Pr = prandtl
    Sc = schmidt
    LeEQ1 = unityLe .ne. 0
    thickFacTR = one

  end subroutine set_ht_visc_common

!----------------------------------------

  subroutine init_typcals_common()&
             bind(C, name="init_typcals_common")

    implicit none

#include <cdwrk.H>
#include <conp.H>

    typVal_Density = zero
    typVal_Temp    = zero
    typVal_RhoH    = zero
    typVal_Trac    = zero
    typVal_Y       = zero
    typVal_YMAX    = one
    typVal_YMIN    = 1.d-6

  end subroutine init_typcals_common
      
!-----------------------------------------------------------------------

  subroutine get_pamb(pambout)bind(C, name="get_pamb")

    implicit none

#include <htdata.H>

    REAL_T pambout

    pambout = pamb

  return

  end subroutine get_pamb

!-----------------------------------------------------------------------

  subroutine get_closed_chamber(closed_chamberout)bind(C, name="get_closed_chamber")

    implicit none

#include <htdata.H>

    integer closed_chamberout

    closed_chamberout = closed_chamber

  return

  end subroutine get_closed_chamber
      
!-----------------------------------------------------------------------

  subroutine get_dpdt(dpdt_factorout)bind(C, name="get_dpdt")

    implicit none

#include <htdata.H>

    REAL_T dpdt_factorout

    dpdt_factorout = dpdt_factor

  return

  end subroutine get_dpdt

!=======================================================================

  subroutine set_common(time1,iteration1)bind(C, name="set_common")

    implicit none

    REAL_T time1
    integer iteration1

#include <timedata.H>

    time = time1
    iteration = iteration1

  return
  
  end subroutine set_common

!------------------------------------------------------

  subroutine active_control(coft,time,dt,myproc,step,restart,usetemp)bind(C, name="active_control")

    implicit none

#include <cdwrk.H>
#include <probdata.H>
#include <bc.H>

!
! Just stuff in the calling sequence.
!

    REAL_T coft,time,dt
    integer myproc,step,restart,usetemp

!
! ACTIVE_CONTROL_IS_USABLE should be defined in your probdata.H
! if you want to call FORT_ACTIVECONTROL.  This is how we enforce
! that all the necessary variables get defined and included in
! the proper problem-specific common blocks.
!

#if !defined(ACTIVE_CONTROL_IS_USABLE)

    call bl_abort('FORT_ACTIVECONTROL is NOT enabled')

#else

    REAL_T slocal,V_new,dVmax,dVmin
    !REAL_T vslope,
    integer ierr
    REAL_T r1,r2,r3,r4,r5,r6,r7
    REAL_T alpha,xsmb,vpmax,exp1 
    REAL_T rhs1,rhs2,vt_tay,vtt_tay,velintegral,sest_test
    integer i1
    integer nfilled,ifill
    logical found_it
    save nfilled

    if(step .eq. 0 )nfilled = navg_pnts+1

!     print *, " entering control ", restart,step,nfilled

    if (restart.ne.0) then
      nfilled = navg_pnts+1
      open(13,file=ac_hist_file,form='formatted', &
             status='old',iostat=ierr)

      found_it = .false.
!        print *, " opening file ", ierr
      if (ierr .eq. 0) then
        if (myproc.eq.0) then
          print*, 'Setting active control from history file ...'
        endif
        rewind(13)

        do

!
!                 This read(13) must correspond to the below write(13)
!
          read(13,1000,iostat=ierr) i1,r1,r2,r3,r4,r5,r6,r7
!              print *," found stuff",i1,r1,r2,r3,r4,r5,r6,r7
!              print *, "ierr ", ierr
          if (ierr.ne.0) goto 100
          if(step-i1.ge.0 .and. step-i1 .le. navg_pnts)then
            nfilled = nfilled -1
            ifill = step-i1
            time_points(ifill) = r1
            vel_points(ifill) = r2
            cntl_points(ifill) = r7
          endif
          if (i1.eq.step) then
            found_it = .true.
            V_in = r2
!             tbase_control = r3
            tbase_control = r1
            zbase_control = r4
            dV_control = r5
            sest = r6
            coft_old = r7
!               print *," found it " ,V_in
          endif
        enddo

      else

        nfilled = navg_pnts+1

        if (myproc.eq.0) then
          open(13,file=ac_hist_file,form='formatted', status='new')
        endif

      endif

 100      if (found_it .eqv. .false.) then

            if (myproc.eq.0) then
              print*, 'Setting active control to defaults ...'
            endif

        end if
        if (myproc.eq.0) then
          do ifill=0,navg_pnts
            print*,'data',ifill,time_points(ifill),vel_points(ifill),cntl_points(ifill)
          enddo
          call flush(6)
        endif
        close(13)
        return

      else

    end if

!     print *, " entering loop ", sest

    if (usetemp.eq.0 .and. coft_old .lt. zero) coft_old = coft

    zbase_control = zbase_control + V_in*dt + dV_control*dt**2.0d0
    V_in_old = V_in
    V_in = V_in + dt*dV_control

    slocal = half*(V_in_old + V_in) - (coft - coft_old)/(dt*scale_control)
 
    do ifill = navg_pnts,1,-1
      time_points(ifill) = time_points(ifill-1)
      vel_points(ifill) = vel_points(ifill-1)
      cntl_points(ifill) = cntl_points(ifill-1)
    enddo

    time_points(0) = time
    vel_points(0) = V_in
    cntl_points(0) = coft
!     sest = (one - corr)*sest + corr*slocal
!     print *," starting main control ", nfilled,sest,slocal
 
    if(nfilled .le. 0)then
!           if (myproc.eq.0) then
!                do ifill=0,navg_pnts
!                   print*,'data',ifill,time_points(ifill),vel_points(ifill),cntl_points(ifill)
!                enddo
!                call flush(6)
!           endif
      velintegral = 0.d0
      do ifill = 1,navg_pnts
        velintegral = velintegral+0.5d0*(vel_points(ifill-1)+vel_points(ifill)) &
                     *(time_points(ifill-1)-time_points(ifill))
      enddo

      sest = (velintegral - (cntl_points(0) - cntl_points(navg_pnts)) &
                                                    /scale_control)   &
                               /(time_points(0)-time_points(navg_pnts))

!        if (myproc.eq.0) then
!           print *,' step, sest_test ', step,sest_test
!           print *,'velintegral ',velintegral,velintegral/(time_points(0)-time_points(navg_pnts))
!           print *,'movement ',(cntl_points(0) - cntl_points(navg_pnts))
!           call flush(6)
!        endif
!     print *," first ", nfilled,sest,slocal

    else
      nfilled = nfilled - 1
      if(step.ne. 0)then
        sest = (one - corr)*sest + corr*slocal
      endif
!     print *," second ", nfilled,sest,slocal
    endif
!     print *," third ", nfilled,sest,slocal

#if 1
!    linear
!     vslope = two*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control

!     V_new = V_in + dt*vslope

!    quadratic 1
!     vslope = 3.d0*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control
!     V_new = V_in + (dt-0.5d0*dt**2/tau_control)*vslope

!    quadratic 2
    rhs2 = 2.d0*((cfix-coft)/(scale_control*tau_control) + sest - V_in)/tau_control
    rhs1 = (sest - V_in)/tau_control

    vt_tay = 3.d0*rhs2 - 2.d0*rhs1
    vtt_tay = 6.d0*(rhs1-rhs2)/tau_control
    V_new = V_in + dt*vt_tay + 0.5d0*dt**2*vtt_tay


    dVmax = changeMax_control * one
    dVmin = changeMax_control * max(one,V_in)

    if(myproc.eq.0)then
      write(6,*)" in contol", coft,cfix,scale_control,tau_control, sest,V_in
      write(6,*)" in contol results",rhs1,rhs2,vt_tay,vtt_tay, V_new
      write(6,*)" in control limits",dVmax,dVmin
    endif

#else

    vpmax = max(V_in,V_in_old,one)*changeMax_control/dt
    xsmb = -(cfix-coft)/scale_control

    alpha = sqrt(abs(xsmb)/vpmax)

    exp1 = dexp(-2.d0*tau_control/alpha)

    V_new = sest+ (xsmb)/alpha* &
        (exp1*dexp(dt/alpha)-dexp(-dt/alpha))/  &
           (1.d0+exp1)


    dVmax = changeMax_control * max(one,V_in)

#endif

    V_new = MAX(V_new,zero)
    V_new = MIN(MAX(V_new,V_in-dVmin),V_in+dVmax)
    V_new = MAX(zero,V_new)
    V_new = Min(V_new,controlVelMax)

    if(myproc.eq.0)then
      write(6,*)" after limiting",controlVelMax,V_new
    endif

    tbase_control = time
    dV_control = (V_new - V_in)/dt

    if (myproc.eq.0) then
      print *
      print *,'****************** control:', scale_control
      print *,'time,dt,coft,cfix,V_new:',time,dt,coft,cfix,V_new
      print *,'changeMax_control: ', changeMax_control
      print *,'V_in,dVmax,dVmin ', V_in,dVmax,dVmin

#if 1
!  vslope not defined yet         print *,'vslope,sest,cfix,slocal:',vslope,sest,cfix,slocal
#else
      print *,'alpha,sest,cfix,slocal:',alpha,sest,cfix,slocal
#endif
      print *,'coft_old,V_in_old:',coft_old,V_in_old
      call flush(6)
    endif

    coft_old = coft

    if (myproc.eq.0) then
      open(13,file=ac_hist_file,form='formatted',position='append')
!        write(13,1000) step,time,V_in,tbase_control,zbase_control,
      write(13,1000) step,time,V_in,slocal,zbase_control, &
             dV_control,sest,coft_old
      close(13)
    endif

1000 format(i7,7g26.18)

#endif /*ACTIVE_CONTROL_IS_USABLE*/

  end subroutine active_control

!-------------------------------------

  subroutine pphys_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)&
                  bind(C, name="pphys_TfromHYpt")

      implicit none

      REAL_T T,Y(*),Hin,errMax

      integer NiterMAX,Niter,NiterDAMP,ihitlo,ihithi
      REAL_T  T0,cp,dH
      REAL_T  res(0:NiterMAX-1),dT, Htarg,HMIN,cpMIN,HMAX,cpMAX
      logical converged, soln_bad, stalled
      REAL_T  H, old_T, old_H, Tsec, Hsec

      integer, parameter :: Discont_NiterMAX = 100
      REAL_T,  parameter :: TMIN = 250.d0, TMAX = 5000.d0

      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
            T0 = T
      else
            T0 = half*(TMIN+TMAX)
            T  = T0
      end if

      NiterDAMP = NiterMAX
      Niter     = 0
      soln_bad  = .FALSE.
      Htarg     = Hin * 1.d4
      ihitlo    = 0
      ihithi    = 0

      CALL CKHBMS(T,Y,H)

      old_T = T
      old_H = H

      dH         = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged  = dH.le.errMAX
      stalled    = .false.

      do while ((.not.converged) .and. (.not.stalled) .and. (.not.soln_bad))

          CALL CKCPBS(T,Y,cp)
          dT = (Htarg - H)/cp
          old_T = T
          if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
                  T = TMAX
                  ihithi = 1
          else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
                  T = TMIN
                  ihitlo = 1
          else
                  T = T + dT
          end if
          soln_bad = (T.lt.TMIN-one) .or. (T.gt.TMAX)
          if (soln_bad) then
                  Niter = -1
                  exit
          else
                  old_H = H
                  CALL CKHBMS(T,Y,H)
                  dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
                  res(Niter) = min(dH,abs(dT))
                  Niter = Niter + 1
          end if
          converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)
          if (Niter .ge. NiterMAX) then
                  if(abs(T-1000.d0).le.1.d-3.and. dH.le.1.d-5)then
                          converged = .true.
                  else
                          Niter = -2
                          exit
                  end if
          end if

          if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
                  T = TMIN
                  CALL CKHBMS(T,Y,HMIN)
                  CALL CKCPBS(T,Y,cpMIN)
                  T=TMIN+(Htarg-HMIN)/cpMIN
                  converged = .true.
          end if
          if ((ihithi.eq.1).and.(H.lt.Htarg)) then 
                  T = TMAX
                  CALL CKHBMS(T,Y,HMAX)
                  CALL CKCPBS(T,Y,cpMAX)
                  T=TMAX+(Htarg-HMAX)/cpMAX
                  converged = .true.
          end if

          if (Niter .ge. NiterMAX) then
              do while (.not. stalled)
                  dT = - (H - Htarg) * (old_T - T)/(old_H - H)
                  Tsec = T + dT
                  soln_bad = (Tsec.lt.TMIN-one) .or. (Tsec.gt.TMAX)
                  if (soln_bad) then
                          Niter = -3
                          exit
                  end if
                  CALL CKHBMS(Tsec,Y,Hsec)
                  if ( (Hsec-Htarg)*(Htarg-H) .gt. 0.d0 ) then
                          old_H = H
                          old_T = T
                  end if
                  H = Hsec
                  T = Tsec
                  stalled = (2*ABS(old_T-T)/(old_T+T).le.errMAX)
                  Niter = Niter + 1
                  if (Niter.gt.NiterMAX+Discont_NiterMAX) then
                          Niter = -2
                          exit
                  endif
              end do
              converged = .true.
          end if
      end do

      if (converged) return

  end subroutine pphys_TfromHYpt

!-------------------------------------

end module PeleLM_F
