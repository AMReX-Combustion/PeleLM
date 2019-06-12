
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
  
  public :: amrex_probinit,setupbc, bcfunction, init_data_new_mech, init_data, &
            den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, chem_fill, press_fill, &
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

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
      use mod_Fvar_def, only : fuelID, domnhi, domnlo, dim
      use mod_Fvar_def, only : ac_hist_file, cfix, changemax_control, &
                               coft_old, controlvelmax, corr, dv_control, &
                               h_control, navg_pnts, scale_control, sest, &
                               tau_control, tbase_control, v_in_old, zbase_control, &
                               pseudo_gravity
      use probdata_module, only : standoff, V_in, pertmag, rho_bc, Y_bc
      use probdata_module, only : flame_dir
      
      
      implicit none
      
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i,istemp
      REAL_T area

      namelist /fortin/ V_in, &
                        standoff, pertmag
      namelist /heattransin/ pamb, dpdt_factor

      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
          zbase_control, pseudo_gravity, istemp,corr,controlVelMax,navg_pnts

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
      dV_control = zero
      tbase_control = zero
      h_control = -one
      pseudo_gravity = 0
      istemp = 0
      navg_pnts = 10

      read(untin,fortin)
      
!     Initialize control variables that depend on fortin variables
      V_in_old = V_in
      
      read(untin,heattransin)
 
      read(untin,control)
      close(unit=untin)

!     Set up boundary functions
      call setupbc()
      
      area = 1.d0
      do i=1,dim
        if (flame_dir /= i) then
         area = area*(domnhi(i)-domnlo(i))
        endif
      enddo
      scale_control = Y_bc(fuelID-1) * rho_bc(1) * area

      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif


      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
         write(6,control)
      end if

  end subroutine amrex_probinit
  
!------------------------------------
  
  subroutine setupbc()bind(C, name="setupbc")

    use network,   only: nspec
    use PeleLM_F, only: pphys_getP1atm_MKS
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, maxspec, maxspnml
    use probdata_module, only : standoff, V_in, Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit
  
    implicit none

    REAL_T Patm, pmf_vals(maxspec+3)
    REAL_T Xt(maxspec), Yt(maxspec), loc
    
    integer n
    integer b(2)
    data  b / 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()
             
  !     Take fuel mixture from pmf file
        loc = (domnlo(2)-standoff)*100.d0
        call pmf(loc,loc,pmf_vals,n)
        if (n.ne.Nspec+3) then
          call bl_pd_abort('setupbc: n(pmf) .ne. Nspec+3')
        endif
              
        do n = 1,Nspec
          Xt(n) = pmf_vals(3+n)
        end do 
              
        CALL CKXTY (Xt, Yt)
  
        do n=1,Nspec
          Y_bc(n-1) = Yt(n)
        end do
        
        T_bc = pmf_vals(1)
        u_bc = zero
        if (V_in .lt. 0) then
          v_bc = pmf_vals(2)*1.d-2
        else
          v_bc = V_in
        endif
              

!     Set density and hmix consistent with data

      call pphys_RHOfromPTY(b, b, &
                           rho_bc(1), DIMARG(b), DIMARG(b), &
                           T_bc(1),   DIMARG(b), DIMARG(b), &
                           Y_bc(0), DIMARG(b), DIMARG(b), Patm)
      call pphys_HMIXfromTY(b, b, &
                           h_bc(1),   DIMARG(b), DIMARG(b), &
                           T_bc(1),   DIMARG(b), DIMARG(b), &
                           Y_bc(0), DIMARG(b), DIMARG(b))

    bcinit = .true.

  end subroutine setupbc

      
!-----------------------

  subroutine bcfunction(x,y,time,u,v,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspec
      use mod_Fvar_def, only : dim
      use mod_Fvar_def, only : dv_control, tbase_control, f_flag_active_control
      use probdata_module, only : V_in, bcinit, rho_bc, Y_bc, T_bc, h_bc, v_bc
      
      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      logical getuv

      integer n

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      rho = rho_bc(1)
      do n = 0, Nspec-1
        Yl(n) = Y_bc(n)
      end do
      T = T_bc(1)
      h = h_bc(1)
         
      if (getuv .eqv. .TRUE.) then
            
        u = zero
        if (f_flag_active_control == 1) then               
          v =  V_in + (time-tbase_control)*dV_control
        else 
          v = v_bc
        endif
      endif
         

  end subroutine bcfunction

! ::: -----------------------------------------------------------
      
  subroutine init_data_new_mech (level,time,lo,hi,nscal, &
          vel,scal,DIMS(state),press,DIMS(press), &
          delta,xlo,xhi)&
          bind(C, name="init_data_new_mech")
          

      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      
      implicit none
      integer  level, nscal
      integer  lo(dim), hi(dim)
      integer  DIMDEC(state)
      integer  DIMDEC(press)
      REAL_T   xlo(dim), xhi(dim)
      REAL_T   time, delta(dim)
      REAL_T   vel(DIMV(state),dim)
      REAL_T   scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
  
      integer i, j, n
      REAL_T Patm
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            scal(i,j,Trac) = zero
         end do
      end do
 
      Patm = pamb / pphys_getP1atm_MKS()
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

  subroutine init_data(level,time,lo,hi,nscal, &
     	 	                   vel,scal,DIMS(state),press,DIMS(press), &
                           delta,xlo,xhi) &
                           bind(C, name="init_data")
                              

      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : bathID, domnhi, domnlo, maxspec, maxspnml 
      use probdata_module, only : standoff, pertmag

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
      integer nPMF

      integer i, j, n
      REAL_T x, y, Yl(maxspec), Xl(maxspec), Patm
      REAL_T pmf_vals(maxspec+3), y1, y2
      REAL_T pert,Lx

!      write(6,*)" made it to initdata"
      if (bathID.lt.1 .or. bathID.gt.Nspec) then
         call bl_pd_abort()
      endif

      do j = lo(2), hi(2)
        y = (float(j)+.5d0)*delta(2)+domnlo(2)
        do i = lo(1), hi(1)
          x = (float(i)+.5d0)*delta(1)+domnlo(1)
               
          pert = 0.d0
          if (pertmag .gt. 0.d0) then
            Lx = domnhi(1) - domnlo(1)
            pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx) &
                      + 1.023 * sin(2*Pi*2*(x-.004598)/Lx) &
                         + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx)  &
                             + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)  &
                                  + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
          endif
                  
          y1 = (y - standoff - 0.5d0*delta(2) + pert)*100.d0
          y2 = (y - standoff + 0.5d0*delta(2) + pert)*100.d0

          call pmf(y1,y2,pmf_vals,nPMF)               
          if (nPMF.ne.Nspec+3) then
            call bl_abort('INITDATA: n .ne. Nspec+3')
          endif
               
          scal(i,j,Temp) = pmf_vals(1)
          do n = 1,Nspec
            Xl(n) = pmf_vals(3+n)
          end do 
               
          CALL CKXTY (Xl, Yl)
               
          do n = 1,Nspec
            scal(i,j,FirstSpec+n-1) = Yl(n)
          end do

          scal(i,j,Trac) = 0.d0

          vel(i,j,1) = 0.d0
          vel(i,j,2) = pmf_vals(2)*1.d-2

        end do
      end do

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

  subroutine den_fill (den,DIMS(den),domlo,domhi,delta, &
                              xlo,time,bc) &
                              bind(C, name="den_fill")
                      
      use mod_Fvar_def, only : domnlo, maxspec, dim
              
      implicit none

      integer DIMDEC(den), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  den(DIMV(den))

      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
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

  subroutine adv_fill (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)&
                           bind(C, name="adv_fill")

      use mod_Fvar_def, only : dim
      
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     adv(DIMV(adv))
      integer    bc(dim,2)

      integer    i,j
      integer lo(dim), hi(dim)

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

      subroutine temp_fill (temp,DIMS(temp),domlo,domhi,delta, &
                              xlo,time,bc)&
                              bind(C, name="temp_fill")

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(temp), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  temp(DIMV(temp))
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
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

  subroutine rhoh_fill (rhoh,DIMS(rhoh),domlo,domhi,delta, &
                              xlo,time,bc)&
                              bind(C, name="rhoh_fill")

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(rhoh), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoh(DIMV(rhoh))
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif

  end subroutine rhoh_fill
  
!
! Fill x & y velocity at once.
!

  subroutine vel_fill (vel,DIMS(vel),domlo,domhi,delta, &
                              xlo,time,bc)&
                              bind(C, name="vel_fill")

      use mod_Fvar_def, only : dim
      
      implicit none
      
      integer DIMDEC(vel), bc(dim,2,dim)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  vel(DIMV(vel),dim)

      call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),1), &
      DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

      call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),2), &
      DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

  end subroutine vel_fill

!
! Fill all chem species at once
!

  subroutine all_chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta, &
                            xlo,time,bc)&
                            bind(C, name="all_chem_fill")

      use network,  only: nspec
      use mod_Fvar_def, only : dim
      
      implicit none

      integer DIMDEC(rhoY), bc(dim,2,Nspec)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoY(DIMV(rhoY),Nspec)

      integer n
      
      do n=1,nspec
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

  subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta, &
                            xlo,time,bc)&
                            bind(C, name="FORT_XVELFILL")
                               
      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none
      
      integer DIMDEC(xvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  xvel(DIMV(xvel))

      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif
      
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

  subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta, &
                            xlo,time,bc)&
                            bind(C, name="FORT_YVELFILL")
                               
      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none
      
      integer DIMDEC(yvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  yvel(DIMV(yvel))
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif
      
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
      
  subroutine chem_fill (rhoY,DIMS(rhoY),domlo,domhi,delta, &
                            xlo,time,bc,id ) &
                            bind(C, name="chem_fill")
                               
      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none
      
      integer DIMDEC(rhoY), bc(dim,2)
      integer domlo(dim), domhi(dim), id
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoY(DIMV(rhoY))
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x, y, time, u, v, rho, Yl, T, h, delta,.false.)
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

  subroutine press_fill (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)&
                            bind(C, name="press_fill")
  
      use mod_Fvar_def, only : dim
      
      implicit none
      
      integer    DIMDEC(p)
      integer    domlo(dim), domhi(dim)
      REAL_T     dx(dim), xlo(dim), time
      REAL_T     p(DIMV(p))
      integer    bc(dim,2)

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

!
!
! ::: -----------------------------------------------------------
!
!     This routine add the forcing terms to the momentum equation
!

  subroutine FORT_MAKEFORCE(time,force,rho, &
                               DIMS(istate),DIMS(state), &
                               dx,xlo,xhi,gravity,scomp,ncomp)&
                               bind(C,name="FORT_MAKEFORCE")

      use mod_Fvar_def, only : dv_control, pseudo_gravity, dim
      
      implicit none

      integer    DIMDEC(state)
      integer    DIMDEC(istate)
      integer    scomp, ncomp
      REAL_T     time, dx(dim)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
      REAL_T     rho    (DIMV(state))
      REAL_T     gravity

      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      REAL_T  hx, hy
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

!     Assumes components are in the following order
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
!     else to zero
         else
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = zero
               enddo
            enddo
         endif
!     Add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nYvel) = force(i,j,nYvel) + dV_control*rho(i,j)
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
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
!     Tracer
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
!     Other scalar
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif

  end subroutine FORT_MAKEFORCE

end module prob_2D_module
