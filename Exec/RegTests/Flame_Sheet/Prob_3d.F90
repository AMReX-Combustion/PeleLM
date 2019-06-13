
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
  
  public :: amrex_probinit,setupbc, bcfunction, init_data_new_mech, init_data, &
            den_fill, adv_fill, &
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

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
      use mod_Fvar_def, only : fuelID, domnhi, domnlo, dim
      use mod_Fvar_def, only : ac_hist_file, cfix, changemax_control, &
                               coft_old, controlvelmax, corr, dv_control, &
                               h_control, navg_pnts, scale_control, sest, &
                               tau_control, tbase_control, V_in, v_in_old, zbase_control, &
                               pseudo_gravity
      use probdata_module, only : standoff, pertmag, rho_bc, Y_bc
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
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, maxspec, maxspnml, V_in
    use probdata_module, only : standoff, Y_bc, T_bc, u_bc, v_bc, w_bc, rho_bc, h_bc
    use probdata_module, only : bcinit
  
    implicit none
      
      REAL_T Patm, pmf_vals(maxspec+3)
      REAL_T Xt(maxspec), Yt(maxspec), loc
      integer n, b(3)
      data  b / 1, 1, 1 /

      Patm = pamb / pphys_getP1atm_MKS()

!     Take fuel mixture from pmf file
      loc = (domnlo(3)-standoff)*100.d0
      call pmf(loc,loc,pmf_vals,n)
      if (n.ne.Nspec+3) then
        call bl_pd_abort('INITDATA: n(pmf) .ne. Nspec+3')
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
      v_bc = zero
      if (V_in .lt. 0) then
        w_bc = pmf_vals(2)*1.d-2
      else
        w_bc = V_in
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

  subroutine bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,dx,getuvw) &
                        bind(C, name="bcfunction")

      use network,   only: nspec
      use mod_Fvar_def, only : dim
      use mod_Fvar_def, only : dv_control, tbase_control, V_in, f_flag_active_control
      use probdata_module, only : bcinit, rho_bc, Y_bc, T_bc, h_bc, w_bc
      
      implicit none

      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(dim)
      logical getuvw

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
         
      if (getuvw .eqv. .TRUE.) then
            
        u = zero
        v = zero
        if (f_flag_active_control == 1) then                
          w =  V_in + (time-tbase_control)*dV_control
        else 
          w = w_bc
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
      use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
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

      integer i, j, k, n
      REAL_T Patm
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               scal(i,j,k,Trac) = zero
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
      use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : domnhi, domnlo, maxspec, maxspnml 
      use probdata_module, only : standoff, pertmag

      
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

      integer i, j, k, n, nPMF
      REAL_T x, y, z, Yl(maxspec), Xl(maxspec), Patm
      REAL_T pmf_vals(maxspec+3), y1, y2
      REAL_T pert,Lx,Ly

         do k = lo(3), hi(3)
            z = (float(k)+.5d0)*delta(3)+domnlo(3)
            do j = lo(2), hi(2)
               y = (float(j)+.5d0)*delta(2)+domnlo(2)
               do i = lo(1), hi(1)
                  x = (float(i)+.5d0)*delta(1)+domnlo(1)
                  
                  pert = 0.d0
                  if (pertmag .gt. 0.d0) then
                     Lx = domnhi(1) - domnlo(1)
                     Ly = domnhi(2) - domnlo(2)
                     pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)             * sin(2*Pi*5*y/Ly) &
                                   + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)   * sin(2*Pi*4*(y-.0053765)/Ly) &
                                   + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) * sin(2*Pi*3*(y-.02137)/Ly)  &
                                   + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)     * sin(2*Pi*6*(y-.018)/Ly)  &
                                               + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
                  endif

                  y1 = (z - standoff - 0.5d0*delta(3) + pert)*100.d0
                  y2 = (z - standoff + 0.5d0*delta(3) + pert)*100.d0
                  
                  call pmf(y1,y2,pmf_vals,nPMF)               
                  if (nPMF.ne.Nspec+3) then
                     call bl_abort('INITDATA: n .ne. Nspec+3')
                  endif
                  
                  scal(i,j,k,Temp) = pmf_vals(1)
                  do n = 1,Nspec
                     Xl(n) = pmf_vals(3+n)
                  end do 
                  
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
                       
       use mod_Fvar_def, only : domnlo, maxspec, dim
              
      implicit none

      integer DIMDEC(den), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  den(DIMV(den))

      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : dim
      
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     adv(DIMV(adv))
      integer    bc(dim,2)

      integer    i,j,k
      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(temp), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  temp(DIMV(temp))
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(rhoh), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoh(DIMV(rhoh))
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(vel), bc(dim,2,dim)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  vel(DIMV(vel),dim)

      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
                  endif

                  if (bc(1,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(1,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
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
                  endif

                  if (bc(1,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(1,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
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
                  endif

                  if (bc(2,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(2,1,3).eq.EXT_DIR) then
                        vel(i,j,k,3) = w
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
                  endif

                  if (bc(2,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(2,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
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
                  endif

                  if (bc(3,1,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(3,1,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
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
                  endif

                  if (bc(3,2,2).eq.EXT_DIR) then
                     vel(i,j,k,2) = v
                  endif

                  if (bc(3,2,3).eq.EXT_DIR) then
                     vel(i,j,k,3) = w
                  endif
               enddo
            enddo
         enddo
      endif

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(xvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  xvel(DIMV(xvel))

      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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

      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  xvel(i,j,k) = u
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
                  xvel(i,j,k) = u
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
                  xvel(i,j,k) = u
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
                  xvel(i,j,k) = u
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
                  xvel(i,j,k) = u
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
                  xvel(i,j,k) = u
               enddo
            enddo
         enddo
      endif
      
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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(yvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  yvel(DIMV(yvel))

      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
      
      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  yvel(i,j,k) = v
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
                  yvel(i,j,k) = v
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
                  yvel(i,j,k) = v
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
                  yvel(i,j,k) = v
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
                  yvel(i,j,k) = v
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
                  yvel(i,j,k) = v
               enddo
            enddo
         enddo
      endif

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(zvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  zvel(DIMV(zvel))
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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
      
      call filcc (zvel,DIMS(zvel),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do k = lo(3),hi(3)
               z = (float(k)+.5)*delta(3)+domnlo(3)
               do j = lo(2), hi(2)
                  y = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x,y,z,time,u,v,w,rho,Yl,T,h,delta,.true.)
                  zvel(i,j,k) = w
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
                  zvel(i,j,k) = w
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
                  zvel(i,j,k) = w
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
                  zvel(i,j,k) = w
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
                  zvel(i,j,k) = w
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
                  zvel(i,j,k) = w
               enddo
            enddo
         enddo
      endif

  end subroutine FORT_ZVELFILL
  
  subroutine all_chem_fill(rhoY,DIMS(rhoY),domlo,domhi,delta, &
                           xlo,time,bc) &
                           bind(C, name="all_chem_fill")

      use network,  only: nspec
      use mod_Fvar_def, only : domnlo, maxspec, dim

      implicit none
      
      integer DIMDEC(rhoY), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoY(DIMV(rhoY),Nspec)

      integer i, j, k, n
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(maxspec), T, h

      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : domnlo, maxspec, dim
      
      implicit none

      integer DIMDEC(rhoY), bc(dim,2)
      integer domlo(dim), domhi(dim), id
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoY(DIMV(rhoY))
      
      integer i, j, k
      integer ilo, ihi, jlo, jhi, klo, khi
      REAL_T  z, y, x
      REAL_T  u, v, w, rho, Yl(0:maxspec-1), T, h

      integer lo(dim), hi(dim)

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

      use mod_Fvar_def, only : dim
      
      implicit none

      integer    DIMDEC(p)
      integer    domlo(dim), domhi(dim)
      REAL_T     dx(dim), xlo(dim), time
      REAL_T     p(DIMV(p))
      integer    bc(dim,2)

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

  subroutine FORT_MAKEFORCE(time,force,rho, &
                            DIMS(istate),DIMS(state), &
                            dx,xlo,xhi,gravity,scomp,ncomp) &
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

!
!     ::::: local variables
!
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      REAL_T  hx, hy, hz
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
            else if (n.eq.nTrac) then
!     Tracer
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
