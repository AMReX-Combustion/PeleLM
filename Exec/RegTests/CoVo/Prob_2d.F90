
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
  
  public :: amrex_probinit, init_data_new_mech, init_data, &
            den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            FORT_XVELFILL, FORT_YVELFILL, chem_fill, press_fill, &
            FORT_MAKEFORCE, zero_visc 

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
      use mod_Fvar_def, only : dim
      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xvort, yvort, rvort, forcevort
      
      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i
 
      namelist /fortin/ meanFlowMag, meanFlowDir, T_mean, P_mean, &
                       xvort, yvort, rvort, forcevort  
      namelist /heattransin/ pamb, dpdt_factor, closed_chamber


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

      meanFlowDir = 1.0
      meanFlowMag = 1.0d0
      T_mean = 298.0d0
      P_mean = pamb
      xvort = 0.5
      yvort = 0.5
      rvort = 0.01945
      forcevort = 1.0

      read(untin,fortin)
      
      read(untin,heattransin)
 
      close(unit=untin)

!     Do some checks on COVO params     
      IF (       meanFlowDir /= 1 .AND. meanFlowDir /= -1 &
           .AND. meanFlowDir /= 2 .AND. meanFlowDir /= -2 &
           .AND. meanFlowDir /= 3 .AND. meanFlowDir /= -3  ) THEN
          WRITE(*,*) " meanFlowDir should be either: "
          WRITE(*,*) " +/-1 for x direction" 
          WRITE(*,*) " +/-2 for y direction" 
          WRITE(*,*) " +/-3 for diagonal direction" 
          WRITE(*,*) " Note: the mean flow direction(s) must be periodic "
          CALL bl_abort('Correct meanFlowDir value !')
      END IF
      
      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if

  end subroutine amrex_probinit
        
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
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, Trac, dim
      use mod_Fvar_def, only : domnlo, maxspec
      use probdata_module, only: meanFlowDir, meanFlowMag, &
                                 T_mean, P_mean, &
                                 xvort, yvort, rvort, forcevort
      
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


      integer i, j, n
      REAL_T x, y, Yl(maxspec), Patm
      REAL_T dx
      REAL_T :: dy, d_sq, r_sq, u_vort, v_vort 

      do j = lo(2), hi(2)
         y = (float(j)+.5d0)*delta(2)+domnlo(2)
         do i = lo(1), hi(1)
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            
            scal(i,j,Temp) = T_mean
            Yl(1) = 0.233
            Yl(2) = 0.767
            
            do n = 1,Nspec
               scal(i,j,FirstSpec+n-1) = Yl(n)
            end do

            scal(i,j,Trac) = 0.d0

            dx = x - xvort
            dy = y - yvort
            d_sq = dx*dx + dy*dy
            r_sq = rvort*rvort

            u_vort = -forcevort*dy/r_sq * exp(-d_sq/r_sq/two)
            v_vort = forcevort*dx/r_sq * exp(-d_sq/r_sq/two)

            SELECT CASE ( meanFlowDir )
               CASE (1)
                  vel(i,j,1) = meanFlowMag + u_vort
                  vel(i,j,2) = v_vort
               CASE (-1)
                  vel(i,j,1) = -meanFlowMag + u_vort
                  vel(i,j,2) = v_vort
               CASE (2)
                  vel(i,j,1) = u_vort
                  vel(i,j,2) = meanFlowMag + v_vort
               CASE (-2)
                  vel(i,j,1) = u_vort
                  vel(i,j,2) = -meanFlowMag + v_vort
               CASE (3)
                  vel(i,j,1) = meanFlowMag + u_vort
                  vel(i,j,2) = meanFlowMag + v_vort
               CASE (-3)
                  vel(i,j,1) = -meanFlowMag + u_vort
                  vel(i,j,2) = -meanFlowMag + v_vort
            END SELECT

         end do
      end do

      Patm = P_mean / pphys_getP1atm_MKS()

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
                              
      use mod_Fvar_def, only : maxspec, dim
      
      implicit none

      integer DIMDEC(den), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  den(DIMV(den))
     
      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)


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

      use mod_Fvar_def, only : maxspec, dim
      
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(dim), domhi(dim)
      REAL_T     delta(dim), xlo(dim), time
      REAL_T     adv(DIMV(adv))
      integer    bc(dim,2)

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)


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

      use mod_Fvar_def, only : maxspec, dim
      
      implicit none

      integer DIMDEC(temp), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  temp(DIMV(temp))
      
      call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)
      
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

      use mod_Fvar_def, only : maxspec, dim
      
      implicit none

      integer DIMDEC(rhoh), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoh(DIMV(rhoh))
      
      call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

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

  subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta, &
                            xlo,time,bc)&
                            bind(C, name="FORT_XVELFILL")
                               
      use mod_Fvar_def, only : maxspec, dim
      
      implicit none
      
      integer DIMDEC(xvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  xvel(DIMV(xvel))

      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)
      
      
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
                               
      use mod_Fvar_def, only : maxspec, dim
      
      implicit none
      
      integer DIMDEC(yvel), bc(dim,2)
      integer domlo(dim), domhi(dim)
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  yvel(DIMV(yvel))
      
      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)
      
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
                               
      use mod_Fvar_def, only :  maxspec, dim
      
      implicit none
      
      integer DIMDEC(rhoY), bc(dim,2)
      integer domlo(dim), domhi(dim), id
      REAL_T  delta(dim), xlo(dim), time
      REAL_T  rhoY(DIMV(rhoY))
      
      call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)
            
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

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                           dx,problo,bc,idir,isrz,id,ncomp) &
                           bind(C, name="zero_visc")   

      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec
      use mod_Fvar_def, only : domnhi, domnlo, dim
      
      implicit none
      integer DIMDEC(diff)
      integer lo(dim), hi(dim)
      integer domlo(dim), domhi(dim)
      integer bc(2*dim)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(dim)
      REAL_T  problo(dim)



! Routine compiled but should be set by the user
! if there is a mix of inflox/wall at a boundary



  end subroutine zero_visc


end module prob_2D_module
