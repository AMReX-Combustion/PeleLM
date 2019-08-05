#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_2d_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private
  
  public :: bcfunction, zero_visc

contains
  


!-----------------------

  subroutine bcfunction(x,y,dir,norm,time,u,v,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspecies
      use mod_Fvar_def, only : dv_control, tbase_control, V_in, f_flag_active_control
      use probdata_module, only : bcinit, rho_bc, Y_bc, T_bc, h_bc, v_bc
      
      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      integer dir, norm  ! This specify the direction and orientation of the face
      logical getuv

      integer n

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      if ((dir == 2).and.(norm == 1)) then
        rho = rho_bc(1)
        do n = 0, nspecies-1
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
      endif  

  end subroutine bcfunction

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
      use mod_Fvar_def, only : domnhi, domnlo
      
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

end module user_defined_fcts_2d_module

