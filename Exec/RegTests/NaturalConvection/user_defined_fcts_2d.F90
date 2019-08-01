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

      use PeleLM_F,  only: pphys_getP1atm_MKS
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : pamb
      use probdata_module, only : Tc, Th 
      
      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      REAL_T rho_temp(1), h_temp(1), T_temp(1), Patm
      integer dir, norm  ! This specify the direction and orientation of the face
      logical getuv
      integer dimloc(2)
      data dimloc / 1,  1 /

      integer n

      if (dir == 1) then
        if (norm == 1) then
          T_temp(1) = Th
        elseif (norm == -1) then
          T_temp(1) = Tc
        else
          call bl_abort('Wring norm in bcfunction')
        endif
        Yl(0) = 0.233d0
        Yl(1) = 0.767d0

        Patm = pamb / pphys_getP1atm_MKS()

        call pphys_RHOfromPTY(dimloc, dimloc, &
                              rho_temp(1), DIMARG(dimloc), DIMARG(dimloc), &
                              T_temp(1),   DIMARG(dimloc), DIMARG(dimloc), &
                              Yl, DIMARG(dimloc), DIMARG(dimloc), Patm)
        call pphys_HMIXfromTY(dimloc, dimloc, &
                              h_temp(1),   DIMARG(dimloc), DIMARG(dimloc), &
                              T_temp(1),   DIMARG(dimloc), DIMARG(dimloc), &
                              Yl, DIMARG(dimloc), DIMARG(dimloc))
        rho = rho_temp(1)
        h = h_temp(1)
        T = T_temp(1)

         
        if (getuv .eqv. .TRUE.) then
            
          u = zero
          v = zero
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

