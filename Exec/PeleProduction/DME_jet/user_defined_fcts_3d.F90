#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_3d_module

implicit none
  
  private
  
  public :: bcfunction, zero_visc, getZone

contains
  
  ! ::: -----------------------------------------------------------
      
      integer function getZone(x, y, z)
      
      use probdata_module, only : splitx, IDX_AMBIENT, IDX_FUELPIPE, IDX_VOLUME, IDX_COFLOW
      use mod_Fvar_def, only : domnlo
      
      implicit none

      REAL_T x, y, z

      getZone = IDX_VOLUME
         
      if (z.gt.domnlo(3)) then
        getZone = IDX_COFLOW
      else
        getZone = IDX_FUELPIPE
      endif
         
      end function getZone
      

! ::: -----------------------------------------------------------
      
  subroutine bcfunction(x,y,z,dir,norm,time,u,v,w,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : pamb, dim, domnlo, domnhi
      use probdata_module, only : blobr, bcinit, xfrontw, splitx, Tfrontw, &
                                  Y_bc, T_bc, u_bc, v_bc, w_bc
      use probdata_module, only : IDX_FUELPIPE, IDX_COFLOW
      
      implicit none

      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(dim)
      REAL_T rho_temp(1), h_temp(1), T_temp(1)
      REAL_T  Patm
      logical getuv
      integer dimloc(3)
      integer dir, norm  ! This specify the direction and orientation of the face

      integer n, airZone,fuelZone
      REAL_T eta, eta1, sigma, pi
      data dimloc / 1,  1,  1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      eta = 0.5d0*(1.d0 - TANH(2.d0*(ABS(y)-blobr)/Tfrontw))
      do n = 0, Nspec-1
        Yl(n) = Y_bc(n,IDX_FUELPIPE)*eta + (1.d0-eta)*Y_bc(n,IDX_COFLOW)
      end do
      T = T_bc(IDX_FUELPIPE)*eta + (1.d0-eta)*T_bc(IDX_COFLOW)
      T_temp(1) = T
         
      if (getuv .eqv. .TRUE.) then
        eta1 = 0.5d0*(1.d0 - TANH(2.d0*(ABS(y)-splitx)/xfrontw))
        u = u_bc(IDX_FUELPIPE)*eta1 + (1.d0-eta1)*u_bc(IDX_COFLOW)
        v = v_bc(IDX_FUELPIPE)*eta1 + (1.d0-eta1)*v_bc(IDX_COFLOW)
        w = w_bc(IDX_FUELPIPE)*eta1 + (1.d0-eta1)*w_bc(IDX_COFLOW)
      endif

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

      use mod_Fvar_def, only : Temp, FirstSpec, RhoH, LastSpec
      use mod_Fvar_def, only : domnlo, dim

      implicit none
      integer DIMDEC(diff)
      integer lo(dim), hi(dim)
      integer domlo(dim), domhi(dim)
      integer bc(2*dim)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(dim)
      REAL_T  problo(dim)
      


      end subroutine zero_visc

end module user_defined_fcts_3d_module

