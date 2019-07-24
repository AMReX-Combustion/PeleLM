#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_3d_module

implicit none
  
  private
  
  public :: bcfunction, zero_visc

contains
  


!!-----------------------

  subroutine bcfunction(x,y,z,dir,norm,time,u,v,w,rho,Yl,T,h,dx,getuvw) &
                        bind(C, name="bcfunction")

      use mod_Fvar_def, only : dim
       
      implicit none

      REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(dim)
      integer dir, norm  ! This specify the direction and orientation of the face
      logical getuvw

      print *,'You are imposing Dirichlet conditions'
      print *,'you need to initialize boundary condition function'
      print *,'in user_defined_ftts_2d.F90'
      call bl_abort('')

  end subroutine bcfunction

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

end module user_defined_fcts_3d_module

