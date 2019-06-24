#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_2d_module

implicit none
  
  private
  
  public :: bcfunction, zero_visc, getZone

contains
  
  ! ::: -----------------------------------------------------------
      
      integer function getZone(x, y)
      
      use probdata_module, only : splitx, IDX_AMBIENT, IDX_FUELPIPE, IDX_VOLUME
      implicit none

      REAL_T x, y

      getZone = IDX_VOLUME
         
         
      if (x .ge. -splitx) then
         getZone = IDX_AMBIENT
      else
         getZone = IDX_FUELPIPE
      endif
         
      end function getZone
      

! ::: -----------------------------------------------------------
      
  subroutine bcfunction(x,y,dir,norm,time,u,v,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : pamb, dim, domnlo, domnhi
      use probdata_module, only : bcinit, xfrontw, splitx, Y_bc, T_bc, u_bc, v_bc

      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      REAL_T rho_temp(1), h_temp(1), T_temp(1)
      REAL_T  Patm
      logical getuv
      integer dimloc(2)
      integer dir, norm  ! This specify the direction and orientation of the face

      integer n, airZone,fuelZone
      REAL_T eta, sigma, pi
      data dimloc / 1,  1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      fuelZone = getZone(0.d0, domnlo(2))
      airZone  = getZone(domnlo(1), domnlo(2))

      sigma = 2.5d0*xfrontw*splitx
      eta = 0.5d0 * ( tanh((x + splitx)/sigma) &
                    - tanh((x - splitx)/sigma))
                  
      do n=1,Nspec
        Yl(n-1) = Y_bc(n-1,airZone)*(1.d0-eta) &
                     + eta*Y_bc(n-1,fuelZone)
      enddo
      T = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)
      T_temp(1) = T 
!  need to define rho and h from rest of state
!           rho = rho_bc(airZone)*(1.d0-eta) + eta*rho_bc(fuelZone)
!           h = h_bc(airZone)*(1.d0-eta) + eta*h_bc(fuelZone)

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

        u = u_bc(airZone)*(1.d0-eta) + eta*u_bc(fuelZone)
        v = v_bc(airZone)*(1.d0-eta) + eta*v_bc(fuelZone)

        ! sinusoidal variation of inflow
        pi = 4.d0*atan(1.d0)
        v = v + 10.d0*eta*sin(2.d0*pi*x/(domnhi(1)-domnlo(1)))*sin(2.d0*pi*time/1.d-5)

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
      
      integer i, j, n, Tid, RHid, YSid, YEid, ys, ye
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, y

         Tid  = Temp      - id + dim
         RHid = RhoH      - id + dim
         YSid = FirstSpec - id + dim
         YEid = LastSpec  - id + dim
         
         do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
         do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
         ys = MAX(YSid,1)
         ye = MIN(YEid,ncomp)
         do_Y = (ye - ys + 1) .GE. 1
!     
!     Do species, Temp, rhoH
!     
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2)) &
                .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
            y = DBLE(j)*dx(2)+domnlo(2)
            j = lo(2)
            do i = lo(1), hi(1)
               
               xl = DBLE(i)*dx(1)+domnlo(1) 
               xr = (DBLE(i)+1.d0)*dx(1)+domnlo(1) 
               xh = 0.5d0*(xl+xr)
                  
                  
!                 if (do_T)  diff(i,j,Tid ) = 0.d0
!                 if (do_RH) diff(i,j,RHid) = 0.d0
                  if (do_Y) then
                     do n=ys,ye
                        diff(i,j,n) = 0.d0
                     enddo
                  endif
                     
            end do
         endif

      end subroutine zero_visc

end module user_defined_fcts_2d_module

