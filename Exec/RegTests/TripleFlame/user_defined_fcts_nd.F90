#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include <PPHYS_CONSTANTS.H>
#include "mechanism.h"

module user_defined_fcts_nd_module

   use amrex_fort_module, only : dim=>amrex_spacedim
   use fuego_chemistry

implicit none

  private
  
  public :: bcfunction, zero_visc, set_Y_from_Phi, set_Y_from_ksi, set_Zst

contains

!=========================================================
!  Dimension agnostic bndy function for Dirichlet
!=========================================================

   subroutine bcfunction( x, dx, dir, norm, time, getuvw, &
                          vel, rho, Yl, T, h)&
                          bind(C, name="bcfunction")

      use mod_Fvar_def, only : dv_control, tbase_control, V_in, f_flag_active_control, pamb
      use probdata_module, only : bcinit, T_bc, u_bc, v_bc, midtanh, widthtanh
      use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      
      implicit none

! In/Out      
      REAL_T, intent(in)  :: x(3)
      REAL_T, intent(in)  :: dx(3)
      integer, intent(in) :: dir, norm  ! This specify the direction and orientation of the face
      REAL_T, intent(in)  :: time
      logical, intent(in) :: getuvw
      REAL_T, intent(out) :: vel(3)
      REAL_T, intent(out) :: rho
      REAL_T, intent(out) :: Yl(0:*)
      REAL_T, intent(out) :: T
      REAL_T, intent(out) :: h

      REAL_T  :: tanhval, vbase, ksi
      REAL_T  :: h_tmp(1), rho_tmp(1)

      integer :: n
      
      integer :: b_lo(3), b_hi(3)
      data  b_lo / 1, 1, 1 /
      data  b_hi / 1, 1, 1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      if ((dir == 2).and.(norm == 1)) then
        Patm = pamb / PP_PA_MKS
        tanhval = 0.5d0*(1.0d0+TANH((x(1)-midtanh)/widthtanh)) 
        ksi = tanhval
        call set_Y_from_Ksi(ksi,Yl(0:NUM_SPECIES-1))
        T = T_bc(1)
        call pphys_RHOfromPTY(b_lo, b_hi, &
                             rho_tmp(1), b_lo, b_hi, &
                             T_bc(1),    b_lo, b_hi, &
                             Yl(0),      b_lo, b_hi, Patm)
        call pphys_HMIXfromTY(b_lo, b_hi, &
                             h_tmp(1),   b_lo, b_hi, &
                             T_bc(1),    b_lo, b_hi, &
                             Yl(0),      b_lo, b_hi)
        rho = rho_tmp(1)
        h = h_tmp(1)

        if (getuvw .eqv. .TRUE.) then
          vel(1) = u_bc
          if (f_flag_active_control == 1 ) then               
            vbase =  V_in + (time-tbase_control)*dV_control
            vel(2) = vbase 
          else 
            vel(2) = V_in
          endif
          vel(3) = 0.0d0
        endif
      endif  

  end subroutine bcfunction

!=========================================================
!  This routine will zero out diffusivity on portions of the
!  boundary that are inflow, allowing that a "wall" block
!  the complement aperture
!
!  INPUTS/OUTPUTS:
!  
!  beta      <=> diffusivity on edges
!  b_lo(hi)   => index extent of beta array
!  lo(hi)     => region of interest, edge-based
!  domlo(hi)  => index extent of problem domain, edge-based
!  dx         => cell spacing
!  problo     => phys loc of lower left corner of prob domain
!  bc         => boundary condition flag (on orient)
!                    in BC_TYPES::physicalBndryTypes
!  idir       => which face, 0=x, 1=y
!  isrz       => 1 if problem is r-z
!  id         => index of state, 0=u
!  ncomp      => components to modify
!=========================================================

   subroutine zero_visc( beta, b_lo, b_hi, lo, hi, domlo, domhi,&
                         dx, problo, bc, idir, isrz, id, ncomp)&
                         bind(C, name="zero_visc")

      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec
      use mod_Fvar_def, only : domnhi, domnlo

      implicit none

      integer :: lo(3), hi(3)
      integer :: b_lo(3), b_hi(3)
      integer :: domlo(3), domhi(3)
      integer :: bc(2*3)
      integer :: idir, isrz, id, ncomp
      REAL_T  :: dx(3)
      REAL_T  :: problo(3)
      REAL_T, dimension(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3),*) :: beta

! Routine compiled but should be set by the user
! if there is a mix of inflox/wall at a boundary

   end subroutine zero_visc

!=========================================================
! Get mass fractions of fuel, H2, oxidizer and bath at the input phi value
!=========================================================

  subroutine set_Y_from_Phi(phi,Yt)bind(C, name="set_Y_from_Phi")
  
      use mod_Fvar_def, only : fuelID, oxidID, bathID, H2ID
      use probdata_module, only : H2_enrich

      implicit none

      REAL_T, INTENT(IN)  :: phi
      REAL_T, INTENT(OUT) :: Yt(NUM_SPECIES)

      REAL_T :: a, a_hyd
      REAL_T :: Xt(NUM_SPECIES)
      REAL_T :: a_fact(NUM_ELEMENTS)
      INTEGER :: ename(NUM_ELEMENTS*L_elem_name)
      INTEGER :: speccomp(NUM_SPECIES*NUM_ELEMENTS)
      INTEGER ::  n, found_CH
      CHARACTER(LEN=16) :: name

      Xt(:) = zero
      
!     Set "a" for computing X from phi
!     hc + a.O2 -> b.CO2 + c.H2O
      
!     ((1-alpha)*hc+alpha*H2) + a_hyd.O2 -> b_hyd.CO2 + c_hyd.H2O
!     a_hyd = (1-alpha)*a_hc + alpha*a_H2

      call cksyme(ename, L_elem_name)
      found_CH = 0
      do n = 1 , NUM_ELEMENTS
          a_fact(n) = 0.0d0
          if ( char(ename((n-1)*L_elem_name + 1)) == trim('C')) then
              a_fact(n) = 1.0d0 
              found_CH = found_CH + 1
          endif
          if ( char(ename((n-1)*L_elem_name + 1)) .eq. trim('H')) then
              a_fact(n) = 0.25d0
              found_CH = found_CH + 1
          endif
      end do
      if (found_CH /= 2) then
          call bl_abort('someting went wrong computing a_fact')
      endif

      call ckncf(speccomp)
      a = 0.0d0
      do n = 1 , NUM_ELEMENTS
          a = a + speccomp((fuelID-1)*NUM_ELEMENTS + n) * a_fact(n) 
      end do

      a_hyd = (1.0d0 - H2_enrich) * a + H2_enrich * 0.5d0

      Xt(oxidID) = 1.d0/(1.d0 + phi/a_hyd  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * (1.0d0 - H2_enrich) * Xt(oxidID) / a_hyd
      Xt(H2ID) = phi * H2_enrich * Xt(oxidID) / a_hyd
      Xt(bathID) = 1.d0 - Xt(fuelID) - Xt(H2ID) - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)
      
  end subroutine set_Y_from_Phi

!=========================================================
! Convert from mixture fraction to equivalence ratio
!=========================================================

  subroutine set_Y_from_Ksi(ksi,Yt)bind(C, name="set_Y_from_ksi")
  
      use probdata_module, only : Zst

      implicit none

      REAL_T, INTENT(IN)  :: ksi
      REAL_T, INTENT(OUT) :: Yt(NUM_SPECIES)

      INTEGER ::  n
      REAL_T :: phi
      CHARACTER(LEN=16) :: name

      phi = ksi / ( 1.0 - min(ksi,0.999999) ) * ( 1.0 - Zst ) / Zst

      call set_Y_from_Phi(phi,Yt)
      
  end subroutine set_Y_from_ksi

!=========================================================
! Compute stoichiometric mixture fraction for the fuel mixture defined by 
! the fuel name in the input file and the H2_enrich set in probin file
!=========================================================

  subroutine set_Zst( )bind(C, name="set_Zst")
  
      use mod_Fvar_def, only : fuelID, oxidID, H2ID
      use probdata_module, only : H2_enrich, Zst

      implicit none

      REAL_T :: a, a_hyd, mwt_fuel
      REAL_T :: MWT(NUM_SPECIES)
      REAL_T :: a_fact(NUM_ELEMENTS)
      INTEGER :: ename(NUM_ELEMENTS*L_elem_name)
      INTEGER :: speccomp(NUM_SPECIES*NUM_ELEMENTS)
      INTEGER ::  n, found_CH
      CHARACTER(LEN=16) :: name


      call cksyme(ename, L_elem_name)
      found_CH = 0
      do n = 1 , NUM_ELEMENTS
          a_fact(n) = 0.0d0
          if ( char(ename((n-1)*L_elem_name + 1)) == trim('C')) then
              a_fact(n) = 1.0d0 
              found_CH = found_CH + 1
          endif
          if ( char(ename((n-1)*L_elem_name + 1)) .eq. trim('H')) then
              a_fact(n) = 0.25d0
              found_CH = found_CH + 1
          endif
      end do
      if (found_CH /= 2) then
          call bl_abort('someting went wrong computing a_fact')
      endif

      call ckncf(speccomp)
      a = 0.0d0
      do n = 1 , NUM_ELEMENTS
          !print *, "elem ", char(ename((n-1)*L_elem_name + 1)), " contribution is ", speccomp((fuelID-1)*NUM_ELEMENTS + n) 
          a = a + speccomp((fuelID-1)*NUM_ELEMENTS + n) * a_fact(n) 
      end do
      print *, "a for fuel is ", a

      call CKWT(MWT)

      a_hyd = (1.0d0 - H2_enrich) * a + H2_enrich * 0.5d0
      mwt_fuel = (1.0d0 - H2_enrich) * MWT(fuelID) + H2_enrich * MWT(H2ID)

      Zst = 1.0d0 / ( 1.0d0 + MWT(OxidID) / mwt_fuel * a_hyd / 0.233 ) 

      write(*,*) " Zst = ", Zst
      IF ( H2_enrich > 0.0 ) THEN
         write(*,*) " Y("//TRIM(name)//") in fuel = ", (1.0d0 - H2_enrich) * MWT(fuelID) / mwt_fuel
         write(*,*) " Y(H2) in fuel = ", H2_enrich * MWT(H2ID) / mwt_fuel
      END IF
      
  end subroutine set_Zst

end module user_defined_fcts_nd_module
