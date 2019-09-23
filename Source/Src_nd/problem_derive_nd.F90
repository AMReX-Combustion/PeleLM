module problem_derive_module

  implicit none

  public

contains

  ! This is a template routine for users to compute their own derived based on the state.
  ! It will be overwritten by having a copy of this file in the user's problem setup.
  
  subroutine derUserDefined(der, der_lo,  der_hi, nder, &
                            dat, dat_lo,  dat_hi, ncomp, &
                            lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                            level,grid_no) &
                            bind(C, name="derUserDefined")

    use chemistry_module, only : nspecies                     
    use amrex_fort_module, only : dim=>amrex_spacedim

    implicit none

    integer :: lo(3), hi(3)
    integer :: der_lo(3), der_hi(3), nder
    integer :: dat_lo(3), dat_hi(3), ncomp
    integer :: domlo(3), domhi(3)
    integer :: bc(3,2,ncomp)
    double precision :: delta(3), xlo(3), time, dt
    double precision, dimension(der_lo(1):der_hi(1),der_lo(2):der_hi(2),der_lo(3):der_hi(3),nder) :: der
    double precision, dimension(dat_lo(1):dat_hi(1),dat_lo(2):dat_hi(2),dat_lo(3):dat_hi(3),ncomp) :: dat
    integer :: level, grid_no

    integer :: i, j, k
    integer :: ivel, irho, irhoY, iTemp

    ! Pointer into state (do not change)
    ivel  = 1
    irho  = dim + 1
    irhoY = irho + 1
    iTemp = irhoY + nspecies

    call bl_abort("derUserDefined: UserDefined derived variable requested derUserDefined has not been overwritten. &
                   Copy Source/Src_nd/problem_derive_nd.F90 in your run folder and whip up your function !")

!  Example: assuming you want to use UserDefined to tag cells within the flame front
!   1) Compute a temperature based progress variable (T_in = 300 K, T_burn = 2000 K): 
!   2) Set your UserDefined to 1 if 0.01<progvar<0.99 and 0 otherwise
!     do k=lo(3),hi(3)
!       do j=lo(2),hi(2)
!         do i=lo(1),hi(1)
!            progvar = ( dat(i,j,k,iTemp) - T_in ) / ( T_burn - T_in )
!            floormin = merge(progvar,0.0d0,progvar>0.01d0)
!            der(i,j,k,1) = merge(floormin,0.0d0,floormin<0.99d0)
!         enddo
!       enddo
!     enddo
!   3) Use UserDefined to refine the flame in your input file
!     amr.refinement_indicators = flamefront
!     amr.flamefront.max_level = 2 
!     amr.flamefront.value_greater = 0.005
!     amr.flamefront.field_name = UserDefined

  end subroutine derUserDefined

end module problem_derive_module
