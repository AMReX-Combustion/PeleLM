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

    double precision :: progvar, floormin

    ! Pointer into state (do not change)
    ivel  = 1
    irho  = dim + 1
    irhoY = irho + 1
    iTemp = irhoY + nspecies

    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           progvar = ( dat(i,j,k,iTemp) - 300.0d0 ) / ( 1560.0d0 - 300.0d0 )
           floormin = merge(progvar,0.0d0,progvar>0.01d0)
           der(i,j,k,1) = merge(floormin,0.0d0,floormin<0.99d0)
        enddo
      enddo
    enddo

  end subroutine derUserDefined

end module problem_derive_module
