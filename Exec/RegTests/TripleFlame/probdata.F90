#include <AMReX_REAL.H>

module probdata_module

  implicit none

  !   Flag checking if setupbc() has been called    
    logical :: bcinit
    
!   BC variables: filled in setupbc(), read in bcfunction()    
    REAL_T :: u_bc, v_bc, T_bc(1)

!   Input data from probin file    
    REAL_T :: midtanh         ! Position of the mixing layer (in x direction)
    REAL_T :: splitx          ! Half of the domain in x direction
    REAL_T :: widthtanh       ! Width of the mixing layer TANH
    REAL_T :: H2_enrich       ! H2 enrichment in volume (or mole fraction)
    REAL_T :: T_in            ! Inflow temperature (uniform)

!   Definition of the stoichiometric mixture fraction (useful here)   
    REAL_T :: Zst             
    
contains

!subroutines here

end module probdata_module
