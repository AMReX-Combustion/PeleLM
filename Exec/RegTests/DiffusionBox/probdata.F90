#include <AMReX_REAL.H>


module probdata_module

  use network, only : nspecies

  implicit none

  ! from probdata.H
    REAL_T :: standoff
    REAL_T :: pertmag
    
    ! from bc.H
    
    logical :: bcinit
    
    REAL_T :: u_bc, v_bc, w_bc
    REAL_T :: Y_bc(0:nspecies-1), T_bc(1), h_bc(1), rho_bc(1)
    
    integer, parameter :: flame_dir = 2
  
contains

!subroutines here

end module probdata_module
