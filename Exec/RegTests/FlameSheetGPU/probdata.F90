#include <AMReX_REAL.H>
#include "mechanism.h"


module probdata_module

  implicit none

  ! from probdata.H
    REAL_T :: standoff
    REAL_T :: pertmag
    
    ! from bc.H
    
    logical :: bcinit
    
    REAL_T :: u_bc, v_bc, w_bc
    REAL_T :: Y_bc(0:NUM_SPECIES-1), T_bc(1), h_bc(1), rho_bc(1)
    
    integer, parameter :: flame_dir = 2
  
contains

!subroutines here

end module probdata_module
