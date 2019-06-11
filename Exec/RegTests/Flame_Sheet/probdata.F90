#include <AMReX_REAL.H>

module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none

  ! from probdata.H
    REAL_T :: standoff
    REAL_T :: V_in
    
    REAL_T :: pertmag
    
    ! from bc.H
    
    logical :: bcinit
    integer, parameter :: Nzones=5
    
    REAL_T :: u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
    REAL_T :: Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
    REAL_T :: h_bc(Nzones)
    
    integer, parameter :: flame_dir = 2
  
contains

!subroutines here

end module probdata_module
