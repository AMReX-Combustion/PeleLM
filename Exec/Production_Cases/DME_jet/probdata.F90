#include <AMReX_REAL.H>


module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none

  ! from probdata.H
    REAL_T :: T_in, V_co, phi_in, T_co
    REAL_T :: splitx, xfrontw
    REAL_T :: blobr, bloby, blobx, Tfrontw, blobT
    
    REAL_T :: fuel_N2_vol_percent
 
    integer, parameter :: IDX_FUELPIPE = 1
    integer, parameter :: IDX_AMBIENT = 2
    integer, parameter :: IDX_VOLUME = 3
    integer, parameter :: IDX_COFLOW = 4
    
 
    ! from bc.H
    
    logical :: bcinit
    integer, parameter :: Nzones=5
    
    REAL_T :: u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
    REAL_T :: Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
    REAL_T :: h_bc(Nzones)
    
contains

!subroutines here

end module probdata_module
