#include <AMReX_REAL.H>


module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none

  ! from probdata.H
  REAL_T :: T_mean, P_mean
  REAL_T :: xgauss, ygauss, rgauss 
  
  logical :: bcinit

  
contains

!subroutines here

end module probdata_module
