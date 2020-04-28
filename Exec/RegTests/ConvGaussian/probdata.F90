#include <AMReX_REAL.H>


module probdata_module

  implicit none

  ! from probdata.H
  INTEGER :: meanFlowDir
  REAL_T  :: meanFlowMag
  REAL_T  :: T_mean, P_mean
  REAL_T  :: xgauss, ygauss, rgauss, ampgauss

  CHARACTER(LEN=256) :: gauss_type
  
  logical :: bcinit
  
contains

!subroutines here

end module probdata_module
