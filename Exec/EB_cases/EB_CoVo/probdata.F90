#include <AMReX_REAL.H>


module probdata_module

  implicit none

  ! from probdata.H
  INTEGER :: meanFlowDir
  REAL_T  :: meanFlowMag
  REAL_T  :: T_mean, P_mean
  REAL_T  :: xvort, yvort, rvort, forcevort
  
  logical :: bcinit

  
contains

!subroutines here

end module probdata_module
