#include <AMReX_REAL.H>


module probdata_module

  implicit none

  ! from probdata.H
  REAL_T  :: T_mean, P_mean
  
  logical :: bcinit

  character(len=255), save :: iname
  logical, save            :: binfmt
  logical, save            :: restart
  double precision, save   :: lambda0, reynolds_lambda0
  integer, save            :: inres
  double precision, save   :: uin_norm
  double precision, save   :: Linput
  double precision, save, dimension(:,:,:), allocatable :: xinput, yinput, zinput, &
       uinput, vinput, winput
  double precision, save, dimension(:), allocatable :: xarray, xdiff
  double precision, save   :: k0, urms0, tau

 
contains

!subroutines here

end module probdata_module
