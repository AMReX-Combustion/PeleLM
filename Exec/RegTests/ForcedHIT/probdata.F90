module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none

  double precision :: T_in
   
  ! From forcedata.H
  integer, save :: spectrum_type, mode_start, nmodes
  double precision, save   :: turb_scale, force_scale, forcing_time_scale_min, forcing_time_scale_max
  double precision, save :: time_offset

  integer, save :: div_free_force
  integer, save :: use_rho_in_forcing

  integer, parameter :: blrandseed = 0
  integer, parameter :: moderate_zero_modes = 0
  double precision, parameter :: forcing_epsilon = 1.0d-4

  double precision, save, dimension(:,:,:), allocatable :: FTX, FTY, FTZ, TAT, TAP, &
                                                FPX, FPY, FPZ, FAX, FAY, FAZ, &
                                                FPXX, FPYX, FPZX, FPXY, FPYY, FPZY, FPXZ, FPYZ, FPZZ
  
contains

!subroutines here

end module probdata_module
