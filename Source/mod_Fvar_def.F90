#include <AMReX_REAL.H>

module mod_Fvar_def

  implicit none
  
  
  ! From visc.H
  
  logical :: use_constant_mu, use_constant_lambda, use_constant_rhoD
  logical :: LeEQ1

  REAL_T  :: constant_mu_val, constant_lambda_val, constant_rhoD_val
  REAL_T  :: Pr, Sc
  REAL_T  :: thickFac
  
  ! From htdata.H
  
  REAL_T  :: pamb, dpdt_factor
  integer :: closed_chamber
  
  integer :: Density, Temp, RhoH, Trac, FirstSpec, LastSpec
  
  ! From timedata.H
  
  integer :: iteration
  REAL_T  :: time
  
end module mod_Fvar_def
