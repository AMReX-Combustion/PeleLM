#include <AMReX_REAL.H>

module mod_Fvar_def

  implicit none
  
  
  ! From visc.H
  
  logical, save :: use_constant_mu, use_constant_lambda, use_constant_rhoD
  logical, save :: LeEQ1

  REAL_T,  save :: constant_mu_val, constant_lambda_val, constant_rhoD_val
  REAL_T,  save :: Pr, Sc
  REAL_T,  save :: thickFacTR
  
  ! From htdata.H
  
  REAL_T,  save :: pamb, dpdt_factor
  integer, save :: closed_chamber
  
  integer, save :: Density, Temp, RhoH, Trac, FirstSpec, LastSpec
  
  ! From timedata.H
  
  integer, save :: iteration
  REAL_T,  save :: time
  
end module mod_Fvar_def