#include <AMReX_REAL.H>

module mod_Fvar_def

  use network,   only: nspecies

  implicit none


  ! From visc.H

  logical :: LeEQ1

  REAL_T  :: Pr, Sc
  REAL_T  :: thickFac = 1.d0

  ! From htdata.H

  REAL_T  :: pamb = 101325.d0

  integer :: Density, Temp, RhoH, Trac, FirstSpec, LastSpec

  ! From timedata.H

  integer :: iteration
  REAL_T  :: time

  ! From probdata.H
  integer :: bathID, fuelID, oxidID, prodID

  integer         , save :: f_flag_active_control
  ! geometry information
  double precision, save :: domnlo(3), domnhi(3)

  ! From bc.H
  integer, parameter :: MAXPNTS = 50
  REAL_T :: time_points(0:MAXPNTS),vel_points(0:MAXPNTS),cntl_points(0:MAXPNTS)

  character(50) :: ac_hist_file
  REAL_T :: tau_control, cfix, coft_old, sest, V_in, V_in_old, corr, &
          changeMax_control, tbase_control, dV_control, scale_control, &
          zbase_control, h_control, controlVelMax
  integer :: navg_pnts, flame_dir, pseudo_gravity

  ! Oxidiser and fuel stream mass fractions. Required to compute mixture
  ! fraction, see init_mixture_fraction() in Src_3d/Derive_PLM_3d.F90.
  ! Needs to be set by the user, e.g. in amrex_probinit(...).
  REAL_T, dimension(nspecies) :: Y_fu, Y_ox

end module mod_Fvar_def
