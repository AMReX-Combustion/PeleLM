#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  amrex::Real xlo,
  amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector)
{
  amrex::Real sum = 0.0, xmid = 0.0;
  int lo_loside = 0, lo_hiside = 0;
  int hi_loside = 0, hi_hiside = 0;
  int loside = 0, hiside = 0;
  amrex::Real ylo = 0.0, yhi = 0.0, x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0,
              dydx = 0.0;

  if (ProbParm::pmf_do_average) {
    lo_loside = 0;
    lo_hiside = 0;
    hi_loside = 0;
    hi_hiside = 0;
    if (xlo < ProbParm::d_pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xhi < ProbParm::d_pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xlo > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      lo_loside = ProbParm::pmf_N - 1;
      lo_hiside = ProbParm::pmf_N - 1;
    }
    if (xhi > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      hi_loside = ProbParm::pmf_N - 1;
      hi_hiside = ProbParm::pmf_N - 1;
    }
    if (lo_loside == 0) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if ((xlo > ProbParm::d_pmf_X[i]) && (xlo < ProbParm::d_pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == 0) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if ((xhi > ProbParm::d_pmf_X[i]) && (xhi < ProbParm::d_pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
for (int j = 0; j < ProbParm::pmf_M; j++) {
      x1 = ProbParm::d_pmf_X[lo_loside];
      y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + lo_loside];
      x2 = ProbParm::d_pmf_X[lo_hiside];
      y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + lo_hiside];
      if (lo_loside == lo_hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      ylo = y1 + dydx * (xlo - x1);
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = ProbParm::d_pmf_X[hi_loside];
        y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hi_loside];
        x2 = ProbParm::d_pmf_X[hi_hiside];
        y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hi_hiside];
        if (hi_loside == hi_hiside)
          dydx = 0.0;
        else
          dydx = (y2 - y1) / (x2 - x1);
        yhi = y1 + dydx * (xhi - x1);
        sum = sum + (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k < hi_loside - 1; k++) {
          sum = sum + (ProbParm::d_pmf_X[k + 1] - ProbParm::d_pmf_X[k]) * 0.5 *
                        (ProbParm::d_pmf_Y[ProbParm::pmf_N * j + k] +
                         ProbParm::d_pmf_Y[ProbParm::pmf_N * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }
  } else {
    xmid = 0.5 * (xlo + xhi);
    loside = -1;
    hiside = -1;
    if (xmid < ProbParm::d_pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      loside = ProbParm::pmf_N - 1;
      hiside = ProbParm::pmf_N - 1;
    }
    if (loside == -1) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if (
          (xmid >= ProbParm::d_pmf_X[i]) &&
          (xmid <= ProbParm::d_pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < ProbParm::pmf_M; j++) {
      x1 = ProbParm::d_pmf_X[loside];
      y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + loside];
      x2 = ProbParm::d_pmf_X[hiside];
      y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hiside];
      if (loside == hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }
}


AMREX_GPU_DEVICE
inline
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real L_y = prob_hi[1] - prob_lo[1];

    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real u[3] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert;
    amrex::Real y1, y2;

    massfrac[N2_ID] = 0.767;
    massfrac[O2_ID] = 0.233;

    state(i,j,k,DEF_Temp) = 300.0;

    if ( x >= ProbParm::splitx ) {
       if ( (x-ProbParm::splitx) < ( std::sqrt((y+0.0009) / 0.01))/100.0) {
          amrex::Real pl1_angle = std::atan(0.01/ProbParm::slot_width);
          y1 = - ( x - ProbParm::splitx - ProbParm::slot_width ) * std::cos(pl1_angle) + y * std::sin(pl1_angle);
          y2 =   ( x - ProbParm::splitx - ProbParm::slot_width ) * std::sin(pl1_angle) + y * std::cos(pl1_angle);
          y1 += ProbParm::standoff;
          y2 += ProbParm::standoff;

          pmf(y2, y2, pmf_vals);
          state(i,j,k,DEF_Temp) = pmf_vals[0];
          for (int n = 0; n < NUM_SPECIES; n++){
            molefrac[n] = pmf_vals[3 + n];
          }
          EOS::X2Y(molefrac, massfrac);
       }
    } else {
       if ( std::abs((x-ProbParm::splitx)) < ( std::sqrt((y+0.0009) / 0.01))/100.0) {
          amrex::Real pl1_angle = std::atan(-0.01/ProbParm::slot_width);
          y1 = ( x - ProbParm::splitx + ProbParm::slot_width ) * std::cos(pl1_angle) - y * std::sin(pl1_angle);
          y2 = ( x - ProbParm::splitx + ProbParm::slot_width ) * std::sin(pl1_angle) + y * std::cos(pl1_angle);
          y1 += ProbParm::standoff;
          y2 += ProbParm::standoff;

          pmf(y2, y2, pmf_vals);
          state(i,j,k,DEF_Temp) = pmf_vals[0];
          for (int n = 0; n < NUM_SPECIES; n++){
            molefrac[n] = pmf_vals[3 + n];
          }
          EOS::X2Y(molefrac, massfrac);
       }
    }

    state(i,j,k,Xvel) = 0;
    state(i,j,k,Yvel) = 0.1;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

#ifdef PLM_USE_EFIELD
    state(i,j,k,DEF_nE) = std::max(0.0,massfrac[H_ID]) * 1.0e22;
    state(i,j,k,DEF_PhiV) = 0.0;
#endif

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};

  massfrac_air[N2_ID] = 0.767;
  massfrac_air[O2_ID] = 0.233;
  massfrac_jet[N2_ID] = 0.725;
  massfrac_jet[O2_ID] = 0.22;
  massfrac_jet[CH4_ID] = 0.055;

  for (int n = 0; n < DEF_NUM_STATE; n++){
     s_ext[n] = 0.0;
  }

// Get zone
  int zone = -1;
  if ( std::abs(x[0]-ProbParm::splitx) < ProbParm::slot_width ) {
     zone = 3; // jet
  } else {
     zone = 1; // coflow
  }

  if (idir == 1 && sgn == 1) {
    s_ext[Xvel] = 0.0;
    amrex::Real rad = 0.0;
    if ( zone == 3 ) {
       rad = std::abs(x[0]-ProbParm::splitx);
       s_ext[Yvel] = ProbParm::Vin * ( 1.0 - (rad/ProbParm::slot_width) * (rad/ProbParm::slot_width) );
    } else if ( zone == 1 ) {
       rad = std::abs(x[0]-ProbParm::splitx);
       s_ext[Yvel] = ProbParm::Vcoflow * ( 1.0 + (rad)/0.025);
       if ( rad <= ProbParm::slot_width*1.4 ) s_ext[Yvel] = 0.0;
    }

    s_ext[DEF_Temp] = 300.0;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = ProbParm::P_mean * 10.0;

    if ( zone == 1 ) {
       EOS::PYT2R(P_cgs, massfrac_air, s_ext[DEF_Temp], rho_cgs);
       s_ext[Density] = rho_cgs * 1.0e3;


       EOS::TY2H(s_ext[DEF_Temp], massfrac_air, RhoH_temp);
       s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

       for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[DEF_first_spec+n] = massfrac_air[n] * s_ext[Density];
       }
    } else if ( zone == 3 ) {
       EOS::PYT2R(P_cgs, massfrac_jet, s_ext[DEF_Temp], rho_cgs);
       s_ext[Density] = rho_cgs * 1.0e3;


       EOS::TY2H(s_ext[DEF_Temp], massfrac_jet, RhoH_temp);
       s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

       for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[DEF_first_spec+n] = massfrac_jet[n] * s_ext[Density];
       }
    }
  }

#ifdef PLM_USE_EFIELD
  if (idir == 0) {  // x-dir
    if ( sgn == 1 ) { // lo
       s_ext[DEF_nE] = 0.0;
       s_ext[DEF_PhiV] = ProbParm::phiV_lox;
    } else {          // hi
       s_ext[DEF_nE] = 0.0;
       s_ext[DEF_PhiV] = ProbParm::phiV_hix;
    }
  }
#endif
}

AMREX_GPU_DEVICE
inline
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

   amrex::Real rad = std::abs(x-ProbParm::splitx);
   if (  rad >= ProbParm::slot_width &&
         rad <= ProbParm::slot_width*1.4 &&
         j <= dBox.smallEnd(1)) {
      for (int n = 0; n < nComp; n++ ) {
         beta(i,j,k,n) = 0.0;
      }
   }
}
#endif
