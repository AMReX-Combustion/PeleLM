#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>

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
                 amrex::Array4<amrex::Real> const& pressure,
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
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert;
    amrex::Real y1, y2;

    pert = 0.0;
    if (ProbParm::pertmag > 0.0) 
    {

#if ( AMREX_SPACEDIM == 2 ) 
       pert = ProbParm::pertmag * 
             (1.0 * std::sin(2 * Pi * 4 * x / L_x) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.4598) / L_x) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.712435) / L_x) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.33) / L_x) +
              0.982 * std::sin(2 * Pi * 5 * (x - 1.4234) / L_x));


    }

    y1 = (y - ProbParm::standoff - 0.5*dx[1] + pert)*100;
    y2 = (y - ProbParm::standoff + 0.5*dx[1] + pert)*100;
#elif ( AMREX_SPACEDIM == 3 ) 
       pert = ProbParm::pertmag *
             (1.0 * std::sin(2 * Pi * 4 * x / L_x) *
                std::sin(2 * Pi * 5 * y / L_y) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.4598) / L_x) *
                std::sin(2 * Pi * 4 * (y - 0.53765) / L_y) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.712435) / L_x) *
                std::sin(2 * Pi * 3 * (y - 2.137) / L_y) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.33) / L_x) *
                std::sin(2 * Pi * 6 * (y - 1.8) / L_y) +
              0.982 * std::sin(2 * Pi * 5 * (x - 1.4234) / L_x));



    }

    y1 = (z - ProbParm::standoff - 0.5*dx[2] + pert)*100;
    y2 = (z - ProbParm::standoff + 0.5*dx[2] + pert)*100;
#endif

    pmf(y1, y2, pmf_vals);

    state(i,j,k,DEF_Temp) = pmf_vals[0];;

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    EOS::X2Y(molefrac, massfrac);


    state(i,j,k,Xvel) = 0;
#if ( AMREX_SPACEDIM == 2 ) 
    state(i,j,k,Yvel) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Yvel) = 0;
    state(i,j,k,Zvel) = pmf_vals[1]*1e-2;
#endif

    pressure(i,j,k,0) = ProbParm::P_mean;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}


AMREX_GPU_DEVICE
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
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if (sgn == 1) {
    pmf(prob_hi[idir], prob_hi[idir], pmf_vals);

    s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    s_ext[Yvel] = pmf_vals[1]*1e-2;
#elif (AMREX_SPACEDIM == 3)
    s_ext[Yvel] = 0.0;
    s_ext[Zvel] = pmf_vals[1]*1e-2;
#endif

    s_ext[DEF_Temp] = pmf_vals[0];

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    EOS::X2Y(molefrac, massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;


    EOS::TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
    s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
    }
  }

}

#endif
