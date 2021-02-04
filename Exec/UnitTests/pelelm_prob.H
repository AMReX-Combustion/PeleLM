#ifndef _PROB_H_
#define _PROB_H_

//#include <AMReX_GpuQualifiers.H>
//#include <AMReX_Extension.H>
//#include <AMReX_REAL.H>

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include "pelelm_prob_parm.H"
#include "pmf.H"
#include "pmf_data.H"
#include "PeleLM_parm.H"

#include "NavierStokesBase.H"
#include "IndexDefines.H"
#include "EOS.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 PmfData const *pmf_data)
{}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/,
  ACParm const& /*lacparm*/,
  PmfData const* /*pmf_data*/)
{
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}

#if 0
struct MyProbTagStruct
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int /*i*/,
    const int /*j*/,
    const int /*k*/,
    amrex::Array4<char> const& /*tag*/,
    amrex::Array4<amrex::Real const> const& /*field*/,
    char /*tagval*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
    const amrex::Real /*time*/,
    const int /*level*/) noexcept
  {
    // could do problem specific tagging here
  }
};

using ProblemTags = MyProbTagStruct;

struct MyProbDeriveStruct
{
  static void
  add(amrex::DeriveList& /*derive_lst*/, amrex::DescriptorList& /*desc_lst*/)
  {
    // Add derives as follows and define the derive function below:
    // derive_lst.add(
    //  "varname", amrex::IndexType::TheCellType(), 1, pc_varname,
    //  the_same_box);
    // derive_lst.addComponent("varname", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_varname(
    const amrex::Box& /*bx*/,
    amrex::FArrayBox& /*derfab*/,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& /*datfab*/,
    const amrex::Geometry& /*geomdata*/,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
    // auto const dat = datfab.array();
    // auto arr = derfab.array();
    // amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    // { do something with arr
    // });
  }
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;
#endif
#endif
