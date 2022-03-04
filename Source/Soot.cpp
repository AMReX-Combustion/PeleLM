
#ifdef SOOT_MODEL
#include "SootModel.H"
#include "SootModel_derive.H"
#include "PeleLM.H"
#include "IndexDefines.H"

void
PeleLM::setSootIndx()
{
  SootComps sootComps;
  // Indices for states we are retrieving, so everything except velocity
  sootComps.qRhoIndx = Density - AMREX_SPACEDIM;
  sootComps.qSpecIndx = first_spec - AMREX_SPACEDIM;
  sootComps.qTempIndx = Temp - AMREX_SPACEDIM;
  sootComps.qSootIndx = first_passive - AMREX_SPACEDIM;
  sootComps.rhoIndx = Density - AMREX_SPACEDIM;
  sootComps.specIndx = first_spec - AMREX_SPACEDIM;
  sootComps.engIndx = RhoH - AMREX_SPACEDIM;
  sootComps.sootIndx = first_passive - AMREX_SPACEDIM;
  soot_model->setIndices(sootComps);
}

void
PeleLM::restartFromNoSoot()
{
  get_old_data(sootsrc_Type).setVal(0.);
  get_new_data(sootsrc_Type).setVal(0.);
  Real moments[NUM_SOOT_MOMENTS+1];
  SootData* const sd = PeleLM::soot_model->getSootData();
  sd->initialSmallMomVals(moments);
  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& S_old = get_old_data(State_Type);
  for (int m = 0; m < NUM_SOOT_MOMENTS+1; ++m) {
    S_new.setVal(moments[m], DEF_first_soot + m, 1, S_new.nGrow());
    S_old.setVal(moments[m], DEF_first_soot + m, 1, S_old.nGrow());
  }
}

void
PeleLM::addSootDerivePlotVars(DeriveList& derive_lst,
                              const DescriptorList& desc_lst)
{
  // Add in soot variables
  Vector<std::string> sootNames = {"rho_soot", "sum_rho_soot"};
  derive_lst.add(
    "soot_vars", IndexType::TheCellType(), sootNames.size(), sootNames,
    soot_genvars, amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("soot_vars", desc_lst, State_Type, Density, 1);
  derive_lst.addComponent(
    "soot_vars", desc_lst, State_Type, DEF_first_soot, NUM_SOOT_MOMENTS + 1);

  // Variables associated with the second mode (large particles)
  Vector<std::string> large_part_names = {"NL", "soot_V_L", "soot_S_L"};
  derive_lst.add(
    "soot_large_particles", IndexType::TheCellType(), large_part_names.size(),
    large_part_names, soot_largeparticledata, amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent(
    "soot_large_particles", desc_lst, State_Type, DEF_first_soot,
    NUM_SOOT_MOMENTS + 1);
}

void
PeleLM::computeSootSrc(Real time, Real dt)
{
  // Get soot source data
  MultiFab& soot_mf = get_new_data(sootsrc_Type);
  soot_mf.setVal(0.);
  // Get viscosity
  const TimeLevel whichTime = which_time(State_Type, time);
  AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
  auto vel_visc_cc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
  MultiFab& mf = (whichTime == AmrOldTime) ? get_old_data(State_Type): get_new_data(State_Type);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  bool pres_term = false;
  for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& gbx = mfi.growntilebox();
    auto const& q_arr = mf.array(mfi, Density); // Only need the scalars
    auto const& mu_arr = vel_visc_cc->const_array(mfi);
    auto const& soot_arr = soot_mf.array(mfi);
    soot_model->computeSootSourceTerm(gbx, q_arr, mu_arr, soot_arr, time, dt, pres_term);
  }
}

void
PeleLM::clipSootMoments()
{
  MultiFab& mf = get_new_data(State_Type);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box gbx = mfi.tilebox();
    auto const& q_arr = mf.array(mfi, DEF_first_soot);
    SootData* sd = soot_model->getSootData_d();
    amrex::ParallelFor(
      gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      GpuArray<Real, NUM_SOOT_MOMENTS + 1> moments;
      for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
        moments[mom] = q_arr(i, j, k, mom);
      }
      sd->momConvClipConv(moments.data());
      for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
        q_arr(i, j, k, mom) = moments[mom];
      }
    });
  }
}

void
PeleLM::estSootTimeStep(Real& est_dt)
{
  MultiFab& S_new = get_new_data(State_Type);
  Real local_dt = 1.E20;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& gbx = mfi.tilebox();
    auto const& q_arr = S_new.const_array(mfi, Density);
    Real sootdt = soot_model->estSootDt(gbx, q_arr);
    local_dt = amrex::min(local_dt, sootdt);
  }
  ParallelDescriptor::ReduceRealMin(local_dt);
  est_dt = amrex::min(est_dt, local_dt);
}
#endif
