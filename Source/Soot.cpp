#include "PeleLM.H"

void
PeleLM::setSootIndx ()
{
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
PeleLM::computeSootSrc (Real time,
                        Real dt)
{
  const int nGrowOp = 0;
  // Get soot source data
  MultiFab& soot_mf = get_new_data(sootsrc_Type);
  soot_mf.setVal(0.);
  // Need all scalar state variables
  // Get fillpatched state
  FillPatchIterator fpi(*this,soot_mf,nGrowOp,time,State_Type,Density,DEF_NUM_SCALARS);
  MultiFab& mf= fpi.get_mf();
  // Get viscosity
  const TimeLevel whichTime = which_time(State_Type,time);
  AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
  auto vel_visc_cc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& gbx = mfi.growntilebox();
    auto const& q_arr    = mf.const_array(mfi);
    auto const& mu_arr   = vel_visc_cc->const_array(mfi);
    auto const& soot_arr = soot_mf.array(mfi);
    soot_model->addSootSourceTerm(gbx, q_arr, mu_arr, soot_arr, time, dt);
  }
  amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;
  amrex::Real soot_dt = std::numeric_limits<amrex::Real>::max();
  SootData const* sd = soot_model->getSootData();
  for (MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const int rhoIndx = sootComps.rhoIndx;
    const int qRhoIndx = sootComps.qRhoIndx;
    const int specIndx = sootComps.specIndx;
    const int qSpecIndx = sootComps.qSpecIndx;
    const int sootIndx = sootComps.sootIndx;
    const int qSootIndx = sootComps.qSootIndx;
    const amrex::Box& bx = mfi.tilebox();
    auto const& s_arr = mf.const_array(mfi);
    auto const& soot_arr = soot_mf.const_array(mfi);
    {
      BL_PROFILE("PeleLM::retrieveSootTimeStep()");
      reduce_op.eval(
        bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          amrex::Real max_rate = std::numeric_limits<amrex::Real>::min();
          amrex::Real rho_rate = soot_arr(i, j, k, rhoIndx) / s_arr(i, j, k, qRhoIndx);
          max_rate = amrex::max(std::abs(rho_rate), max_rate);
          for (int n = 0; n < NUM_SOOT_GS; ++n) {
            int scomp = sd->refIndx[n];
            int indx = specIndx + scomp;
            int qindx = qSpecIndx + scomp;
            amrex::Real cur_rate =
              soot_arr(i, j, k, indx) / (s_arr(i, j, k, qindx) + 1.E-12);
            max_rate = amrex::max(max_rate, std::abs(cur_rate));
          }
          // Limit time step based only on positive moment sources
          for (int n = 0; n < DEF_NUM_SOOT_VARS; ++n) {
            int qindx = qSootIndx + n;
            int indx = sootIndx + n;
            max_rate =
              amrex::max(max_rate, -soot_arr(i, j, k, indx) / s_arr(i, j, k, qindx));
          }
          return 1. / max_rate;
        });
      ReduceTuple hv = reduce_data.value();
      Real ldt_cpu = amrex::get<0>(hv);
      soot_dt = amrex::min(soot_dt, ldt_cpu);
    }
  }
  ParallelDescriptor::ReduceRealMin(soot_dt);
  soot_model->setTimeStep(soot_dt);
}
