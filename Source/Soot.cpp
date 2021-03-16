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
  // TODO: Determine the optimal way to limit time step
  // Currently, we rely on subcycling which seems to work better
  return;
}
