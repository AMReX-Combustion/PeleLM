
#include "IndexDefines.H"
#include "PeleLM.H"
#include <AMReX_Extrapolater.H>

using namespace amrex;

#ifdef PELELM_USE_SPRAY
#include "SprayParticles.H"

namespace {
bool virtual_particles_set = false;
//
// Containers for the real "active" Particles
//
SprayParticleContainer* SprayPC = nullptr;
//
// Container for temporary, virtual Particles
//
SprayParticleContainer* VirtPC = nullptr;
//
// Container for temporary, ghost Particles
//
SprayParticleContainer* GhostPC = nullptr;

SprayData sprayData;
// Indices for spray source MultiFab
int sprayMomSrcIndx = Xvel;
int sprayRhoSrcIndx = Density;
int spraySpecSrcIndx = DEF_first_spec;
int sprayEngSrcIndx = DEF_first_spec + SPRAY_FUEL_NUM;
SprayComps scomps;

void
RemoveParticlesOnExit()
{
  delete SprayPC;
  SprayPC = nullptr;
  delete GhostPC;
  GhostPC = nullptr;
  delete VirtPC;
  VirtPC = nullptr;
}
std::string init_file;
int init_function = 1;
int particle_verbose = 0;
Real particle_cfl = 0.5;
Real wall_temp = 300.;
int mass_trans = 1;
int mom_trans = 1;
} // namespace

int PeleLM::do_spray_particles = 1;
// momentum + density + fuel species + enthalpy
int PeleLM::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

int PeleLM::write_spray_ascii_files = 0;
int PeleLM::plot_spray_src = 0;
std::string PeleLM::sprayFuelNames[SPRAY_FUEL_NUM];

SprayParticleContainer*
PeleLM::theSprayPC()
{
  return SprayPC;
}

SprayParticleContainer*
PeleLM::theVirtPC()
{
  return VirtPC;
}

SprayParticleContainer*
PeleLM::theGhostPC()
{
  return GhostPC;
}

void
PeleLM::resetSprayMF(amrex::Real time)
{
  int finest_level = parent->finestLevel();
  int ncycle = parent->nCycle(level);
  int spray_state_ghosts =
    SprayParticleContainer::getStateGhostCells(level, finest_level, ncycle);
  tmp_spray_source.setVal(0.);
  FillPatch(*this, Sborder, spray_state_ghosts, time, State_Type, 0, NUM_STATE);
}

void
PeleLM::particleEstTimeStep(Real& est_dt)
{
  if (do_spray_particles != 1 || theSprayPC() == nullptr) {
    return;
  }
  BL_PROFILE("PeleLM::particleEstTimeStep()");
  Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

  if (est_dt_particle > 0) {
    est_dt = amrex::min(est_dt, est_dt_particle);
  }

  if (verbose && ParallelDescriptor::IOProcessor()) {
    if (est_dt_particle > 0) {
      amrex::Print() << "...estdt from particles at level " << level << ": "
                     << est_dt_particle << '\n';
    } else {
      amrex::Print() << "...there are no particles at level " << level << '\n';
    }
  }
}

void
PeleLM::readSprayParams()
{
  ParmParse pp("peleLM");

  pp.query("do_spray_particles", do_spray_particles);
  if (do_spray_particles != 1) {
    do_spray_particles = 0;
    return;
  }
  SprayParticleContainer::readSprayParams(
    particle_verbose, particle_cfl, wall_temp, mass_trans, mom_trans,
    write_spray_ascii_files, plot_spray_src, init_function, init_file,
    sprayData, sprayFuelNames);
}

std::string
PeleLM::spraySrcName(const int i)
{
  if (i >= spraySpecSrcIndx && i <= spraySpecSrcIndx + SPRAY_FUEL_NUM - 1) {
    const int sp = i - spraySpecSrcIndx;
    return "I_R_spray_" + sprayFuelNames[sp];
  } else if (i <= AMREX_SPACEDIM) {
    return "I_R_spray_" + desc_lst[State_Type].name(i);
  } else if (i == sprayEngSrcIndx) {
    return "I_R_spray_" + desc_lst[State_Type].name(DEF_RhoH);
  } else {
    amrex::Abort("Should not be here");
  }
  return "";
}

// Define gas phase state MF for computing spray source terms
void
PeleLM::defineSprayStateMF()
{
  int finest_level = parent->finestLevel();
  int ncycle = parent->nCycle(level);
  int spray_state_ghosts =
    SprayParticleContainer::getStateGhostCells(level, finest_level, ncycle);
  int spray_source_ghosts =
    SprayParticleContainer::getSourceGhostCells(level, finest_level, ncycle);
  Sborder.define(
    grids, dmap, NUM_STATE, spray_state_ghosts, amrex::MFInfo(), Factory());
  tmp_spray_source.define(
    grids, dmap, num_spray_src, spray_source_ghosts, amrex::MFInfo(),
    Factory());
}

void
PeleLM::defineParticles()
{
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
#if NUM_SPECIES > 1
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == sprayFuelNames[i]) {
        sprayData.indx[i] = ns;
      }
    }
    if (sprayData.indx[i] < 0) {
      amrex::Print() << "Fuel " << sprayFuelNames[i]
                     << " not found in species list" << std::endl;
      amrex::Abort();
    }
  }
#else
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    sprayData.indx[ns] = 0;
  }
#endif
  amrex::Vector<Real> fuelEnth(NUM_SPECIES);
  auto eos = pele::physics::PhysicsType::eos();
  eos.T2Hi(sprayData.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = sprayData.indx[ns];
    sprayData.latent[ns] -= fuelEnth[fspec] * 1.E-4;
  }
  scomps.mass_tran = mass_trans;
  scomps.mom_tran = mom_trans;
  // Component indices for conservative variables
  scomps.rhoIndx = Density;
  scomps.momIndx = Xvel;
  scomps.engIndx = DEF_RhoH;
  scomps.utempIndx = DEF_Temp;
  scomps.specIndx = DEF_first_spec;
  // Component indices for spray source MultiFab
  scomps.rhoSrcIndx = sprayRhoSrcIndx;
  scomps.momSrcIndx = sprayMomSrcIndx;
  scomps.specSrcIndx = spraySpecSrcIndx;
  scomps.engSrcIndx = sprayEngSrcIndx;
}

void
PeleLM::setupVirtualParticles()
{
  BL_PROFILE("PeleLM::setupVirtualParticles()");
  if (theSprayPC() != nullptr && !virtual_particles_set) {
    if (level < parent->finestLevel()) {
      SprayParticleContainer::AoS virts;
      ((PeleLM*)&parent->getLevel(level + 1))->setupVirtualParticles();
      theVirtPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);

      theSprayPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleLM::removeVirtualParticles()
{
  BL_PROFILE("PeleLM::removeVirtualParticles()");
  if (VirtPC != nullptr) {
    VirtPC->RemoveParticlesAtLevel(level);
  }
  virtual_particles_set = false;
}

void
PeleLM::setupGhostParticles(int ngrow)
{
  BL_PROFILE("PeleLM::setupGhostParticles()");
  AMREX_ASSERT(level < parent->finestLevel());
  if (theSprayPC() != nullptr) {
    SprayParticleContainer::AoS ghosts;
    theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleLM::removeGhostParticles()
{
  BL_PROFILE("PeleLM::removeGhostParticles()");
  if (GhostPC != nullptr) {
    GhostPC->RemoveParticlesAtLevel(level);
  }
}

/**
 * Create new particle data
 **/
void
PeleLM::createParticleData()
{
  SprayPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
  theSprayPC()->SetVerbose(particle_verbose);
  VirtPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
  GhostPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
}

/**
 * Initialize the particles on the grid at level 0
 **/
void
PeleLM::initParticles()
{
  BL_PROFILE("PeleLM::initParticles()");

  if (level > 0) {
    return;
  }

  //
  // Make sure to call RemoveParticlesOnExit() on exit.
  //
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (do_spray_particles) {
    if (theSprayPC() == nullptr) {
      createParticleData();
    }

    if (!init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(init_file, NSR_SPR + NAR_SPR);
    } else if (init_function > 0) {
      ProbParm const* lprobparm = prob_parm;
      theSprayPC()->InitSprayParticles(*lprobparm);
    } else {
      Abort("Must initialize spray particles with particles.init_function or "
            "particles.init_file");
    }
    if (particle_verbose > 1)
      amrex::Print() << "Total number of initial particles "
                     << theSprayPC()->TotalNumberOfParticles(false, false)
                     << std::endl;
  }
}

void
PeleLM::particlePostRestart(bool is_checkpoint)
{
  if (level > 0) {
    return;
  }

  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == nullptr);
    createParticleData();

    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    amrex::ExecOnFinalize(RemoveParticlesOnExit);
    {
      theSprayPC()->Restart(
        parent->theRestartFile(), "particles", is_checkpoint);
      amrex::Gpu::Device::streamSynchronize();
    }
  }
}

void
PeleLM::particleMKD(const Real time, const Real dt)
{
  //
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  //
  int finest_level = parent->finestLevel();
  int ncycle = parent->nCycle(level);
  int spray_state_ghosts =
    SprayParticleContainer::getStateGhostCells(level, finest_level, ncycle);
  int spray_source_ghosts =
    SprayParticleContainer::getSourceGhostCells(level, finest_level, ncycle);
  //
  // Setup the virtual particles that represent particles on finer levels
  // and make a copy of the particles on this level into ghost particles
  // for the finer level
  //
  if (level < finest_level) {
    setupVirtualParticles();
    int finer_ref = parent->MaxRefRatio(level);
    int ghost_width = SprayParticleContainer::getGhostPartCells(
      level + 1, finest_level, finer_ref);
    setupGhostParticles(ghost_width);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (particle_verbose) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }
  auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();
  // Do the valid particles themselves
  bool isVirt = false;
  bool isGhost = false;
  bool doMove = true;
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
    spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);
  // Only need the coarsest virtual particles here.
  if (level < finest_level) {
    isVirt = true;
    isGhost = false;
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);
  }

  // Miiiight need all Ghosts
  if (theGhostPC() != nullptr && level != 0) {
    isVirt = false;
    isGhost = true;
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);
  }
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    spray_source_ghosts, level, tmp_spray_source, spraydot);
  spraydot.FillBoundary(geom.periodicity());
  Extrapolater::FirstOrderExtrap(spraydot, geom, 0, spraydot.nComp());
}

void
PeleLM::particleMK(const Real time, const Real dt)
{
  auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();
  int finest_level = parent->finestLevel();
  int ncycle = parent->nCycle(level);
  int spray_state_ghosts =
    SprayParticleContainer::getStateGhostCells(level, finest_level, ncycle);
  int spray_source_ghosts =
    SprayParticleContainer::getSourceGhostCells(level, finest_level, ncycle);
  bool isVirt = false;
  bool isGhost = false;
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
    spray_state_ghosts, spray_source_ghosts, ltransparm);

  if (level < parent->finestLevel() && theVirtPC() != nullptr) {
    isVirt = true;
    isGhost = false;
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }
  // Ghost particles need to be kicked except during the final iteration.
  if (theGhostPC() != nullptr && level != 0) {
    isVirt = false;
    isGhost = true;
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  theSprayPC()->transferSource(
    spray_source_ghosts, level, tmp_spray_source, spraydot);
  spraydot.FillBoundary(geom.periodicity());
  Extrapolater::FirstOrderExtrap(spraydot, geom, 0, spraydot.nComp());
}

// TODO: This has not been checked or updated, use with caution
#if 0
std::unique_ptr<MultiFab>
PeleLM::particleDerive(const std::string& name, Real time, int ngrow)
{
  BL_PROFILE("PeleLM::particleDerive()");

  if (theSprayPC() && name == "particle_count") {
    amrex::Abort("Should not be called until it is updated");
    std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
    MultiFab temp_dat(grids, dmap, 1, 0);
    temp_dat.setVal(0);
    theSprayPC()->Increment(temp_dat, level);
    MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
    return derive_dat;
  } else if (theSprayPC() && name == "total_particle_count") {
    amrex::Abort("Should not be called until it is updated");
    //
    // We want the total particle count at this level or higher.
    //
    std::unique_ptr<MultiFab> derive_dat =
      particleDerive("particle_count", time, ngrow);

    IntVect trr(AMREX_D_DECL(1, 1, 1));

    for (int lev = level + 1; lev <= parent->finestLevel(); lev++) {
      auto ba = parent->boxArray(lev);
      const auto& dm = parent->DistributionMap(lev);
      MultiFab temp_dat(ba, dm, 1, 0);

      trr *= parent->refRatio(lev - 1);

      ba.coarsen(trr);
      MultiFab ctemp_dat(ba, dm, 1, 0);

      temp_dat.setVal(0);
      ctemp_dat.setVal(0);

      theSprayPC()->Increment(temp_dat, lev);

      for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi) {
        const FArrayBox& ffab = temp_dat[mfi];
        FArrayBox& cfab = ctemp_dat[mfi];
        const Box& fbx = ffab.box();

        AMREX_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

        for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p)) {
          const Real val = ffab(p);
          if (val > 0)
            cfab(amrex::coarsen(p, trr)) += val;
        }
      }

      temp_dat.clear();

      MultiFab dat(grids, dmap, 1, 0);
      dat.setVal(0);
      dat.copy(ctemp_dat);

      MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
    }

    return derive_dat;
  } else {
    return AmrLevel::derive(name, time, ngrow);
  }
}
#endif

void
PeleLM::particle_redistribute(int lbase, bool init_part)
{
  if (do_spray_particles != 1) {
    return;
  }
  BL_PROFILE("PeleLM::particle_redistribute()");
  int flev = parent->finestLevel();
  if (theSprayPC()) {
    //
    // If we are calling with init_part = true, then we want to force the
    // redistribute without checking whether the grids have changed.
    //
    if (init_part) {
      theSprayPC()->Redistribute(lbase);
      return;
    }

    //
    // These are usually the BoxArray and DMap from the last regridding.
    //
    static Vector<BoxArray> ba;
    static Vector<DistributionMapping> dm;

    bool changed = false;

    while (parent->getAmrLevels()[flev] == nullptr) {
      flev--;
    }

    if (ba.size() != flev + 1) {
      ba.resize(flev + 1);
      dm.resize(flev + 1);
      changed = true;
    } else {
      for (int i = 0; i <= flev && !changed; i++) {
        // Check if BoxArrays have changed during regridding
        if (ba[i] != parent->boxArray(i))
          changed = true;

        if (!changed) {
          // Check DistributionMaps have changed during regridding
          if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
            changed = true;
        }
      }
    }

    if (changed) {
      //
      // We only need to call Redistribute if the BoxArrays or DistMaps have
      // changed.
      // We also only call it for particles >= lbase. This is
      // because if we called redistribute during a subcycle, there may be
      // particles not in the proper position on coarser levels.
      //
      if (verbose && ParallelDescriptor::IOProcessor())
        amrex::Print() << "Calling redistribute because grid has changed "
                       << '\n';
      if (flev == 0) {
        // Do a local redistribute
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 1);
      } else {
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 1);
      }
      //
      // Use the new BoxArray and DistMap to define ba and dm for next time.
      //
      for (int i = 0; i <= flev; i++) {
        ba[i] = parent->boxArray(i);
        dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
      }
    } else {
      if (verbose && ParallelDescriptor::IOProcessor())
        amrex::Print()
          << "NOT calling redistribute because grid has NOT changed " << '\n';
    }
  }
}

void
PeleLM::particleTimestamp(int ngrow)
{
  amrex::Abort("Need to fill in PeleLM::TimestampParticles");
}

#endif
