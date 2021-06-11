
#include "IndexDefines.H"
#include "PeleLM.H"
#include <AMReX_Extrapolater.H>
#include "SprayParticles.H"

using namespace amrex;

#ifdef AMREX_PARTICLES

namespace {
bool virtual_particles_set = false;
//
// Containers for the real "active" Particles
//
SprayParticleContainer* SprayPC = 0;
//
// Container for temporary, virtual Particles
//
SprayParticleContainer* VirtPC = 0;
//
// Container for temporary, ghost Particles
//
SprayParticleContainer* GhostPC = 0;

SprayData sprayData;
amrex::Real sprayRefT = 300.;
amrex::Real parcelSize = 1.;
amrex::Real spraySigma = -1.; // Surface tension
amrex::Real wall_temp = -1.;
// Indices for spray source MultiFab
int sprayMomSrcIndx = Xvel;
int sprayRhoSrcIndx = Density;
int spraySpecSrcIndx = DEF_first_spec;
int sprayEngSrcIndx = DEF_first_spec + SPRAY_FUEL_NUM;
SprayComps scomps;
bool splash_model = true;

void
RemoveParticlesOnExit()
{
  delete SprayPC;
  SprayPC = 0;
  delete GhostPC;
  GhostPC = 0;
  delete VirtPC;
  VirtPC = 0;
}
std::string particle_init_file;
int particle_init_function = 1;
} // namespace

int PeleLM::do_spray_particles = 1;
int PeleLM::particle_verbose = 0;
Real PeleLM::particle_cfl = 0.5;
// momentum + density + fuel species + enthalpy
int PeleLM::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

int PeleLM::write_spray_ascii_files = 0;
int PeleLM::plot_spray_src = 0;
int PeleLM::particle_mass_tran = 1;
int PeleLM::particle_mom_tran = 1;
Vector<std::string> PeleLM::sprayFuelNames;

void
getPSatCoef(Real* psat_coef,
            ParmParse& ppp,
            std::string fuel_name,
            const int spf)
{
  std::string psat_read = fuel_name + "_psat";
  std::vector<Real> inp_coef(4);
  ppp.getarr(psat_read.c_str(), inp_coef);
  for (int i = 0; i < 4; ++i) {
    psat_coef[4 * spf + i] = inp_coef[i];
  }
}

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
PeleLM::particleEstTimeStep(Real& est_dt)
{
  if (do_spray_particles != 1 || theSprayPC() == 0) return;
  BL_PROFILE("PeleLM::particleEstTimeStep()");
  Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

  if (est_dt_particle > 0)
    est_dt = amrex::min(est_dt, est_dt_particle);

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
PeleLM::readParticleParams()
{
  ParmParse pp("peleLM");

  pp.query("do_spray_particles", do_spray_particles);
  if (do_spray_particles != 1) {
    do_spray_particles = 0;
    return;
  }
  amrex::ParmParse ppp("particles");
  //
  // Control the verbosity of the Particle class
  ppp.query("v", particle_verbose);

  ppp.get("mass_transfer", particle_mass_tran);
  ppp.get("mom_transfer", particle_mom_tran);
  ppp.query("cfl", particle_cfl);
  if (particle_cfl > 0.5)
    amrex::Abort("particles.cfl must be <= 0.5");
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = ppp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM)
    amrex::Abort("Number of fuel species in input file must match SPRAY_FUEL_NUM");

  if (particle_mass_tran) {
    sprayFuelNames.assign(nfuel, "");
    std::vector<std::string> fuel_names;
    std::vector<Real> crit_T;
    std::vector<Real> boil_T;
    std::vector<Real> spraycp;
    std::vector<Real> latent;
    std::vector<Real> sprayrho;
    std::vector<Real> mu(nfuel, 0.);
    std::vector<Real> lambda(nfuel, 0.);
    ppp.getarr("fuel_species", fuel_names);
    ppp.getarr("fuel_crit_temp", crit_T);
    ppp.getarr("fuel_boil_temp", boil_T);
    ppp.getarr("fuel_cp", spraycp);
    ppp.getarr("fuel_latent", latent);
    ppp.getarr("fuel_rho", sprayrho);
    ppp.queryarr("fuel_mu", mu);
    ppp.queryarr("fuel_lambda", lambda);
    for (int i = 0; i < nfuel; ++i) {
      sprayFuelNames[i] = fuel_names[i];
      sprayData.critT[i] = crit_T[i];
      sprayData.boilT[i] = boil_T[i];
      sprayData.cp[i] = spraycp[i];
      sprayData.latent[i] = latent[i];
      sprayData.ref_latent[i] = latent[i];
      sprayData.rho[i] = sprayrho[i];
      sprayData.mu[i] = mu[i];
      sprayData.lambda[i] = lambda[i];
      getPSatCoef(sprayData.psat_coef.data(), ppp, fuel_names[i], i);
    }
  }

  //
  // Set the number of particles per parcel
  //
  ppp.query("parcel_size", parcelSize);
  ppp.query("use_splash_model", splash_model);
  if (splash_model) {
    if (!ppp.contains("fuel_sigma") || !ppp.contains("wall_temp")) {
      Print() << "fuel_sigma and wall_temp must be set for splash model. "
              << "Set use_splash_model = false to turn off splash model" << std::endl;
      Abort();
    }
    // Set the fuel surface tension and contact angle
    ppp.get("fuel_sigma", spraySigma);
    // TODO: Have this retrieved from proper boundary data
    ppp.get("wall_temp", wall_temp);
  }

  // Must use same reference temperature for all fuels
  // TODO: This means the reference temperature must be the same for all fuel
  // species
  ppp.get("fuel_ref_temp", sprayRefT);
  //
  // Set if spray ascii files should be written
  //
  ppp.query("write_spray_ascii_files", write_spray_ascii_files);
  //
  // Set if gas phase spray source term should be written
  //
  ppp.query("plot_src", plot_spray_src);
  //
  // Used in initData() on startup to read in a file of particles.
  //
  ppp.query("init_file", particle_init_file);
  //
  // Used in initData() on startup to set the particle field using the
  // SprayParticlesInitInsert.cpp problem specific function
  //
  ppp.query("init_function", particle_init_function);

  sprayData.num_ppp = parcelSize;
  sprayData.ref_T = sprayRefT;
  sprayData.sigma = spraySigma;

  if (verbose && ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Spray fuel species " << sprayFuelNames[0];
    for (int i = 1; i < SPRAY_FUEL_NUM; ++i)
      amrex::Print() << ", " << sprayFuelNames[i];
    amrex::Print() << std::endl;
    amrex::Print() << "Number of particles per parcel " << parcelSize << std::endl;
  }
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();
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

//Define gas phase state MF for computing spray source terms
void
PeleLM::defineSprayStateMF()
{
  int nGrowS = 4;
  if (level > 0)
    {
      int cRefRatio = parent->MaxRefRatio(level - 1);
      if (cRefRatio > 4)
        amrex::Abort("Spray particles not supported for ref_ratio > 4");
      else if (cRefRatio > 2)
        nGrowS += 3;
    }
  Sborder.define(grids, dmap, NUM_STATE, nGrowS, amrex::MFInfo(), Factory());
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
      amrex::Print() << "Fuel " << sprayFuelNames[i] << " not found in species list"
                     << std::endl;
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
  scomps.mass_tran = PeleLM::particle_mass_tran;
  scomps.mom_tran = PeleLM::particle_mom_tran;
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
  if (PeleLM::theSprayPC() != 0 && !virtual_particles_set) {
    if (level < parent->finestLevel()) {
#ifdef USE_SPRAY_SOA
      SprayParticleContainer::ParticleTileType virts;
#else
      SprayParticleContainer::AoS virts;
#endif
      ((PeleLM*)&parent->getLevel(level + 1))->setupVirtualParticles();
      PeleLM::theVirtPC()->CreateVirtualParticles(level + 1, virts);
      PeleLM::theVirtPC()->AddParticlesAtLevel(virts, level);

      PeleLM::theSprayPC()->CreateVirtualParticles(level + 1, virts);
      PeleLM::theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleLM::removeVirtualParticles()
{
  BL_PROFILE("PeleLM::removeVirtualParticles()");
  if (VirtPC != 0)
    VirtPC->RemoveParticlesAtLevel(level);
  virtual_particles_set = false;
}

void
PeleLM::setupGhostParticles(int ngrow)
{
  BL_PROFILE("PeleLM::setupGhostParticles()");
  AMREX_ASSERT(level < parent->finestLevel());
  if (PeleLM::theSprayPC() != 0) {
#ifdef USE_SPRAY_SOA
    SprayParticleContainer::ParticleTileType ghosts;
#else
    SprayParticleContainer::AoS ghosts;
#endif
    PeleLM::theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    PeleLM::theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleLM::removeGhostParticles()
{
  BL_PROFILE("PeleLM::removeGhostParticles()");
  if (GhostPC != 0)
    GhostPC->RemoveParticlesAtLevel(level);
}

/**
 * Create new particle data
 **/
void
PeleLM::createParticleData()
{
  SprayPC = new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, parcelSize, wall_temp);
  theSprayPC()->SetVerbose(particle_verbose);
  VirtPC = new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, parcelSize, wall_temp);
  GhostPC = new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, parcelSize, wall_temp);
}

/**
 * Initialize the particles on the grid at level 0
 **/
void
PeleLM::initParticles()
{
  BL_PROFILE("PeleLM::initParticles()");

  if (level > 0)
    return;

  //
  // Make sure to call RemoveParticlesOnExit() on exit.
  //
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (do_spray_particles) {
    if (theSprayPC() == 0) {
      createParticleData();
    }

    if (!particle_init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(particle_init_file, NSR_SPR + NAR_SPR);
    } else if (particle_init_function > 0) {
      ProbParm const* lprobparm = prob_parm.get();
      theSprayPC()->InitSprayParticles(*lprobparm);
    } else {
      Abort("Must initialize spray particles with particles.init_function or particles.init_file");
    }
    if (particle_verbose > 1)
      amrex::Print() << "Total number of initial particles " <<
        theSprayPC()->TotalNumberOfParticles(false, false) << std::endl;
  }
}

void
PeleLM::particlePostRestart(const std::string& restart_file, bool is_checkpoint)
{
  if (level > 0)
    return;

  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == 0);
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
PeleLM::particleMKD (const Real       time,
                     const Real       dt,
                     const int        ghost_width,
                     const int        spray_n_grow,
                     const int        tmp_src_width,
                     const int        where_width,
                     amrex::MultiFab& tmp_spray_source)
{
  //
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  //
  int finest_level = parent->finestLevel();

  //
  // Check if I need to insert new particles
  //
  int nstep = parent->levelSteps(0);

  //
  // Setup the virtual particles that represent particles on finer levels
  //
  if (level < finest_level)
    setupVirtualParticles();

  //
  // Make a copy of the particles on this level into ghost particles
  // for the finer level
  //
  if (level < finest_level)
    setupGhostParticles(ghost_width);

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (particle_verbose)
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";

  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source,
    level, dt, time,
    false, // not virtual particles
    false, // not ghost particles
    spray_n_grow,
    tmp_src_width,
    true, // Move the particles
    where_width);
  // Only need the coarsest virtual particles here.
  if (level < finest_level)
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source,
      level, dt, time, true, false,
      spray_n_grow, tmp_src_width, true, where_width);

  // Miiiight need all Ghosts
  if (theGhostPC() != 0 && level != 0)
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source,
      level, dt, time, false, true,
      spray_n_grow, tmp_src_width, true, where_width);
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, spraydot);
  spraydot.FillBoundary(geom.periodicity());
  Extrapolater::FirstOrderExtrap(spraydot, geom, 0, spraydot.nComp());
}

void
PeleLM::particleMK (const Real       time,
                    const Real       dt,
                    const int        spray_n_grow,
                    const int        tmp_src_width,
                    const int        amr_iteration,
                    const int        amr_ncycle,
                    amrex::MultiFab& tmp_spray_source)
{
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source,
    level, dt, time, false, false,
    spray_n_grow, tmp_src_width);

  if (level < parent->finestLevel() && theVirtPC() != 0)
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source,
      level, dt, time, true, false,
      spray_n_grow, tmp_src_width);
  // Ghost particles need to be kicked except during the final iteration.
  if (theGhostPC() != 0 && level != 0)
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source,
      level, dt, time, false, true,
      spray_n_grow, tmp_src_width);
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, spraydot);
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
  if (do_spray_particles != 1) return;
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
