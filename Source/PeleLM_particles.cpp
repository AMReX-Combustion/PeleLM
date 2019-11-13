#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "PeleLM.H"
#include "PeleLM_F.H"
#include "Spray_F.H"
#include "SprayParticles.H"

using namespace amrex;

#ifdef AMREX_PARTICLES

///#include <AMReX_Particles.H> // mr: test including this
///#include <AMReX_ParticleInit.H> // mr: test including this

#define MAX_NUM_FUELS 1

namespace
{
    bool virtual_particles_set = false;

    std::string ascii_particle_file;
    std::string binary_particle_file;

    const std::string chk_particle_file("Spray");

    //
    // We want to call this routine when on exit to clean up particles.
    //

    //
    // Containers for the real "active" Particles
    //
    SprayParticleContainer* SprayPC = 0;
    //
    // Container for temporary, virtual Particles
    //
    SprayParticleContainer* VirtPC  = 0;
    //
    // Container for temporary, ghost Particles
    //
    SprayParticleContainer* GhostPC  = 0;

    void RemoveParticlesOnExit ()
    {
        delete SprayPC;
        delete GhostPC;
        delete VirtPC;
    }
}

int PeleLM::do_spray_particles          =  0;
int PeleLM::particle_verbose            =  0;
Real PeleLM::particle_cfl               = 0.05;

namespace {
    std::string       particle_init_file;
    int               particle_init_uniform = 0;
    int               is_mass_tran = 0;
    int               is_heat_tran = 0;
    int               is_mom_tran = 0;
    int               nfuel_species = 1;
    Vector<Real>      fuel_mass_frac(MAX_NUM_FUELS);
    Vector<Real>      fuel_density(MAX_NUM_FUELS);
    Vector<Real>      fuel_crit_temp(MAX_NUM_FUELS);
    Vector<Real>      fuel_boil_temp(MAX_NUM_FUELS);
    Vector<Real>      fuel_latent(MAX_NUM_FUELS);
    Vector<Real>      fuel_cp(MAX_NUM_FUELS);
    Vector<Real>      fuel_molwt(MAX_NUM_FUELS);
    Vector<int>       fuel_indx(MAX_NUM_FUELS);
    std::string       particle_output_file;
    std::string       timestamp_dir;
    std::vector<int>  timestamp_indices;
}

SprayParticleContainer*
PeleLM::theSprayPC()
{
    return SprayPC;
}

SprayParticleContainer*
PeleLM::theVirtPC ()
{
    return VirtPC;
}

SprayParticleContainer*
PeleLM::theGhostPC ()
{
    return GhostPC;
}

void
PeleLM::particle_est_time_step (Real& est_dt)
{
    BL_PROFILE("PeleLM::particle_est_time_step()");
    Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

    if (est_dt_particle > 0)
        est_dt = std::min(est_dt, est_dt_particle);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        if (est_dt_particle > 0)
        {
            std::cout << "...estdt from particles at level "
                      << level << ": " << est_dt_particle << '\n';
        }
        else
        {
            std::cout << "...there are no particles at level "
                      << level << '\n';
        }
    }
}

void 
PeleLM::read_particle_params ()
{
    ParmParse pp("pelelm");
    
    pp.query("do_spray_particles",do_spray_particles);

    ParmParse ppp("particles");
    //
    // Control the verbosity of the Particle class
    ppp.query("v",particle_verbose);
    //
    ppp.get("mass_transfer",is_mass_tran);
    ppp.get("heat_transfer",is_heat_tran);
    ppp.get("mom_transfer",is_mom_tran);
    //
    // Used in import_fuel_properties()
    //
    ppp.get("fuel_species", nfuel_species);
    ppp.getarr("fuel_mass_frac", fuel_mass_frac, 0, nfuel_species);
    ppp.getarr("fuel_density", fuel_density, 0, nfuel_species);
    ppp.getarr("fuel_crit_temp", fuel_crit_temp, 0, nfuel_species);
    ppp.getarr("fuel_boil_temp", fuel_boil_temp, 0, nfuel_species);
    ppp.getarr("fuel_latent", fuel_latent, 0, nfuel_species);
    ppp.getarr("fuel_cp", fuel_cp, 0, nfuel_species);
    ppp.getarr("fuel_molwt", fuel_molwt, 0, nfuel_species);
    ppp.getarr("fuel_indx", fuel_indx, 0, nfuel_species);
    // mr: test input reading
    //amrex::Print() << "Read fuel molwt " << fuel_molwt.data() << '\n';
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in initData() on startup to set a uniform particle field
    //
    ppp.query("particle_init_uniform", particle_init_uniform);
    //
    // Used in post_restart() to read in a file of particles.
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    // ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    ppp.query("particle_output_file", particle_output_file);
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(timestamp_dir, 0755))
            amrex::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
}

void
PeleLM::setup_virtual_particles()
{
    BL_PROFILE("PeleLM::setup_virtual_particles()");

    if(PeleLM::theSprayPC() != 0 && !virtual_particles_set)
    {
        SprayParticleContainer::AoS virts;
        if (level < parent->finestLevel())
        {
            ((PeleLM*) &parent->getLevel(level+1))->setup_virtual_particles();
            PeleLM::theVirtPC()->CreateVirtualParticles(level+1, virts);
            PeleLM::theVirtPC()->AddParticlesAtLevel(virts, level);

            PeleLM::theSprayPC()->CreateVirtualParticles(level+1, virts);
            PeleLM::theVirtPC()->AddParticlesAtLevel(virts, level);
        }
        virtual_particles_set = true;
    }
}

void
PeleLM::remove_virtual_particles()
{
    BL_PROFILE("PeleLM::remove_virtual_particles()");
    if (VirtPC != 0)
        VirtPC->RemoveParticlesAtLevel(level);
    virtual_particles_set = false;
}

void
PeleLM::setup_ghost_particles(int ngrow)
{
    BL_PROFILE("PeleLM::setup_ghost_particles()");
    BL_ASSERT(level < parent->finestLevel());
    if(PeleLM::theSprayPC() != 0)
    {
        SprayParticleContainer::AoS ghosts;
        PeleLM::theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
        PeleLM::theGhostPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    }
}

void
PeleLM::remove_ghost_particles()
{
    BL_PROFILE("PeleLM::remove_ghost_particles()");
    if (GhostPC != 0)
        GhostPC->RemoveParticlesAtLevel(level);
}


/**
 * Initialize the particles on the grid at level 0
 **/
void
PeleLM::init_particles ()
{
    BL_PROFILE("PeleLM::init_particles()");

    if (level > 0)
        return;

    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    amrex::ExecOnFinalize(RemoveParticlesOnExit);

    if (do_spray_particles)
    {
  	BL_ASSERT(theSprayPC() == 0);
	
	SprayPC = new SprayParticleContainer(parent,&phys_bc);
	theSprayPC()->SetVerbose(particle_verbose);

        if (parent->subCycle())
        {
            VirtPC = new SprayParticleContainer(parent,&phys_bc);
            GhostPC = new SprayParticleContainer(parent,&phys_bc);
        }
	
        // fuel properties from user input file
        import_fuel_properties(&nfuel_species, fuel_mass_frac.data(),fuel_density.data(),
                               fuel_crit_temp.data(), fuel_latent.data(), fuel_boil_temp.data(), 
                               fuel_cp.data(),fuel_molwt.data(),fuel_indx.data());

        import_control_parameters(&is_heat_tran, &is_mass_tran, &is_mom_tran);

	if (! particle_init_file.empty())
	{
	    theSprayPC()->InitFromAsciiFile(particle_init_file,NSR_SPR);
	} 
        else if (particle_init_uniform > 0) 
        {
            // Initialize uniform particle distribution
            //  {vel(DIM), temperature, diameter, density}
            //ParticleInitType<5, 0, 0, 0> pdata = {{0., 0., 300, 1e-2, 1.}, {}, {}, {}}; //2D
            ParticleInitType<6, 0, 0, 0> pdata = {{0., 0., 0., 298, 1e-6, 681.41}, {}, {}, {}}; //3D
            theSprayPC()->InitOnePerCell(Real(0.5),Real(0.5),Real(0.5),pdata);
          
            if (!particle_output_file.empty())
            {
              long cnt = SprayPC->TotalNumberOfParticles();// (bool only_valid=true, bool only_local=false) const;
              std::cout << "Writing " << cnt << "to " << particle_output_file.c_str() << std::endl;
              //theSprayPC()->WriteAsciiFile(particle_output_file,time);
	      // mr: have to use this below
	      //theSprayPC()->WriteAsciiFile(particle_output_file);
            }
        }
    }
}


void
PeleLM::particle_post_restart (const std::string& restart_file, bool is_checkpoint)
{
  amrex::Print() << "Do particle_post_restart." << '\n';

    if (level > 0)
       return;

    if (do_spray_particles)
    {
        BL_ASSERT(SprayPC == 0);

        SprayPC = new SprayParticleContainer(parent,&phys_bc);
        theSprayPC()->SetVerbose(particle_verbose);

        if (parent->subCycle())
        {
            VirtPC = new SprayParticleContainer(parent,&phys_bc);
            GhostPC = new SprayParticleContainer(parent,&phys_bc);
        }

        // fuel properties from user input file
        import_fuel_properties(&nfuel_species, fuel_mass_frac.data(),fuel_density.data(),
                               fuel_crit_temp.data(), fuel_latent.data(), fuel_boil_temp.data(), 
                               fuel_cp.data(),fuel_molwt.data(),fuel_indx.data());

        import_control_parameters(&is_heat_tran, &is_mass_tran, &is_mom_tran);

        // Make sure to call RemoveParticlesOnExit() on exit.
        //
        amrex::ExecOnFinalize(RemoveParticlesOnExit);

        theSprayPC()->Restart(parent->theRestartFile(), chk_particle_file, is_checkpoint);

        if (!particle_output_file.empty())
           theSprayPC()->WriteAsciiFile(particle_output_file);
    }
}

std::unique_ptr<MultiFab>
PeleLM::particle_derive(const std::string& name,
		       Real               time,
		       int                ngrow)
{
    BL_PROFILE("PeleLM::particle_derive()");

    if (theSprayPC() && name == "particle_count")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
	MultiFab    temp_dat(grids,dmap,1,0);
	temp_dat.setVal(0);
	theSprayPC()->Increment(temp_dat,level);
	MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
	return derive_dat;
    }
    else if (theSprayPC() && name == "total_particle_count")
    {
	//
	// We want the total particle count at this level or higher.
	//
	std::unique_ptr<MultiFab> derive_dat = particle_derive("particle_count",time,ngrow);

	IntVect trr(D_DECL(1,1,1));

	for (int lev = level+1; lev <= parent->finestLevel(); lev++)
	{
            auto ba = parent->boxArray(lev);
            const auto& dm = parent->DistributionMap(lev);
            MultiFab temp_dat(ba, dm, 1, 0);

            trr *= parent->refRatio(lev - 1);

            ba.coarsen(trr);
            MultiFab ctemp_dat(ba, dm, 1, 0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            theSprayPC()->Increment(temp_dat, lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab = temp_dat[mfi];
                FArrayBox& cfab = ctemp_dat[mfi];
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
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
    }
    else
    {
	return AmrLevel::derive(name,time,ngrow);
    }
}

void
PeleLM::particle_redistribute (int lbase, bool init_part)
{
    BL_PROFILE("PeleLM::particle_redistribute()");
    if (theSprayPC())
    {
        //  
        // If we are calling with init_part = true, then we want to force the redistribute
        //    without checking whether the grids have changed.
        //  
        if (init_part)
        {
            theSprayPC()->Redistribute(lbase);
            return;
        }

        //
        // These are usually the BoxArray and DMap from the last regridding.
        //
        static Vector<BoxArray>            ba;
        static Vector<DistributionMapping> dm;

        bool changed = false;

        int flev = parent->finestLevel();
	
        while ( parent->getAmrLevels()[flev] == nullptr ) {
            flev--;
	}
 
        if (ba.size() != flev+1)
        {
            ba.resize(flev+1);
            dm.resize(flev+1);
            changed = true;
        }
        else
        {
            for (int i = 0; i <= flev && !changed; i++)
            {
                if (ba[i] != parent->boxArray(i))
                    //
                    // The BoxArrays have changed in the regridding.
                    //
                    changed = true;

                if ( ! changed)
                {
                    if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
                        //
                        // The DistributionMaps have changed in the regridding.
                        //
                        changed = true;
                }
            }
        }

        if (changed)
        {
            //
            // We only need to call Redistribute if the BoxArrays or DistMaps have changed.
	    // We also only call it for particles >= lbase. This is
	    // because of we called redistribute during a subcycle, there may be particles not in
	    // the proper position on coarser levels.
            //
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "Calling redistribute because changed " << '\n';

            theSprayPC()->Redistribute(lbase);
            //
            // Use the new BoxArray and DistMap to define ba and dm for next time.
            //
            for (int i = 0; i <= flev; i++)
            {
                ba[i] = parent->boxArray(i);
                dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
            }
        }
        else
        {
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "NOT calling redistribute because NOT changed " << '\n';
        }
    }
}


void
PeleLM::TimestampParticles (int ngrow)
{
#if 1
    return;
    amrex::Abort("Need to fill in PeleLM::TimestampParticles");
#else
    static bool first = true;
    static int imax = -1;
    if (first)
    {
	first = false;

	ParmParse ppp("particles");

	// have to do it here, not in read_particle_params, because Density, ..., are set after
	// read_particle_params is called.

	int timestamp_density = 1;
	ppp.query("timestamp_density", timestamp_density);
	if (timestamp_density) {
	    timestamp_indices.push_back(Density);
	    std::cout << "Density = " << Density << std::endl;
	}

	int timestamp_temperature = 0;
	ppp.query("timestamp_temperature", timestamp_temperature);
	if (timestamp_temperature) {
	    timestamp_indices.push_back(Temp);
	    std::cout << "Temp = " << Temp << std::endl;
	}
	
	if (!timestamp_indices.empty()) {
	    imax = *(std::max_element(timestamp_indices.begin(), timestamp_indices.end()));
	}
    }

    if ( SprayPC && !timestamp_dir.empty())
    {
	std::string basename = timestamp_dir;
		
	if (basename[basename.length()-1] != '/') basename += '/';
	
	basename += "Timestamp";
	
	int finest_level = parent->finestLevel();
	Real time        = state[State_Type].curTime();

	for (int lev = level; lev <= finest_level; lev++)
	{
	    if (SprayPC->NumberOfParticlesAtLevel(lev) <= 0) continue;
	    
	    MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

	    if (imax >= 0) {  // FillPatchIterator will fail otherwise
		int ng = (lev == level) ? ngrow : 1;
		FillPatchIterator fpi(parent->getLevel(lev), S_new, 
				      ng, time, State_Type, 0, imax+1);
		const MultiFab& S = fpi.get_mf();
		SprayPC->Timestamp(basename, S    , lev, time, timestamp_indices);
	    } else {
		SprayPC->Timestamp(basename, S_new, lev, time, timestamp_indices);
	    }
	}
    }	
#endif
}

#endif
