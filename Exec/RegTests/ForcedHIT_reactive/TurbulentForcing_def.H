#include <TurbulentForcing_params.H>
#include <AMReX_BaseFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>


// factor by which to reduce sampling for faster performance
AMREX_GPU_MANAGED int TurbulentForcing::ff_factor;
// make the forcing divergence free?
AMREX_GPU_MANAGED bool TurbulentForcing::div_free_force;
// how many modes to use
AMREX_GPU_MANAGED int TurbulentForcing::nmodes;
// don't use any modes below mode_start. We probably don't need this
AMREX_GPU_MANAGED int TurbulentForcing::mode_start;

amrex::Real* TurbulentForcing::forcedata;

void
TurbulentForcing::read_turbulent_forcing_params ()
{
    using namespace amrex;

    // read in parameters for turbulent forcing
    ParmParse pp("turb");

    // Make user set nmodes. This must match up with nmodes used in making MagicFile
    nmodes = -1;
    pp.get("nmodes", nmodes);

    div_free_force = true;
    pp.query("div_free_force", div_free_force);

    ff_factor = 4;
    pp.query("ff_factor", ff_factor);

    mode_start = 0;
    pp.query("mode_start", mode_start);

    //    
    // Read magic forcing data from file
    //

    // tmp CPU storage that holds everything in one flat arrray
    const int num_elmts=array_size*array_size*array_size;
    const int tmp_size = num_fdarray*num_elmts;
    Real tmp[tmp_size];

    std::string forceFile = "forcedata.dat";
    pp.query("force_file", forceFile);
    FILE* binfile = fopen(forceFile.c_str(),"rb") ;
    if (binfile==NULL)
    {
	fputs ("Unable to open forceFile 'forcedata.dat'\n",stderr);
	exit (EXIT_FAILURE);
    }

    int result;
    for (int j=0; j < num_fdarray; j++) {
        //printf("Reading force array %i : ...\n", j);
        result=fread(&tmp[j*num_elmts],sizeof(Real),num_elmts,binfile);
	if (result != num_elmts)
	{
	    fputs ("Error reading forceFile\n",stderr);
	    exit (EXIT_FAILURE);
	}
    }
    fclose(binfile);

// Now allocate forcedata and copy in tmp array.
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
	forcedata = static_cast<Real*>(The_Arena()->alloc(tmp_size*sizeof(Real)));
	Gpu::htod_memcpy_async(forcedata, tmp, tmp_size*sizeof(Real));
    }
    else
#endif
    {
        forcedata = static_cast<Real*>(The_Pinned_Arena()->alloc(tmp_size*sizeof(Real)));
        std::memcpy(forcedata, tmp, tmp_size*sizeof(Real));
    }
}
