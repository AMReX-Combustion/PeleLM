# 1.0.1
Release version 1.0.1 (In preparation)

### Minor changes
   * Revert some of HIP changes (AMREX_GPU_DEVICE -> AMREX_GPU_HOST_DEVICE)


# 1.0.0
Release version 1.0.0 (January, 2021)
### First official GitHub release (Tagged v1.0.0)
### Major changes:
   * Implementation of Embedded Boundary (EB) method to represent geometry
   * The enthalpy equation now involves deltaT iterative solve for the diffusion
     and the advection term uses rhoH face states reconstructed from extrapolated
     species and tempeture.
   * Code kernels moved from F90 to C++ in order to port PeleLM to MPI+CUDA
   
### Minor changes (probably incomplete for this release):
   * Bug fix on the FluxRegisters and MacSync
   * Re-working the wbar term as a runtime option and bug fix
   * Some LES capabilities are inherited from changes in IAMR (for velocity only)
   * Corrections and improvements of the documentation
   * Addition of two complete tutorials: Flow Past Cylinder (EB) and Triple Flame
   * Implementation of GitHub Actions CI testing compilation/execution on the development branch

### Dependencies stable \#
   * PelePhysics : 8730cebf03
   * IAMR        : 09c3848b20
   * AMReX       : 889af2f8a8

# 0.1.0
Release version 0.1.0 (June, 2019)
### Major cleaning of PeleLM:
  * F77 code and .H include files have been removed and all data are defined in F90 modules.
  * Prob_$(dim)d.F90 now only requires amrex_probinit and init_data.
    (See note 1 below, and in documention)
  * Boundary Conditions are now imposed via keywords, rather than integer.
  * All error tagging functions have been removed and replaced by a more generic procedure.
    (See note 2 below, and in documention)
  * PeleLM no longer depends on the DERIVE functions from IAMR. Common functions have been duplicated within PeleLM
    This helps to remove confusing dependencies of user defined routines such as FORT_XVELFILL.
  * User-specific routines "bcfunction" and "zero_visc" are now in user_defined_fcts_$(dim).F90
    These routines may not be needed but are nonetheless compiled, so an empty version is present in
    Src_$(dim) folder, and will be linked into the build if case-specific versions do not exist
    (see note 1 below, and in documention)
  * The Exec folder has been reorganized. All previous cases have been moved to Obsolete and no longer compile.
    The FIAB case has been cleaned and renamed as FlameSheet in RegTests. Production cases have been removed.
    Production cases thar are part of a collaboration with the CCSE combustion team have been moved to the
    GitHub repo PeleProduction.
  * Several derive functions have been added to plot mass fractions, mole fractions and molecular weight. 
    Old flags have been removed.
  * IAMR and PeleLM have now their Initialize_specific() routine. They still share the Initialize() routine.

### New regression test cases in RegTests dir:
    1D PMF in X and Y directions, with and without control,
    as well as a Convected Vortex (CoVo) case and the Natural Convection case.

### Regression test cases and convergence tests are managed with case-specific python scripts (in Exec folders)

### PeleLM now makes full use of PelePhysics for EOS, transport properties and chemistry.  ChemDriver has been removed.
    PeleLM now works with Sutherland and Constant modules from PelePhysics, however note that reference values have to be put in CGS format, not MKS.

### Several bugs where fixed in the MPI/OpenMP tiling implementation

### Wall Boundary condition are expanded to more completely specify the energy boundary conditions (isothermal, adiabatic).

Notes
1. In earlier PeleLM versions, the Prob*.F files contained all the case-specific functions to tag cells for refinement, to implement boundary conditions,
   and to initialize the solution. A usage pattern arose that resulted in substantial code duplication. That pattern has been formalized and
   lifted into the Source folder as "bc_fill_${dim}d.F90".  All boundary fill functions there call "bcfunction" with sufficient location/time
   information to completely specify most common boundary conditions.  If this is sufficient for your case, you need only make a local copy of
   the "user_defined_fcts_$(dim).F90" and edit "bcfunction".  If either of these files are not provided, the default versions will be built.

2. Error tagging functions are now created on the fly, based on data in the inputs file.  Currently, the user can specify a min or max
   threshold value of a specific derivable quantity, and a time period over which to act.  Additionally, local gradient (actually, adjacent
   value difference) and vorticity taggers are also available.  See the example setups in the Exec folder, and the documentation, for details.

