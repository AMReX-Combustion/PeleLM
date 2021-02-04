.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:newcase:

Setting up a new `PeleLM` Case
==============================

In order to set up and run a new case in `PeleLM`, the user must provide problem-specific code for two main tasks

- Initial conditions
- Boundary conditions

These functions are typically collected into a single subfolder in ``${PELELM_HOME}/Exec``, such as ``FlameSheet``.
The user can organize these tasks in any way that is convenient - the examples distributed with ``PeleLM``
represent a certain style for managing this with some level of flexibility, but the basic requirement is
simply that source be linked into the build for the functions ``pelelm_initdata`` for initial conditions
and ``bcnormal`` for boundary conditions.

Initial Conditions
------------------

At the beginning of a `PeleLM` run, for each level, after grids are generated, the cell-centered values of
the state must be initialized.  In the code, this is done in an ``MFIter`` loop over grids, and a call to
the user's initialization function, ``pelelm_initdata``, that must be provided:

::

  for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
      const Box& box = mfi.validbox();
      auto sfab = S_new.array(mfi);

      amrex::ParallelFor(box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {    
        pelelm_initdata(i, j, k, sfab, geomdata, *lprobparm, lpmfdata);
      });  
  }

where ``(i,j,k)`` are the cell indices, ``sfab`` is a light data pointer to the initial state MultiFab and ``geomdata``, ``lprobparm``
and ``lpmfdata`` are container for the geometry, user-define input and PMF data.
The associated user function (in pelelm_prob.H) will provide a value for each entry of the state (velocity, density, mass fraction, ...) :

::

   AMREX_GPU_DEVICE
   AMREX_FORCE_INLINE
   void
   pelelm_initdata (int i, int j, int k, 
                    amrex::Array4<amrex::Real> const& state,
                    amrex::GeometryData const& geomdata,
                    ProbParm const& prob_parm,
                    PmfData const *pmf_data)
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


      ...


      state(i,j,k,DEF_Temp) = prob_parm.Temp;

      for (int n = 0; n < NUM_SPECIES; n++){
        massfrac[n] = 1.0/NUM_SPECIES;
      }


      state(i,j,k,Xvel) = 0.0;
      state(i,j,k,Yvel) =  prob_parm.Vel;
   #elif ( AMREX_SPACEDIM == 3 ) 
      state(i,j,k,Zvel) = 0.0;
   #endif

      amrex::Real rho_cgs, P_cgs;
      P_cgs = prob_parm.P_mean * 10.0;

      EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
      state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

      EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
      state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
      }
   }

Note home the geometry data are retrived from the ``geomdata`` object to obtain the coordinate of each cell center.
The state data indices (``Xvel``, ``Yvel``, ``Density``, ... ) are prescribed in the ``$(PELELM_HOME)/Source/IndexDefines.H`` file. 
Note that the conserved states are stored for species and
enthalpy (i.e., :math:`\rho Y_i` and :math:`\rho h`); these are the variables that the user must fill
in the initial and boundary condition routines.  Typically, however, the primitive state
(i.e., :math:`Y_i` and :math:`T`) is known directly.  If that is the case, the user can make use of
the compiled-in model-specific equation-of-state routines (``EOS::``) to translate primitive to conserved state
values. Consult the example setups provided to see how to call these routines, and how to load the
final values required for initial data.

The runtime option (such as initial temperature, inlet velocity, ...) are gathered in the ``ProbParm`` C++ structure defined in ``pelelm_prob_parm.H`` and filled from the input file in ``pelelm_prob.cpp`` using `AMReX` parser. This structure can be modified by the user to hold any data necessary for initial or boundary conditions.

Boundary Conditions
-------------------

In `PeleLM`, a single function is used to fill all the state
component at physical boundaries. The function ``bcnormal``
is in the ``pelelm_prob.H`` file. The main objective of this
function is to fill the ``s_ext`` array fill boundary state
data.
The function ``bcnormal`` is called on each side (`lo` or `hi`)
for each spatial dimension and will be used to fill the ghost cells
of the state variables for which the `PeleLM` internal boundary condition is 
``EXT_DIR`` (external Dirichlet) on that face. For example, specifying
a `PeleLM` ``Inflow`` boundary condition on the lower face in the y direction in the
input file leads to an ``EXT_DIR`` for species mass fraction, which then need 
to be provided in ``bcnormal``. An example of the ``bcnormal`` of the FlameSheet
is presented here:

::

   AMREX_GPU_DEVICE
   AMREX_FORCE_INLINE
   void
   bcnormal(
     const amrex::Real x[AMREX_SPACEDIM],
     amrex::Real s_ext[DEF_NUM_STATE],
     const int idir,
     const int sgn,
     const amrex::Real time,
     amrex::GeometryData const& geomdata,
     ProbParm const& prob_parm,
     ACParm const& ac_parm,
     PmfData const *pmf_data)
   {
     const amrex::Real* prob_lo = geomdata.ProbLo();
     const amrex::Real* prob_hi = geomdata.ProbHi();
     amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
     amrex::Real molefrac[NUM_SPECIES] = {0.0};
     amrex::Real massfrac[NUM_SPECIES] = {0.0};

     if (sgn == 1) {
       PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);
   
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
       P_cgs = prob_parm.P_mean * 10.0;
   
       EOS::PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
       s_ext[Density] = rho_cgs * 1.0e3;

       EOS::TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
       s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion
   
       for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
       }
     }
   }


The ``sgn`` input takes a value of 1 on the low face and -1 on the high face,
while ``ìdir`` provide the spatial direction (0, 1 or 2 corresponding to  X, Y or Z, respectively).
This allow to differentiate between the various boundary conditions when more than 1 ``ÈXT_DIR``
is needed. In this example, the boundary conditions are extracted from a pre-computed premixed flame
which data are stored in the ``pmf_data`` structure.

Here, we've made use of a local convenience function,
``bcnormal`` endowed with the knowledge of all boundary values, and
extract the appropriate quantity from the results of that call.  This
was done to localize all boundary condition calculations to a single
routine in the code, and helps to preserve consistency.  This is only
one style though, and as long as appropriate Dirichlet values are set
for this state, it makes no difference how the work is organized.
For example, data may be provided by interpolating "live data" being
actively generated by a co-running separate code, by interpolating data
files, evaluating functional forms, etc.

Note that although the array structure to be filled contains valid cell-centered state data where it
overlaps the valid domain, the values set in the grow cells of the container will be applied on the
boundary face of the corresponding cells.  Internally, all `PeleLM` code understands to apply
Dirichlet conditions on the boundary faces.

.. _sec:refcrit:pelelm:

Refinement Criteria
-------------------

The dynamic creation and destruction of grid levels is a fundamental part of `PeleLM`'s capabilities. The
process for this is described in some detail in the `AMReX` documentation, but we summarize the key points
here.

At regular intervals (set by the user), each Amr level that is not the finest allowed for the run
will invoke a "regrid" operation.  When invoked, a list of error tagging functions is traversed. For each,
a field specific to that function is derived from the state over the level, and passed through a kernel
that "set"'s or "clear"'s a flag on each cell.  The field and function for each error tagging quantity is
identified in the setup phase of the code where the state descriptors are defined (i.e., in `PeleLM_setup.cpp`).
Each function in the list adds or removes to the list of cells tagged for refinement. This final list of tagged
cells is sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm to create rectangular grids
which will define a new finer level (or set of levels).  State data is filled over these new grids, copying where
possible, and interpolating from coarser level when no fine data is available.  Once this process is complete,
the existing Amr level(s) is removed, the new one is inserted into the hierarchy, and the time integration
continues.

The traditional `AMReX` approach to setting up and controlling the regrid process involves explicitly
creating ("hard coding") a number of functions directly into `PeleLM`'s setup code. (Consult the source code
and `AMReX` documentation for precisely how this is done).  `PeleLM` provides a limited capability to augment
the standard set of error functions that is based entirely on runtime data specified in the inputs (ParmParse)
data.  The following example portion of a ParmParse'd input file demonstrates the usage of this feature:

::

      amr.refinement_indicators = flame_tracer lo_temp gradT

      amr.flame_tracer.max_level = 3
      amr.flame_tracer.value_greater = 1.e-6
      amr.flame_tracer.field_name = Y(H)

      amr.lo_temp.max_level = 1
      amr.lo_temp.value_less = 450
      amr.lo_temp.field_name = temp

      amr.gradT.max_level = 2
      amr.gradT.adjacent_difference_greater = 20
      amr.gradT.field_name = temp
      amr.gradT.start_time = 0.001
      amr.gradT.end_name = 0.002

Here, we have added three new custom-named criteria -- ``flame_tracer``: cells with the mass fraction of H greater than 1 ppm;
``lo_temp``: cells with T less than 450K, and ``gradT``: cells having a temperature difference of 20K from that of their
immediate neighbor.  The first will trigger up to Amr level 3, the second only to level 1, and the third to level 2.
The third will be active only when the problem time is between 0.001 and 0.002 seconds.

Note that these additional user-created criteria operate in addition to those defined as defaults.  Also note that
these can be modified between restarts of the code.  By default, the new criteria will take effect at the next
scheduled regrid operation.  Alternatively, the user may restart with ``amr.regrid_on_restart = 1`` in order to
do a full (all-levels) regrid after reading the checkpoint data and before advancing any cells.
