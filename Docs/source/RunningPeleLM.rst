.. highlight:: rst

.. _sec:control:

`PeleLM` control
================

Physical Units
^^^^^^^^^^^^^^

`PeleLM` currently supports only MKS units.  All inputs and problem initialization should be
specified in MKS; output is in MKS unless otherwise specified.


Control parameters
^^^^^^^^^^^^^^^^^^

The `PeleLM` executable uses two inputs files at runtime to set and alter the
behavior of the algorithm and initial conditions.

The main inputs file, typically named ``inputs`` is used to
set `AMReX` parameters and the control flow in the C++ portions of
the `PeleLM` code.  Each parameter here has a namespace (like ``amr`` in
a parameter listed as ``amr.max_grid``).  Parameters set here are read using
the ``ParmParse`` class in `AMReX`.  The namespaces are typically used to group
control parameters by source code class or overall functionality.  There are,
for example, a large set of parameters that control the generation of the
solution-adaptive meshes during the run, as well as the location and content of
output files and logging information.  There are also a set of parameters that
control the details of the ``PeleLM`` time-stepping strategy, such as the
number of SDC iterations taken per time step, solver types and tolerances,
and algorithmic variations.  These latter control parameters are detailed
separately, in :ref:`sec:control:pelelm`.

The second inputs file, typically named ``probin`` is used by the
Fortran code that initializes the problem setup.  It is read at
problem initialization (via a Fortran ``namelist``) and the
problem-specific quantities are stored in a Fortran module ``probdata_module`` defined in the problem's ``probdata.f90`` file.

The ``inputs`` file is specified on the command-line directly as an argument to the executable.  The
associated ``probin`` file is specified inside the ``inputs`` file using the ``amr.probin_file`` parameter, e.g.,::

    amr.probin_file = my_special_probin

Here the Fortran code will read a file called ``my_special_probin``.

Working with ``probin`` files
-----------------------------

There are three different Fortran namelist's that can be defined in the
``probin`` file:

- ``&fortin`` is the main namelist read by the problem's ``probinit`` subroutine in the ``Prob_2d.f90`` or ``Prob_3d.f90`` file.

- ``&heattransin`` is used to set the ambient pressure and various algorithmic options.

- ``&flctin`` is used to read in turbulence files to set the space and time-dependent inflow.

- ``&control`` is used dynamically control the inflow velocity.

Consult the problem-specific files for how parameters set in these namelists are used.  When defining your
own problem, you should feel free to add any additional controls here that are useful for your case.


Working with ``inputs`` files
-----------------------------

**Important**: because the ``inputs`` file is handled by the `C++` portion of
the code, any quantities you specify in scientific notation, must take the
form ``1.e5`` and not ``1.d5``---the "d" specifier is not recognized.


Problem Geometry
----------------

The ``geometry.`` namespace is used by `AMReX` to define the
computational domain.  The main parameters here are:

1. ``geometry.prob_lo``: physical location of low corner of the
domain (type: ``Real``; must be set. A number is needed for each dimension in the problem
  
2. ``geometry.prob_hi``: physical location of high corner of the
domain (type: ``Real``; must be set. A number is needed for each dimension in the problem
  
3. ``geometry.coord_sys``: coordinate system, 0 = Cartesian,
1 = :math:`rz` (2D only), 2 = spherical (1D only); must be set.

4. ``geometry.is_periodic``: is the domain periodic in this direction?  ``0`` if false, ``1`` if true  (default: ``0 0 0``). An integer is needed for each dimension in the problem

As an example, the following::

    geometry.prob_lo = -0.1 -0.1 0.0
    geometry.prob_hi = +0.1 +0.1 0.2 
    geometry.coord_sys = 0 
    geometry.is_periodic = 0 1 0 

defines the domain to span the region from (-10,-10,0) cm at the lower left to
(10, 10, 20) cm at the upper right in physical coordinates, specifies a
Cartesian geometry, and makes the domain periodic in the :math:`y`-direction
only.

Domain Boundary Conditions
--------------------------

Boundary conditions are specified using integer keys that are interpreted
by `AMReX`.  The runtime parameters that we use are:

- ``ns.lo_bc``: boundary type of each low face  (must be set)
- ``ns.hi_bc``: boundary type of each high face (must be set)

The valid boundary types are: ::

    Interior
    Inflow
    Outflow
    Symmetry
    SlipWallAdiab
    NoSlipWallAdiab
    SlipWallIsotherm
    NoSlipWallIsotherm

Note: ``ns.lo_bc`` and ``ns.hi_bc`` must be consistent with 
``geometry.is_periodic``---if the domain is periodic in a particular
direction then the low and high bc's must be set to ``Interior`` for that direction.

As an example, the following: ::

    ns.lo_bc = Inflow SlipWallAdiab Interior 
    ns.hi_bc = Outflow SlipWallAdiab Interior

    geometry.is_periodic = 0 0 1

would define a problem with inflow in the low-:math:`x` direction,
outflow in the high-:math:`x` direction, adiabatic slip wall on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.

Resolution
----------

The grid resolution is specified by defining the resolution at the
coarsest level (level 0) and the number of refinement levels and
factor of refinement between levels.  The relevant parameters are:

- ``amr.n_cell``:  number of cells in each direction at the coarsest level (Integer > 0; must be set)

- ``amr.max_level``:  number of levels of refinement above the coarsest level (Integer >= 0; must be set)

- ``amr.ref_ratio``: ratio of coarse to fine grid spacing between subsequent levels (2 or 4; must be set)

- ``amr.regrid_int``: how often (in terms of number of steps) to regrid (Integer; must be set)

- ``amr.regrid_on_restart``: should we regrid immediately after restarting? (0 or 1; default: 0)

Note: if ``amr.max_level = 0`` then you do not need to set ``amr.ref_ratio`` or ``amr.regrid_int``.

Some examples: ::

    amr.n_cell = 32 64 64

would define the domain to have 32 cells in the :math:`x`-direction, 64 cells
in the :math:`y`-direction, and 64 cells in the :math:`z`-direction *at the
coarsest level*.  (If this line appears in a 2D inputs file then the
final number will be ignored.) ::

    amr.max_level = 2 

would allow a maximum of 2 refined levels in addition to the coarse
level.  Note that these additional levels will only be created only if
the tagging criteria are such that cells are flagged as needing
refinement.  The number of refined levels in a calculation must be
less than or equal to ``amr.max_level``, but can change in time and need not
always be equal to ``amr.max_level``. ::
 
    amr.ref_ratio = 2 4 

would set factor of 2 refinement between levels 0 and 1, and factor of 4
refinement between levels 1 and 2.  Note that you must have at least
``amr.max_level`` values of ``amr.ref_ratio`` (Additional values
may appear in that line and they will be ignored). Ratio values must be either or 2 or 4. ::

    amr.regrid_int = 2 2

tells the code to regrid every 2 steps.  Thus in this example, new
level 1 grids will be created every 2 level-0 time steps, and new
level 2 grids will be created every 2 level-1 time steps. If ``amr.regrid_int`` is less than 0 for any level, then regridding starting at that level will be disabled. If ``amr.regrid_int`` = -1 only, then we
never regrid for any level. Note that this is not compatible with ``amr.regrid_on_restart = 1``.


Regridding
----------

The details of the regridding strategy are described elsewhere; here we 
cover how the input parameters can control the gridding. The user defines functions which tag individual
cells at a given level if they need refinement (this is discussed in :ref:`sec:refcrit:pelelm`).
This list of tagged cells is
sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm
to create rectangular grids that contain the tagged cells. The relevant runtime parameters are:

- ``amr.regrid_file``: name of file from which to read the grids (text; default: no file)

If set to a filename, e.g.\ ``fixed_grids``, then list of grids
at each fine level are read in from this file during the gridding
procedure. These grids must not violate the ``amr.max_grid_size`` criterion.  The rest of the gridding procedure
described below will not occur if ``amr.regrid_file`` is set.

- ``amr.grid_eff``: grid efficiency (Real >0 and <1; default: 0.7)

- ``amr.n_error_buf``: radius of additional tagging around already tagged cells (Integer >= 0; default: 1)

- ``amr.max_grid_size``: maximum size of a grid in any direction (Integer > 0; default: 128 (2D), 32 (3D))

Note: ``amr.max_grid_size`` must be even, and a multiple of ``amr.blocking_factor`` at every level.
   
- ``amr.blocking_factor``:  all generated grid dimensions will be a multiple of this (Integer > 0; default: 2)

Note: ``amr.blocking_factor`` at every level must be a power of
2 and the domain size must be a multiple of ``amr.blocking_factor`` at level 0.
   
- ``amr.refine_grid_layout``: refine grids more if the number of processors is greater than the number of grids
  (0 if false, 1 if true; default: 1) 

Note also that ``amr.n_error_buf``, ``amr.max_grid_size`` and
``amr.blocking_factor`` can be read in as a single value which is
assigned to every level, or as multiple values, one for each level.

As an example, consider: ::

    amr.grid_eff = 0.9
    amr.max_grid_size = 64 
    amr.blocking_factor = 32

The grid efficiency, ``amr.grid_eff``, here means that during the grid
creation process, at least 90% of the cells in each grid at the level
at which the grid creation occurs must be tagged cells.  A higher
grid efficiency means fewer cells at higher levels, but may result
in the production of lots of small grids, which have inefficient cache
and OpenMP performance and higher communication costs.

The ``amr.max_grid_size`` parameter means that each of the final grids
will be no longer than 64 cells on a side at every level.
Alternately, we could specify a value for each level of refinement as:
``amr.max_grid_size = 64 32 16``, in which case our final grids
will be no longer than 64 cells on a side at level 0, 32 cells on a
side at level 1, and 16 cells on a side at level 2.  The ``amr.blocking_factor``
means that all of the final grids will be multiples of 32 at all levels.
Again, this can be specified on a level-by-level basis, like
``amr.blocking_factor = 32 16 8``, in which case the 
dimensions of all the final grids will be multiples of 32
at level 0, multiples of 16 at level 1, and multiples of 8 at level 2.


Getting good performance
------------------------

These parameters can have a large impact on the performance
of `PeleLM`, so taking the time to experiment with is worth the effort.
For example, having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance. The gridding algorithm proceeds in this order:

1. Grids are created using the Berger-Rigoutsos clustering algorithm, modified to ensure that all new fine grids are divisible by ``amr.blocking_factor``.

2. Next, the grid list is chopped up if any grids are larger than ``max_grid_size``. Note that because ``amr.max_grid_size`` is a multiple of ``amr.blocking_factor`` the ``amr.blocking_factor`` criterion is still satisfied.

3. Next, if ``amr.refine_grid_layout = 1`` and there are more processors than grids, and if ``amr.max_grid_size`` / 2 is a multiple of ``amr.blocking_factor``, then the grids will be redefined, at each level independently, so that the maximum length of a grid at level :math:`\ell`, in any dimension, is ``amr.max_grid_size``:math:`[\ell]` / 2.

4. Finally, if ``amr.refine_grid_layout = 1``,  and there are still more processors than grids, and if ``amr.max_grid_size`` / 4 is a multiple of ``amr.blocking_factor``, then the grids will be redefined, at each level independently, so that the maximum length of a grid at level :math:`\ell`, in any dimension, is ``amr.max_grid_size``:math:`[\ell]` / 4.


Simulation Time
---------------

There are two parameters that can define when a simulation ends:

- ``max_step``: maximum number of level 0 time steps (Integer greater than 0; default: -1)
- ``stop_time``: final simulation time (Real greater than 0;  default: -1.0)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (``max_step``) or by the final
simulation time (``stop_time``), or both. The code will stop at
whichever criterion comes first. Note that if the code reaches ``stop_time`` then the final time
step will be shortened so as to end exactly at ``stop_time``, not
past it.

As an example: ::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or 
the number of level 0 steps taken equals 1000, whichever comes first.


Time Step
---------

The following parameters affect the timestep choice:

- ``ns.cfl``: CFL number (Real > 0 and <= 1; default: 0.8)

- ``ns.init_shrink``: factor by which to shrink the initial time step (Real > 0 and <= 1; default: 1.0)

- ``ns.change_max``: factor by which the time step can grow in subsequent steps (Real >= 1; default: 1.1)

- ``ns.fixed_dt``: level 0 time step regardless of cfl or other settings (Real > 0; unused if not set)

- ``ns.dt_cutoff``: time step below which calculation will abort (Real > 0; default: 0.0)

As an example, consider: ::

    ns.cfl = 0.9 
    ns.init_shrink = 0.01 
    ns.change_max = 1.1
    ns.dt_cutoff = 1.e-20

This defines the ``cfl`` parameter to be 0.9,
but sets (via ``init_shrink``) the first timestep we take
to be 1% of what it would be otherwise.  This allows us to
ramp up to the numerical timestep at the start of a simulation.
The ``change_max`` parameter restricts the timestep from increasing
by more than 10\% over a coarse timestep.    Note that the time step
can shrink by any factor; this only controls the extent to which it can grow.
The ``dt_cutoff`` parameter will force the code to abort if the
timestep ever drops below :math:`10^{-20}`.  This is a safety feature---if the
code hits such a small value, then something likely went wrong in the
simulation, and by aborting, you won't burn through your entire allocation
before noticing that there is an issue.

Occasionally, the user will want to set the timestep explicitly, using ::

    ns.fixed_dt = 1.e-4

If ``ns.init_shrink`` not equal 1 then the first time step will in fact be
``ns.init_shrink`` * ``ns.fixed_dt``.


Restart
-------

`PeleLM` has a standard sort of checkpointing and restarting capability. 
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

- ``amr.check_file``: prefix for restart files (text; default: ``chk``) 

- ``amr.check_int``: how often (by level 0 time steps) to write restart files (Integer > 0; default: -1)

- ``amr.check_per``: how often (by simulation time) to write restart files (Real > 0; default: -1.0) Note that ``amr.check_per`` will write a checkpoint at the first timestep whose ending time is past an integer multiple of this interval. In particular, the timestep is not modified to match this interval, so you won't get a checkpoint at exactly the time you requested.

- ``amr.restart``: name of the file (directory) from which to restart
  (Text; not used if not set)

- ``amr.checkpoint_files_output``: should we write checkpoint files? (0 or 1; default: 1).  If you are doing a scaling study then set ``amr.checkpoint_files_output = 0`` so you can test scaling of the algorithm without I/O.

- ``amr.check_nfiles``: how parallel is the writing of the checkpoint files? (Integer $\geq 1$; default: 64). See the Software Section for more details on parallel I/O and the ``amr.check_nfiles`` parameter.

- ``amr.checkpoint_on_restart``: should we write a checkpoint immediately after restarting? (0 or 1; default: 0)


Note:

- You can specify both ``amr.check_int`` or ``amr.check_per``, if you so desire; the code will print a warning in case you did this unintentionally. It will work as you would expect -- you will get checkpoints at integer multiples of ``amr.check_int`` timesteps and at integer multiples of ``amr.check_per`` simulation time intervals.

- ``amr.plotfile_on_restart`` and ``amr.checkpoint_on_restart`` only take effect if ``amr.regrid_on_restart`` is in effect.

As an example,::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix ``chk_run`` will be generated every 10 level-0 time steps.  The directory names will be ``chk_run00000``, ``chk_run00010``, ``chk_run00020``, etc.  If instead you specify::

    amr.check_file = chk_run
    amr.check_per = 0.5

then restart files (really directories) starting with the prefix ``chk_run`` will be generated every 0.1 units of simulation time.  The directory names will be ``chk_run00000``, ``chk_run00043``, ``chk_run00061``, etc, where t = 0.1 after 43 level-0 steps, t = 0.2 after 61 level-0 steps, etc. To restart from ``chk_run00061``, for example, then set ::

    amr.restart = chk_run00061


Controlling Plotfile Generation
-------------------------------

The main output from `PeleLM` is in the form of plotfiles (which are
really directories).  The following options in the inputs file control
the generation of plotfiles:

- ``amr.plot_file``: prefix for plotfiles (text; default:
  ``plt``)

- ``amr.plot_int``: how often (by level-0 time steps) to write
  plot files (Integer > 0; default: -1)

- ``amr.plot_per``: how often (by simulation time) to write
  plot files (Real > 0; default: -1.0)

Note that ``amr.plot_per`` will write a plotfile at the first
timestep whose ending time is past an integer multiple of this interval.
In particular, the timestep is not modified to match this interval, so
you won't get a checkpoint at exactly the time you requested.

- ``amr.plot_vars``: name of state variables to include in plotfiles (valid options: ``ALL``, ``NONE`` or a list; default: ``ALL``)

- ``amr.derive_plot_vars``: name of derived variables to include in plotfiles (valid options: ``ALL``, ``NONE`` or a list; default: ``NONE``)

- ``amr.plot_files_output``: should we write plot files? (0 or 1; default: 1)

If you are doing a scaling study then set ``amr.plot_files_output = 0`` so you can test scaling of the algorithm without I/O.

- ``amr.plotfile_on_restart``: should we write a plotfile immediately after restarting?  (0 or 1; default: 0)
  
- ``amr.plot_nfiles``: how parallel is the writing of the plotfiles?  (Integer >= 1; default: 64)

All the options for ``amr.derive_plot_vars`` are kept in ``derive_lst`` in ``Pelelm_setup.cpp``.  Feel free to look at
it and see what's there. Also, you can specify both ``amr.plot_int`` or ``amr.plot_per``, if you so desire; the code will print a warning in case you did this unintentionally. It will work as you would expect -- you will get plotfiles at integer multiples of ``amr.plot_int`` timesteps and at integer multiples of ``amr.plot_per`` simulation time intervals. As an example: ::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plot files (really directories) starting with the prefix
``plt_run`` will be generated every 10 level-0 time steps.  The
directory names will be ``plt_run00000``, ``plt_run00010``, ``plt_run00020``, etc.


If instead you specify::

    amr.plot_file = plt_run
    amr.plot_per = 0.5

then restart files (really directories) starting with the prefix
``plt_run`` will be generated every 0.1 units of simulation time.  The
directory names will be ``plt_run00000``, ``plt_run00043``, ``plt_run00061``, etc, where t = 0.1 after 43 level-0 steps, t = 0.2 after 61 level-0 steps, etc.



Screen Output
-------------

There are several options that set how much output is written to the
screen as `PeleLM` runs:

- ``amr.v``: verbosity of ``Amr.cpp`` (0 or 1; default: 0)
- ``ns.v``: verbosity of ``NavierStokesBase.cpp`` (0 or 1; default: 0)
- ``diffusion.v``: verbosity of ``Diffusion.cpp`` (0 or 1; default: 0)
- ``mg.v``: verbosity of multigrid solver (allow values: 0,1,2,3,4; default: 0)  
- ``amr.grid_log``: name of the file to which the grids are written (text; not used if not set)  
- ``amr.run_log``: name of the file to which certain output is written (text; not used if not set)  
- ``amr.run_log_terse``: name of the file to which certain (terser) output is written (text; not used if not set)  
- ``amr.sum_interval``:  if > 0, how often (in level-0 time steps) to compute and print integral quantities (Integer; default: -1)

The integral quantities include total mass, momentum and energy in
the domain every ``ns.sum_interval`` level-0 steps.
The print statements have the form::

    TIME= 1.91717746 MASS= 1.792410279e+34

for example.  If this line is commented out then it will not compute and print these quanitities.


As an example: ::

    amr.grid_log = grdlog
    amr.run_log = runlog 

Every time the code regrids it prints a list of grids at all relevant
levels.  Here the code will write these grids lists into the file ``grdlog``.  Additionally, every time step the code prints certain statements to the screen (if ``amr.v = 1``), such as: ::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746 
    PLOTFILE: file = plt00001 

The ``run_log`` option will output these statements into ``runlog`` as well.

Terser output can be obtained via: ::

    amr.run_log_terse = runlogterse

This file, ``runlogterse`` differs from ``runlog``, in that it only contains lines of the form ::

    10  0.2  0.005

in which 10 is the number of steps taken, 0.2 is the
simulation time, and 0.005 is the level-0 time step.  This file
can be plotted very easily to monitor the time step.


.. _sec:control:pelelm:

`PeleLM` algorithm controls
---------------------------



Here, we document `PeleLM`-specific controls. --TODO
