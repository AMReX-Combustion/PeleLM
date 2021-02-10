.. role:: cpp(code)
   :language: c++

.. _sec:build:make:

Building with GNU Make
======================

The build of `PeleLM` is managed with GNUmake.  For a specific case setup and
run configuration, you write your own make input files that define a number of
variables and build rules, and then invoke ``make`` to initiate the build process.
This will result in an executable upon successful completion. Temporary
files generated in the building process (such as object files) are stored in a
directory named  ``tmp_build_dir`` in such a way that allows multiple build
configurations to co-exist.

Dissecting a Simple Make File
-----------------------------

An example input file for GNU Make can be found in any of the example setup,
such as ``$(PELELM_HOME)/Exec/RegTests/FlameSheet``. Table :numref:`tab:makevars`
below shows a list of important variables.

.. raw:: latex

   \begin{center}

.. _tab:makevars:

.. table:: Important make variables

   +--------------------+-------------------------------------+-------------+
   | Variable           | Value                               | Default     |
   +====================+=====================================+=============+
   | AMREX_HOME         | Path to amrex                       | environment |
   +--------------------+-------------------------------------+-------------+
   | IAMR_HOME          | Path to IAMR                        | environment |
   +--------------------+-------------------------------------+-------------+
   | PELELM_HOME        | Path to PeleLM                      | environment |
   +--------------------+-------------------------------------+-------------+
   | PELE_PHYSICS_HOME  | Path to PelePhysics                 | environment |
   +--------------------+-------------------------------------+-------------+
   | COMP               | gnu, cray, ibm, intel, llvm, or pgi | none        |
   +--------------------+-------------------------------------+-------------+
   | DEBUG              | TRUE or FALSE                       | TRUE        |
   +--------------------+-------------------------------------+-------------+
   | DIM                | 2 or 3                              | none        |
   +--------------------+-------------------------------------+-------------+
   | USE_MPI            | TRUE or FALSE                       | FALSE       |
   +--------------------+-------------------------------------+-------------+
   | USE_OMP            | TRUE or FALSE                       | FALSE       |
   +--------------------+-------------------------------------+-------------+
   | USE_CUDA           | TRUE or FALSE                       | FALSE       |
   +--------------------+-------------------------------------+-------------+
   | USE_HIP            | TRUE or FALSE                       | FALSE       |
   +--------------------+-------------------------------------+-------------+

.. raw:: latex

   \end{center}

At the beginning of ``$(PELELM_HOME)/Exec/FlameSheet/GNUmakefile``, the make
variable ``AMREX_HOME`` is set to the path to the top directory of AMReX.  Note that in
the example :cpp:`?=` is a conditional variable assignment operator that only
has an effect if ``AMREX_HOME`` has not been defined (including in the
environment). The make variable can also take it value from the corresponding
environment variable, ``AMREX_HOME``, if it exists.  For
example in bash, one can set

.. highlight:: bash

::

    export AMREX_HOME=/path/to/amrex

prior to running ``make``.  alternatively, in tcsh one can set

.. highlight:: bash

::

    setenv AMREX_HOME /path/to/amrex

Path to `IAMR` (``IAMR_HOME``), `PelePhysics` (``PELE_PHYSICS_HOME``) and `PeleLM` (``PELELM_HOME``)
should also be provided in the same manner.

One must set the ``COMP`` variable to choose a compiler suite (for C, C++, f90).
Currently the list of supported compiler suites includes gnu, cray, ibm, intel, llvm,
and pgi. One must also set the ``DIM`` variable to either 1, 2, or 3, depending
on the dimensionality of the problem.

Variables ``DEBUG``, ``USE_MPI``, ``USE_OMP``, ``USE_CUDA`` and ``USE_HIP`` are optional with default set
to TRUE, FALSE, FALSE, FALSE and FALSE, respectively. Note that the last three entries are mutually exclusive.
The meaning of these variables should be obvious.  When ``DEBUG=TRUE``, aggressive compiler optimization flags are
turned off and assertions in source code are turned on. For production runs, ``DEBUG`` should be set to FALSE.

After defining these make variables, an application code may also wish to
include its own Make.package file (e.g., ``./Make.package``) or otherwise
directly append source files to the build system, using operator ``+=``.
Variables for various source file types are shown below.

    CEXE_sources
        C++ source files. Note that C++ source files are assumed to have a .cpp
        extension.

    CEXE_headers
        C++ headers with .h or .H extension.

    cEXE_sources
        C source files with .c extension.

    cEXE_headers
        C headers with .h extension.

    f90EXE_sources
        Free format Fortran source with .f90 extension.

    F90EXE_sources
        Free format Fortran source with .F90 extension.  Note that these
        Fortran files will go through preprocessing.

In the ``FlameSheet`` example, the extra source file, ``drm19Soln_seed_0.50.f`` is in a
directory that is already in the build system's search path.  Additional files,
that are local to this setup, such as ``pele_prob.cpp`` need to be added to the appropriate
file list explicitly as well.  If this case included files in a separate folder
(e.g., ``mysrcdir``), you will then need to add the following:

::

        VPATH_LOCATIONS += mysrcdir
        INCLUDE_LOCATIONS += mysrcdir

Here ``VPATH_LOCATIONS`` and ``INCLUDE_LOCATIONS`` are the search path for
source and header files, respectively.

Finally, `PeleLM` requires a number of defines and setup for every case that must be processed
into final filelists for building, and various defines for complilation -- these are managed
in the make include file ``$(PELELM_HOME)/Tools/Make/Make.PeleLM``.  In particular, this
file contains macros to find the chemistry mechanism/model files associated with the string
value of the ``Chemistry_Model`` variable.  Look in ``$(PELELM_HOME)/Tools/Make/Make.PeleLM``
for a list of currently recognized models, and to see which folder that the string maps to
in ``$(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models`` folder.  That folder will contain
a ``Make.package`` that appends the model-specific source files to the build list (typically
a C-source file generated by `FUEGO` from a CHEMKIN-compatible set of specification files -- see
the file ``$(PELE_PHYSICS_HOME)/README.rst`` for more information on model generation.

Tweaking the Make System
------------------------

The GNU Make build system is located in the `AMReX` source code distribution in
``$(AMREX_HOME)/Tools/GNUMake``.  You can read ``README.md`` and the make files there for more information.
Here we will give a brief overview.

Besides building executable, other common make commands include:

    ``make clean``
        This removes the executable, .o files, and the temporarily generated
        files. Note that one can add additional targets to this rule using the
        double colon (::)

    ``make realclean``
        This removes all files generated by make.

    ``make help``
        This shows the rules for compilation.

    ``make print-xxx``
        This shows the value of variable xxx. This is very useful for debugging
        and tweaking the make system.

Compiler flags are set in ``$(AMREX_HOME)/Tools/GNUMake/comps/``. Note that variables
like ``CC`` and ``CFLAGS`` are reset in that directory and their values in
environment variables are disregarded.  Site-specific setups (e.g., the MPI
installation) are in ``$(AMREX_HOME)/Tools/GNUMake/sites/``, which includes a generic
setup in ``Make.unknown``. You can override the setup by having your own
``sites/Make.$(host_name)`` file, where variable ``host_name`` is your host
name in the make system and can be found via ``make print-host_name``.  You can
also have an ``$(AMREX_HOME)/Tools/GNUMake/Make.local`` file to override various
variables. See ``$(AMREX_HOME)/Tools/GNUMake/Make.local.template`` for an example.


.. _sec:build:local:

Specifying your own compiler / GCC on macOS
-------------------------------------------

The ``$(AMREX_HOME)/Tools/GNUMake/Make.local`` file can also be used to specify your
own compile commands by setting the valiables ``CXX``, ``CC``, ``FC``, and
``F90``. This might be neccarry if your systems contains non-standard names for
compiler commands.

For example, mac OSX Xcode ships with its own (woefully outdated) version of GCC
(4.2.1). It is therefore recommended to install GCC using the `homebrew
<https://brew.sh>`_ package manager. Running ``brew install gcc`` installs gcc
with names reflecting the version number. If GCC 8.2 is installed, homebrew
installs it as ``gcc-8``. AMReX can be built using ``gcc-8`` without MPI by
using the following ``$(AMREX_HOME)/Tools/GNUMake/Make.local``:

::

    ifeq ($(USE_MPI),TRUE)
      CXX = mpicxx
      CC  = mpicc
      FC  = mpif90
      F90 = mpif90
    else
      CXX = g++-8
      CC  = gcc-8
      FC  = gfortran-8
      F90 = gfortran-8
    endif

For building with MPI, we assume ``mpicxx``, ``mpif90``, etc. provide access to
the correct underlying compilers.

Note that if you are building `PeleLM` using homebrew's gcc, it is recommended
that you use homebrew's mpich. Normally is it fine to simply install its
binaries: ``brew install mpich``. But if you are experiencing problems, we
suggest building mpich usinging homebrew's gcc: ``brew install mpich
--cc=gcc-8``.
