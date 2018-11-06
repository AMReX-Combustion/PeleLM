.. highlight:: rst

`PeleLM` Quickstart
===================

For the impatient, the following summarizes how to obtain `PeleLM` and all the supporting software
required, and how to build and run a simple case in order to obtain a first set of results.
A thorough discussion of the model equations, and time stepping algorithms in `PeleLM` is
given in :ref:`sec:model`.  More details about the make system are given in :ref:`sec:build:make`.
Parameters provided for runtime control of `PeleLM` are discussed in :ref:`sec:control`.  Visualization
of the results from `PeleLM` is discussed in :ref:`sec:visualization`.

Obtaining `PeleLM`
------------------

First, make sure that "git" is installed on your machine---we recommend version 1.7.x or higher. Then...

1. Download the `AMReX` repository by typing: ::

    git clone https://github.com/AMReX-Codes/amrex.git

This will create a folder called ``amrex/`` on your machine. Set the environment variable, ``AMREX_HOME``, on your
machine to point to the path name where you have put `AMReX`::

        export AMREX_HOME=/path/to/amrex/
        
2. Download the `IAMR` repository by typing: ::

    git clone https://github.com/AMReX-Codes/IAMR.git
    
This will create a folder called ``IAMR/`` on your machine.
Set the environment variable, ``IAMR_HOME``.

3. Clone the `Pele` repositories: ::

    git clone git@github.com:AMReX-Combustion/PeleLM.git
    git clone git@github.com:AMReX-Combustion/PelePhysics.git

This will create folders called ``PeleLM`` and ``PelePhysics`` on your machine.
Set the environment variables, ``PELELM_HOME`` and ``PELE_PHYSICS_HOME``, respectively to where you put these.

4. Periodically update each of these repositories by typing ``git pull`` within each repository.


Building `PeleLM`
-----------------

In `PeleLM` each different problem setup is stored in its own
sub-folder under ``$(PELELM_HOME)/Exec/``, and a local version of the 
`PeleLM` executable is built directly in that folder (object libraries are not used to manage `AMReX`
and the application code).  In the following, we step through building a representative `PeleLM` executable.

1. We will work in the folder containing setup for the `FlameInABox` problem in 2D
(``$(PELELM_HOME)/Exec/FlameInABox``).
In this setup, cold fuel enters the domain bottom and passes through a flame sheet.
Hot products exit the domain at the top.  The sides of the domain are periodic, and the coordinates are
cartesian. Go to the problem-specific source folder::

    cd $(PELELM_HOME)/Exec/FlameInABox

2. Edit the ``GNUmakefile`` to ensure that the following are set::

    DIM = 2
    COMP = gnu (or your favorite C++/F90 compiler suite)
    DEBUG = TRUE
    USE_MPI = FALSE
    USE_OMP = FALSE

If you want to try compilers other than those in the GNU suite, and you find that they don't
work, please let us know.  Note that for centers managing their enviroments with "modules", the
programming environment determining your available compiler should agree with your choice of ``COMP``
in the ``GNUmakefile`` (e.g., ``PrgEnv-gnu`` module requires ``COMP=gnu``).
If successful, the resulting executable name will look something like ``PeleLM2d.gnu.ex``.


Running `PeleLM`
----------------

1. `PeleLM` takes an input file as its first command-line argument.  The file
contains a set of parameter definitions that will override defaults set in the code.
To run `PeleLM` in serial with an example inputs file, type::

    ./PeleLM2d.gnu.ex inputs.2d-regt

2. While running, `PeleLM` typically generates subfolders in the current folder that are named ``plt00000/``, ``plt00020/``, etc, and ``chk00000/``, ``chk00020/``, etc. These are "plotfiles" and "checkpoint" files. The plotfiles are used for visualization of derived fields; the checkpoint files are used for restarting the code.


The output folders contain a collection of ASCII and binary files.  The field data is generally written in a self-describing binary format; the ASCII header files provide additional metadata to give the AMReX-compatible readers context to the field data.


Visualization of the results
----------------------------

There are several options for visualizing the data.  The popular
packages `Vis-It` and `Paraview` support the `AMReX` file format natively (currently called ``BoxLib`` format),
as does the `yt` python package.  The standard tool used within the
`AMReX`-community is `Amrvis`, a package developed and supported 
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data, however there are limited visualization
tools available in `Amrvis`, so most users make use of multiple tools depending on their needs.

For more information on how to use `Amrvis` and `VisIt`, refer to the `AMReX`
User's guide in the `AMReX` git repository for download/build/usage instructions.

