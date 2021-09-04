.. highlight:: rst

..  _sec:QUICKSTART:

`PeleLM` Quickstart
===================
`PeleLM` was created in 2017 by renaming `LMC`, the low Mach code from `CCSE <https://ccse.lbl.gov>`_, 
and is built on the `AMReX` library, the `AMReX-Hydro` set of advection schemes, the `IAMR` code and the `PelePhysics` chemistry and thermodynamics library.
For the impatient, the following summarizes how to obtain `PeleLM` and all the supporting software
required, and how to build and run a simple case in order to obtain a first set of results.
A thorough discussion of the model equations, and time stepping algorithms in `PeleLM` is
given in :ref:`sec:model`.  More details about the make system are given in :ref:`sec:build:make`.
Parameters provided for runtime control of `PeleLM` are discussed in :ref:`sec:control`.  Visualization
of the results from `PeleLM` is discussed in :ref:`sec:visualization`.

Obtaining `PeleLM`
------------------

First, make sure that "git" is installed on your machine---we recommend version 1.7.x or higher.

Then, there are two options to obtain `PeleLM` and its dependencies:

1. PeleProduction
^^^^^^^^^^^^^^^^^

`PeleProduction` enables the user to obtain a consistent version of `PeleLM` and all its dependencies
 with a single git clone (from the user). This is the prefered option when one wants to use `PeleLM` 
but do not intend to make development into the code. More information on `PeleProduction` can be found 
on the `GitHub page <https://github.com/AMReX-Combustion/PeleProduction.git>`_.

   a. Download the `PeleProduction` repository and : ::

        git clone https://github.com/AMReX-Combustion/PeleProduction.git 

        cd PeleProduction 

   b. The first time you do this, you will need to tell git that there are submodules. Git will look at the ``.gitmodules`` file in this branch and use that : ::

        cd Submodules
        git submodule init
        git submodule update 

   c. Finally, get into the FlameSheet folder of the `PeleLM` submodule: ::

        cd PeleLM/Exec/RegTests/FlameSheet

2. Individual repositories
^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, all the individual dependencies of `PeleLM` can be obtained independently.
The user then needs to provide environment variables for each of `AMReX`, `IAMR`, `AMReX-Hydro`, `PelePhysics` and `PeleLM` installation path.
This method is intended for users wanting to modify the `PeleLM` source code and who are more comfortable with maintaining up-to-date the four repositories.

   a. Download the `AMReX` repository by typing: ::

        git clone https://github.com/AMReX-Codes/amrex.git

     This will create a folder called ``amrex/`` on your machine. Set the environment variable, ``AMREX_HOME``, on your
     machine to point to the path name where you have put `AMReX`::

        export AMREX_HOME=/path/to/amrex/
        
   b. Download the `IAMR` repository by typing: ::

        git clone https://github.com/AMReX-Codes/IAMR.git
    
     This will create a folder called ``IAMR/`` on your machine.
     Set the environment variable, ``IAMR_HOME``.
     Then switch to the ``development`` branch of IAMR: ::
     
        cd IAMR
        git checkout -b development origin/development

   c. Download the `AMReX-Hydro` repository by typing: ::

        git clone https://github.com/AMReX-Codes/AMReX-Hydro.git
    
     This will create a folder called ``AMReX-Hydro/`` on your machine.
     Set the environment variable, ``AMREX_HYDRO_HOME``.

   d. Clone the `PeleLM` and `PelePhysics` repositories: ::

        git clone git@github.com:AMReX-Combustion/PeleLM.git
        git clone git@github.com:AMReX-Combustion/PelePhysics.git

     This will create folders called ``PeleLM`` and ``PelePhysics`` on your machine.
     Set the environment variables, ``PELELM_HOME`` and ``PELE_PHYSICS_HOME``, respectively to where you put these.

   d. Periodically update each of these repositories by typing ``git pull`` within each repository.

   e. Finally, get into the ``FlameSheet`` folder of `PeleLM` : ::

        cd PeleLM/Exec/RegTests/FlameSheet

Building `PeleLM`
-----------------

In `PeleLM` each different problem setup is stored in its own
sub-folder under ``$(PELELM_HOME)/Exec/``, and a local version of the 
`PeleLM` executable is built directly in that folder (object libraries are not used to manage `AMReX`
and the application code).  In the following, we step through building a representative `PeleLM` executable.

1. Regardless of which path you decided to choose in order to get the `PeleLM` code and its dependencies, you should be now be in the ``FlameSheet`` folder.
If you have chosen Option 2 to get the `PeleLM` sources, you have already set the environement variable necessary to compile the executable.
If you have chosen the first option, you now have to modify the ``GNUmakefile`` to ensure that the variable ``TOP`` define on the first line
points to the ``Submodules`` folder of `PeleProduction` : ::

    TOP = /path/to/PeleProduction/Submodules

such that the following lines provide path to `PeleLM` and its dependencies. Note that an absolute path in needed.

2. Edit the ``GNUmakefile`` to ensure that the following are set: ::

    DIM = 2
    COMP = gnu (or your favorite C++/F90 compiler suite)
    DEBUG = FALSE
    USE_MPI = FALSE
    USE_OMP = FALSE

   If you want to try compilers other than those in the GNU suite, and you find that they don't
   work, please let us know.  Note that for centers managing their enviroments with "modules", the
   programming environment determining your available compiler should agree with your choice of ``COMP``
   in the ``GNUmakefile`` (e.g., ``PrgEnv-gnu`` module requires ``COMP=gnu``).

3. Start by building the Sundials Third Party Library used to integrate the chemistry: ::
   
    make TPL

   and finally build `PeleLM` executable: ::

    make

If successful, the resulting executable name will look something like ``PeleLM2d.gnu.ex``. Depending on your
compilation option the actual name of the executable might vary (including ``MPI``, or ``DEBUG``, ...).

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
packages `Vis-It` and `Paraview` `support the AMReX file format natively <https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html>`_,
as does the `yt` python package.  The standard tool used within the
`AMReX`-community is `Amrvis <https://github.com/AMReX-Codes/Amrvis>`_, a package developed and supported 
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data, however there are limited visualization
tools available in `Amrvis`, so most users make use of multiple tools depending on their needs.

For more information on how to use `Amrvis` and `VisIt`, refer to the `AMReX`
User's guide in the `AMReX` git repository for download/build/usage instructions.
