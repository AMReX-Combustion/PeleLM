`PeleLM` is an adaptive-mesh low Mach number hydrodynamics code for reacting flows.  `PeleLM` has a project
`homepage <https://amrex-combustion.github.io/PeleLM/>`_, and can be obtained via
`GitHub <https://github.com/AMReX-Combustion/PeleLM>`_.

.. image:: https://readthedocs.org/projects/pelelm/badge/?version=latest:target: https://pelelm.readthedocs.io/en/latest/?badge=latest:alt: Documentation Status


Getting Started
---------------

* To compile and run the `Pele` suite of codes, one needs a C++ compiler that supports the C++11 standard and a Fortran compiler that supports the 2003 standard.  A hierarchical strategy for parallelism is supported, based MPI + OpenMP.  The codes work with all major MPI and OpenMP implementations.  The codes should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 4.8.4 and above.


To build `PeleLM` and run a sample 2D flame problem:

1. Set the environment variable, AMREX_HOME, and clone a copy of `AMReX` there ::

    export AMREX_HOME=<location for AMReX>    
    git clone git@github.com:AMReX-Codes/amrex.git ${AMREX_HOME}

2. Set the environment variable, IAMR_HOME, and clone a copy of `IAMR` there ::

    export IAMR_HOME=<location for IAMR>
    git clone git@github.com:AMReX-Codes/IAMR.git ${IAMR_HOME}

3. Set the environment variable, PELELM_HOME, and clone a copy of `PeleLM` there ::

    export PELELM_HOME=<location for PeleLM>
    git clone git@github.com:AMReX-Combustion/PeleLM.git ${PELELM_HOME}

4. Set the environment variable, PELE_PHYSICS_HOME, and clone a copy of `PelePhysics` there ::

    export PELE_PHYSICS_HOME=<location for PelePhysics>
    git clone git@github.com:AMReX-Combustion/PelePhysics.git ${PELEPHYSICS_HOME}

5. Move to an example build folder, build an executable ::

    cd ${PELELM_HOME}/Exec/FlameInABox
    make

Dependencies
------------

`PeleLM` was created in 2017 by renaming `LMC`, the low Mach code from
`CCSE <https://ccse.lbl.gov>`_, and is built on the `AMReX` library
and the `IAMR` code (see above).

Development model
-----------------

To add a new feature to PeleLM, the procedure is:

1. Create a branch for the new feature (locally) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleLM repository ::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4.  Submit a merge request through the github project page - be sure you are requesting to merge your branch to the development branch.

Documentation
-------------
Documentation for PeleLM is under development in the Docs directory.  To build ::

    cd ${PELELM_DIR}/Docs/UsersGuide
    make


