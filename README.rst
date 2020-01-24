.. image:: https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg
  :target: https://amrex-codes.github.io/amrex/

PeleLM - README
===============

`PeleLM` is an adaptive-mesh low Mach number hydrodynamics code for reacting flows.  `PeleLM` has a project
`homepage <https://amrex-combustion.github.io/PeleLM/>`_, and can be obtained via
`GitHub <https://github.com/AMReX-Combustion/PeleLM>`_.  
Use `this link <https://groups.google.com/forum/#!forum/pelelmusers/join>`_ to 
to sign up for the PeleLM user forum, where
updates and significant changes will be posted.  The forum is also where general questions can be posted about
building and running the code, processing code output, and details about the algorithm and its implementation.

Current project build status
----------------------------

.. image:: https://readthedocs.org/projects/pelelm/badge
  :target: https://pelelm.readthedocs.io/en/latest/index.html

Getting started with PeleLM
---------------------------

* To compile and run the `Pele` suite of codes, one needs a C++ compiler that supports the C++11 standard and a Fortran compiler that supports the 2003 standard.  A hierarchical strategy for parallelism is supported, based MPI + OpenMP.  The codes work with all major MPI and OpenMP implementations.  The codes should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 4.8.4 and above.


* PeleLM depends on several separate GitHub repositories, each under active development. This can significantly complicate the required effort to keep all the required software up to date and internally compatible. Recently, we moved the Pele codes to a new software management style based on git "submodules", which dramatically simplifies the initial install/build/run procedure: ::

    git clone --recursive https://github.com/AMReX-Combustion/PeleProduction.git
    cd PeleLMruns/FlameSheet2D
    make -j 12
    mpiexec -np 8 ./PeleLM2d.gnu.MPI.ex inputs.2d-regt
            
* Notes

   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least 4.8.4, and MPI should be at least version 3).
   B. The example is 2D premixed flame, flowing vertically upward through the domain with no gravity. The lateral boundaries are periodic.  A detailed methane model is used.  The solution is initialized with a wrinkled (perturbed) 1D steady flame solution computed using the PREMIX code.  Two levels of solution-adaptive refinement are automatically triggered by the presence of the flame intermediate, H.
   C. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with CCSE's Amrvis (<https://github.com/AMReX-Codes/Amrvis>) for example. Please look at the AMReX documentation for further options about visualization (<https://amrex-codes.github.io/amrex/docs_html/Visualization.html>).


Dependencies
------------

`PeleLM` was created in 2017 by renaming `LMC`, the low Mach code from
`CCSE <https://ccse.lbl.gov>`_, and is built on the `AMReX` library
and the `IAMR` code (see above).

Development model
-----------------

To add a new feature to PeleLM (or the other subrepositories), the procedure is:

1. Create a branch for the new feature (locally, within the appropriate submodule folder) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleLM repository ::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4.  Submit a merge request through the github project page - be sure you are requesting to merge your branch to the development branch of any of the repositories (`master` is only updated via offline merges from `development`).

Documentation
-------------
`Documentation <https://pelelm.readthedocs.io/en/latest/index.html>`_ for PeleLM is under continuous development.  Please forward any suggestions as a GitHub "issue" at the PeleLM code page


Acknowledgment
--------------
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
