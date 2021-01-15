# PeleLM

[![AMReX Badge](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/)
[![Documentation Status](https://readthedocs.org/projects/pelelm/badge/?version=latest)](https://pelelm.readthedocs.io/en/latest/?badge=latest)
[![License PeleLM](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)
![ConvergenceTesting](https://github.com/AMReX-Combustion/PeleLM/workflows/ConvergenceTesting/badge.svg)

## Overview

*PeleLM* is an adaptive-mesh low Mach number hydrodynamics code for reacting flows.
*PeleLM* supports Embedded Boundary method to represent complex geometries and is parallelized
with MPI + OpenMP for CPUs and MPI + CUDA or MPI + HIP for GPUs.

*PeleLM* is part of the Pele combustion Suite and *PeleLM* has a project [homepage](https://amrex-combustion.github.io/PeleLM/).
Use [this link](https://groups.google.com/forum/#!forum/pelelmusers/join) to sign up for the PeleLM user forum, where
updates and significant changes will be posted.  The forum is also where general questions can be posted about
building and running the code, processing code output, and details about the algorithm and its implementation.

## Documentation

*PeleLM* complete documentation is available on [ReadTheDoc](https://pelelm.readthedocs.io/en/latest/index.html).
It is also possible to build a local version of the documentation once you have obtained the source code using :

        cd ${PELELM_HOME}/Docs
        make html

### Getting started

A first simple 2D flame problem is available in the *PeleLM* QuickStart section:

https://pelelm.readthedocs.io/en/latest/GettingStarted.html

### Core Algorithm

The *PeleLM* governing equations and core algorithm are described in:

https://pelelm.readthedocs.io/en/latest/Model.html

### Tutorials

A set of self-contained tutorials describing more complex problems is also provided:

https://pelelm.readthedocs.io/en/latest/Tutorials.html

## Contributing

New contributions to *PeleLM* are welcome !

The *PeleLM* contributions workflow follows these steps:
1. Fork the main repository
2. Create an `AmazingNewFeature` branch implementing your changes 
3. Open a Pull Request from `AmazingNewFeature` on your fork to branch `development` of the main *PeleLM* repository

Follow [GitHub directions](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) 
to fork *PeleLM* main repo on your GitHub account, then clone the *PeleLM* dependencies 
([PelePhysics](https://github.com/AMReX-Combustion/PelePhysics),
[IAMR](https://github.com/AMReX-Codes/IAMR),
[AMReX](https://github.com/AMReX-Codes/amrex)) along with your own *PeleLM* fork on your local machine.

Then step into the *PeleLM* folder and add the main *PeleLM* repository as the `upstream` remote in order to keep track of the main repo :

       git add remote upstream https://github.com/AMReX-Combustion/PeleLM

At any point, you can update the `developement` branch of your local repository with changes implemented in the main *PeleLM* repo by pulling from `upstream` : 

        git checkout development
        git pull upstream development

You are now free to modify your own fork of *PeleLM*. To add a new feature to *PeleLM*, the procedure is:

1. Create a branch for the new feature from the `development` branch (locally) :

        git checkout development 
        git checkout -b AmazingNewFeature

2. and commit your changes to your local repo : 

        git commit -m "Developed AmazingNewFeature"

3. Alongside your development, regularly merge changes from the main repo `development` branch into your `AmazingNewFeature` branch,
fix any conficts, and push your changes to your GitHub fork :
   
        git pull upstream development     [Fix arising conflicts]
        git push -u origin AmazingNewFeature 

4. When you are ready to propose your new feature/improvement/bug fix to the main *PeleLM* repo, reiterate Step 3 and submit a Pull Request through the GitHub page from your fork onto the `development` branch of the main repo:

 - Click on the ``compare & pull request`` button to start your PR.
 - Provide a title and a short description for your PR:
   * what feature/fix do you propose
   * how did you test it
   * any other information deemed useful : does it modify the default *PeleLM* behavior ? ...
 - Press ``Create pull request``.

Please DO NOT write large Pull Requests, as they are very difficult and time-consuming to review.
As much as possible, split them into small targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.

## Acknowledgment

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
