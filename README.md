# PeleLM

[![AMReX Badge](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/)
[![Documentation Status](https://readthedocs.org/projects/pelelm/badge/?version=latest)](https://pelelm.readthedocs.io/en/latest/?badge=latest)
[![License PeleLM](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)

## Overview

`PeleLM` is an adaptive-mesh low Mach number hydrodynamics code for reacting flows.  `PeleLM` has a project
`homepage <https://amrex-combustion.github.io/PeleLM/>`_, and can be obtained via
`GitHub <https://github.com/AMReX-Combustion/PeleLM>`_.  
Use `this link <https://groups.google.com/forum/#!forum/pelelmusers/join>`_ to 
to sign up for the PeleLM user forum, where
updates and significant changes will be posted.  The forum is also where general questions can be posted about
building and running the code, processing code output, and details about the algorithm and its implementation.

## Documentation

`PeleLM` complete documentation is available on `ReadTheDoc <https://pelelm.readthedocs.io/en/latest/index.html>`_.
It is also possible to build a local version of the documentation using: ::
    cd ${PELELM_HOME}/Docs
    make html

Getting started
---------------

A first sample 2D flame problem is available in the `PeleLM` QuickStart section:

`https://pelelm.readthedocs.io/en/latest/GettingStarted.html <https://pelelm.readthedocs.io/en/latest/GettingStarted.html>`_

Core Algorithm
--------------

The `PeleLM` governing equations and core algorithm are described in:

`https://pelelm.readthedocs.io/en/latest/Model.html <https://pelelm.readthedocs.io/en/latest/Model.html>`_

Tutorials
---------

A set of self-contained tutorials describing more complex problems is also provided:

`https://pelelm.readthedocs.io/en/latest/Tutorials.html <https://pelelm.readthedocs.io/en/latest/Tutorials.html>`_

## Contributing

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

4. Submit a merge request through the github project page - be sure you are requesting to merge your branch to the development branch.

## Acknowledgment

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
