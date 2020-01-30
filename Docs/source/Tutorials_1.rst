.. role:: cpp(code)
   :language: c++

.. _sec:tutorial1:

Tutorial - A simple triple flame
===============================

.. _sec:TUTO1::Intro:

Introduction
------------------------------

Laminar flames have the potential to reveal the fundamental structure of combustion 
without the added complexities of turbulence. 
They also aid in our understanding of the more complex turbulent flames. 
Depending on the fuel involved and the flow configuration, the laminar flames can take on a number of interesting geometries. 
For example, as practical combustion systems often operate in partially premixed mode,
with one or more fuel injections, a wide range of fresh gas compositions can be observed; 
and these conditions favor the appearance of edge flames, see Fig. :numref:`fig:TripleFlameIntro`. 

.. |a| image:: ./Visualization/TripleFlame_C2H4300.png
     :width: 100%

.. _fig:TripleFlameIntro:

.. table:: Normalized heat release rate (top) and temperature (bottom) contours of two-dimensional (2D) laminar lifted flames of ethylene.
     :align: center

     +-----+
     | |a| |
     +-----+

Edge flames are composed of lean and rich premixed flame wings usually surrounding a central
anchoring diffusion flame extending from a single point [PCI2007]_. Edge flames play
an important role in flame stabilization, re-ignition and propagation.
Simple fuels can exhibit up to three burning branches while diesel fuel, with a low temperature combustion mode, 
can exhibit up to 5 branches.

The goal of this tutorial is to setup a simple 2D laminar triple edge flame configuration with `PeleLM`. 
This document provides step by step instructions to properly set-up the domain and boundary conditions, 
construct an initial solution, and provides guidance on how to monitor and influence the initial transient to reach
a final steady-state solution. 
In a final Section, post-processing tools available in `PeleAnalysis` are used to extract information about 
the triple flame.

..  _sec:TUTO1::PrepStep:

Preparation step
-----------------------

Setting-up your environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As explained in section :ref:`sec:QUICKSTART`, `PeleLM` relies on a number of supporting softwares: 

- `AMREX` is a software frameworks that provides the data structure and enable massive parallelization.
- `IAMR` is a parallel, adaptive mesh refinement (AMR) code that solves the variable-density incompressible Navier-Stokes equations.
- `PelePhysics` is a repository of physics databases and implementation code. In particular, the choice of chemistry and transport models as well as associated functions and capabilities are managed in `PelePhysics`.

All of these codes have their own development cycle, and it can make the setup of a `PeleLM` run a bit tricky. i
To simplify the process, `PeleProduction` will be employed. `PeleProduction` is a collection of run folders for 
various `Pele` codes and processing. It includes git submodules for the dependent codes 
(such as `PeleLM`, `PelePhysics`, `AMReX`, etc), that can be frozen to a particular commit. 
This organizational strategy enables to manage the interactions between the various dependent repositories 
(to keep them all compatible with each other).

Numerical setup
-----------------------

Test case
^^^^^^^^^^^^^^^^^^^^^
Direct Numerical Simulations (DNS) are performed on a 2x4:math:`cm^2` 2D computational domain 
using a 64x128 base grid and 4 levels of refinement. The refinement ratio is set to 2, corresponding to 
a minimum grid size inside the reaction layer just below 20:math:`μm`. 
The maximum box size is fixed at 32, and the base grid is composed of 8 boxes, 
as shown in Fig :numref:`fig:NumSetup`.
The edge flame is stabilized against an incoming mixing layer with a uniform velocity profile. The mixing
layer is prescribed using an hyperbolic tangent of mixture fraction :math:`z` between 0 and 1:

.. math::

    z(x) = 0.5(1 + tanh(\frac{x - 0.6(x_{hi} + x_{lo})}{0.05(x_{hi} - x_{lo})}))

where :math:`z` is based on the classical elemental composition:

.. math::

    z = ...

Fig :numref:`fig:NumSetup` 

.. |b| image:: ./Visualization/SetupSketch.png
     :width: 50%

.. _fig:NumSetup:

.. table:: Sketch of the computational domain with level 0 box decomposition and input mixture fraction profile.
     :align: center

     +-----+
     | |b| |
     +-----+

Numerical scheme
^^^^^^^^^^^^^^^^^^^^^

Boundary conditions
^^^^^^^^^^^^^^^^^^^^^

Initialization and transient phase
----------------------------------

Initial solution
^^^^^^^^^^^^^^^^^^^^^

Initial transient
^^^^^^^^^^^^^^^^^^^^^

Refinement of the computation
-----------------------------

Analysis
-----------------------

.. [PCI2007] S. Chung, Stabilization, propagation and instability of tribrachial triple flames, Proceedings of the Combustion Institute 31 (2007) 877–892

