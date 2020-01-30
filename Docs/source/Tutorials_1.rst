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

Practical combustion systems often operate in partially premixed mode,
with one or more fuel injections, resulting in a wide range of fresh gas compositions.  
These conditions favor the appearance edge flames, which are
composed of a lean and a rich premixed flame wing surrounding a central
anchoring diffusion flame extending from a single point [PCI2007]_. Edge flames play
an important role in flame stabilization, re-ignition and propagation.
Simple fuels can exhibit up to three burning branches while diesel fuel, with a low temperature combustion mode, 
can exhibit up to 5 branches.

The goal of this tutorial is to setup a simple 2D laminar triple flame configuration in `PeleLM`. 
This document provides step by step instructions to properly set-up the domain and boundary conditions, 
construct an initial solution, and provides guidance on how to monitor and influence the initial transient to reach
a final steady-state solution. 
In a final Section, post-processing tools available in PeleAnalysis are used to extract information about 
the triple flame.

..  _sec:TUTO1::PrepStep:

Preparation step
-----------------------

Setting-up your environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

    f(z) = ...

where :math:`z` is based on the classical elemental composition:
.. math::

    z = ...

Fig :numref:`fig:NumSetup` 

.. |a| image:: ./Visualization/tmp.png
     :width: 100%

.. _fig:NumSetup:

.. table:: Sketch of the computational domain with intial box decomposition.
     :align: center

     +-----+
     | |a| |
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

