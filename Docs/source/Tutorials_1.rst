.. highlight:: rst

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

Setting-up your environment
---------------------------

PeleProduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As explained in section :ref:`sec:QUICKSTART`, `PeleLM` relies on a number of supporting softwares: 

- `AMREX` is a software frameworks that provides the data structure and enable massive parallelization.
- `IAMR` is a parallel, adaptive mesh refinement (AMR) code that solves the variable-density incompressible Navier-Stokes equations.
- `PelePhysics` is a repository of physics databases and implementation code. In particular, the choice of chemistry and transport models as well as associated functions and capabilities are managed in `PelePhysics`.

All of these codes have their own development cycle, and it can make the setup of a `PeleLM` run a bit tricky.
To simplify the process, `PeleProduction <https://github.com/AMReX-Combustion/PeleProduction>`_ will be employed. `PeleProduction` is a collection of run folders for various `Pele` codes and processing. It includes git submodules for the dependent codes 
(such as `PeleLM`, `PelePhysics`, `AMReX`, etc), that can be frozen to a particular commit. 
This organizational strategy enables to manage the interactions between the various dependent repositories 
(to keep them all compatible with each other).

Step by step instructions 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First, make sure that "git" is installed on your machine---we recommend version 1.7.x or higher.
Then, follow these few steps to setup your run environment:

1. Download the `PeleProduction` repository and : ::

    git clone https://github.com/AMReX-Combustion/PeleProduction.git 

    cd PeleProduction 

2. Switch to the TripleFlame branch : ::

    git checkout -b Tutorial_TripleFlame remotes/origin/Tutorial_TripleFlame 

3. The first time you do this, you will need to tell git that there are submodules. Git will look at the ``.gitmodules`` file in this branch and use that : ::

    git submodule int 

4. Finally, get the correct commits of the sub-repos set up for this branch: ::

    git submodule update

You are now ready to build the ``TripleFlame`` case associated with this branch. To do so: ::

   cd PeleLMruns/TripleFlame

And follow the next steps !


Numerical setup
-----------------------

In this section we review the content of the various input files for the Triple Flame test case. To get additional information about the keywords discussed, the user is referred to section :ref:`sec:control`.

Test case and boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Direct Numerical Simulations (DNS) are performed on a 2x4 :math:`cm^2` 2D computational domain 
using a 64x128 base grid and up to 4 levels of refinement (although we will start with a lower number of levels). 
The refinement ratio between each level is set to 2. With 4 levels, this means that the minimum grid size inside the reaction layer will be just below 20 :math:`μm`. 
The maximum box size is fixed at 32, and the base (level 0) grid is composed of 8 boxes, 
as shown in Fig :numref:`fig:NumSetup`.

Symmetric boundary conditions are used in the transversal (:math:`x`) direction, while ``Inflow`` (dirichlet) and ``Outflow`` (neumann) boundary conditions are used in the main flow direction (:math:`y`). The flow goes from the bottom to the top of the domain. The specificities of the ``Inflow`` boundary condition are explained in subsection :ref:`sec:TUTO1::InflowSpec`

.. |b| image:: ./Visualization/SetupSketch.png
     :width: 100%

.. _fig:NumSetup:

.. table:: Sketch of the computational domain with level 0 box decomposition (left) and input mixture fraction profile (right).
     :align: center

     +-----+
     | |b| |
     +-----+

The geometrical problem is specified in the first block of the ``inputs.2d-regt``: ::

   #----------------------DOMAIN DEFINITION------------------------                                                                        
   geometry.is_periodic = 0 0       # Periodicity in each direction: 0 => no, 1 => yes
   geometry.coord_sys   = 0         # 0 => cart, 1 => RZ
   geometry.prob_lo     = 0. 0.     # x_lo y_lo
   geometry.prob_hi     = 0.02 0.04 # x_hi y_hi

The second block determines the boundary conditions. Refer to Fig :numref:`fig:NumSetup`: ::

   # >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
   # Interior, Inflow, Outflow, Symmetry,
   # SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm
   peleLM.lo_bc = Symmetry  Inflow
   peleLM.hi_bc = Symmetry  Outflow

The number of levels, refinement ratio, maximium grid size as well as other related refinement parameters are set under the third block  : ::

   #-------------------------AMR CONTROL----------------------------
   amr.n_cell          = 64 128     # Level 0 number of cells in each direction
   amr.v               = 1          # amr verbosity level
   amr.max_level       = 1          # maximum level number allowed
   amr.ref_ratio       = 2 2 2 2    # refinement ratio
   amr.regrid_int      = 2          # how often to regrid
   amr.n_error_buf     = 1 1 1 2    # number of buffer cells in error est
   amr.grid_eff        = 0.9        # what constitutes an efficient grid
   amr.grid_eff        = 0.7        # what constitutes an efficient grid
   amr.blocking_factor = 16         # block factor in grid generation
   amr.max_grid_size   = 32         # maximum box size


..  _sec:TUTO1::InflowSpec:

Inflow specification
^^^^^^^^^^^^^^^^^^^^^

The edge flame is stabilized against an incoming mixing layer with a uniform velocity profile. The mixing
layer is prescribed using an hyperbolic tangent of mixture fraction :math:`z` between 0 and 1, as can be seen in Fig :numref:`fig:NumSetup`:

.. math::

    z(x) = 0.5 \Big(1 + tanh \Big( \frac{x - 0.6(x_{hi} + x_{lo})}{0.05(x_{hi} - x_{lo})} \Big) \Big)

where :math:`z` is based on the classical elemental composition [CF1990]_:

.. math::

    z =  \frac{\beta - \beta_{ox}}{\beta_{fu} - \beta_{ox}}}
    
where :math:`\beta` is Bilger's coupling function, and subscript :math:`ox` and :math:`fu` correspond to oxidizer and fuel streams respectively.

Specifying dirichlet ``Inflow`` conditions in `PeleLM` can seem daunting at first. But it is actually a very flexible process. We walk the user through the details of it for the Triple Flame case just described. The files involved are:

- ``probdata.F90``, where the input variables (``*_in`` and ``*_bc``) are defined (they are part of the ``probdata_module`` module)
- ``Prob_nd.F90``, where these input variables are filled -- only once at the beginning of the program:

  - ``*_in`` are filled based on what is defined under the ``&fortin`` namelist in the ``probin.2d.test``
  - ``*_bc`` are filled in the routine ``setupbc()``. They are usually a function of the ``*_in`` variables. In our case, a simple copy for the velocity and temperature.
  
- ``user_defined_fcts_nd.F90``, where the ``*_bc`` variables are used in the routine ``bcfunction()`` which is called every time step to prescribe the dirichlet inflow conditions.

Note that in our specific case, we compute the input value of the mass fractions (Y) *directly* in ``bcfunction()``, using the ``probdata_module`` variable ``H2_enrich``. We do not need any additional information, because we hard coded the hyperbolic tangent profile of :math:`z` (see previous formula) and there is a direct relation with the mass fraction profiles. The interested reader can look at the function ``set_Y_from_Ksi`` and ``set_Y_from_Phi`` in ``user_defined_fcts_nd.F90``.


Initial solution
^^^^^^^^^^^^^^^^^^^^^

Numerical scheme
^^^^^^^^^^^^^^^^^^^^^

The ``NUMERICS CONTROL`` block can be modified by the user to increase the number of SDC iterations. Note that there are many other parameters controlling the numerical algorithm that the advanced user can tweak, but we will not talk about them in the present Tutorial. The interested user can refer to section :ref:`sec:control:pelelm`.

Initial transient phase
----------------------------------

Build the executable
^^^^^^^^^^^^^^^^^^^^^

Refinement of the computation
-----------------------------

Analysis
-----------------------

.. [PCI2007] S. Chung, Stabilization, propagation and instability of tribrachial triple flames, Proceedings of the Combustion Institute 31 (2007) 877–892
.. [CF1990] R. Bilger, S. Starner, R. Kee, On reduced mechanisms for methane-air combustion in nonpremixed flames, Combustion and Flames 80 (1990) 135-149

