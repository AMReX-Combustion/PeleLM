---
title: 'PeleLM: an AMR Low Mach Number Reactive Flow Simulation Code'
tags:
  - C++
  - adaptive mesh refinement
  - hydrodynamics
  - combustion
  - reactions
  - CFD
  - low-Mach number
authors:
  - name: Lucas Esclapez
    orcid: 0000-0002-2438-7292
    affiliation: 1, 2
  - name: Marc Day
    orcid: 0000-0002-1711-3963
    affiliation: 1, 2
  - name: Emmanuel Motheau
    orcid: 0000-0003-1968-1611
    affiliation: 1
  - name: Candace Gilet
  - name: Andy Nonaka
    affiliation: 1
  - name: Anne Felden
    affiliation: 1
  - name: Weiqun Zhang
    affiliation: 1
  - name: John Bell
    orcid: 0000-0002-5749-334X
    affiliation: 1
  - name: Valentina Ricciuti
  - name: Jon Rood
    affiliation: 2
  - name: Brian Friesen
    affiliation: 3
  - name: Deepak Dalakoti
    affiliation: 4
  - name: Armin Wehrfritz
    affiliation: 4
  - name: Nicolas Wimer
    affiliation: 2
  - name: Ray Grout
    affiliation: 2
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory, USA
    index: 1
  - name: High Performance Algorithms and Complex Fluids, National Renewable Energy Laboratory, USA
    index: 2
  - name: National Energy Research Scientific Computing Center, Lawrence Berkeley National Laboratory, USA
    index: 3
  - name: University of New South Wales, Australia
    index: 4
date: 01 February 2022
bibliography: paper.bib
---

# Summary

PeleLM evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). 
The code depends upon the AMReX [@AMReX] library to provide the underlying data structures and tools to manage 
and operate on them across massively parallel computing architectures. Together with its compressible flow counterpart 
PeleC, the thermo-chemistry library PelePhysics and the multi-physics library PeleMP, it forms the Pele suite of 
open-source simulation codes.

PeleLM uses a finite volume appraoch to solve for the multi-species reacting Navier-Stokes equations in 
their low Mach number limit [@Day:2000], where the characteristic fluid velocity is small compared to the sound speed, 
and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly, 
acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time 
step based on an advective CFL condition.
This low Mach number limit mathematically translates into a constraint on the divergence of the velocity field [@Majda:1986]. The 
momemtum equation is then solved for using a predictor/corrector method initially developed for incompressible flows [@Almgren1998]
and later extended to reactive, variable-density flows [@Pember:1998]. A constrained time-centered 
face-centered velocity field is first obtained using a second-order Godunov procedure followed by a discrete projection 
using time-centered source terms. At this point, the velocity field is used to advance the species and energy equations. 
Then a semi-implicit Crank-Nicholson method is used to obtain a cell-centered velocity upon which the low Mach contrain
is subsequently enforced using an approximate projection.
The species and energy conservation equations include advection, diffusion and reactions processes as well as
any external forcing (gravity, ...). The advection term is discretized using a second-order Godunov scheme, the 
diffusion term is obtained using a semi-implicit Crank-Nicholson scheme while the much stiffer reaction term
is obtained using a fully implicit Backward Differentiation Formulas scheme (specifically, the CVODE integrator 
[@cvodeDocumentation] of the Sundials suite [@SUNDIALS]).
Because of the wide scale separation between the slow hydronamics and the fast reations, this operator splitting 
appraoch can lead to significant splitting error if not carefully integrated on time. PeleLM uses an iterative 
Spectral Defered Correction (SDC) scheme [@Nonaka12,@Nonaka18] to ensure a tight coupling of all the processes while
iteratively enforcing the low Mach number contrain.
On an AMR hierarchy, each level evolves at its own time scale using a sub-cycling approach where each next finer level
performs two steps to reach the same physical time at which a synchronization operation enforce conservation across levels.

In addition, PeleLM uses an Embedded Boundary (EB) approach to represent complex geometries: an arbitrary surface can 
be intersected with the Cartesian matrix of uniform cells, and the numerical stencils are modified near cells that are cut 
by the EB. Redistribution schemes [@Berger:2021] are then used for the explicit advection and diffusion updates in order to alleviate the 
constrain associated with small cut cells. Through its dependency to the multi-physics library PeleMP, PeleLM also inherits
the ability to include Lagrangian sprays as well as soot and radiation models.

PeleLM is written in C++ and is built upon the AMReX [@AMReX] library from which it inherits its parallel paradigm.
It uses a MPI+X approach where MPI is used to distribute AMR grids across nodes and each grid can be further divided into 
logical tiles spread across threads using OpenMP for multi-core CPU machines, or spread across GPU threads using CUDA/HIP 
on GPU-based machines.

# Statement of Need

While there exist several reactive flow Direct Numerical Simulation codes, PeleLM presents a unique set of features. 
From its creation, under the name LMC in the early 2000, the motivation was to combine an AMR approach with a low Mach number 
formulation to achieve high performances from a small desktop stations to the world largest supercomputer, and to this day
it remains the only publicly available code to offers these features.

PeleLM is predominantly used to study the fine scale interactions between turbulence and chemical reactions occuring in many
combustion applications. A better understanding of these interactions is the basis for developing accurate modeling approaches
that can be used to design the next generation of low-emission combustion devices. To this end, PeleLM has been used to study
diesel jet flame [@Dalakoti:2020], premixed dodecane/air flame chemical pathways [@Dasgupta:2019] or the effect of increasingly 
intense turbulence [@Aspden:2019]. Notable applications also includes large scale helium plumes enabled by AMR [@Wimer:2021].

Although PeleLM can be used on small desktop stations, our main focus is enabling massively parallel simulations at scale on 
high-performance computer architectures to tackle the challenging requirements of fundamental and applied combustion research.

# Acknowledgments

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE 
organizations -- the Office of Science and the National Nuclear Security Administration -- responsible for the planning and 
preparation of a capable exascale ecosystem -- including software, applications, hardware, advanced system engineering, and 
early testbed platforms -- to support the nation's exascale computing imperative.

# References
