---
title: 'PeleLM: an AMR Low Mach Number Reactive Flow Simulation Code'
tags:
  - C++
  - hydrodynamics
  - combustion
  - reactions
  - CFD
  - low-Mach number
authors:
  - name: John Bell
    orcid: 0000-0002-5749-334X
    affiliation: 1
  - name: Marc Day
    orcid: 0000-0002-1711-3963
    affiliation: 1
  - name: Lucas Esclapez
    orcid: 0000-0002-2438-7292
    affiliation: 1
  - name: Anne Felden
    affiliation: 1
  - name: Emmanuel Motheau
    orcid: 0000-0003-1968-1611
    affiliation: 1
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
date: 16 January 2021
bibliography: paper.bib
---

# Summary

PeleLM evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). 
The code depends upon the AMReX [@AMReX] library to provide the underlying data structures, and tools to manage 
and operate on them across massively parallel computing architectures. Together with its compressible flow counterpart 
PeleC, the thermo-chemistry library PelePhysics and the multi-physics library PeleMP, it forms the Pele suite of 
open-source simulation codes.

PeleLM uses a finite volume appraoch to solve for the multi-species reacting Navier-Stokes equations in 
their low Mach number limit [@LMC], where the characteristic fluid velocity is small compared to the sound speed, 
and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly, 
acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time 
step based on an advective CFL condition.
The low Mach number limit mathematically translates into a constraint on the divergence of the velocity field [@Majda]. The 
momemtum equation is then solved for using a predictor/corrector method [@Almgren]. A constrained time-centered 
face-centered velocity field is first obtained using a second-order Godunov procedure followed by a discrete projection 
using time-centered source terms. At this point, the velocity field is used to advance the species and energy equations. 
Then a semi-implicit Crank-Nicholson method is used to obtain a cell-centered velocity upon which the low Mach contrain
is subsequently enforced using an approximate projection.
The species and energy conservation equations include advection, diffusion and reactions processes as well as
any external forcing (gravity, ...). The advection term is discretized using a second-order Godunov scheme, the 
diffusion term is obtained using a semi-implicit Crank-Nicholson scheme while the much stiffer reaction term
is obtained using a fully implicit Backward Differentiation Formulas scheme (such as implemented in Sundials [@SUNDIALS]).
Because of the wide scale separation between the slow hydronamics and the fast reations, this operator splitting 
appraoch can lead to significant splitting error if not carefully integrated on time. PeleLM uses an iterative 
Spectral Defered Correction (SDC) scheme [@Nonaka12,@Nonaka18] to ensure a tight coupling of all the processes while
iteratively enforcing the low Mach number contrain.
On an AMR hierarchy, each level evolves at its own time scale  

In addition, PeleLM uses an Embedded Boundary (EB) approach to represent complex geometries: an arbitrary surface can 
be intersected with the Cartesian matrix of uniform cells, and the numerical stencils are modified near cells that are cut 
by the EB.

PeleLM is written in C++ and is built upon the AMReX [@AMReX] library from which it inherits its parallel paradigm.
MPI is used to distribute AMR grids across nodes. Each grid can be further divided into logical tiles spread
across threads using OpenMP for multi-core CPU machines, or spread across GPU threads using CUDA/HIP on GPU-based machines.

# Statement of Need

While there exist several reactive flow Direct Numerical Simulation codes, PeleLM present a unique set of features. 
From its creation, under the name LMC in the early 2000, the motivation was to combine an AMR approach with a low Mach number 
formulation to achieve high performances from a small desktop stations to the world largest supercomputer.

# Acknowledgments

# References
