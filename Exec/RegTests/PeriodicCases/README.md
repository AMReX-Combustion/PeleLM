Periodic cases - README
================================

Description
-----------

This folder contains a set of fully periodic cases used to test PeleLM advection and diffusion algorithm. A specific case can be selected between the three following options:
 - `ConvectedVortex` : is a classical CFD test case consisting in convecting an isentropic vortex across a periodic box in absence of molecular viscosity (Euler equation). The analytical solution to this problem after `n` rotations is given by the initial solution, and comparing the final solution to the initial one enable to evalute the error on the advection term in the momemtum equation.
 - `ConvectedGaussian` : is a counterpart of the "ConvectedVortex" case for scalar advection. In this case, the solution is initialized with a Gaussian perturbation of species or temperature from a background state and, in the absence of molecular diffusion, convecting this pertubation `n` times in a periodic domain should left the perturbation unchanged. This allows to evaluate the error on the scalar advection term in PeleLM.
 - `DiffusedGaussian` : in contrast with the previous two cases, this problem is initially velocity free and the same species or temperature perturbation used in "ConvectedGaussian" is evolve in time under the effect of molecular diffusion. At the present time, this case is not rigged up to provide an analytical solution (but it could be ...) but it allows to test basic PeleLM bulding block in presence of a divergence of U /= 0 constraint.


Convergence testing
----------------------
The "ConvectedVortex" and "ConvectedGaussian" cases allow to measure the convergence rate of error as the grid is refined. To this end, two scripts are provided:
 - `multiRuns.py` starts a batch of simulation at increasing resolution.
 - `pprocConvOrder.py` will run the fcompare AMReX tool to extract the error of each simulation and draw convergence plots.
Details on the tunable parameters of these are provided in the header part of each one. For instance, to use these scripts to evaluate the convergence order of the momentum equation advection term, use the following:

`./multiRuns.py --test_name Convergence_COVO --input_file inputs.2d_CoVo_Conv_pos45d`

`./pprocConvOrder.py <path_to_fcompare.exe> --test_name Convergence_COVO`
Note that you will need to provide a path the AMReX fcompare tool executable in order to extract the error. 
