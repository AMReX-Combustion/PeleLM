.. highlight:: rst

.. _sec:model:

The `PeleLM` Model
==================

In this section, we present the actual model that is evolved numerically by `PeleLM`, and the numerical algorithms
to do it.  There are many control parameters to customize the solution strategy and process, and in order to actually
set up and run specific problems with `PeleLM`, the user must specific the chemical model, and provide routines
that implement initial and boundary data and refinement criteria for the adaptive mesh refinement.  We discuss
setup and control of `PeleLM` in later sections.

Overview of `PeleLM`
--------------------

`PeleLM` evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). The code depends upon the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library to provide the underlying data structures, and tools to manage and operate on them across massively parallel computing architectures. `PeleLM` also borrows heavily from the source code and algorithmic infrastructure of the `IAMR <https://github.com/AMReX-Codes/IAMR>`_. `IAMR` implements an AMR integration for the variable-density incompressible Navier-Stokes equations. `PeleLM` extends `IAMR` to include complex coupled models for generalized thermodynamic relationships, multi-species transport and chemical reactions.  The core algorithms in `PeleLM` (and `IAMR`) are described in the following papers:

* *A conservative, thermodynamically consistent numerical approach for low Mach number combustion. I. Single-level integration*, A. Nonaka, J. B. Bell, and M. S. Day, *Combust. Theor. Model.*, **22** (1) 156-184 (2018)

* *A Deferred Correction Coupling Strategy for Low Mach Number Flow with Complex Chemistry*, A. Nonaka, J. B. Bell, M. S. Day, C. Gilet, A. S. Almgren, and M. L. Minion, *Combust. Theory and Model*, **16** (6) 1053-1088 (2012)

* *Numerical Simulation of Laminar Reacting Flows with Complex Chemistry*, M. S. Day and J. B. Bell, *Combust. Theory Model* **4** (4) 535-556 (2000)

* *An Adaptive Projection Method for Unsteady, Low-Mach Number Combustion*, R. B. Pember, L. H. Howell, J. B. Bell, P. Colella, W. Y. Crutchfield, W. A. Fiveland, and J. P. Jessee, *Comb. Sci. Tech.*, **140** 123-168 (1998)

* *A Conservative Adaptive Projection Method for the Variable Density Incompressible Navier-Stokes Equations,* A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, and M. L. Welcome, *J. Comp. Phys.*, **142** 1-46 (1998)

The low Mach number flow equations
----------------------------------

`PeleLM` solves the reacting Navier-Stokes flow equations in the *low Mach number* regime, where the characteristic fluid velocity is small compared to the sound speed, and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly, acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time step based on an advective CFL condition, and this leads to an increase in the allowable time step of order :math:`1/M` over an explicit, fully compressible method (:math:`M` is the Mach number).  In this mathematical framework, the total pressure is decomposed into the sum of a spatially constant (ambient) thermodynamic pressure :math:`P_0` and a perturbational pressure, :math:`\pi({\vec x})` that drives the flow.  Under suitable conditions, :math:`\pi/P_0 = \mathcal{O} (M^2)`. 

The set of conservation equations specialized to the low Mach number regime is a system of PDEs with advection, diffusion and reaction (ADR) processes that are constrained to evolve on the manifold of a spatially constant :math:`P_0`:

.. math::

    &\frac{\partial (\rho \boldsymbol{u})}{\partial t} + 
    \nabla \cdot \left(\rho  \boldsymbol{u} \boldsymbol{u} + \tau \right)
    = -\nabla \pi + \rho \boldsymbol{F},\\
    &\frac{\partial (\rho Y_m)}{\partial t} +
    \nabla \cdot \left( \rho Y_m \boldsymbol{u}
    + \boldsymbol{\mathcal{F}}_{m} \right)
    = \rho \dot{\omega}_m,\\
    &\frac{ \partial (\rho h)}{ \partial t} +
    \nabla \cdot \left( \rho h \boldsymbol{u}
    + \boldsymbol{\mathcal{Q}} \right) = 0 ,

where :math:`\rho` is the density, :math:`\boldsymbol{u}` is the velocity, :math:`h` is the mass-weighted enthalpy, :math:`T` is temperature and :math:`Y_m` is the mass fraction of species :math:`m`. :math:`\dot{\omega}_m` is the molar production rate for species :math:`m`, the modeling of which will be described later in this section. :math:`\tau` is the stress tensor, :math:`\boldsymbol{\mathcal{Q}}` is the heat flux and :math:`\boldsymbol{\mathcal{F}}_m` are the species diffusion fluxes. These transport fluxes require the evaluation of transport coefficients (e.g., the viscosity :math:`\mu`, the conductivity :math:`\lambda` and the diffusivity matrix :math:`D`) which are computed using the library EGLIB, as will be described in more depth in the diffusion section. The momentum source, :math:`\boldsymbol{F}`, is an external forcing term.  For example, we have used :math:`\boldsymbol{F}` to implement a long-wavelength time-dependent force to establish and maintain quasi-stationary turbulence.

These evolution equations are supplemented by an equation of state for the thermodynamic pressure.  For example, the ideal gas law,

.. math::

    P_0(\rho,Y_m,T)=\frac{\rho \mathcal{R} T}{W}=\rho \mathcal{R} T
    \sum_m \frac{Y_m}{W_m}

can be used, although `PeleLM` will soon support other more general expressions, such as Soave-Redlich-Kwong.  In the above, :math:`W_m` and :math:`W` are the species :math:`m`, and mean molecular weights, respectively.  To close the system we also require a relationship between enthalpy, species and temperature.  We adopt the definition used in the CHEMKIN standard,

.. math::

    h=\sum_m Y_m h_m(T)

where :math:`h_m` is the species :math:`m` enthalpy.  Note that expressions for :math:`h_m(T)` see <section on thermo properties> incorporate the heat of formation for each species.


Neither species diffusion nor reactions redistribute the total mass, hence we have :math:`\sum_m \boldsymbol{\mathcal{F}}_m = 0` and :math:`\sum_m \dot{\omega}_m = 0`. Thus, summing the species equations and using the definition :math:`\sum_m Y_m = 1` we obtain the continuity equation:

.. math::

    \frac{\partial \rho}{\partial t} + \nabla \cdot \rho \boldsymbol{u} = 0

This, together with the conservation equations form a differential-algebraic equation (DAE) system that describes an evolution subject to a constraint.  A standard approach to attacking such a system computationally is to differentiate the constraint until it can be recast as an initial value problem.  Following this procedure, we set the thermodynamic pressure constant in the frame of the fluid,

.. math::

    \frac{DP_0}{Dt} = 0

and observe that if the initial conditions satisfy the constraint, an evolution satisfying the above will continue to satisfy the constraint over all time.  Expanding this expression via the chain rule and continuity:

.. math::

    \nabla \cdot \boldsymbol{u} = \frac{1}{T}\frac{DT}{Dt}
    + W \sum_m \frac{1}{W_m} \frac{DY_m}{Dt} = S

The constraint here take the form of a condition on the divergence of the flow.  Note that the actual expressions to use here will depend upon the chosen models for evaluating the transport fluxes.


Transport Fluxes
^^^^^^^^^^^^^^^^

Expressions for the transport fluxes appearing above can be approximated in the Enskog-Chapman expansion as:

.. math::

    &&\boldsymbol{\mathcal{F}}_{m} = \rho Y_m \boldsymbol{V_m} \\
    &&\tau_{i,j} = - \Big(\kappa - \frac{2}{3} \mu \Big) \delta_{i,j}
    \frac{\partial {u_k}}{\partial x_k}
    - \mu \Big(\frac{\partial u_i}{\partial x_j} +
    \frac{\partial u_j}{\partial x_i}\Big) \\
    &&\boldsymbol{\mathcal{Q}} =  \sum_m h_m \boldsymbol{\mathcal{F}}_{m}
    - \lambda' \nabla T - P_0 \sum_m \theta_m \boldsymbol{d_m}

where :math:`\mu` is the shear viscosity, :math:`\kappa` is the bulk viscosity, and :math:`\lambda'` is the partial thermal conductivity. In the *full matrix diffusion model*, the vector of :math:`m` species diffusion velocities, :math:`\boldsymbol{V_m}`, is given by:

.. math::

    \boldsymbol{V_m} = - \sum_j  {D}_{m,j} \boldsymbol{d_j}
    - \theta_m \nabla ln(T)

where :math:`{D}_{m,j}` is the diffusion matrix, and :math:`\boldsymbol{\theta}` are thermal diffusion coefficients associated with the Soret (mass concentration flux due to an energy gradient) and Dufour (the energy flux due to a mass concentration gradient) effects. The :math:`m` species transport driving force due to composition gradients, :math:`\boldsymbol{d_m}`, is given by:

.. math::

    \boldsymbol{d_m} = \nabla X_m + (X_m -Y_m) \frac{\nabla P_0}{P_0}

Alternatively (as in the transport library, EGLIB) the thermal diffusion *ratios* :math:`\boldsymbol{\chi}` may be preferred and the diffusion velocities and energy flux recast as:

.. math::

    \boldsymbol{V_m} = - \sum_j  {D}_{m,j} ( \boldsymbol{d_j}
    + \chi_j \nabla ln(T))\\
    \boldsymbol{\mathcal{Q}} =  \sum_m h_m \boldsymbol{\mathcal{F}}_{m}
    - \lambda \nabla T + P_0 \sum_m \chi_m \boldsymbol{V_m}

where  :math:`{D} \boldsymbol{\chi} = \boldsymbol{\theta}`.

As can be seen, the expression for these fluxes relies upon several transport coefficients that need to be evaluated. However, in the present framework several effects are neglected, thus simplifying the fluxes evaluation.

The `PeleLM` Equation Set
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _sec:model:EqSets:

The full diffusion model couples together the advance of all thermodynamics fields, including a dense matrix transport operator that is cumbersome to deal with computationally, while also being generally viewed as an overkill for most practical combustion applications -- particularly those involving turbulent fluid dynamics.  For `PeleLM`, we make the following simplifying assumptions:

1. The bulk viscosity, :math:`\kappa`, is usually negligible, compared to the shear viscosity,

2. The low Mach limit implies that there are no spatial gradients in the thermodynamic pressure,

3. The *mixture-averaged* diffusion model is assumed,

4. Dufour and Soret effects are negligible

With these assumptions, the conservation equations take the following form:

.. math::

    &&\frac{\partial (\rho \boldsymbol{u})}{\partial t} +
    \nabla \cdot \left(\rho  \boldsymbol{u} \boldsymbol{u} + \tau \right)
    = -\nabla \pi + \rho \boldsymbol{F}, \\
    &&\frac{\partial (\rho Y_m)}{\partial t} +
    \nabla \cdot \left( \rho Y_m \boldsymbol{u} + \boldsymbol{\mathcal{F}}_{m} \right)
    = \rho \dot{\omega}_m \\
    &&\frac{ \partial (\rho h)}{ \partial t} +
    \nabla \cdot \left( \rho h \boldsymbol{u} + \boldsymbol{\mathcal{Q}} \right) = 0,

with

.. math::

    &&\boldsymbol{\mathcal{F}}_{m} = \rho Y_m \boldsymbol{V_m} = - \rho \sum_k \widetilde{D_{m,k}} \nabla X_m \\
    &&\tau_{i,j} = \frac{2}{3} \mu \delta_{i,j} \frac{\partial {u_k}}{\partial x_k} - \mu \Big(
    \frac{\partial  u_i}{\partial x_j} + \frac{\partial  u_j}{\partial x_i}\Big) \\
    &&\boldsymbol{\mathcal{Q}} =  \sum_m h_m \boldsymbol{\mathcal{F}}_{m}  - \lambda \nabla T

where :math:`\boldsymbol{d_m} = \nabla X_m` and :math:`\widetilde{D_{m,k}} = Y_m D_{m,k}`.
Using these expressions, we can write an equation for :math:`T` that is needed in order to evaluate the right-hand side of the divergence constraint:

.. math::

    \rho C_p \frac{DT}{Dt} = \nabla \cdot \lambda \nabla T + \sum_m \Big( h_m \nabla \cdot \boldsymbol{\mathcal{F}}_{m} - \nabla \cdot h_m \boldsymbol{\mathcal{F}}_{m} - h_m \rho \dot\omega_m \Big)

where :math:`C_p = \partial h/\partial T` is the specific heat of the mixture at constant pressure. For an ideal gas, the constraint then becomes:

.. math::

    \nabla \cdot \boldsymbol{u} &=&\frac{1}{\rho C_p T}\Big[ \nabla \cdot \lambda \nabla T
    + \sum_m \Big( h_m \nabla \cdot \boldsymbol{\mathcal{F}}_{m}
    - \nabla \cdot h_m \boldsymbol{\mathcal{F}}_{m}\Big) \Big] \\
    &&- \frac{W}{\rho} \sum_m \frac{1}{W_m} \nabla \cdot \boldsymbol{\mathcal{F}}_{m}
    + \sum_m \Big( \frac{W}{W_m} -\frac{h_m(T)}{c_{p} T} \Big)\dot{\omega}_m

The mixture-averaged transport coefficients discussed above (:math:`\mu`, :math:`\lambda` and :math:`D_{m,mix}`) can be evaluated from transport properties of the pure species. We follow the treatment used in the EGLib library, based on the theory/approximations developed by Ern and Givangigli (however, `PeleLM` uses a recoded version of these routines that are thread safe and vectorize well on suitable processors).


The following choices are currently implemented in `PeleLM`

* The viscosity, :math:`\mu`, is estimated based on one step of the conjugate gradient method, using temperature dependent ratios of collisions integrals (EGZE3).

* The conductivity, :math:`\lambda`, is based on an empirical mixture formula (EGZL1):

.. math::

    \lambda = \mathcal{A}_{0.25}

with

.. math::

    \mathcal{A}_{\alpha}= \Big( \sum_m X_m (\lambda_m)^{\alpha} \Big)^{1/\alpha}

* The flux diffusion matrix is approximated using the diagonal of the flux diffusion vector :math:`\rho \widetilde{ \Upsilon}`, where:

.. math::

    \rho \widetilde{ \Upsilon}_m =  \rho \frac{W_m}{W}D_{m,mix}, \;\;\;\mbox{where} \;\;\;
    D_{m,mix} = \frac{1-Y_m}{ \sum_{j \neq m} X_j / \mathcal{D}_{m,j}}

and the :math:`\mathcal{D}_{m,j}` are the binary diffusion coefficients of the pair (m,j). This leads to a mixture-averaged approximation that is similar to that of Hirschfelder-Curtiss (EGZVR1):

.. math::

    \rho Y_m \boldsymbol{V_m} = - \rho D_{m,mix} \frac{W_m}{W} \nabla X_m 

Note that with these definitions, there is no guarantee that :math:`\sum \boldsymbol{\mathcal{F}}_{m} = 0`, as required for mass conservation. An arbitrary *correction flux,* consistent with the mixture-averaged diffusion approximation, is added in `PeleLM` to enforce conservation.

The pure species and mixture transport properties are evaluated with (thread-safe, vectorized) EGLib functions, which require as input polynomial fits of the logarithm of each quantity versus the logarithm of the temperature.

.. math::

    ln(q_m) = \sum_{n=1}^4 a_{q,m,n} ln(T)^{(n-1)} 

:math:`q_m` represents :math:`\eta_m`, :math:`\lambda_m` or :math:`D_{m,j}`. These fits are generated as part of a preprocessing step managed by the tool `FUEGO` based on the formula (and input data) discussed above. The role of `FUEGO` is to preprocess the model parameters for transport as well as chemical kinetics and thermodynamics.


Chemical kinetics and the reaction source term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Chemistry in combustion systems involves the :math:`N_s` species interacting through a set of :math:`M_r` elementary reaction steps, expressed as

.. math::

    \sum_{m=1}^{N_s} \nu_{m,j}'[X_m] \rightleftharpoons \sum_{m=1}^{N_s} \nu_{m,j}''[X_m],\quad for \quad j \in [1,M_r] 

where :math:`[X_m]` is the molar concentration of species :math:`m`, and :math:`\nu_{m,j}'`, :math:`\nu_{m,j}''` are the stoichiometric coefficients on the reactant and product sides of reaction :math:`j`, associated with :math:`m`. For such a system, the rate of reaction :math:`j` (:math:`R_j`) can be expressed in terms of the the forward (:math:`k_{f,j}`) and backward (:math:`k_{r,j}`) rate coefficients,

.. math::

    R_{j} = k_{f,j}\prod_{m=1}^{N_s}  [X_{m}]^{\nu_{m,j}'}-k_{r,j}\prod_{m=1}^{N_s} [X_{m}]^{\nu_{m,j}''}

The net molar production rate, :math:` \dot{\omega}_m` of species :math:`m` is obtained by
collating the rate of creation and destruction over reactions:

.. math::

    \dot{\omega}_m = \sum_{j=1}^{M_r} \nu_{m,j} R_j 

where :math:`\nu_{m,j} =\nu_{m,j}'' - \nu_{m,j}'`. Expressions for the reaction rates coefficients :math:`k_{(f,r),j}` depend on the type of reaction considered. We use the CHEMKIN modified Arrhenius reaction format:

.. math::

    k_f = AT^{\beta} exp \left( \frac{-E_a}{RT}\right)

where :math:`A` is the pre-exponential (frequency) factor, :math:`\beta` is the temperature exponent and :math:`E_a` is the activation energy. The CHEMKIN format additionally allows for a number of specializations of this format to represent pressure dependencies and third-body enhancements -- see the CHEMKIN Manual or Cantera website for additional information.

Most fundamental Arrhenius reactions are bidirectional, and typically only the forward rates are specified. In this case, the balance of forward and reverse rates are dictacted by equilibrium thermodynamics, via the equilibrium constant, :math:`K_{c,j}`.  In a low Mach system, :math:`K_{c,j}` is a function only of temperature and the thermodynamic properties of the reactants and products of reaction :math:`j`,

.. math::

    &&k_{r,j} = \frac{k_{f,j}}{K_{c,j}(T)} \;\;\; \mbox{where} \;\;\; K_{c,j}=K_{p,j} \left( \frac{P_{0}}{RT} \right)^{\sum_{k=1}^{N_s} \nu_{k,j}}\\
    &&\mbox{and} \;\;\; K_{p,j}=\exp \left( \frac{\Delta {S_j}^{0}}{R} - \frac{\Delta {H_j}^{0}}{RT} \right)

:math:`\Delta H_j` and :math:`\Delta S_j` are the change in enthalpy and entropy of the reaction :math:`j`, and :math:`P_0` is the ambient thermodynamic pressure.

Species production rates are evaluated via functions that are generated as part of a preprocessing step managed by the tool `FUEGO`.

Thermodynamic properties
^^^^^^^^^^^^^^^^^^^^^^^^

Currently, expressions for the thermodynamic properties in `PeleLM` follow those of CHEMKIN, which assume a mixture of ideal gases. Species enthalpies and entropies are thus functions of only temperature (for perfect gases, they are independent of pressure) and are given in terms of polynomial fits to the species molar heat capacities (:math:`C_{p,\cdot}`),

.. math::

    \frac{C_{p,m}(T)}{\mathcal{R}} = \sum_{k=1}^{N_s} a_{k,m}T^{k-1}

where, in the standard CHEMKIN framework (the 7-coefficients NASA format), :math:`N =5` and

.. math::

    \frac{C_{p,m}(T)}{\mathcal{R}} = a_{1,m} + a_{2,m} T + a_{3,m} T^2 + a_{4,m} T^3 + a_{5,m} T^4

Accordingly, the standard-state molar enthalpy of species :math:`m` is given by:

.. math::

    \frac{H_{m}(T)}{\mathcal{R}T} = a_{1,m} +\frac{a_{2,m}}{2} T   + \frac{a_{3,m}}{3} T^2 +  \frac{a_{4,m}}{4} T^3 + \frac{ a_{5,m}}{5} T^4 + a_{6,m}/T

Note that the standard specifies that the heat of formation for the molecule is included in this expression.
Similarly, the standard-state molar entropy is written as:

.. math::

    \frac{S_{m}(T)}{\mathcal{R}} = a_{1,m}ln(T) + {a_{2,m}} T   + \frac{a_{3,m}}{2} T^2 +  \frac{a_{4,m}}{3} T^3 + \frac{ a_{5,m}}{4} T^4 + a_{7,m}

For each species :math:`m`, in the model the user must specify the 7 :math:`k` coefficients :math:`a_{k,m}`. All other required thermodynamic properties are then determined (see, e.g., the CHEMKIN manual for additional details). Thermodynamic properties of the species, and those of the mixture, are evaluated via functions that are generated as part of a preprocessing step managed by the tool `FUEGO`.


`FUEGO` chemistry preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical model for `PeleLM` contains all the information associated with the CHEMKIN parameterization of the Arrhenius reaction set, as well as fitting coefficients for the thermodynamic relationships, and the specification of the species including data required to compute pure-species transport properties. In the combustion community, this information is communicated for each complete model --or *mechanism*, through multiple text files that conform to the CHEMKIN standards. The CHEMKIN driver code (or equivalent) can then be used to ingest the large number of parameters contained in these files and provide a set of functions for evaluating all the properties and rates required.  Earlier versions of `PeleLM` linked to the CHEMKIN codes directly (and thereby assumed that all problems consisted of a mixture of ideal gases).  However, evaluations were not very efficient because the functions stepped through generic expressions that included a large number of conditional statements and unused generality.  Direct evaluation of these complex expressions allows for a much more efficient code that optimizes well with modern compilers. This is important because an appreciable fraction of `PeleLM` runtime is spent in these functions. Performance issues notwithstanding, customized evaluators will be necessary to extend `PeleLM` to a larger class of (*real*) gas models outside the CHEMKIN standard, such as SRK, that are already part of the `PeleC` code capabilities (`PeleC` shares use of `PelePhysics` for combustion model specification).

For these reasons, `PeleLM` no longer uses CHEMKIN functions directly, but instead relies on a preprocessing tool, `FUEGO`, to generate highly efficient C code implementations of the necessary thermodynamic, transport and kinetics evaluations.  The source code generated from `FUEGO` is linked into the `PeleLM` executable, customizing each executable for a specific model at compile time.  The implementation source code files can also be linked conveniently to post-processing analysis tools. The `FUEGO` processing tool, and the functions necessary to interface the generated functions to `PeleLM` are distributed in the auxiliary code package, `PelePhysics`.  Included in the `PelePhysics` distribution is a broad set of models for the combustion of hydrogen, carbon-monoxide, methane, heptane, :math:`n`-dodecane, dimethyl ether, and others, as well as instructions for users to extend this set using `FUEGO`, based on their own CHEMKIN-compliant inputs. `PelePhysics` also provides support for simpler *gama-law* equations-of-state, and simple/constant transport properties.


The `PeleLM` temporal integration
---------------------------------

The temporal discretization in `PeleLM` combines a modified spectral deferred correction (SDC) coupling of chemistry and transport with a density-weighted approximate projection method for low Mach number flow.  The projection method enforces a constrained evolution of the velocity field, and is implemented iteratively in such a way as to ensure that the update simultaneously satisfies the  equation of state and discrete conservation of mass and total enthalpy.  A time-explicit approach is used for advection; faster diffusion and chemistry processes are treated time-implicitly, and iteratively coupled together within the deferred corrections strategy. The integration algorithm, discussed in the following sections, is second-order accurate in space and time, and is implemented in the context of a subcycled approach for a nested hierarchy of mesh levels, where each level consists of logically rectangular patches of rectangular cells.  All cells at a level have the same size in all coordinates.

Due to the complexity of the `PeleLM` algorithm, it is best presented in a number of passes.  Focusing first on the single-level advance, we begin with a general discussion of the SDC-based time step iteration, which is designed to couple together the various physics processes.  We then describe the projection steps used to enforce the constraint in the context of this iterative update.  Next, we dive a little deeper into precisely how the advance of the thermodynamic components of the state is sequenced.  There are a few crucial nuances to the formulation/sequencing of the energy advection, energy diffusion, conservative corrections to the species diffusion fluxes, and of the projection that can then be discussed in the context of overall single-level time step.  Finally, with all these aspects defined, we give an overview of the modifications necessary to support the AMR subcycling strategy.

SDC preliminaries
^^^^^^^^^^^^^^^^^

The basic idea of SDC is to write the solution of an ODE

.. math::

    &&\phi_t = F(t,\phi(t)), \qquad t\in[t^n,t^{n+1}];\\
    &&\phi(t^n) = \phi^n,

as an integral,

.. math::

    \phi(t) = \phi^n + \int_{t^n}^{t} F(\phi)~d\tau,

where we suppress explicit dependence of :math:`F` and :math:`\phi` on :math:`t` for notational simplicity.
Given an approximation :math:`\phi^{(k)}(t)` to :math:`\phi(t)`, one can then define a residual,

.. math::

    E(t,\phi^{(k)}) = \phi^n + \int_{t^n}^t F(\phi^{(k)})~d\tau - \phi^{(k)}(t).\label{eq:residual}

Defining the error as :math:`\delta^{(k)}(t) = \phi(t) - \phi^{(k)}(t)`, one can then show that

.. math::

    \delta^{(k)}(t) = \int_{t^n}^t \left[F(\phi^{(k)}+ \delta^{(k)}) - F(\phi^{(k)})\right]d\tau + E(t,\phi^{(k)}).

In SDC algorithms, the integral in the above equation
is evaluated with a higher-order quadrature rule.
By using a low-order discretization of the integral one can construct
an iterative scheme that improves the overall order of accuracy of the approximation by one per
iteration, up to the order of accuracy of the underlying quadrature rule 
used to evaluate the integral.
Specifically, if we let :math:`\phi^{(k)}` represent the current approximation and define 
:math:`\phi^{(k+1)} = \phi^{(k)} + \delta^{(k)}` to be the iterative update, 
then arrive at the update equation,

.. math::

    \phi^{(k+1)}(t) = \phi^n + \int_{t^n}^t \left[F(\phi^{(k+1)}) - F(\phi^{(k)})\right]d\tau +
    \int_{t^n}^t F(\phi^{(k)})~d\tau,\label{eq:update}

where a low-order discretization (e.g., forward or backward Euler) is used for the first integral 
and a higher-order quadrature is used to evaluate the second integral.  For our reacting flow model,
the underlying projection methodology for the time-advancement of velocity is second-order,
so we require the use of second-order (or higher) numerical quadrature for the second integral.

MISDC Correction Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^

We based the time advance here on a variant of SDC, referred to as MISDC, in which :math:`F` is decomposed into distinct
processes, each treated separately with methods appropriate to its own time scale.  Here, we write

.. math::

    \phi_t = F \equiv A(\phi) + D(\phi) + R(\phi),\label{eq:multi}

to refer to advection, diffusion, and reaction processes.
For this construction we assume that we are given an approximate solution :math:`\phi^{(k)}` that
we want to improve. 
A series of correction equations is developed to update :math:`\phi^{(k)}` that uses relatively
simple second-order discretizations of :math:`A(\phi)` and :math:`D(\phi)` but a high-accuracy 
treatment of :math:`R(\phi)`.  In our approach, :math:`A(\phi^{(k)})` is piecewise-constant over 
each time step, and is evaluated using a second-order Godunov procedure.
The Godunov procedure computes a time-centered 
advection term at :math:`t^{n+1/2}`, and incorporates an explicit diffusion source term and an 
iteratively lagged reaction source term, i.e.,

.. math::

    A(\phi^{(k)}) \equiv A^{n+1/2,(k)} = A\left(\phi^n,D(\phi^n),I_R^{(k-1)}\right),

where :math:`I_R^{(k-1)}` is the effective contribution due to reactions from the previous iteration, i.e.,

.. math::

    I_R^{(k-1)} = \frac{1}{\Delta t^n}\int_{t^n}^{t^{n+1}} R({\phi}^{(k-1)})~d\tau.\label{eq:IR}

where :math:`\Delta t^n = t^{n+1} - t^n`.  Here :math:`I_R^{(k-1)}` is computed from a high-accuracy
integration of the reaction kinetics equations,
augmented with piecewise constant-in-time representation of advection and diffusion.
Details of this procedure are given below. 

We also represent :math:`D(\phi^{(k)})` as piecewise constant 
over the time step, using a mid-point rule:

.. math::

    D(\phi^k) =  \frac{1}{2} (D(\phi^n) +  D(\phi^{(n+1,k)}))

In the spirit of MISDC, we solve correction equations for the individual processes
sequentially.  We begin by discretizing the update equation, but only
including the advection and diffusion terms in the correction integral,

.. math::

    \phi_{\rm AD}^{(k+1)}(t) = \phi^n + \int_{t^n}^t
    \left[A^{(k+1)} - A^{(k)} + D^{(k+1)} - D^{(k)}\right]d\tau
    + \int_{t^n}^t F^{(k)}~d\tau.

Thus, :math:`\phi_{\rm AD}^{(k+1)}(t)` represents an updated approximation of the solution after correcting the
advection and diffusion terms only.  For the first integral, we use an explicit update for the advection term and a 
backward Euler discretization for the diffusion term.
For the second integral, we represent :math:`F` in terms of :math:`A`, :math:`D`, and :math:`R` and
use the definition
of :math:`A^{(k)}`, :math:`D^{(k)}`, and :math:`I_R^{(k-1)}` to obtain
a discrete update for 
:math:`\phi_{\rm AD}^{n+1,(k+1)}`:

.. math::

    \phi_{\rm AD}^{n+1,(k+1)} &=& \phi^n + \Delta t
    \left[A^{n+1/2,(k+1)} - A^{n+1/2,(k)} + D_{\rm AD}^{n+1,(k+1)} - D^{n+1,(k)}\right] \\
    &&\hspace{0.5cm}+ \Delta t\left[A^{n+1/2,(k)} + \frac{1}{2}\left(D^n + D^{n+1,(k)}\right) + I_R^{(k)}\right],

This equation simplifies to the following backward Euler type linear system, with the
right-hand-side consisting of known quantities:

.. math:: \phi_{\rm AD}^{n+1,(k+1)} - \Delta t D_{\rm AD}^{n+1,(k+1)} = \phi^n + \Delta t \left[A^{n+1/2,(k+1)} + \frac{1}{2}\left(D^n - D^{n+1,(k)}\right) + I_R^{(k)}\right],
    :label: ADimplicit

After computing :math:`\phi_{\rm AD}^{n+1,(k+1)}`, we complete the update by solving a correction equation for
the reaction term.  Standard MISDC approaches would formulate the reaction correction equation as

.. math::

    {\phi}^{(k+1)}(t) = \phi^n &+& \int_{t^n}^t \left[ A^{(k+1)} - A^{(k)}
    + D_{\rm AD}^{(k+1)} - D^{(k)} \right]~d\tau\\
    &+& \int_{t^n}^t \left[R^{(k+1)} - R^{(k)}\right]d\tau + \int_{t^n}^t F^{(k)}~d\tau,

and use a backward Euler type discretization for the integral of the reaction terms.
Here, to address stiffness issues with detailed chemical kinetics, we will instead
formulate the correction equation for the 
reaction as an ODE, which is treated separately with an ODE integrator package.
In particular, by differentiating the SDC update we obtain

.. math::

    {\phi}^{(k+1)}_t &=& \left[ A^{n+1/2,(k+1)} - A^{n+1/2,(k)} + D_{\rm AD}^{n+1,(k+1)} - D^{n+1,(k)} \right]\\
    &&\hspace{-0.5cm}+ \left[R^{(k+1)} - R^{(k)}\right] + \left[A^{n+1/2,(k)} +
    \frac{1}{2}\left(D^n + D^{n+1,(k)}\right) + R^{(k)}\right]\\
    &=& R^{(k+1)} + \underbrace{A^{n+1/2,(k+1)} + D_{\rm AD}^{n+1,(k+1)} +
    \frac{1}{2}\left[D^n - D^{n+1,(k)}\right]}_{F_{\rm AD}^{n+1,(k+1)}},

which we then advance with the ODE integrator over :math:`\Delta t` to obtain :math:`\phi^{n+1,(k+1)}`.
After the integration, we can evaluate :math:`I_R^{(k+1)}`, which is required for the next iteration

.. math::

    I_R^{(k+1)} = \frac{\phi^{n+1,(k+1)} - \phi^n}{\Delta t} - F_{\rm AD}^{n+1,(k+1)}.

Summarizing, the variant of SDC used in the single-level time-step of `PeleLM` integrates the :math:`A`, :math:`D` and :math:`R` components of the discretization scheme in an iterative fashion, and each process incorporates a source term that is constructed using a lagged approximation of the other processes. In the case of the implicit diffusion, an additional source term arises from the SDC formulation.  If the SDC iterations were allowed to fully converge, all the process advanced implicitly would be implicitly coupled to all others.  Moreover, each process is discretized using methods that are tailored specifically to the needs of that operator. In the next section, we give more details for each of the components, including how and where the *velocity projections* play a role.

Data centering, :math:`A`-:math:`D`-:math:`R`, and the projections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`PeleLM` implements a finite-volume, Cartesian grid discretization approach with constant grid spacing, where
:math:`U`, :math:`\rho`, :math:`\rho Y_m`, :math:`\rho h`, and :math:`T` represent cell averages, and the pressure field, :math:`\pi`, is defined on the nodes
of the grid, and is temporally constant on the intervals over the time step. There are three major steps in the algorithm:\

**Step 1**: (*Compute advection velocities*) Use a second-order Godunov procedure to predict a time-centered
velocity, :math:`U^{{\rm ADV},*}`, on cell faces using the cell-centered data (plus sources due to any auxiliary forcing) at :math:`t^n`,
and the lagged pressure gradient from the previous time interval, which we denote as :math:`\nabla \pi^{n-1/2}`.  
The provisional field, :math:`U^{{\rm ADV},*}`, fails to 
satisfy the divergence constraint.  We apply a discrete projection by solving the elliptic equation
with a time-centered source term:

.. math::

    D^{{\rm FC}\rightarrow{\rm CC}}\frac{1}{\rho^n}G^{{\rm CC}\rightarrow{\rm FC}}\phi
    = D^{{\rm FC}\rightarrow{\rm CC}}U^{{\rm ADV},*} - \left(\widehat S^n
    + \frac{\Delta t^n}{2}\frac{\widehat S^n - \widehat S^{n-1}}{\Delta t^{n-1}}\right),

for :math:`\phi` at cell-centers, where :math:`D^{{\rm FC}\rightarrow{\rm CC}}` represents a cell-centered divergence of face-centered data,
and :math:`G^{{\rm CC}\rightarrow{\rm FC}}` represents a face-centered gradient of cell-centered data, and :math:`\rho^n` is computed on
cell faces using arithmetic averaging from neighboring cell centers.  Also, :math:`\widehat S` refers to the RHS of the constraint
equation, with adjustments to be discussed in the next section -- these adjustments are computed to ensure that the final update satisfied the equation of state. The solution, :math:`\phi`, is then used to define

.. math::

    U^{\rm ADV} = U^{{\rm ADV},*} - \frac{1}{\rho^n}G^{{\rm CC}\rightarrow{\rm FC}}\phi,

After the *MAC*-projection, :math:`U^{\rm ADV}` is a second-order accurate, staggered grid vector
field at :math:`t^{n+1/2}` that discretely satisfies the constraint.  This field is the advection velocity used for computing
the time-explicit advective fluxes for :math:`U`, :math:`\rho h`, and :math:`\rho Y_m`.

**Step 2**: (*Advance thermodynamic variables*) Integrate :math:`(\rho Y_m,\rho h)` over the full time step.  The details of this are presented in the next subsection.

**Step 3**: (*Advance the velocity*) Compute an intermediate cell-centered velocity field, 
:math:`U^{n+1,*}` using the lagged pressure gradient, by solving

.. math::

    \rho^{n+1/2}\frac{U^{n+1,*}-U^n}{\Delta t}
    + \rho^{n+1/2}\left(U^{\rm ADV}\cdot\nabla U\right)^{n+1/2} = \\
    \frac{1}{2}\left(\nabla\cdot\tau^n
    + \nabla\cdot\tau^{n+1,*}\right) - \nabla\pi^{n-1/2} + \frac{1}{2}(F^n + F^{n+1}),

where :math:`\tau^{n+1,*} = \mu^{n+1}[\nabla U^{n+1,*} +(\nabla U^{n+1,*})^T - 2\mathcal{I}\widehat S^{n+1}/3]` and 
:math:`\rho^{n+1/2} = (\rho^n + \rho^{n+1})/2`, and :math:`F` is the velocity forcing.  This is a semi-implicit discretization for :math:`U`, requiring
a linear solve that couples together all velocity components.  The time-centered velocity in the advective derivative,
:math:`U^{n+1/2}`, is computed in the same way 
as :math:`U^{{\rm ADV},*}`, but also includes the viscous stress tensor evaluated at :math:`t^n` as a source term
in the Godunov integrator.  At 
this point, the intermediate velocity field :math:`U^{n+1,*}` does not satisfy the constraint.  Hence, we apply an 
approximate projection to update the pressure and to project :math:`U^{n+1,*}` onto the constraint surface.  
In particular, we compute :math:`\widehat S^{n+1}` from the new-time 
thermodynamic variables and an estimate of :math:`\dot\omega_m^{n+1}`, which is evaluated
directly from the new-time thermodynamic variables. We project the new-time velocity by solving the elliptic equation,

.. math::

    L^{{\rm N}\rightarrow{\rm N}}\phi = D^{{\rm CC}\rightarrow{\rm N}}\left(U^{n+1,*}
    + \frac{\Delta t}{\rho^{n+1/2}}G^{{\rm N}\rightarrow{\rm CC}}\pi^{n-1/2}\right) - \widehat S^{n+1}

for nodal values of :math:`\phi`.  Here, :math:`L^{{\rm N}\rightarrow{\rm N}}` represents a nodal Laplacian of nodal data, computed
using the standard bilinear finite-element approximation to :math:`\nabla\cdot(1/\rho^{n+1/2})\nabla`.
Also, :math:`D^{{\rm CC}\rightarrow{\rm N}}` is a discrete
second-order operator that approximates the divergence at nodes from cell-centered data 
and :math:`G^{{\rm N}\rightarrow{\rm CC}}` approximates a cell-centered gradient from nodal data.  Nodal 
values for :math:`\widehat S^{n+1}` required for this equation are obtained by interpolating the cell-centered values.  Finally, we 
determine the new-time cell-centered velocity field using

.. math::

    U^{n+1} = U^{n+1,*} - \frac{\Delta t}{\rho^{n+1/2}}G^{{\rm N}\rightarrow{\rm CC}}(\phi-\pi^{n-1/2}),

and the new time-centered pressure using :math:`\pi^{n+1/2} = \phi`.

Thus, there are three different types of linear solves required to advance the velocity field.  The first is the *MAC* solve in order to obtain *face-centered* velocities used to compute advective fluxes.  The second is the multi-component *cell-centered* solver is used to obtain the provisional new-time velocities.  Finally, a *nodal* solver is used to project the provisional new-time velocities so that they satisfy the constraint.

Thermodynamic Advance
^^^^^^^^^^^^^^^^^^^^^

Here we describe the details of **Step 2** above, in
which we iteratively advance :math:`(\rho Y_m,\rho h)` over the full time step.
We begin by computing the diffusion
operators at :math:`t^n` that will be needed throughout the iteration.  Specifically, we evaluate the transport coefficients
:math:`(\lambda,C_p,\mathcal D_m,h_m)^n` from :math:`(Y_m,T)^n`, and the provisional diffusion
fluxes, :math:`\widetilde{\boldsymbol{\cal F}}_m^n`.  These fluxes are conservatively
corrected (i.e., adjusted to sum to zero by adding a mass-weighted "correction velocity") to obtain :math:`{\boldsymbol{\cal F}}_m^n` such that :math:`\sum {\boldsymbol{\cal F}}_m^n = 0`.
Finally, we copy the transport coefficients, diffusion fluxes and the thermodynamic state from :math:`t^n` as starting values for
:math:`t^{n+1}`, and initialize the reaction terms, :math:`I_R` from the values used in the previous step.
The following sequence is then repeated for each iteration, :math:`k<k_{max}`

**Step 2-I:** Use a second-order Godunov integrator to predict
species time-centered edge states, :math:`(\rho Y_m)^{n+1/2,(k+1)}` and their advection terms at :math:`t^{n+1/2}`,
:math:`(A_m^{n+1/2,(k+1)})`. Source terms for this prediction include
explicit diffusion forcing, :math:`D^{n}`, and an iteration-lagged reaction term, :math:`I_R^{(k)}`.
Since the remaining steps of the algorithm for this iteration (including diffusion and chemistry advances) 
will not affect the new-time density for this iteration, we can already compute :math:`\rho^{n+1,(k+1)}`.  
This will be needed in the trapezoidal-in-time diffusion solves.

.. math::

    \frac{\rho^{n+1,(k+1)} - \rho^n}{\Delta t} = A_{\rho}^{n+1/2,(k+1)} = \sum A_{m}^{n+1/2,(k+1)}
    = -\sum_m\nabla\cdot\left(U^{\rm ADV}\rho Y_m\right)^{n+1/2,(k+1)}.

In addition to predicting :math:`\rho` and :math:`\rho Y_m` to the faces to compute advective fluxes, we also 
need :math:`\rho h` there. We could use a Godunov scheme as well, however, because :math:`h` contains the heat of formation
scaled to an arbitrary reference state, it is not generally monotonic through flames. 
Also, because the equation of state is generally nonlinear, this will often lead to numerically-generated n
on-mononoticity in the temperature field. 
An analytically equivalent approach, based on the fact that temperature should be smoother and monotonic 
through the flame, is to instead predict temperature with the Godunov scheme to the cell faces directly. 
Then, using :math:`T`, :math:`\rho = \sum (\rho Y_m)` and :math:`Y_m = (\rho Y_m)/\rho` on the cell faces directly, 
we can evaluate :math:`h` instead of extrapolating. 
We can then evaluate the enthalpy advective flux divergence, :math:`A_h^{n+1/2,(k+1)}`, for :math:`\rho h`. 


**Step 2-II:** Update the transport coefficients (if necessary) with the most current cell-centered thermodynamic
state, then interpolate those values to the cell faces.
We now compute provisional, time-advanced species mass fractions, :math:`\widetilde Y_{m,{\rm AD}}^{n+1,(k+1)}`,
by solving a backward Euler type correction equation for the Crank-Nicolson update, using Eq. :eq:`ADimplicit`. 
Note that the provisional species diffusion fluxes reads 
:math:`\widetilde{\boldsymbol{\cal F}}_{m,{\rm AD}}^{(0)} = -\rho^n D_{m,mix}^n \nabla \widetilde X_{m,{\rm AD}}^{(0)}`. 
This expression couples together all of the species mass fractions (:math:`Y_m`) in the update of each, 
even for the mixture-averaged model. Computationally, it is much more tractable to write this as a diagonal matrix update 
with a lagged correction by noting that :math:`X_m = (W/W_m)Y_m`.  
Using the chain rule, :math:`\widetilde{\boldsymbol{\cal F}}_{m,{\rm AD}}^{(0)}` then has components proportional to :math:`\nabla Y_m` and :math:`\nabla W`. The latter is lagged in the iterations, and is typically very small. In the limit of sufficient iterations, diffusion is driven by the true form of the the driving force, :math:`d_m`, but in this form, each iteration involves decoupled diagonal solves (following the SDC formalism used above):

.. math::

    \frac{\rho^{n+1,(k+1)}\widetilde Y_{m,{\rm AD}}^{n+1,(k+1)} - (\rho Y_m)^n}{\Delta t}
    = A_m^{{n+1/2,(k+1)}} + \widetilde D_{m,AD}^{n+1,(k+1)} + \frac{1}{2}(D_m^n - D_m^{n+1,(k)}) + I_{R,m}^{(k)}

where

.. math::

    &D_m^n = - \nabla \cdot {\boldsymbol{\cal F}}_m^n\\
    &D_m^{n+1,(k)} = - \nabla \cdot {\boldsymbol{\cal F}}_m^{n+1,(k)}\\
    &\widetilde D_{m,AD}^{n+1,(k+1)} = - \nabla \cdot \widetilde {\boldsymbol{\cal F}}_{m,AD}^{n+1,(k+1)}\\
    &\widetilde D_{m,AD}^{n+1,(k+1)} = \nabla \cdot \Big[ \rho^{n+1,(k+1)} D_{m,mix}^{n+1,(k)}\frac{W^{n+1,(k)}}{W_m}\nabla\widetilde Y_{m,{\rm AD}}^{n+1,(k+1)}
    \; + \; \rho^{n+1,(k+1)} D_{m,mix}^{n+1,(k)} \frac{Y_m^{n+1,(k)}}{W_m} \nabla W^{n+1,(k)} \Big]

By iteratively lagging the :math:`\nabla W` term (and :math:`D_{m,mix}`), this equation is a scalar, time-implicit, 
parabolic and linear for the updated :math:`\widetilde Y_{m,{\rm AD}}^{n+1,(k+1)}` (and requires a linear solve).  
The form of this solve, from a software perspective, is identical to that of the *MAC* projection discussed above.

Once all the species equations are updated, we compute :math:`{\boldsymbol{\cal F}}_{m,{\rm AD}}^{n+1,(k+1)}`,
which are conservatively corrected versions of :math:`\widetilde{\boldsymbol{\cal F}}_{m,{\rm AD}}^{n+1,(k+1)}`,
and then the species mass fractions are updated too, using

.. math:: \frac{\rho^{n+1,(k+1)}Y_{m,{\rm AD}}^{n+1,(k+1)} - (\rho Y_m)^n}{\Delta t}
    = A_m^{{n+1/2,(k+1)}} + D_{m,AD}^{n+1,(k+1)} + \frac{1}{2}(D_m^n - D_m^{n+1,(k)}) + I_{R,m}^{n+1,(k)}
    :label: SDCspec

where

.. math::

    D_{m,AD}^{n+1,(k+1)} = - \nabla \cdot {\boldsymbol{\cal F}}_{m,{\rm AD}}^{n+1,(k+1)}

Next, we compute the time-advanced enthalpy, :math:`h_{\rm AD}^{n+1,(k+1)}`.  Much like for the diffusion of the species densities,
:math:`\rho Y_m`, where a :math:`\nabla X_m` driving force leads to a nonlinear, coupled Crank-Nicolson update; the
enthalpy diffuses with a :math:`\nabla T` driving force. We define an alternative linearized strategy.
We begin by following the same SDC-correction formalism used for the species, and write
the nonlinear update for :math:`\rho h` (noting that there is no reaction source term here):

.. math:: \frac{\rho^{n+1,(k+1)} h_{{\rm AD}}^{n+1,(k+1)} - (\rho h)^n}{\Delta t} = A_h^{n+1/2,(k+1)} + D_{T,AD}^{n+1,(k+1)} + H_{AD}^{n+1,(k+1)}\\
    + \frac{1}{2} \Big( D_T^n - D_T^{n+1,(k)} + H^n - H^{n+1,(k)} \Big)
    :label: SDCrhoH

where

.. math::

    &D_T^n = \nabla \cdot \lambda^n \nabla T^n,   \hspace{2cm}
    &H^n = - \nabla \cdot \sum h_m(T^n) \; {\boldsymbol{\cal F}}_m^n\\
    &D_T^{n+1,(k)} = \nabla \cdot \lambda^{n+1,(k)} \nabla T^{{k}},
    &H^{n+1,(k)} = - \nabla \cdot  \sum h_m(T^{n+1,(k)}) \; {\boldsymbol{\cal F}}_m^{n+1,(k)}\\
    &D_{T,AD}^{n+1,(k+1)} = \nabla \cdot \lambda_{AD}^{n+1,(k+1)} \nabla T_{AD}^{n+1,(k+1)},
    &H_{AD}^{n+1,(k+1)} = - \nabla \cdot \sum h_m(T_{AD}^{n+1,(k+1)}) \; {\boldsymbol{\cal F}}_{m,AD}^{n+1,(k+1)}

However, since we cannot compute :math:`h_{{\rm AD}}^{n+1,(k+1)}` directly, we solve this iteratively based on the approximation
:math:`h_{{\rm AD}}^{(k+1),\ell+1} \approx h_{{\rm AD}}^{(k+1),\ell} + C_{p}^{(k+1),\ell} \delta T^{(k+1),\ell+1}`, with
:math:`\delta T^{(k+1),\ell+1} = T_{{\rm AD}}^{(k+1),\ell+1} - T_{{\rm AD}}^{(k+1),\ell}`, and iteration index, 
:math:`\ell` = 1::math:`\,\ell_{MAX}`. 
The enthalpy update equation is thus recast into a linear equation for :math:`\delta T^{(k+1);\ell+1}`

.. math::

    \rho^{n+1,(k+1)} C_p^{(k+1),\ell} \delta T^{(k+1),\ell+1}
    &-& \Delta t \, \nabla \cdot \lambda^{(k)} \nabla (\delta T^{(k+1),\ell +1}) \nonumber  \\
    &=& \rho^n h^n - \rho^{n+1,(k+1)} h_{AD}^{(k+1),\ell} + \Delta t \Big( A_h^{n+1/2,(k+1)} + D_{T,AD}^{(k+1),\ell}
    + H_{AD}^{(k+1),\ell} \Big) \\
    &&+ \; \frac{\Delta t}{2} \Big( D_T^n - D_T^{n+1,(k)} + H^n - H^{n+1,(k)} \Big) \nonumber

where :math:`H_{AD}^{(k+1),\ell} = - \nabla \cdot \sum h_m(T_{AD}^{(k+1),\ell}) \, {\boldsymbol{\cal F}}_{m,AD}^{n+1,(k+1)}`
and :math:`D_{T,AD}^{(k+1),\ell} = \nabla \cdot \lambda^{(k)} \, \nabla T_{AD}^{(k+1),\ell}`.
Note that again the solve for this
Crank-Nicolson update has a form that is identical to that of
the *MAC* projection discussed above.  After each 
iteration, update :math:`T_{{\rm AD}}^{(k+1),\ell+1} = T_{{\rm AD}}^{(k+1),\ell} + \delta T^{(k+1),\ell+1}` and 
re-evaluate :math:`(C_p ,h_m)^{(k+1),\ell+1}` using :math:`(T_{{\rm AD}}^{(k+1),\ell+1}, Y_{m,{\rm AD}}^{n+1,(k+1)}`).

**Step 2-III:**
Based on the updates above, we define an effective contribution of advection and diffusion to the
update of :math:`\rho Y_m` and :math:`\rho h`:

.. math::

    &&Q_{m}^{n+1,(k+1)} = A_m^{n+1/2,(k+1)} + D_{m,AD}^{(n+1,k+1)} + \frac{1}{2}(D_m^n - D_m^{n+1,(k)}) \\
    &&Q_{h}^{n+1,(k+1)} = A_h^{n+1/2,(k+1)} + D_{T,AD}^{n+1,(k+1)} + H_{AD}^{n+1,(k+1)} + \frac{1}{2}(D_T^n - D_T^{n+1,(k)} + H^n - H^{n+1,(k)} )

that we treat as piecewise-constant source terms to advance :math:`(\rho Y_m,\rho h)^n` to :math:`(\rho Y_m,\rho h)^{n+1,(k+1)}`.
The ODE system for the reaction part over :math:`\Delta t^n` then takes the following form:

.. math::

    \frac{\partial(\rho Y_m)}{\partial t} &=& Q_{m}^{n+1,(k+1)} + \rho\dot\omega_m(Y_m,T),\label{eq:MISDC VODE 3}\\
    \frac{\partial(\rho h)}{\partial t} &=& Q_{h}^{n+1,(k+1)}.\label{eq:MISDC VODE 4}

After the integration is complete, we make one final call to the equation of state
to compute :math:`T^{n+1,(k+1)}` from :math:`(Y_m,h)^{n+1,(k+1)}`.  We also can compute the effect of reactions
in the evolution of :math:`\rho Y_m` using,

.. math::

    I_{R,m}^{(k+1)} = \frac{(\rho Y_m)^{n+1, (k+1)} - (\rho Y_m)^n}{\Delta t} - Q_{m}^{n+1,(k+1)}.

If :math:`k<k_{\rm max}-1`, set :math:`k=k+1` and return to **Step 2-I**.  Otherwise, the 
time-advancement of the thermodynamic variables is complete, and set 
:math:`(\rho Y_m,\rho h)^{n+1} = (\rho Y_m,\rho h)^{n+1,(k+1)}`.
If :math:`k+1 = k_{max}`, **Step 2** of our algorithm is complete.


Modifications for AMR
^^^^^^^^^^^^^^^^^^^^^

The framework to manage adaptive mesh refinement (AMR) used in `PeleLM` borrows heavily from the `AMReX` library,
and the `IAMR` code; the reader is referred to documentation of both of these components in order to understand the
distributed, logically rectangular data structures used, and the recursive time-stepping strategy for
advancing a hierarchy of nested grid levels.  

Summarizing, there is a bulk-synchronous advance of each level over its
respective time step, :math:`dt`, followed recursively by a number of (sub-)steps of the next-finer AMR level.
Each fine level advanced is over an interval :math:`(1/R) dt`, if the fine cells are a factor of :math:`R`
smaller, and in this scenario, the coarser level provides Dirichlet boundary condition data
for the fine-level advances.  Note that the levels are properly nested so that the finer level is fully contained within
the coarser level, except perhaps at physical boundaries, where their edges can be coincident - thus, the fine level has
sufficient boundary data for a well-posed advance.

After two adjacent levels in the hierarchy reach the same physical time, a *synchronization* operation is performed
to ensure that the coarse data is consistent with the volume integral of the fine data that covers it, and the
fluxes across of the coarse-fine interface are those of the fine solution.  The latter of these two operations can be
quite complex, as it must correct coarse-grid errors committed by each of the operators used to perform the original
advance.  It may also be non-local, in that cells far away from the coarse-fine interface may need to
incorporate flux increments due to the mismatched coarse and fine solutions. Formally, the synchronzation is a bilevel correction
that should be computed as a sequence of two-level solves. However, this would lead to the same amound of work
that was required to create the original (pre-sync) data. We assume that the corrections computed for the synchronization are smooth 
enough to be well represented by an increment on the coaser of the two-levels, and interpolated to the finer grid. Note that the transport
coefficients are not updated to account for the state changes during the synchronization.

Generically, the synchronization procedure in `PeleLM` follows that described for the `IAMR` code, but with
modifications to explicitly enforce that the sum of the species diffusion correction fluxes is zero, that the
nonlinear enthalpy update is solved (similar to described above for the single-level advance), and the
corection for the advection velocity is adjusted iteratively so that the final synchronized state satisfies the EOS.

There are several components in `PeleLM` that contribute to the flux mismatch at the coarse-fine interface.
The first component arise from the face-centered velocity, :math:`U^{ADV,\ell}`, used to advect the scalars
at each AMR level :math:`\ell`, since the field satifies a divergence constraint on the coarse and fine levels
separately.  We compute a velocity mismatch

.. math::

    \delta U^{ADV,\ell} = -U^{ADV,\ell,n+1/2} + \frac{1}{R^{d-1}}\sum_{k=0}^{R-1} \sum_{edges} U^{ADV,\ell+1,n+k+1/2}

(where :math:`d` is the number of spatial dimensions) along the coarse-fine boundary.  We then solve the elliptic
projection equation

.. math::

    D^{MAC} \frac{1}{\rho} \delta e^{\ell} = D^{MAC} \delta U^{ADV,\ell} + \delta_{\chi}^{\ell}

where :math:`\delta_{\chi}^{\ell}` is incremented iteratively to enfoce the final state to satisfy the EOS
and compute the correction velocity

.. math::

    U^{ADV,\ell,corr} = -\frac{1}{\rho} G^{MAC} \delta e^{\ell}

which is the increment of velocity required to carry advection fluxes needed to correct the errors made by advancing the
coarse state with the wrong velocities.

The second part of the mismatch arises because the advective and diffusive fluxes on the coarse grid were computed 
without explicitly accounting for the fine grid, while on the fine grid the fluxes were computed using coarse-grid 
Dirichlet boundary data.  We define the flux discrepancies on the coarser level :math:`\ell` of the pair of levels considered:

.. math::

    \delta \boldsymbol{\cal F^{\ell}} = \Delta t^{\ell} \Big(
    -\boldsymbol{\cal F}^{\ell,n+1/2} + \frac{1}{R^{d-1}} \sum_{k=0}^{R-1} \sum_{edges}
    \boldsymbol{\cal F}^{\ell+1,n+k+1/2} \Big)

where :math:`\boldsymbol{\cal F}` is the total (advective + diffusive) flux through a face on the coarse-fine
interface prior the synchronization operations. Since all operations are performed on the coarse level we will drop 
the :math:`\ell` in the following.

Since mass is conserved, corrections to density, :math:`\delta \rho^{sync}` on the coarse grid associated with
mismatched advection fluxes may be computed explicitly

.. math::

    \delta \rho^{sync} = -D^{MAC} \Big( \sum_m U^{ADV,corr} \rho Y_m \Big)^{n+1/2} 
    + \sum_m \nabla \cdot \delta \boldsymbol{\cal F}_{m}

We can compute the post-sync new-time value of density, :math:`\rho^{n+1} = \rho^{n+1,p} + \delta \rho^{sync}`,
where :math:`^p` denotes *pre*-sync quantities.
The synchronization correction of a state variable :math:`\delta ( \rho \phi )^{sync}`
(where :math:`\phi \in (Y_m, h)`) is obtained by subtracting the pre-sync state value :math:`( \rho \phi )^{n+1,p}` from the corrected
one :math:`( \rho \phi )^{n+1}`, both of which expressed from an SDC iteration update (see Eq. :eq:`SDCspec` and :eq:`SDCrhoH`)
but with the divergence of the correction velocity fluxes (:math:`U^{ADV,corr}`) 
and fluxes mismatch (:math:`\delta \boldsymbol{\cal F}`) included in the advection and diffusion corrected operators.

Given the corrected density we can decompose the sync corrections 
:math:`\delta ( \rho \phi )^{sync} = \phi^{n+1,p} \delta \rho^{sync} + \rho^{n+1} \delta \phi^{sync}` 
and obtain the linear system for :math:`\delta \phi^{sync}` since the fluxes mismatch contain implicit 
diffusion fluxes from the Crank-Nicolson discretization. For species :math:`m` the implicit system reads:

.. math::
    \rho^{n+1} \widetilde{\delta Y_m^{sync}} - \Delta t \nabla \cdot 
    \widetilde{\cal F_{m}}(\widetilde{\delta Y_m^{sync}})
    = -D^{MAC} (U^{ADV,corr} \rho Y_m)^{n+1/2} + \nabla \cdot \delta \boldsymbol{\cal F}_{m} - Y_m^{n+1,p} \delta \rho^{sync}
    :label: specSyncEq

where :math:`\widetilde{\cal F_{m}}` is the species correction flux due to the sync correction, 
:math:`\widetilde{\delta Y_m^{sync}}`. However, as in the single-level algorithm, the species
fluxes must be corrected to sum to zero.  These adjusted fluxes are then used to recompute a :math:`\delta Y_m^{sync}`,
which is then used via the expression above to compute :math:`\delta (\rho Y_m)^{sync}`, the increment to the
species mass densities.

In order to get the equation for the enthalpy sync correction, we operate as for species mass fractions.
We will present the details of the method.
The SDC advection-diffusion udpate for *pre*-sync enthalpy is :eq:`SDCrhoH` (now including the superscript :math:`p`) 
and its corrected counterpart reads:

.. math:: \frac{\rho^{n+1} h_{{\rm AD}}^{n+1} - (\rho h)^n}{dt}
    = A_h^{n+1/2,(k+1),*} + D_T^{n+1,(k+1),*} + H^{n+1,(k+1),*} \\
    + \frac{1}{2} \Big( D_T^{n,*} - D_T^{n+1,(k),*} + H^{n,*} - H^{n+1,(k),*} \Big)
    :label: SDCrhoHcorr

where

.. math::
    A_h^{n+1/2,(k+1),*} = - \nabla \cdot ( \rho h (U^{ADV} + U^{ADV,corr})^{n+1/2,(k+1)} + \delta \boldsymbol{\cal F}_{Adv,h} ) \\
                        = A_h^{n+1/2,(k+1),p} - \nabla \cdot ( \rho h ( U^{ADV,corr})^{n+1/2,(k+1)} +  \delta \boldsymbol{\cal F}_{Adv,h} ) \\
    \\
    D_T^{n,*}           = \nabla \cdot ( \lambda^{n} \nabla T^{n} + \delta \boldsymbol{\cal F}_{DT,h}^{n} ) = D_T^{n,p} + \nabla \cdot ( \delta \boldsymbol{\cal F}_{DT,h}^{n} ) \\
    \\
    D_T^{n+1,(k),*}     = \nabla \cdot ( \lambda^{n+1,(k)} \nabla T^{n+1,(k)} + \delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k)} ) \\
                        = D_T^{n+1,(k),p} + \nabla \cdot (\delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k)} ) \\
    \\
    D_T^{n+1,(k+1),*}   = \nabla \cdot ( \lambda^{n+1,(k+1),p} \nabla (T^{n+1,(k+1),p}+\delta T^{sync}) + \delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k+1)} ) \\    
                        = D_T^{n+1,(k+1),p} + \nabla \cdot ( \lambda^{n+1,(k+1),p} \nabla \delta T^{sync} + \delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k+1)}) \\
    \\
    H^{n,*}            = - \nabla \cdot (\sum_m h_m(T^{n}) \; {\boldsymbol{\cal F}}_m^{n} + \delta \boldsymbol{\cal F}_{DH,h}^{n} ) = H^{n,p} - \nabla \cdot (\delta \boldsymbol{\cal F}_{DH,h}^{n}) \\
    \\
    H^{n+1,(k),*}      = - \nabla \cdot (\sum_m h_m(T^{n+1,(k)}) \; {\boldsymbol{\cal F}}_m^{n+1,(k)} + \delta \boldsymbol{\cal F}_{DH,h}^{n+1,(k)} )\\
                       = H^{n+1,(k),p} - \nabla \cdot ( \delta \boldsymbol{\cal F}_{DH,h}^{n+1,(k)} )\\
    \\
    H^{n+1,(k+1),*} = - \nabla \cdot (\sum_m h_m(T^{n+1,(k+1),p}+\delta T^{sync}) \; {\boldsymbol{\cal F}}_m^{n+1,(k+1)} + \delta \boldsymbol{\cal F}_{DH,h}^{n+1,(k+1)} )
   
and with :math:`\delta T^{sync} = T^{n+1,(k+1)} - T^{n+1,(k+1),p}`. 
Subtracting the *pre*-sync eq. :eq:`SDCrhoH` (with the superscript :math:`p`) from the above 
equation :eq:`SDCrhoHcorr` and gathering the flux mismatches and correction velocity fluxes in :math:`S_h^{sync}` we obtain: 

.. math:: \frac{\delta(\rho h)^{sync}}{\Delta t} = S_h^{sync} + D_T^{n+1,(k+1)} - D_T^{n+1,(k+1),p} + H^{n+1,(k+1)} - H^{n+1,(k+1),p}
    :label: rhoHsyncEq

where

.. math::

    S_h^{sync} = - \nabla \cdot ( \rho h ( U^{ADV,corr})^{n+1/2,(k+1)} +  \delta \boldsymbol{\cal F}_{Adv,h} ) \\
                 + \frac{1}{2} \nabla \cdot ( \delta \boldsymbol{\cal F}_{DT,h}^{n} - \delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k)} )
                 - \frac{1}{2} \nabla \cdot ( \delta \boldsymbol{\cal F}_{DH,h}^{n} - \delta \boldsymbol{\cal F}_{DH,h}^{n+1,(k)} ) \\
                 + \nabla \cdot (\delta \boldsymbol{\cal F}_{DT,h}^{n+1,(k+1)}) - \nabla \cdot (\delta \boldsymbol{\cal F}_{DH,h}^{n+1,(k+1)} )

and the updated (:math:`n+1`) fluxes (without the :math:`*`) only now contain the implicit contribution, i.e:

.. math::

    D_T^{n+1,(k+1)}   =  D_T^{n+1,(k+1),p} + \nabla \cdot ( \lambda^{n+1,(k+1),p} \nabla \delta T^{sync} ) \\
    H^{n+1,(k+1)}     = - \nabla \cdot (\sum_m h_m(T^{n+1,(k+1),p}+\delta T^{sync}) \; {\boldsymbol{\cal F}}_m^{n+1,(k+1)}  

To go further we note that:

.. math::

    D_T^{n+1,(k+1)} - D_T^{n+1,(k+1),p} = \nabla \cdot ( \lambda^{n+1,(k+1),p} \nabla \delta T^{sync} ) \\ 
    H^{n+1,(k+1)}   - H^{n+1,(k+1),p}   = - \nabla \cdot \sum_m (h_m (T^{n+1,(k+1),p} + \delta T^{sync}) \; \delta {\boldsymbol{\cal F}}_m^{sync}  \\
                                          + \delta {h}_m^{sync} \;  {\boldsymbol{\cal F}}_m^{n+1,(k+1),p}) 

where :math:`\delta h_m^{sync} = h_m(T^{n+1,(k+1)}) - h_m(T^{n+1,(k+1),p})` and :math:`\delta {\boldsymbol{\cal F}}_m^{sync}` is 
the species flux increment due to the species sync correction appearing on the LHS of eq. :eq:`specSyncEq`. 
Eq. :eq:`rhoHsyncEq` is the equation for the sync correction. At this point, we can drop the SDC iteration index :math:`k+1` 
for simplicity (all :math:`k` related quantities are contained in :math:`S_h^{sync}`). 
Note that the evaluation of the transport properties is relatively expensive, such that we don't want to update 
the conductivity in :math:`D_T^{n+1}` since a lagged (pre-sync) version is sufficient for second-order accuracy. 
However we do want to use an updated version of :math:`h_m`. 

Just as in the level advance, we cannot compute :math:`h^{n+1}` directly, so we solve this iteratively based on the approximation
:math:`h^{n+1,\eta+1} \approx h^{n+1,\eta} + C_{p}^{n+1,\eta} \Delta T^{\eta+1}`, with
:math:`\Delta T^{\eta+1} = T^{n+1,\eta+1} - T^{n+1,\eta}`, and iteration index, :math:`\eta` = 1::math:`\,\eta_{MAX}`.
The sync equation is thus recast into a linear equation for :math:`\Delta T^{\eta+1}`, and we
lag the :math:`H` terms in iteration :math:`\eta`,

.. math::

    \rho^{n+1} C_p^{n+1,\eta} \Delta T^{\eta +1}
    - dt \, \nabla \cdot \lambda^{(n+1,p)} \nabla (\Delta T^{\eta +1}) \hspace{10em}\\
    = \rho^{(n+1,p)} h^{(n+1,p)} - \rho^{n+1} h^{n+1,\eta} + dt \Big( S_h^{sync} + 
    \nabla \cdot \lambda^{(n+1,p)} \nabla (\delta T^{sync,\eta}) \\
    - \nabla \cdot \sum_m \Big( h_m^{n+1} \delta {\boldsymbol{\cal F}}_{m}^{sync}
    + \delta h_m^{sync} {\boldsymbol{\cal F}}_{m}^{(m+1)} \Big) \Big)

After each :math:\eta iteration, update :math:`T^{n+1,\eta+1} = T^{n+1,\eta} + \Delta T^{\eta+1}`, 
:math:`\delta T^{sync,\eta+1} = T^{n+1,\eta+1} - T^{(n+1,p)}`, and 
re-evaluate :math:`(C_p,h_m)^{n+1,\eta+1}` using :math:`(T^{n+1,\eta+1}, Y_m^{n+1}`).  Iterations are continued 
until the norm of :math:`\Delta T^{\eta+1}` drops below a tolerance threshold. 
Then set :math:`T^{n+1} = T^{(n+1,p)} + \delta T^{sync,\eta_{MAX}}`, and compute :math:`h^{n+1} = h(T^{n+1},Y_m^{n+1})`.
