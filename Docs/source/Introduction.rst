.. highlight:: rst

Introduction
============

`PeleLM` evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). The code depends upon the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library to provide the underlying data structures, and tools to manage and operate on them across massively parallel computing architectures. `PeleLM` also borrows heavily from the source code and algorithmic infrastructure of the `IAMR <https://github.com/AMReX-Codes/IAMR>`_. `IAMR` implements an AMR integration for the variable-density incompressible Navier-Stokes equations. `PeleLM` extends `IAMR` to include complex coupled models for generalized thermodynamic relationships, multi-species transport and chemical reactions.  The core algorithms in `PeleLM` (and `IAMR`) are described in the following papers:

* *A conservative, thermodynamically consistent numerical approach for low Mach number combustion. I. Single-level integration*, A. Nonaka, J. B. Bell, and M. S. Day, *Combust. Theor. Model.*, **22** (1) 156-184 (2018)

* *A Deferred Correction Coupling Strategy for Low Mach Number Flow with Complex Chemistry*, A. Nonaka, J. B. Bell, M. S. Day, C. Gilet, A. S. Almgren, and M. L. Minion, *Combust. Theory and Model*, **16** (6) 1053-1088 (2012)

* *Numerical Simulation of Laminar Reacting Flows with Complex Chemistry*, M. S. Day and J. B. Bell, *Combust. Theory Model* **4** (4) 535-556 (2000)

* *An Adaptive Projection Method for Unsteady, Low-Mach Number Combustion*, R. B. Pember, L. H. Howell, J. B. Bell, P. Colella, W. Y. Crutchfield, W. A. Fiveland, and J. P. Jessee, *Comb. Sci. Tech.*, **140** 123-168 (1998)

* *A Conservative Adaptive Projection Method for the Variable Density Incompressible Navier-Stokes Equations,* A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, and M. L. Welcome, *J. Comp. Phys.*, **142** 1-46 (1998)

The low Mach number flow equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

where :math:`\rho` is the density, :math:`\boldsymbol{u}` is the velocity, :math:`h` is the mass-weighted enthalpy, :math:`T` is temperature and :math:`Y_m` is the mass fraction of species :math:`m`. :math:`\dot{\omega}_m` is the molar production rate for species :math:`m`, the modeling of which will be described later in this section. :math:`\tau` is the stress tensor, :math:`\boldsymbol{\mathcal{Q}}` is the heat flux and \:math:`\boldsymbol{\mathcal{F}}_m` are the species diffusion fluxes. These transport fluxes require the evaluation of transport coefficients (e.g., the viscosity :math:`\mu`, the conductivity :math:`\lambda` and the diffusivity matrix :math:`D`) which are computed using the library EGLIB, as will be described in more depth in the diffusion section. The momentum source, :math:`\boldsymbol{F}`, is an external forcing term.  For example, we have used :math:`\boldsymbol{F}` to implement a long-wavelength time-dependent force to establish and maintain quasi-stationary turbulence.

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

    \nabla \cdot \boldsymbol{u} = \frac{1}{\rho T}\frac{DT}{Dt}
    + \frac{W}{\rho} \sum_m \frac{1}{W_m} \frac{DY_m}{Dt} = S

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

The full diffusion model couples together the advance of all thermodynamics fields, including a dense matrix transport operator that is cumbersome to deal with computationally, while also being generally viewed as an overkill for most practical combustion applications -- particularly those involving turbulent fluid dynamics.  For `PeleLM`, we make the following simplifying assumptions:

1. The bulk viscosity, :math:`\kappa`, is negligible, compared to the shear viscosity,

2. The low Mach limit implies that there are no spatial gradients in the thermodynamic pressure,

3. The *mixture-averaged* diffusion model is assumed,

4. Dufour and Soret effects are negligible

With these assumptions, the conservation equations take the following form:

.. math::

    &&\frac{\partial (\rho \boldsymbol{u})}{\partial t} +
    \nabla \cdot \left(\rho  \boldsymbol{u} \boldsymbol{u} + \tau \right)
    = -\nabla \pi + \rho \boldsymbol{F}, \\
    &&\frac{\partial (\rho Y_m)}{\partial t} +
    \nabla \cdot \left( \rho Y_m \boldsymbol{u} + \boldsymbol{\mathcal{F}}_{m} \right) \\
    &&\frac{ \partial (\rho h)}{ \partial t} +
    \nabla \cdot \left( \rho h \boldsymbol{u} + \boldsymbol{\mathcal{Q}} \right) = 0,

with

.. math::

    &&\boldsymbol{\mathcal{F}}_{m} = \rho Y_m \boldsymbol{V_m} = - \rho D_{m,mix} \nabla X_m \\
    &&\tau_{i,j} = \frac{2}{3} \mu \delta_{i,j} \frac{\partial {u_k}}{\partial x_k} - \mu \Big(
    \frac{\partial  u_i}{\partial x_j} + \frac{\partial  u_j}{\partial x_i}\Big) \\
    &&\boldsymbol{\mathcal{Q}} =  \sum_m h_m \boldsymbol{\mathcal{F}}_{m}  - \lambda \nabla T

Using these expressions, we can write an equation for :math:`T` that is needed in order to evaluate the right-hand side of the divergence constraint:

.. math::

    \rho C_p \frac{DT}{Dt} = \nabla \cdot \lambda \nabla T + \sum_m \Big( h_m \nabla \cdot \boldsymbol{\mathcal{F}}_{m} - \nabla \cdot h_m \boldsymbol{\mathcal{F}}_{m} - h_m \rho \dot\omega_m \Big)

where :math:`C_p = \partial h/\partial T` is the specific heat of the mixture at constant pressure. The constraint then becomes:

.. math::

    \nabla \cdot \boldsymbol{u} &=&\frac{1}{\rho C_p T}\Big[ \nabla \cdot \lambda \nabla T
    + \sum_m \Big( h_m \nabla \cdot \boldsymbol{\mathcal{F}}_{m}
    - \nabla \cdot h_m \boldsymbol{\mathcal{F}}_{m}\Big) \Big] \\
    &&- \frac{W}{\rho} \sum_m \frac{1}{W_m} \nabla \cdot \boldsymbol{\mathcal{F}}_{m}
    + \frac{1}{\rho} \sum_m \Big( \frac{W}{W_m} -\frac{h_m(T)}{c_{p} T} \Big)\dot{\omega}_m

The mixture-averaged transport coefficients discussed above (:math:`\mu`, :math:`\lambda` and :math:`D_{m,mix}`) can be evaluated from transport properties of the pure species. We follow the treatment used in the EGLib library, based on the theory/approximations developed by Ern and Givangigli.


The following choices are currently implemented in `PeleLM`

* The viscosity, :math:`\mu`, is estimated based \textcolor{red}{FIXME}

* The conductivity, :math:`\lambda`, is based on an empirical mixture formula:

.. math::

    \lambda = \frac{1}{2} (\mathcal{A}_{-1} + \mathcal{A}_{1})

with

.. math::

    \mathcal{A}_{\alpha}= \Big( \sum_m X_m (\lambda_m)^{\alpha} \Big)^{1/\alpha}

* The diffusion flux is approximated using the diagonal matrix :math:`diag(\widetilde{ \Upsilon})`, where:

.. math::

    \widetilde{ \Upsilon}_m =  D_{m,mix}, \;\;\;\mbox{where} \;\;\;
    D_{m,mix} = \frac{1-Y_m}{ \sum_{j \neq m} X_j / \mathcal{D}_{m,j}}

This leads to a mixture-averaged approximation that is similar to that of Hirschfelder-Curtiss:

.. math::

    \rho Y_m \boldsymbol{V_m} = - \rho D_{m,mix} \nabla X_m 

Note that with these definitions, there is no guarantee that :math:`\sum \boldsymbol{\mathcal{F}}_{m} = 0`, as required for mass conservation. An arbitrary *correction flux,* consistent with the mixture-averaged diffusion approximation, is added in \pelelm\ to enforce conservation.

Pure species transport properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mixture-averaged transport coefficients require expressions for the pure species binary transport coefficients.  These, in turn, depend upon the forces of interaction between colliding molecules, which are complex functions of the shape and properties of each binary pair of species involved, as well as of their environment, intermolecular distance, etc. In practice, these interactions are usually described by a Lennard-Jones 6-12 potential (for non polar molecules, Stockmayer potential otherwise) that relates the evolution of the potential energy of the pair of species to their intermolecular distance. Here, the single component viscosities and binary diffusion coefficients are given by Hirschfelder:1954:

.. math::

    \eta_m = \frac{5}{16} \frac{\sqrt{\pi m_m k_B T}}{\pi \sigma^2_m \Omega^{(2,2)*}},
    \hspace{4mm}
    \mathcal{D}_{m,j} = \frac{3}{16}\frac{\sqrt{2 \pi k^3_B T^3/m_{m,j}}}{P_0 \pi \sigma^2_{m,j} \Omega^{(1,1)*}}

where :math:`k_B` is the Boltzmann constant, :math:`\sigma_m` is the Lennard-Jones collision diameter and :math:`m_m (= W_k/\mathcal{A})` is the molecular mass of species :math:`m`. :math:`m_{m,j}` is the reduced molecular mass and :math:`\sigma_{m,j}` is the reduced collision diameter of the :math:`(m,j)` pair, given by:

.. math::

    m_{m,j} = \frac{m_m m_j }{ (m_m + m_j)},
    \sigma_{m,j} = \frac{1}{2} \zeta^{-\frac{1}{6}}(\sigma_m + \sigma_j)

where :math:`\zeta=1` if the partners are either both polar or both nonpolar, but in the case of a polar molecule (:math:`p`) interacting with a nonpolar (:math:`n`) molecule:

.. math::

    \zeta=1 + \frac{1}{4} \alpha^*_n (\mu^*_p)^2 \sqrt{\frac{\epsilon_p}{\epsilon_n}}

with :math:` \alpha^*_n = \alpha_n / \sigma^3_n` the reduced polarizability of the nonpolar molecule and  :math:`\mu^*_p = \mu_p/\sqrt{\epsilon_p \sigma^3_p}` the reduced dipole moment of the polar molecule, expressed in function of the Lennard-Jones potential :math:`\epsilon_p` of the :math:`p` molecule.

Both quantities rely upon the evaluation of *collision integrals* :math:`\Omega^{(\cdot,\cdot)*}`, which account for inter-molecular interactions, and are usually tabulated in function of reduced variables:

* :math:`\Omega^{(2,2)*}` is tabulated in function of a reduced temperature, :math:`T^*_m` and a reduced dipole moment, :math:`\delta^*_m`, given by:

.. math::

    T^*_m = \frac{k_BT}{\epsilon_m},
    \delta^*_m = \frac{1}{2} \frac{\mu^2_m}{\epsilon_m \sigma^3_m}

%where :math:`\epsilon_m` is the Lennard-Jones potential well depth and :math:`\mu_m` is the dipole moment of species :math:`m`. 

* :math:`\Omega^{(1,1)*}` is tabulated in function of a reduced temperature, :math:`T^*_{m,j}` and a reduced dipole moment, :math:`\delta^*_{m,j}`, given by:

.. math::

    T^*_{m,j} = \frac{k_BT}{\epsilon_{m,j}},
    \delta^*_{m,j} = \frac{1}{2} \frac{\mu^2_{m,j}}{\epsilon_{m,j} \sigma^3_{m,j}}

where the reduced collision diameter of the pair (:math:`\sigma_{m,j}`) is given by <redCollision>; and the Lennard-Jones potential :math:`\epsilon_{m,j}` and dipole moment :math:`\mu_{m,j}` of the :math:`(m,j)` pair are given by:

.. math::

    \frac{\epsilon_{m,j}}{k_B} = \zeta^2 \sqrt{\frac{\epsilon_m}{k_B} \frac{\epsilon_j}{k_B}},
    \mu^2_{m,j} = \xi \mu_m \mu_j 

with :math:`\xi = 1` if :math:`\zeta = 1` and :math:`\xi = 0` otherwise.

The expression for the pure species thermal conductivities are more complex. They are assumed to be composed of translational, rotational and vibrational contributions:

.. math::

    \lambda_m = \frac{\eta_m}{W_m} (f_{tr}C_{v,tr} + f_{rot}C_{v,rot} + f_{vib}C_{v,vib})

where

.. math::

    &&f_{tr} = \frac{5}{2}\Big(1-\frac{2}{\pi} \frac{C_{v,rot}}{C_{v,tr}} \frac{A}{B} \Big)\\
    &&f_{rot} = \frac{\rho \mathcal{D}_{m,m}}{\eta_m} \Big( 1 + \frac{2}{\pi} \frac{A}{B}  \Big)\\
    &&f_{vib} = \frac{\rho \mathcal{D}_{m,m}}{\eta_m}

and

.. math::

    A = \frac{5}{2} - \frac{\rho \mathcal{D}_{m,m}}{\eta_m},
    B = Z_{rot} + \frac{2}{\pi} \Big( \frac{5}{3} \frac{C_{v,rot}}{\mathcal{R}} + \frac{\rho \mathcal{D}_{m,m}}{\eta_m} \Big)

The molar heat capacities :math:`C_{v,\cdot}` depend on the molecule shape. In the case of a linear molecule:

.. math::

    \frac{C_{v,tr}}{\mathcal{R}} = \frac{3}{2},
    \hspace{1.5em}
    \frac{C_{v,rot}}{\mathcal{R}} = 1,
    \hspace{1.5em} 
    {C_{v,vib}} = C_v - \frac{5}{2} \mathcal{R}

In the case of a nonlinear molecule, the expressions are

.. math::

    \frac{C_{v,tr}}{\mathcal{R}} = \frac{3}{2},
    \hspace{1.5em} 
    \frac{C_{v,rot}}{\mathcal{R}} =  \frac{3}{2},
    \hspace{1.5em} 
    {C_{v,vib}} = C_v - 3 \mathcal{R}

For single-atom molecules the thermal conductivity reduces to:

.. math::

    \lambda_m = \frac{\eta_m}{W_m} (f_{tr}C_{v,tr} ) = \frac{15 \eta_m \mathcal{R}}{4 W_m}

Finally, :math:`Z_{rot}` is the rotational relaxation number, a parameter given by:

.. math::

    Z_{rot}(T) = Z_{rot} (298) \frac{F(298)}{F(T)}

with 

.. math::

    F(T) = 1 + \frac{\pi^{(3/2)}}{2} \sqrt{\frac{\epsilon/k_B}{T} } + \Big( \frac{\pi^2}{4} +2 \Big) \Big( \frac{\epsilon/k_B}{T} \Big) + \pi^{(3/2)}\Big( \frac{\epsilon/k_B}{T} \Big)^{(3/2)} 

The pure species and mixture transport properties are evaluated with EGLib functions, which are linked directly into \pelelm.  EGLib requires as input polynomial fits of the logarithm of each quantity versus the logarithm of the temperature.

.. math::

    ln(q_m) = \sum_{n=1}^4 a_{q,m,n} ln(T)^{(n-1)} 

where :math:`q_m` represents :math:`\eta_m`, :math:`\lambda_m` or :math:`D_{m,j}`. These fits are generated as part of a preprocessing step managed by the tool `FUEGO` based on the formula (and input data) discussed above. The role of `FUEGO` to preprocess the model parameters for transport as well as chemical kinetics and thermodynamics, is discussed in some detail in <Section FuegoDescr>.


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

where :math:`\nu_{m,j} =\nu_{m,j}'' - \nu_{m,j}'`. Expressions for the reaction rates coefficients :math:`k_{(f,r),j}` depend on the type of reaction considered. \pelelm \; relies on the CHEMKIN Arrhenius reaction format:

.. math::

    k_f = AT^{\beta} exp \left( \frac{-E_a}{RT}\right)

where :math:`A` is the pre-exponential (frequency) factor, :math:`\beta` is the temperature exponent and :math:`E_a` is the activation energy. The CHEMKIN format additionally allows for a number of specializations of this format to represent pressure dependencies and third-body enhancements -- see the CHEMKIN Manual or Cantera website for additional information.

Most fundamental Arrhenius reactions are bidirectional, and typically only the forward rates are specified. In this case, the balance of forward and reverse rates are dictacted by equilibrium thermodynamics, via the equilibrium constant, :math:`K_{c,j}`.  In a low Mach system, :math:`K_{c,j}` is a function only of temperature and the thermodynamic properties of the reactants and products of reaction :math:`j`,

.. math::

    &&k_{r,j} = \frac{k_{f,j}}{K_{c,j}(T)} \;\;\; \mbox{where} \;\;\; K_{c,j}=K_{p,j} \left( \frac{P_{0}}{RT} \right)^{\sum_{k=1}^{N_s} \nu_{k,j}}\\
    &&\mbox{and} \;\;\; K_{p,j}=\exp \left( \frac{\Delta {S_j}^{0}}{R} - \frac{\Delta {H_j}^{0}}{RT} \right)

:math:`\Delta H_j` and :math:`\Delta S_j` are the change in enthalpy and entropy of the reaction :math:`j`, and :math:`P_0` is the ambient thermodynamic pressure.

Species production rates are evaluated via functions that are generated as part of a preprocessing step managed by the tool `FUEGO` (see <Section FuegoDescr>).

Thermodynamic properties
^^^^^^^^^^^^^^^^^^^^^^^^

Currently, expressions for the thermodynamic properties in \pelelm\ follow those of CHEMKIN, which assume a mixture of ideal gases. Species enthalpies and entropies are thus functions of only temperature (for perfect gases, they are independent of pressure) and are given in terms of polynomial fits to the species molar heat capacities (:math:`C_{p,\cdot}`),

.. math::

    \frac{C_{p,m}(T)}{\mathcal{R}} = \sum_{k=1}^{N_s} a_{k,m}T^{k-1}

where, in the standard CHEMKIN framework (the 7-coefficients NASA format), :math:`N =5`,

.. math::

    \frac{C_{p,m}(T)}{\mathcal{R}} = a_{1,m} + a_{2,m} T + a_{3,m} T^2 + a_{4,m} T^3 + a_{5,m} T^4

Accordingly, the standard-state molar enthalpy of species :math:`m` is given by:

.. math::

    \frac{H_{m}(T)}{\mathcal{R}T} = a_{1,m} +\frac{a_{2,m}}{2} T   + \frac{a_{3,m}}{3} T^2 +  \frac{a_{4,m}}{4} T^3 + \frac{ a_{5,m}}{5} T^4 + a_{6,m}/T

Note that the standard specifies that the heat of formation for the molecule is included in this expression.
Similarly, the standard-state molar entropy is written as:

.. math::

    \frac{S_{m}(T)}{\mathcal{R}} = a_{1,m}ln(T) + {a_{2,m}} T   + \frac{a_{3,m}}{2} T^2 +  \frac{a_{4,m}}{3} T^3 + \frac{ a_{5,m}}{4} T^4 + a_{7,m}

For each species, :math:`m`, in the model the user must specify the coefficients :math:`a_{k,m}`. All other required thermodynamic properties are then determined (see, e.g., the CHEMKIN manual for additional details. Thermodynamic properties of the species, and those of the mixture, are evaluated via functions that are generated as part of a preprocessing step managed by the tool `FUEGO` (see next <Section FuegoDescr>).


`FUEGO` chemistry preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical model for `PeleLM` contains all the information associated with the CHEMKIN parameterization of the Arrhenius reaction set, as well as fitting coefficients for the thermodynamic relationships, and the specification of the species including data required to compute pure-species transport properties. In the combustion community, this information is communicated for each complete model --or *mechanism*, through multiple text files that conform to the CHEMKIN standards. The CHEMKIN driver code (or equivalent) can then be used to ingest the large number of parameters contained in these files and provide a set of functions for evaluating all the properties and rates required.  Earlier versions of `PeleLM` linked to the CHEMKIN codes directly (and thereby assumed that all problems consisted of a mixture of ideal gases).  However, evaluations were not very efficient because the functions stepped through generic expressions that included a large number of conditional statements and unused generality.  Direct evaluation of these complex expressions allows for a much more efficient code that optimizes well with modern compilers. This is important because an appreciable fraction of `PeleLM` runtime is spent in these functions. Performance issues notwithstanding, customized evaluators will be necessary to extend `PeleLM` to a larger class of (*real*) gas models outside the CHEMKIN standard, such as SRK, that are already part of the `PeleC` code capabilities (`PeleC` shares use of `PelePhysics` for combustion model specification).

For these reasons, `PeleLM` no longer uses CHEMKIN functions directly, but instead relies on a preprocessing tool, `FUEGO`, to generate highly efficient C code implementations of the necessary thermodynamic, transport and kinetics evaluations.  The source code generated from `FUEGO` is linked into the `PeleLM` executable, customizing each executable for a specific model at compile time.  The implementation source code files can also be linked conveniently to post-processing analysis tools. The `FUEGO` processing tool, and the functions necessary to interface the generated functions to `PeleLM` are distributed in the auxiliary code package, `PelePhysics`.  Included in the `PelePhysics` distribution is a broad set of models for the combustion of hydrogen, carbon-monoxide, methane, heptane, :math:`n`-dodecane, dimethyl ether, and others, as well as instructions for users to extend this set using `FUEGO`, based on their own CHEMKIN-compliant inputs. `PelePhysics` also provides support for simpler *gama-law* equations-of-state, and simple/constant transport properties.


