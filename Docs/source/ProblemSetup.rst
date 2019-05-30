.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Setting up a new `PeleLM` Case
==============================

In order to set up and run a new case in `PeleLM`, the user must provide problem-specific code for three main tasks

- Initial conditions
- Boundary conditions
- Grid refinement criteria (optional, provided default functions are discussed below)

These functions are typically collected into a single subfolder in ``${PELELM_HOME}/Exec``, such as ``FlameInABox``.
The user can organize these tasks in any way that is convenient - the examples distributed with ``PeleLM``
represent a certain style for managing this with some level of flexibility, but the basic requirement is
simply that source be linked into the build for the functions ``FORT_INITDATA`` (for initial conditions)
and ``FORT_XXXFILL`` (where ``XXX`` is ``DEN``, ``TEMP``, ``ADV``, ``RHOH``, ``VEL``, ``CHEM``, ``PRES``)
for boundary conditions.

Initial Conditions
------------------

At the beginning of a `PeleLM` run, for each level, after grids are generated, the cell-centered values of
the state must be initialized.  In the code, this is done in an ``MFIter`` loop over grids, and a call to
the user's initialization function, ``FORT_INITDATA``, that must be provided:

::

  // Initialize phi_new by calling a Fortran routine.
  // MFIter = MultiFab Iterator
  for (MFIter snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
  {
    const int  i       = snewmfi.index();
    RealBox    gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
    const Box& vbx     = snewmfi.validbox();
    const int* lo      = vbx.loVect();
    const int* hi      = vbx.hiVect();
    const int* s_lo    = S_new[snewmfi].loVect();
    const int* s_hi    = S_new[snewmfi].hiVect();
    const int* p_lo    = P_new[snewmfi].loVect();
    const int* p_hi    = P_new[snewmfi].hiVect();

    FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                   S_new[snewmfi].dataPtr(Xvel),
                   S_new[snewmfi].dataPtr(BL_SPACEDIM),
                   ARLIM(s_lo), ARLIM(s_hi),
                   P_new[snewmfi].dataPtr(),
                   ARLIM(p_lo), ARLIM(p_hi),
                   dx,gridloc.lo(),gridloc.hi() );
  }

The associated fortran routine must shape the arrays accordingly, and fill the data:

::

      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
     &                         vel,scal,DIMS(state),press,DIMS(press),
     &                         delta,xlo,xhi)
      implicit none
      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     time, delta(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
      integer i, j, n

      press(:,:) = 0.d0
      do j = lo(2), hi(2)
         y = (float(j)+.5d0)*delta(2)+domnlo(2)
         do i = lo(1), hi(1)
            x = (float(i)+.5d0)*delta(1)+domnlo(1)
            scal(i,j,1:n) = ...
            vel(i,j,1:2) = ...
         enddo
      enddo
      end

Note that in this f77 example, the Fortran is preprocessed to define ``SDIM``, ``REAL_T``, ``DIMS``, ``DIMV``
convenience macros (TODO: clean this up...). Also, the user will need to know the order and definition of
each state component.  The variables, and their order, are set for all problems setups in the function:

::

  void PeleLM::variableSetUp ()
  {
    ...
    FirstSpec = ++counter;
    nspecies  = getChemSolve().numSpecies();
    counter  += nspecies - 1;
    RhoH = ++counter;
    Trac = ++counter;
    Temp = ++counter;
    RhoRT = ++counter;
    NUM_STATE = ++counter;
    NUM_SCALARS = NUM_STATE - Density;
    ...
    
in the ``$(PELELM_HOME)/Source`` folder.  Note that the conserved states are stored for species and
enthalpy (i.e., :math:`\rho Y_i` and :math:`\rho h`); these are the variables that the user must fill
in the initial and boundary condition routines.  Typically, however, the primitive state
(i.e., :math:`Y_i` and :math:`T`) is known directly.  If that is the case, the user can make use of
the compiled-in model-specific equation-of-state routines to translate primitive to conserved state
values. Consult the example setups provided to see how to call these routines, and how to load the
final values required for initial data.

Boundary Conditions
-------------------

In `PeleLM` separate functions are provided to fill each state
component at physical boundaries.  The ``variableSetup`` routine
discussed above sets, for each state, which function will be called,
but all of them have the same form/aguments.  A box of data will be
provided with some overlap of the valid domain, and some overlap of
the grow region outside the domain.  The region of index space
defining the domain is level-specific, and so is passed directly to
the boundary function, as is the time, the grid spacing, and an 2 x D
array indicating the numerical boundary condition to apply (adapted
from the ``inputs`` file parameters of the run). The task of this
routine is to set values in the grow cells of the input array
accordingly.  Generally, this is done by first calling a utility
function, ``filcc``, that can fill grow cells for all of the boundary
condition types, **except** ``EXT_DIR`` (external Dirichlet) --
Dirichlet values must be set directly by the user.  Below, we include
an example of typical logic for carrying this out.  First ``filcc`` is
called, and then each boundary orientation is checked for whether the
Dirichlet conditions need to be applied.  If so, corresponding values
are set.  Here, we've made use of a local convenience function,
``bcfunction`` endowed with the knowledge of all boundary values, and
extract the appropriate quantity from the results of that call.  This
was done to localize all boundary condition calculations to a single
routine in the code, and helps to preserve consistency.  This is only
one style though, and as long as appropriate Dirichlet values are set
for this state, it makes no difference how the work is organized.
For example, data may be provided by interpolating "live data" being
actively generated by a co-running separate code, by interpolating data
files, evaluating functional forms, etc.

::

      subroutine FORT_DENFILL (den,DIMS(den),domlo,domhi,delta,xlo,time,bc)
      implicit none

      integer DIMDEC(den), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  den(DIMV(den))

      include 'cdwrk.H'
      include 'bc.H'
      include 'probdata.H`
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      lo(2) = ARG_L2(den)
      hi(1) = ARG_H1(den)
      hi(2) = ARG_H2(den)

      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+domnlo(2)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (float(j)+.5)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (float(i)+.5)*delta(1)+domnlo(1)
               call bcfunction(x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif

      end

Note that although the array structure to be filled contains valid cell-centered state data where it
overlaps the valid domain, the values set in the grow cells of the container will be applied on the
boundary face of the corresponding cells.  Internally, all `PeleLM` code understands to apply
Dirichlet conditions on the boundary faces.

.. _sec:refcrit:pelelm:

Refinement Criteria
-------------------

The dynamic creation and destruction of grid levels is a fundamental part of `PeleLM`'s capabilities. The
process for this is described in some detail in the `AMReX` documentation, but we summarize the key points
here.

At regular intervals (set by the user), each Amr level that is not the finest allowed for the run
will invoke a "regrid" operation.  When invoked, a list of error tagging functions is traversed. For each,
a field specific to that function is derived from the state over the level, and passed through a kernel
that "set"'s or "clear"'s a flag on each cell.  The field and function for each error tagging quantity is
identified in the setup phase of the code where the state descriptors are defined (i.e., in `PeleLM_setup.cpp`).
Each function in the list adds or removes to the list of cells tagged for refinement. This final list of tagged
cells is sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm to create rectangular grids
which will define a new finer level (or set of levels).  State data is filled over these new grids, copying where
possible, and interpolating from coarser level when no fine data is available.  Once this process is complete,
the existing Amr level(s) is removed, the new one is inserted into the hierarchy, and the time integration
continues.

The traditional `AMReX` approach to setting up and controlling the regrid process involves explicitly
creating ("hard coding") a number of functions directly into `PeleLM`'s setup code. (Consult the source code
and `AMReX` documentation for precisely how this is done).  `PeleLM` provides a limited capability to augment
the standard set of error functions that is based entirely on runtime data specified in the inputs (ParmParse)
data.  The following example portion of a ParmParse'd input file demonstrates the usage of this feature:

::

      amr.refinement_indicators = flame_tracer lo_temp gradT

      amr.flame_tracer.max_level = 3
      amr.flame_tracer.value_greater = 1.e-6
      amr.flame_tracer.field_name = Y(H)

      amr.lo_temp.max_level = 1
      amr.lo_temp.value_less = 450
      amr.lo_temp.field_name = temp

      amr.gradT.max_level = 2
      amr.gradT.adjacent_difference_greater = 20
      amr.gradT.field_name = temp
      amr.gradT.start_time = 0.001
      amr.gradT.end_name = 0.002

Here, we have added three new custom-named criteria -- ``flame_tracer``: cells with the mass fraction of H greater than 1 ppm;
``lo_temp``: cells with T less than 450K, and ``gradT``: cells having a temperature difference of 20K from that of their
immediate neighbor.  The first will trigger up to Amr level 3, the second only to level 1, and the third to level 2.
The third will be active only when the problem time is between 0.001 and 0.002 seconds.

Note that these additional user-created criteria operate in addition to those defined as defaults.  Also note that
these can be modified between restarts of the code.  By default, the new criteria will take effect at the next
schedule regrid operation.  Alternatively, the user may restart with ``amr.regrid_on_restart = 1`` in order to
do a full (all-levels) regrid after reading the checkpoint data and before advancing any cells.


