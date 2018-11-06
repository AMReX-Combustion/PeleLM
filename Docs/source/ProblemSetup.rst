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
the user's initialization function:

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

The associated fortran routines must shape the data accordingly:

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

Note that in this example, the Fortran is preprocessed to define ``SDIM``, ``REAL_T``, ``DIMS``, ``DIMV``
convenience macros (TODO: clean this up...). Also, the user will need to know the order and definition of
each state component.  The variables, and their order, are set in the function:

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
in the initial and boundary condition routines.  Typically, however, the primitive state is known instead
(i.e., :math:`Y_i` and :math:`T`).  If that is the case, the user can make use of the compiled-in
equation-of-state routines to translate primitive to conserved state values. Consult the example setups
provided to see how to call these routines, and load the final values of the initial data.

Boundary Conditions
-------------------

In `PeleLM` separate functions are provided to fill each state component at physical boundaries.
The ``variableSetup`` routine discussed above sets, for each state, which function will be called,
but all of them have the same form/aguments.  A box of data will be provided with some overlap
of the valid domain, and some overlap of the grow region outside the domain.  The region of
index space defining the domain is level-specific, and so is passed directly to the boundary
function, as is the time, the grid spacing, and an 2 x D array indicating the numerical boundary
condition to apply (adapted from the ``inputs`` file parameters of the run). The task of this
routine is to set values in the grow cells of the input array accordingly.  Generally, this is
done by first calling a utility function, ``filcc`` that can fill grow cells for all of the boundary
condition types, **except** ``EXT_DIR`` (external Dirichlet) -- those must be set directly by the
user.  (so, ``filcc`` handles reflecting even/odd, extrapolation, etc).  Below, we include
an example of typical logic for carrying this out.  First ``filcc`` is called, and then each
boundary orientation is checked for whether the Dirichlet conditions need to be applied.  If so,
corresponding values are set.  Here, we've made use of a local convenience function, ``bcfunction``
endowed with the knowledge of all boundary values, and extract the appropriate quantity from the
results of that call.  This was done to localize all boundary condition calculations to a single
routine in the code, and helps to preserve consistency.  This is only one style though, and
as long as appropriate Dirichlet values are set for this state, it makes no difference how the
work is organized.

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

Refinement Criteria
-------------------

TODO: Add junk here for refinement criteria.
