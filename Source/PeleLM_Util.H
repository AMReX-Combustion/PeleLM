#ifndef PeleLM_Util_H 
#define PeleLM_Util_H

#include <AMReX_MultiFab.H>

namespace amrex
{
  amrex::Vector<amrex::Real>
  VectorMax(const amrex::Vector<const amrex::MultiFab *>& mf,
            const amrex::IntVect&                         tilesize,
            int                                           sComp,
            int                                           nComp);

  amrex::Vector<amrex::Real>
  VectorMin(const amrex::Vector<const amrex::MultiFab *>& mf,
            const amrex::IntVect&                         tilesize,
            int                                           sComp,
            int                                           nComp);

  amrex::Vector<amrex::Real>
  VectorMaxAbs(const amrex::Vector<const amrex::MultiFab *>& mf,
               const amrex::IntVect&                         tilesize,
               int                                           sComp,
               int                                           nComp);
}

#endif
