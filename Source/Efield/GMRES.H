#ifndef PLM_GMRES_H_
#define PLM_GMRES_H_
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

typedef void (*JtimesVFunc) (const amrex::MultiFab& a_x, const amrex::MultiFab& a_rhs, amrex::MultiFab& a_Ax);
typedef void (*PrecondFunc) (const amrex::MultiFab& v, amrex::MultiFab& Pv);
typedef void (*NormFunc)    (const amrex::MultiFab& MFin, amrex::Real norm);

class GMRESSolver
{
   public:
     GMRESSolver();
     virtual ~GMRESSolver();

     GMRESSolver (const GMRESSolver&) = delete;
     GMRESSolver& operator= (const GMRESSolver&) = delete;
     GMRESSolver& operator= (GMRESSolver&&) = delete;

     void define(const amrex::Geometry& a_geom,
                 const amrex::BoxArray& a_grids,
                 const amrex::DistributionMapping& a_dmap,
                 const int a_KrylovSize,
                 const int a_nComp,
                 const int a_nGrow);

     void setJtimesV(JtimesVFunc a_jtv);
     void setPrecond(PrecondFunc a_prec);
     void setNorm(NormFunc a_norm);

     JtimesVFunc jtimesv() const noexcept;
     PrecondFunc precond() const noexcept;
     NormFunc    norm() const noexcept;

     int solve(amrex::MultiFab& sol,
               const amrex::MultiFab& rhs,
               amrex::Real a_abs_tol,
               amrex::Real a_rel_tol);

   private:

     int one_restart(amrex::MultiFab& a_x,
                     const amrex::MultiFab& a_rhs);

     void computeResidual(const amrex::MultiFab& a_x,
                          const amrex::MultiFab& a_rhs,
                          amrex::MultiFab& a_res);

     amrex::Real computeResidualNorm(const amrex::MultiFab& a_x,
                                     const amrex::MultiFab& a_rhs);

     int m_verbose;

//   GMRES attributes
     bool m_converged;
     int m_krylovSize;
     int m_restart;
     int m_nComp;

//   GMRES data
     amrex::Vector<amrex::MultiFab> KspBase;
     amrex::MultiFab Ax;
     amrex::MultiFab res;

//   AMReX attributes
     int m_nGrow;
     amrex::Geometry m_geom;
     amrex::BoxArray m_grids;
     amrex::DistributionMapping m_dmap;

//   GMRES user-provided functions
     JtimesVFunc m_jtv = nullptr;     // JtimesV function. Always needed
     PrecondFunc m_prec = nullptr;    // Preconditioner
     NormFunc m_norm = nullptr;       // User-prefered MF norm. L2 used if not provided
};

#endif
