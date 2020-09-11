#ifndef PLM_GMRES_H_
#define PLM_GMRES_H_
#include <AMReX_REAL.H>
#include <PeleLM.H>
#include <AMReX_MultiFab.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

#define MEMBER_FUNC_PTR(object,ptrToMember)  ((object).*(ptrToMember))

typedef void (PeleLM::*JtimesVFunc) (const amrex::MultiFab& a_x, amrex::MultiFab& a_Ax);
typedef void (PeleLM::*PrecondFunc) (const amrex::MultiFab& v, amrex::MultiFab& Pv);
typedef void (PeleLM::*NormFunc)    (const amrex::MultiFab& MFin, amrex::Real& norm);

class GMRESSolver
{
   public:
     GMRESSolver();
     virtual ~GMRESSolver();

     GMRESSolver (const GMRESSolver&) = delete;
     GMRESSolver& operator= (const GMRESSolver&) = delete;
     GMRESSolver& operator= (GMRESSolver&&) = delete;

     void define(PeleLM* a_level,
                 const int a_KrylovSize,
                 const int a_nComp,
                 const int a_nGrow);

     void setJtimesV(JtimesVFunc a_jtv);
     void setPrecond(PrecondFunc a_prec);
     void setNorm(NormFunc a_norm);

     void setVerbose(int a_v){ m_verbose = a_v;};
     void setMaxRestart(int a_maxRestart){ m_restart = a_maxRestart;};

     JtimesVFunc jtimesv() const noexcept;
     PrecondFunc precond() const noexcept;
     NormFunc    norm() const noexcept;

     int solve(amrex::MultiFab& sol,
               const amrex::MultiFab& rhs,
               amrex::Real a_abs_tol,
               amrex::Real a_rel_tol);

   private:

     void one_restart(amrex::MultiFab& a_x,
                      const amrex::MultiFab& a_rhs);

     void one_iter(const int iter, amrex::Real &resNorm);

     void appendBasisVector(const int iter, amrex::Vector<amrex::MultiFab>& Base);

     void gramSchmidtOrtho(const int iter, amrex::Vector<amrex::MultiFab>& Base);

     amrex::Real givensRotation(const int iter);

     void prepareForSolve();

     void computeResidual(const amrex::MultiFab& a_x,
                          const amrex::MultiFab& a_rhs,
                          amrex::MultiFab& a_res);

     amrex::Real computeResidualNorm(const amrex::MultiFab& a_x,
                                     const amrex::MultiFab& a_rhs);

     amrex::Real computeNorm(const amrex::MultiFab& a_vec);

     int m_verbose;

//   GMRES attributes
     int check_GramSchmidtOrtho = 1;
     bool m_converged;
     int m_krylovSize;
     int m_restart;
     int m_nComp;
     amrex::Real initResNorm;
     amrex::Real targetResNorm;
     int iter_count;
     int restart_count;

//   GMRES data
//   MultiFab data
     amrex::Vector<amrex::MultiFab> KspBase;  // Krylov basis
     amrex::MultiFab Ax;                      // A container for A*x
     amrex::MultiFab res;                     // A container for residual

//   Real data
     amrex::Vector<amrex::Vector<amrex::Real>> H;        // Hessenberg matrix
     amrex::Vector<amrex::Real>                y;        // Solution vector
     amrex::Vector<amrex::Real>                g;        // Residual
     amrex::Vector<amrex::Vector<amrex::Real>> givens;   // Givens rotation

//   AMReX attributes
     PeleLM* m_level;
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
