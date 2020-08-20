#include <GMRES.H>

using namespace amrex;

GMRESSolver::GMRESSolver() {}
GMRESSolver::~GMRESSolver() {}

void
GMRESSolver::define(const Geometry& a_geom,
                    const BoxArray& a_grids,
                    const DistributionMapping& a_dmap,
                    const int a_KrylovSize, 
                    const int a_nComp,
                    const int a_nGrow)
{
   BL_PROFILE("GMRESSolver::define()");

   m_geom = a_geom;
   m_grids = a_grids;
   m_dmap = a_dmap;
   m_nComp = a_nComp;
   m_nGrow = a_nGrow;
   m_krylovSize = a_KrylovSize;

// Build krylov base memory
   KspBase.resize(m_krylovSize+1);
   for (int n = 0; n <= m_krylovSize ; ++n) {
      KspBase[n].define(m_grids,m_dmap,m_nComp,m_nGrow);
   }

// Work MultiFabs
   Ax.define(m_grids,m_dmap,m_nComp,m_nGrow);
   res.define(m_grids,m_dmap,m_nComp,m_nGrow);
}

JtimesVFunc
GMRESSolver::jtimesv() const noexcept
{
   return m_jtv;
}

PrecondFunc
GMRESSolver::precond() const noexcept
{
   return m_prec;
}

NormFunc
GMRESSolver::norm() const noexcept
{
   return m_norm;
}

void
GMRESSolver::setJtimesV(JtimesVFunc a_jtv)
{
   m_jtv = a_jtv;
}

void
GMRESSolver::setPrecond(PrecondFunc a_prec)
{
   m_prec = a_prec;
}

void
GMRESSolver::setNorm(NormFunc a_norm)
{
   m_norm = a_norm;
}

int
GMRESSolver::solve(MultiFab& a_sol,
                   const MultiFab& a_rhs,
                   amrex::Real a_abs_tol,
                   amrex::Real a_rel_tol)
{
   BL_PROFILE("GMRESSolver::solve()");

// Checks
   AMREX_ALWAYS_ASSERT(m_jtv != nullptr);
   AMREX_ALWAYS_ASSERT(a_sol.nComp() == m_nComp && a_rhs.nComp() == m_nComp);

   Real beta = computeResidualNorm(a_sol,a_rhs);                      // Initial resisual
   Real rel_tol = beta * a_rel_tol;                                   // Target relative tolerance

   int GMRES_count = 0;
   int GMRES_restart = 0;
   do {
      GMRES_count += one_restart(a_sol,a_rhs);
      GMRES_restart++;
   } while( !m_converged && GMRES_restart < m_restart );
   return GMRES_count;
}

int
GMRESSolver::one_restart(MultiFab& a_x, const MultiFab& a_rhs)
{
   Real H[m_krylovSize+1][m_krylovSize] = {0.0};     // Hessenberg matrix
   Real y[m_krylovSize+1] = {0.0};                   // Solution vector
   Real g[m_krylovSize+1] = {0.0};                   // Residual
   Real givens[m_krylovSize+1][2] = {0.0};           // Givens rotations

   return m_krylovSize;
}

Real
GMRESSolver::computeResidualNorm(const MultiFab& a_x, const MultiFab& a_rhs)
{
   computeResidual(a_x,a_rhs,res);
   Real resNorm = 0.0;
   if ( m_norm != nullptr ) {
      m_norm(res,resNorm);
   } else {
      for ( int comp = 0; comp < m_nComp; comp++ ) {
         resNorm += MultiFab::Dot(res,comp,res,comp,1,0);
      }
      resNorm = std::sqrt(resNorm);
   }
   return resNorm;
}

void
GMRESSolver::computeResidual(const MultiFab& a_x, const MultiFab& a_rhs, MultiFab& a_res)
{
   BL_PROFILE("GMRESSolver::computeResidual()");
   m_jtv(a_x,a_rhs,Ax);
   MultiFab::Xpay(Ax, -1.0, a_rhs, 0, 0, m_nComp, 0);
   if ( m_prec != nullptr ) m_prec(Ax,a_res);
}
