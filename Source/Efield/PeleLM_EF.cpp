namespace EFConst{
    AMREX_GPU_DEVICE_MANAGED amrex::Real eps0 = 8.854187817e-12;          //Free space permittivity (C/(V.m))
    AMREX_GPU_DEVICE_MANAGED amrex::Real epsr = 1.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real elemCharge = 1.60217662e-19;     //Coulomb per charge
    AMREX_GPU_DEVICE_MANAGED amrex::Real Na = 6.022e23;                   //Avogadro's number
}
  
void PeleLM::ef_init() {
    amrex::Print() << " Init EFIELD solve options \n";  

    // Params defaults
    PeleLM::ef_verbose                = 0;
    PeleLM::nE                        = -1;
    PeleLM::PhiV                      = -1;
    PeleLM::ef_PoissonTol             = 1.0e-12;
    PeleLM::ef_PoissonMaxIter         = 100;
    PeleLM::ef_PoissonVerbose         = 0;
    PeleLM::ef_PoissonMaxOrder        = 2;

    for (int n = 0; n  < NUM_SPECIES; n++) 
    {
      PeleLM::zk[n] = 0.0;
    }

//  Abort if EB + Efield    
#ifdef AMREX_USE_EB    
    Abort("No support for EB + Efield yet");
#endif

    ParmParse pp("ef");

   // Get the phiV bc
   Vector<std::string> lo_bc(AMREX_SPACEDIM);
   Vector<std::string> hi_bc(AMREX_SPACEDIM);
   pp.getarr("phiV_lo_bc",lo_bc,0,AMREX_SPACEDIM);
   pp.getarr("phiV_hi_bc",hi_bc,0,AMREX_SPACEDIM);
   for (int i = 0; i < AMREX_SPACEDIM; i++) 
   {    
      if (!lo_bc[i].compare("Interior")){
         phiV_bc.setLo(i,0);
      } else if (!lo_bc[i].compare("Dirichlet")) { 
         phiV_bc.setLo(i,1);
      } else if (!lo_bc[i].compare("Neumann")) { 
         phiV_bc.setLo(i,2);
      } else {
         amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
      }
      if (!hi_bc[i].compare("Interior")){
         phiV_bc.setHi(i,0);
      } else if (!hi_bc[i].compare("Dirichlet")) { 
         phiV_bc.setHi(i,1);
      } else if (!hi_bc[i].compare("Neumann")) { 
         phiV_bc.setHi(i,2);
      } else {
         amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
      }
   }

   pp.query("verbose",ef_verbose);
   pp.query("Poisson_tol",ef_PoissonTol);
   pp.query("Poisson_maxiter",ef_PoissonMaxIter);
   pp.query("Poisson_verbose",ef_PoissonVerbose);
   pp.query("Poisson_maxorder",ef_PoissonMaxOrder);
 }

void PeleLM::ef_solve_phiv(const Real &time) {
   BL_PROFILE("PLM_EF::ef_solve_phiv()");

   if ( ef_verbose ) amrex::Print() << " Solve for electro-static field on level " << level << "\n";

   MultiFab& S = (time == state[State_Type].prevTime() ) ? get_old_data(State_Type) 
                                                         : get_new_data(State_Type);

   MultiFab rhs_poisson(grids,dmap,1,1);

// Get alias to PhiV in S     
   MultiFab PhiV_alias(S,amrex::make_alias,PhiV,1);

// Get RHS: charge distribution
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhs    = rhs_poisson.array(mfi);
      auto const& rhoY   = S.const_array(mfi,first_spec);
      auto const& nE_arr = S.const_array(mfi,nE);
      Real        factor = -1.0 / ( EFConst::eps0  * EFConst::epsr);
      amrex::ParallelFor(bx, [rhs, rhoY, nE_arr, factor]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         rhs(i,j,k) = - nE_arr(i,j,k) * EFConst::elemCharge * factor;
         for (int n = 0; n < NUM_SPECIES; n++) {
            rhs(i,j,k) += zk[n] * rhoY(i,j,k,n) * factor;
         }
      });
   }

// Set-up solver tolerances
   const Real S_tol     = ef_PoissonTol; 
   const Real S_tol_abs = std::max(rhs_poisson.norm0(),PhiV_alias.norm0()) * ef_PoissonTol;

// Set-up Poisson solver
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);

   MLPoisson phiV_poisson({geom}, {grids}, {dmap}, info);   
   phiV_poisson.setMaxOrder(ef_PoissonMaxOrder);

// Set-up BC's
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
   ef_set_PoissonBC(mlmg_lobc, mlmg_hibc);
   phiV_poisson.setDomainBC(mlmg_lobc, mlmg_hibc);

   MultiFab phiV_crse;
   if (level > 0) {
      auto& crselev = getLevel(level-1);
      phiV_crse.define(crselev.grids, crselev.dmap, 1, 0, MFInfo(), crselev.Factory());
      FillPatch(crselev,phiV_crse,0,time, State_Type, PhiV, 1, 0);   
      phiV_poisson.setCoarseFineBC(&phiV_crse, crse_ratio[0]);
   }
   MultiFab PhiVmf(grids,dmap,1,1,MFInfo(),Factory());
   FillPatch(*this,PhiVmf,1,time,State_Type,PhiV,1,0);
   phiV_poisson.setLevelBC(0, &PhiVmf);

// LinearSolver options
   MLMG mlmg(phiV_poisson);
   mlmg.setMaxIter(ef_PoissonMaxIter);
   mlmg.setVerbose(ef_PoissonVerbose);

   PhiV_alias.mult(0.0,0,1,0);
   mlmg.solve({&PhiV_alias}, {&rhs_poisson}, S_tol, S_tol_abs);

// TODO: will need the flux later
}

void PeleLM::ef_define_data() {
   BL_PROFILE("PLM_EF::ef_define_data()");

   Ke_cc.define(grids,dmap,1,1);
   De_cc.define(grids,dmap,1,1);
}

void PeleLM::ef_advance_setup(const Real &time) {
   // Solve for phiV_time and get gradPhiV_tn
   ef_solve_phiv(time);

   // Calc De_time, Ke_time, Kspec_time
   ef_calc_transport(time);

   // Calc Udrift: here phiV_tnp1 = phiV_tn
}

void PeleLM::ef_calc_transport(const Real &time) {
   BL_PROFILE("PLM_EF::ef_calc_transport()");

   if ( ef_verbose ) amrex::Print() << " Compute EF transport prop. \n";

   const TimeLevel whichTime = which_time(State_Type, time);

   BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);

   // Fillpatch the state
   FillPatchIterator fpi(*this,S,Ke_cc.nGrow(),time,State_Type,first_spec,NUM_SPECIES+3);
   MultiFab& S_cc = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& gbx = mfi.growntilebox();
      auto const& rhoY = S_cc.array(mfi,0);
      auto const& T    = S_cc.array(mfi,NUM_SPECIES+1);
      auto const& Ke   = Ke_cc.array(mfi);
      auto const& De   = De_cc.array(mfi);
      Real factor = PP_RU_MKS / ( EFConst::Na * EFConst::elemCharge );
      amrex::ParallelFor(gbx, [rhoY, T, factor, Ke, De]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         getKappaE(i,j,k,Ke);
         getDiffE(i,j,k,factor,T,Ke,De);
      });
   }
}

// Setup BC conditions for linear Poisson solve on PhiV. Directly copied from the diffusion one ...
void PeleLM::ef_set_PoissonBC(std::array<LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                              std::array<LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc) {

    const BCRec& bc = get_desc_lst()[State_Type].getBC(PhiV);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geom().isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}
