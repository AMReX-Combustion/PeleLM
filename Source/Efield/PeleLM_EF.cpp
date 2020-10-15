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
    PeleLM::ef_debug                  = 0;
    PeleLM::ef_substep                = 1;
    PeleLM::nE                        = -1;
    PeleLM::PhiV                      = -1;
    PeleLM::ef_PoissonTol             = 1.0e-12;
    PeleLM::ef_PoissonMaxIter         = 100;
    PeleLM::ef_PoissonVerbose         = 0;
    PeleLM::ef_PoissonMaxOrder        = 2;
    PeleLM::ef_use_PETSC_direct       = 0;
    PeleLM::ef_lambda_jfnk            = 1.0e-7;
    PeleLM::ef_diffT_jfnk             = 1;
    PeleLM::ef_maxNewtonIter          = 10;
    PeleLM::ef_newtonTol              = std::pow(1.0e-13,2.0/3.0);
    PeleLM::ef_GMRES_size             = 20;
    PeleLM::ef_GMRES_reltol           = 1.0e-10;
    PeleLM::ef_GMRES_maxRst           = 1;
    PeleLM::ef_GMRES_verbose          = 0;
    PeleLM::ef_PC_fixedIter           = -1;
    PeleLM::ef_PC_MG_Tol              = 1.0e-6;

    // Get the charge per unit mass
    Real zk[NUM_SPECIES];
    EOS::charge_mass(zk);
    for (int k = 0; k < NUM_SPECIES; k++) {
       PeleLM::zk[k] = zk[k];
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
   pp.query("debug",ef_debug);
   pp.query("Poisson_tol",ef_PoissonTol);
   pp.query("Poisson_maxiter",ef_PoissonMaxIter);
   pp.query("Poisson_verbose",ef_PoissonVerbose);
   pp.query("Poisson_maxorder",ef_PoissonMaxOrder);
   pp.query("JFNK_Substeps",ef_substep);
   pp.query("JFNK_lambda",ef_lambda_jfnk);
   pp.query("JFNK_difftype",ef_diffT_jfnk);
   pp.query("JFNK_maxNewton",ef_maxNewtonIter);
   pp.query("JFNK_newtonTol",ef_newtonTol);
   pp.query("GMRES_restart_size",ef_GMRES_size);
   pp.query("GMRES_rel_tol",ef_GMRES_reltol);
   pp.query("GMRES_max_restart",ef_GMRES_maxRst);
   pp.query("GMRES_verbose",ef_GMRES_verbose);
   pp.query("Precond_MG_tol",ef_PC_MG_Tol);
   pp.query("Precond_fixedIter",ef_PC_fixedIter);
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

   diff_e.define(this);
   De_ec = diff_e.get();

   mob_e.define(this);
   Ke_ec = mob_e.get();

   KSpec_old.define(grids,dmap,NUM_SPECIES,1);
   KSpec_new.define(grids,dmap,NUM_SPECIES,1);

   ef_state_old.define(grids,dmap,2,1);
   ef_state_refGhostCell.define(grids,dmap,2,2);
   bg_charge.define(grids,dmap,1,1);
   nl_state.define(grids,dmap,2,2);
   nl_resid.define(grids,dmap,2,2);

   elec_Ueff = new MultiFab[AMREX_SPACEDIM];
   for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      const BoxArray& edgeba = getEdgeBoxArray(d);
      elec_Ueff[d].define(edgeba, dmap, 1, 1);
   }

   gphiV_fb.define(this);
   gphiV_old = gphiV_fb.get();

   Udrift_spec = new MultiFab[AMREX_SPACEDIM];
   for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      const BoxArray& edgeba = getEdgeBoxArray(d);
      Udrift_spec[d].define(edgeba, dmap, NUM_SPECIES, 1);
      Udrift_spec[d].setVal(0.0);
   }
}

void PeleLM::ef_advance_setup(const Real &time) {
   // Solve for phiV_time and get gradPhiV_tn
   ef_solve_phiv(time);
}

void PeleLM::ef_calc_transport(const Real &time) {
   BL_PROFILE("PLM_EF::ef_calc_transport()");

   if ( ef_verbose ) amrex::Print() << " Compute EF transport prop.\n";

   const TimeLevel whichTime = which_time(State_Type, time);

   BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab& S     = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
   MultiFab& diff  = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
   MultiFab& Kspec = (whichTime == AmrOldTime) ? KSpec_old : KSpec_new;

   // Get the cc transport coeffs. These are temporary.
   MultiFab Ke_cc(grids,dmap,1,1);
   MultiFab De_cc(grids,dmap,1,1);

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
      auto const& rhoD = diff.array(mfi);
      auto const& Ke   = Ke_cc.array(mfi);
      auto const& De   = De_cc.array(mfi);
      auto const& Ks   = Kspec.array(mfi);
      Real factor = PP_RU_MKS / ( EFConst::Na * EFConst::elemCharge );
      amrex::ParallelFor(gbx, [rhoY, T, factor, Ke, De]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         getKappaE(i,j,k,Ke);
         getDiffE(i,j,k,factor,T,Ke,De);
      });
      Real mwt[NUM_SPECIES];
      EOS::molecular_weight(mwt);
      amrex::ParallelFor(gbx, [rhoY, rhoD, T, Ks, mwt]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         getKappaSp(i,j,k, mwt, zk, rhoY, rhoD, T, Ks);
      });
   }
   if ( ef_debug ) {
      std::string timetag = (whichTime == AmrOldTime) ? "old" : "new";
      VisMF::Write(Kspec,"KappaSpec"+timetag+"_Lvl"+std::to_string(level));
   }

   // CC -> EC transport coeffs. These are PeleLM class object used in the non-linear residual.
   auto math_bc = fetchBCArray(State_Type,nE,1);
   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(De_cc,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         const Box ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         const auto& diff_c  = De_cc.array(mfi);
         const auto& diff_ed = De_ec[dir]->array(mfi);
         const auto& kappa_c  = Ke_cc.array(mfi);
         const auto& kappa_ed = Ke_ec[dir]->array(mfi);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, use_harmonic_avg, diff_c, diff_ed,
                                  kappa_c, kappa_ed, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, 1, use_harmonic_avg, on_lo, on_hi, diff_c, diff_ed);
            cen2edg_cpp( i, j, k, dir, 1, use_harmonic_avg, on_lo, on_hi, kappa_c, kappa_ed);
         });
      }
   }
}

void PeleLM::jtimesv(const MultiFab &v,
                           MultiFab &Jv)
{
   Real vNorm;
   ef_normMF(v,vNorm);

   // x is zero, Ax is zero and return
   if ( vNorm == 0.0 ) {
      Jv.setVal(0.0);
      return;
   }

   Real delta_pert = ef_lambda_jfnk * ( ef_lambda_jfnk + nl_stateNorm / vNorm );

   if ( ef_diffT_jfnk == 1 ) {
      // Create perturbed state
      MultiFab S_pert(grids,dmap,2,2);
      MultiFab::Copy(S_pert, nl_state, 0, 0, 2, 2);
      MultiFab::Saxpy(S_pert,delta_pert,v, 0, 0, 2 ,0);

      // Get perturbed residual
      MultiFab res_pert(grids,dmap,2,1);
      ef_nlResidual(dtsub,S_pert,res_pert);
      res_pert.mult(-1.0);

      // Get Ax by finite differece
      MultiFab::LinComb(Jv,1.0,res_pert,0,-1.0,nl_resid,0,0,2,0);
      Jv.mult(-1.0/delta_pert);
   } else if ( ef_diffT_jfnk == 2 ) {
      // Create perturbed states
      MultiFab S_pertm(grids,dmap,2,2);
      MultiFab S_pertp(grids,dmap,2,2);
      MultiFab::Copy(S_pertp, nl_state, 0, 0, 2, 2);
      MultiFab::Copy(S_pertm, nl_state, 0, 0, 2, 2);
      MultiFab::Saxpy(S_pertp,delta_pert,v, 0, 0, 2 ,0);
      MultiFab::Saxpy(S_pertm,-delta_pert,v, 0, 0, 2 ,0);

      // Get perturbed residuals
      MultiFab res_pertp(grids,dmap,2,1);
      MultiFab res_pertm(grids,dmap,2,1);
      ef_nlResidual(dtsub,S_pertp,res_pertp);
      ef_nlResidual(dtsub,S_pertm,res_pertm);
      res_pertm.mult(-1.0);
      res_pertp.mult(-1.0);

      // Get Ax by finite differece
      MultiFab::LinComb(Jv,1.0,res_pertp,0,-1.0,res_pertm,0,0,2,0);
      Jv.mult(-0.5/delta_pert);
   } else {
      Abort(" Unrecognized ef_diffT_jfnk. Should be either 1 (one-sided) or 2 (centered)");
   }

}

void PeleLM::ef_solve_PNP(      int      misdc,
                          const Real     &dt,
                          const Real     &time,
                          const MultiFab &Dn,
                          const MultiFab &Dnp1,
                          const MultiFab &Dhat,
                                MultiFab &ForcingnE)
{
   BL_PROFILE("PLM_EF::ef_solve_PNP()");

   const Real strt_time = ParallelDescriptor::second();

   // Substepping of non-linear solve
   dtsub = dt/ef_substep;

   // Get a copy of the old nE and PhiV into ef_state_old
   Real prev_time = state[State_Type].prevTime();
   FillPatchIterator nEfpi(*this,ef_state_refGhostCell,2,prev_time,State_Type,nE,2);
   MultiFab& nEfpi_mf = nEfpi.get_mf();
   MultiFab::Copy(ef_state_old,nEfpi_mf,0,0,2,1);                            // Copy into a storage to hold the old state
   MultiFab::Copy(ef_state_refGhostCell,nEfpi_mf,0,0,2,2);  // Copy into a storage to hold the reference ghost cell value

   // Need the gradient of old phiV
   {
      MultiFab phi_a(ef_state_old,amrex::make_alias,1,1);
      ef_calcGradPhiV(prev_time,phi_a,gphiV_old);
   }

   // Non-linear state & residual from FPI
   MultiFab::Copy(nl_state, ef_state_refGhostCell, 0, 0, 2, 2);
   if ( ef_debug ) VisMF::Write(nl_state,"PrePNPState_"+std::to_string(level));

   // GMRES
   GMRESSolver gmres;
   int GMRES_tot_count = 0;
   if ( !ef_use_PETSC_direct ) {
      gmres.define(this,ef_GMRES_size,2,1);        // 2 component in GMRES, 1 GC (needed ?)
      JtimesVFunc jtv = &PeleLM::jtimesv;
      gmres.setJtimesV(jtv);
      NormFunc normF = &PeleLM::ef_normMF;         // Right now, same norm func as default in GMRES.
      gmres.setNorm(normF);
      PrecondFunc prec = &PeleLM::ef_applyPrecond;
      gmres.setPrecond(prec);
      gmres.setVerbose(ef_GMRES_verbose);
      gmres.setMaxRestart(ef_GMRES_maxRst);
   }

   // Need to create the preconditioner LinOp
   PCLinOp_needUpdate = 1;
   PCMLMG_needUpdate = 1;

   int NK_tot_count = 0;
   for (int sstep = 0; sstep < ef_substep; sstep++) {
      curtime = prev_time + (sstep+1) * dtsub;
      // -----------------
      // Pre Newton
      // Set up the NL state scaling
      nE_scale = (nl_state.norm0(0) > 1.0e-12) ? nl_state.norm0(0) : 1.0;
      phiV_scale = (nl_state.norm0(1) > 1.0e-6 ) ? nl_state.norm0(1) : 1.0;
      nl_state.mult(1.0/nE_scale,0,1);
      nl_state.mult(1.0/phiV_scale,1,1);
      if ( ef_verbose ) {
         amrex::Print() << "(" << sstep << ") ne scaling: " << nE_scale << "\n";
         amrex::Print() << "(" << sstep << ") phiV scaling: " << phiV_scale << "\n";
      }

      // Compute the background charge distribution
      ef_bg_chrg((sstep+1)*dtsub, Dn, Dnp1, Dhat);
      if ( ef_debug ) VisMF::Write(bg_charge,"NL_BGcharge_"+std::to_string(level));

      // Newton initial guess
      newtonInitialGuess(dtsub,nl_state);
      ef_normMF(nl_state,nl_stateNorm);

      // Initial NL residual: update residual scaling and preconditioner
      ef_nlResidual( dtsub, nl_state, nl_resid, true, true );
      nl_resid.mult(-1.0);
      ef_normMF(nl_resid,nl_residNorm);
      if ( ef_debug ) VisMF::Write(nl_resid,"NLResInit_Lvl"+std::to_string(level));

      // Check for convergence
      Real max_nlres = std::max(nl_resid.norm0(0),nl_resid.norm0(1));
      if ( max_nlres <= ef_newtonTol ) {
         ForcingnE.setVal(0.0);
         if ( ef_verbose ) {
            amrex::Print() << "No Newton iteration needed, exiting. \n";
         }
         return;
      }

      // -----------------
      // Newton iteration
      int exit_newton = 0;
      int NK_ite = 0;
      do {
         NK_ite += 1;

         // Verbose
         if ( ef_verbose ) {
            amrex::Print() << " Newton it: " << NK_ite << " L2**2 residual: " << 0.5*nl_residNorm*nl_residNorm
                                                       << ". Linf residual: " << max_nlres << "\n";
         }

         // Solve for Newton direction
         MultiFab newtonDir(grids,dmap,2,1);
         newtonDir.setVal(0.0,0,2,1);
         if ( !ef_use_PETSC_direct ) {
            const Real S_tol     = ef_GMRES_reltol;
            const Real S_tol_abs = max_nlres * ef_GMRES_reltol;
            GMRES_tot_count += gmres.solve(newtonDir,nl_resid,S_tol_abs,S_tol);
            if ( ef_debug ) VisMF::Write(newtonDir,"NLDir_NewtIte"+std::to_string(NK_ite)+"_Lvl"+std::to_string(level));
         } else {
         }

         // Linesearch & update state: TODO
         nl_state.plus(newtonDir,0,2,0);
         ef_normMF(nl_state,nl_stateNorm);
         ef_nlResidual( dtsub, nl_state, nl_resid, false, true );
         nl_resid.mult(-1.0);
         if ( ef_debug ) VisMF::Write(nl_resid,"NLRes_NewtIte"+std::to_string(NK_ite)+"_Lvl"+std::to_string(level));
         ef_normMF(nl_resid,nl_residNorm);
         max_nlres = std::max(nl_resid.norm0(0),nl_resid.norm0(1));

         // Exit condition
         exit_newton = testExitNewton(nl_resid,NK_ite);

      } while( !exit_newton );
      NK_tot_count += NK_ite;

      // -----------------
      // Post Newton

      // Increment the forcing term
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ForcingnE,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& old_nE   = ef_state_old.const_array(mfi);
         auto const& new_nE   = nl_state.const_array(mfi);
         auto const& I_R_nE   = get_new_data(RhoYdot_Type).const_array(mfi,NUM_SPECIES);
         auto const& force    = ForcingnE.array(mfi);
         Real scaling         = nE_scale;
         Real dtinv           = 1.0 / dtsub;
         amrex::ParallelFor(bx, [old_nE, new_nE, I_R_nE, force, dtinv, scaling]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            force(i,j,k) += (new_nE(i,j,k) * scaling - old_nE(i,j,k)) * dtinv - I_R_nE(i,j,k);
         });
      }

      // Unscale nl_state and if not last subcycle update 'old' state
      nl_state.mult(nE_scale,0,1);
      nl_state.mult(phiV_scale,1,1);
      if ( sstep != ef_substep - 1 ) {
         MultiFab::Copy(ef_state_old,nl_state,0,0,2,1);
         // Need to update the gradient of old phiV
         {
            MultiFab phi_a(ef_state_old,amrex::make_alias,1,1);
            ef_calcGradPhiV(curtime,phi_a,gphiV_old);
         }
      }

   }

   // Update the state
   MultiFab&  S = get_new_data(State_Type);
   MultiFab::Copy(S,nl_state,0,nE,2,0);
   if ( ef_debug ) VisMF::Write(nl_state,"PostPNPState_"+std::to_string(level));
   //ef_solve_phiv(time+dt);
   //if ( time == state[State_Type].curTime() ) Print() << " Updating new PhiV\n";

   if ( ef_verbose )
   {
     const int IOProc = ParallelDescriptor::IOProcessorNumber();

     Real mx = ParallelDescriptor::second() - strt_time, mn = mx;

     ParallelDescriptor::ReduceRealMin(mn,IOProc);
     ParallelDescriptor::ReduceRealMax(mx,IOProc);

     if ( !ef_use_PETSC_direct ) {
        Real avgGMRES = (float)GMRES_tot_count/(float)NK_tot_count;
     //   Real avgMG = (float)MGitcount/(float)GMRES_tot_count;
        amrex::Print() << "Avg GMRES/Newton: " << avgGMRES << "\n";
     //   amrex::Print() << "Avg MGVcycl/GMRES: " << avgMG << "\n";
     }
     amrex::Print() << "PeleLM_EF::ef_solve_PNP(): lev: " << level << ", time: ["
                    << mn << " ... " << mx << "]\n";
   }

//   if ( level == 2 ) Abort("Stop right there");

}

int PeleLM::testExitNewton(const MultiFab  &res,
                                 int       newtonIter){

   int exit = 0;
   Real max_res = std::max(res.norm0(0),res.norm0(1));
   if ( max_res <= ef_newtonTol ) {
      exit = 1;
      if ( ef_verbose ) {
         amrex::Print() << " Newton iterations converged: \n";
         amrex::Print() << " Final Newton L2**2 res norm : " << 0.5*nl_residNorm*nl_residNorm << "\n";
         amrex::Print() << " Final Newton Linf res norm : " << max_res << "\n";
      }
   }

   if ( newtonIter >= ef_maxNewtonIter && exit == 0 ) {
      exit = 1;
      amrex::Print() << " Max Newton iteration reached without convergence !!! \n";
   }

   return exit;
}

void PeleLM::newtonInitialGuess(const Real      &dt_lcl,
                                      MultiFab  &a_nl_state) {

   // Get the unscaled non-linear state
   MultiFab nl_state_usc(grids,dmap,2,2);
   MultiFab::Copy(nl_state_usc, a_nl_state, 0, 0, 2, 2);
   nl_state_usc.mult(nE_scale,0,1);
   nl_state_usc.mult(phiV_scale,1,1);

   // Get aliases to make it easier
   MultiFab nE_usc(nl_state_usc,amrex::make_alias,0,1);
   MultiFab phi_usc(nl_state_usc,amrex::make_alias,1,1);
   MultiFab a_nE(a_nl_state,amrex::make_alias,0,1);
   MultiFab a_phi(a_nl_state,amrex::make_alias,1,1);

   // Lap(PhiV) and grad(PhiV)
   FluxBoxes gphi_fb(this, 1, 0);
   MultiFab** gphiV = gphi_fb.get();
   MultiFab laplacian_term(grids, dmap, 1, 0);
   compPhiVLap(phi_usc,laplacian_term,gphiV);

   // Diffusion term nE
   MultiFab diffnE(grids, dmap, 1, 0);
   compElecDiffusion(nE_usc,diffnE);

   // Advective term nE
   MultiFab advnE(grids, dmap, 1, 0);
   compElecAdvection(nE_usc,phi_usc,gphiV,advnE);

   // Get an estimate of the electron advective step size
   const Real* dx = geom.CellSize();
   Real nE_adv_dt = 1.0e20;
   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      Real umax = elec_Ueff[dir].norm0(0);
      nE_adv_dt = amrex::min(nE_adv_dt,dx[0]/umax);
   }

   ParallelDescriptor::ReduceRealMin(nE_adv_dt);

   nE_adv_dt = amrex::min(dt_lcl,nE_adv_dt);

   // Advance explicitly electron AD at nE_adv_dt
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(a_nl_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& ne_diff  = diffnE.const_array(mfi);
      auto const& ne_adv   = advnE.const_array(mfi);
      auto const& ne_curr  = a_nE.array(mfi);
      Real nE_scale_inv = 1.0 / nE_scale;
      amrex::ParallelFor(bx, [ne_curr,ne_diff,ne_adv,nE_adv_dt,nE_scale_inv]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         ne_curr(i,j,k) += 0.0 * nE_adv_dt * nE_scale_inv * ( ne_diff(i,j,k) + ne_adv(i,j,k) );
      });
   }
}

void PeleLM::ef_nlResidual(const Real      &dt_lcl,
                           const MultiFab  &a_nl_state,
                                 MultiFab  &a_nl_resid,
                                 int       update_res_scaling,
                                 int       update_precond){
   BL_PROFILE("PLM_EF::ef_nlResidual()");

   // Get the unscaled non-linear state
   MultiFab nl_state_usc(grids,dmap,2,2);
   MultiFab::Copy(nl_state_usc, a_nl_state, 0, 0, 2, 2);
   nl_state_usc.mult(nE_scale,0,1);
   nl_state_usc.mult(phiV_scale,1,1);

   // Get aliases to make it easier
   MultiFab nE_a(nl_state_usc,amrex::make_alias,0,1);
   MultiFab phi_a(nl_state_usc,amrex::make_alias,1,1);

   // Lap(PhiV) and grad(PhiV)
   FluxBoxes gphi_fb(this, 1, 0);
   MultiFab** gphiV = gphi_fb.get();
   MultiFab laplacian_term(grids, dmap, 1, 0);
   compPhiVLap(phi_a,laplacian_term,gphiV);
   if ( ef_debug ) VisMF::Write(laplacian_term,"NLRes_phiVLap_"+std::to_string(level));

   // Diffusion term nE
   MultiFab diffnE(grids, dmap, 1, 0);
   compElecDiffusion(nE_a,diffnE);
   if ( ef_debug ) VisMF::Write(diffnE,"NLRes_ElecDiff_"+std::to_string(level));

   // Advective term nE
   MultiFab advnE(grids, dmap, 1, 0);
   compElecAdvection(nE_a,phi_a,gphiV,advnE);
   if ( ef_debug ) VisMF::Write(advnE,"NLRes_ElecAdv_"+std::to_string(level));

   if ( ef_debug ) VisMF::Write(get_new_data(RhoYdot_Type),"NLres_IRnE_"+std::to_string(level));

   // Assemble the non-linear residual
   // res(ne(:)) = dt * ( diff(:) + conv(:) + I_R(:) ) - ( ne(:) - ne_old(:) )
   // res(phiv(:)) = \Sum z_k * \tilde Y_k / q_e - ne + Lapl_PhiV
   a_nl_resid.setVal(0.0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(a_nl_resid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& I_R_nE   = get_new_data(RhoYdot_Type).const_array(mfi,NUM_SPECIES);
      auto const& lapPhiV  = laplacian_term.const_array(mfi);
      auto const& ne_diff  = diffnE.const_array(mfi);
      auto const& ne_adv   = advnE.const_array(mfi);
      auto const& ne_curr  = nE_a.const_array(mfi);
      auto const& ne_old   = ef_state_old.const_array(mfi,0);
      auto const& charge   = bg_charge.const_array(mfi);
      auto const& res_nE   = a_nl_resid.array(mfi,0);
      auto const& res_phiV = a_nl_resid.array(mfi,1);
      Real scalLap         = EFConst::eps0 * EFConst::epsr / EFConst::elemCharge;
      amrex::ParallelFor(bx, [ne_curr,ne_old,lapPhiV,I_R_nE,ne_diff,ne_adv,charge,res_nE,res_phiV,dt_lcl,scalLap]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         res_nE(i,j,k) = ne_old(i,j,k) - ne_curr(i,j,k) + dt_lcl * ( ne_diff(i,j,k) + ne_adv(i,j,k) + I_R_nE(i,j,k) );
         res_phiV(i,j,k) = lapPhiV(i,j,k) * scalLap - ne_curr(i,j,k) + charge(i,j,k);
      });
   }

   // Deal with scaling
   if ( update_res_scaling ) {
      FnE_scale = (a_nl_resid.norm0(0) > 1.0e-12) ? a_nl_resid.norm0(0) : 1.0 ;
      FphiV_scale = (a_nl_resid.norm0(1) > 1.0e-12) ? a_nl_resid.norm0(1) : 1.0 ;
      if ( ef_verbose ) {
         amrex::Print() << " F(ne) scaling: " << FnE_scale << "\n";
         amrex::Print() << " F(PhiV) scaling: " << FphiV_scale << "\n";
      }
   }

   a_nl_resid.mult(1.0/FnE_scale,0,1);
   a_nl_resid.mult(1.0/FphiV_scale,1,1);

   // Update the preconditioner
   if ( update_precond && !ef_use_PETSC_direct ) {
      ef_setUpPrecond(dt_lcl, nl_state_usc);
   }
}

void PeleLM::compPhiVLap(MultiFab& phi,
                         MultiFab& phiLap,
                         MultiFab** gPhiV){

// Set-up Poisson operator
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   info.setMaxCoarseningLevel(0);
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
      FillPatch(crselev, phiV_crse, 0, curtime, State_Type, PhiV, 1, 0);
      phiV_poisson.setCoarseFineBC(&phiV_crse, crse_ratio[0]);
   }
   phiV_poisson.setLevelBC(0, &phi);

// LinearSolver to get divergence
   MLMG solver(phiV_poisson);
   solver.apply({&phiLap},{&phi});

// Need the flux (grad(phi))
   Array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(gPhiV[0],gPhiV[1],gPhiV[2])};
   solver.getFluxes({fp},{&phi});
}

void PeleLM::compElecDiffusion(MultiFab& a_ne,
                               MultiFab& elecDiff)
{
   // Set-up Poisson operator
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   info.setMaxCoarseningLevel(0);
   MLABecLaplacian ne_lapl({geom}, {grids}, {dmap}, info);
   ne_lapl.setMaxOrder(ef_PoissonMaxOrder);

   // Set-up BC's
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
   ef_set_neBC(mlmg_lobc, mlmg_hibc);
   ne_lapl.setDomainBC(mlmg_lobc, mlmg_hibc);

   MultiFab nE_crse;
   if (level > 0) {
      auto& crselev = getLevel(level-1);
      nE_crse.define(crselev.grids, crselev.dmap, 1, 0, MFInfo(), crselev.Factory());
      FillPatch(crselev, nE_crse, 0, curtime, State_Type, nE, 1, 0);
      ne_lapl.setCoarseFineBC(&nE_crse, crse_ratio[0]);
   }
   ne_lapl.setLevelBC(0, &a_ne);

   // Coeffs
   ne_lapl.setScalars(0.0, 1.0);
   Array<const MultiFab*,AMREX_SPACEDIM> bcoeffs{D_DECL(De_ec[0],De_ec[1],De_ec[2])};
   ne_lapl.setBCoeffs(0, bcoeffs);

   // LinearSolver to get divergence
   MLMG solver(ne_lapl);
   solver.apply({&elecDiff},{&a_ne});

   elecDiff.mult(-1.0);

// TODO: will need the flux
   FluxBoxes fluxb(this, 1, 0);
   MultiFab** ne_DiffFlux = fluxb.get();
   Array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(ne_DiffFlux[0],ne_DiffFlux[1],ne_DiffFlux[2])};
   solver.getFluxes({fp},{&a_ne});
//   if ( ef_debug ) VisMF::Write(*ne_DiffFlux[0],"NLRes_nEDiffFlux_"+std::to_string(level));
//   if ( ef_debug ) VisMF::Write(*ne_DiffFlux[1],"NLRes_nEDiffFlux_"+std::to_string(level));
}

void PeleLM::compElecAdvection(MultiFab &a_ne,
                               MultiFab &a_phiV,
                               MultiFab *gphiV[AMREX_SPACEDIM],
                               MultiFab &elecAdv)
{

   // Get the face effective velocity
   // effVel = Umac - \mu_e * 0.5 * ( gradPhiV_old + gradPhiVcurr)
   for (int d = 0; d < AMREX_SPACEDIM; ++d) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(elec_Ueff[d],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& ueff    = elec_Ueff[d].array(mfi);
         auto const& gphi_c  = gphiV[d]->const_array(mfi);
         auto const& gphi_o  = gphiV_old[d]->const_array(mfi);
         auto const& umac    = u_mac[d].const_array(mfi);
         auto const& kappa_e = Ke_ec[d]->const_array(mfi);
         amrex::ParallelFor(bx, [ueff, gphi_c, gphi_o, umac, kappa_e]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            ueff(i,j,k) = umac(i,j,k) - kappa_e(i,j,k) * 0.5 * ( gphi_c(i,j,k) + gphi_o(i,j,k) );
         });
      }
   }
   if ( ef_debug ) VisMF::Write(elec_Ueff[0],"NLRes_ElecUeffX_"+std::to_string(level));
   if ( ef_debug ) VisMF::Write(elec_Ueff[1],"NLRes_ElecUeffY_"+std::to_string(level));

   MultiFab nE_new(grids,dmap,1,1);
   nE_new.setVal(0.0);
   MultiFab::Copy(nE_new,a_ne,0,0,1,1);
   nE_new.FillBoundary();

   {
      FArrayBox cflux[AMREX_SPACEDIM];
      FArrayBox edgstate[AMREX_SPACEDIM];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(elecAdv,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         const Box& bx  = mfi.tilebox();
         const Box& gbx = mfi.growntilebox(1);
         AMREX_D_TERM( const Box& xbx = mfi.grownnodaltilebox(0,0);,
                       const Box& ybx = mfi.grownnodaltilebox(1,0);,
                       const Box& zbx = mfi.grownnodaltilebox(2,0););

         // data arrays
         auto const& ne_arr = nE_new.array(mfi);
         auto const& ne_adv = elecAdv.array(mfi);
         AMREX_D_TERM( Array4<Real const> u = elec_Ueff[0].const_array(mfi);,
                       Array4<Real const> v = elec_Ueff[1].const_array(mfi);,
                       Array4<Real const> w = elec_Ueff[2].const_array(mfi););

         // Set temporary edge FABs
         AMREX_D_TERM( cflux[0].resize(xbx,1);,
                       cflux[1].resize(ybx,1);,
                       cflux[2].resize(zbx,1););
         AMREX_D_TERM( edgstate[0].resize(xbx,1);,
                       edgstate[1].resize(ybx,1);,
                       edgstate[2].resize(zbx,1););
         AMREX_D_TERM( Array4<Real> xstate = edgstate[0].array();,
                       Array4<Real> ystate = edgstate[1].array();,
                       Array4<Real> zstate = edgstate[2].array(););
         AMREX_D_TERM( Array4<Real> xflux = cflux[0].array();,
                       Array4<Real> yflux = cflux[1].array();,
                       Array4<Real> zflux = cflux[2].array(););

         // Predict edge states (TODO: ignoring bcs for now !)
         amrex::ParallelFor(xbx, [ne_arr,u,xstate]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            xstate(i,j,k) = ef_edge_state(i,j,k,0,ne_arr,u);
         });
         amrex::ParallelFor(ybx, [ne_arr,v,ystate]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            ystate(i,j,k) = ef_edge_state(i,j,k,1,ne_arr,v);
         });
#if ( AMREX_SPACEDIM ==3 )
         amrex::ParallelFor(zbx, [ne_arr,w,zstate]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            zstate(i,j,k) = ef_edge_state(i,j,k,2,ne_arr,w);
         });
#endif

         // Computing fluxes
         amrex::ParallelFor(xbx, [u,xstate,xflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            xflux(i,j,k) = u(i,j,k) * xstate(i,j,k);
         });
         amrex::ParallelFor(ybx, [v,ystate,yflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            yflux(i,j,k) = v(i,j,k) * ystate(i,j,k);
         });
#if ( AMREX_SPACEDIM ==3 )
         amrex::ParallelFor(zbx, [w,zstate,zflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            zflux(i,j,k) = w(i,j,k) * zstate(i,j,k);
         });
#endif

         // Compute divergence
         const auto dxinv = geom.InvCellSizeArray();
         amrex::ParallelFor(bx, [ ne_adv, D_DECL(xflux,yflux,zflux), dxinv]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            ne_adv(i,j,k) =   dxinv[0] * (xflux(i+1,j,k) - xflux(i,j,k))
                            + dxinv[1] * (yflux(i,j+1,k) - yflux(i,j,k))
#if ( AMREX_SPACEDIM ==3 )
                            + dxinv[2] * (zflux(i,j,k+1) - zflux(i,j,k))
#endif
                            ;
         });
      }

   }

   elecAdv.mult(-1.0);
}
void PeleLM::ef_setUpPrecond (const Real &dt_lcl,
                              const MultiFab& a_nl_state) {
   BL_PROFILE("PLM_EF::ef_setUpPrecond()");

   if ( PCLinOp_needUpdate ) {
      LPInfo info;
      info.setAgglomeration(1);
      info.setConsolidation(1);
      info.setMetricTerm(false);

      if ( pnp_pc_drift != nullptr ) {
         delete pnp_pc_drift;
         delete pnp_pc_Stilda;
         delete pnp_pc_diff;
      }

      pnp_pc_drift = new MLABecLaplacian({geom}, {grids}, {dmap}, info);
      pnp_pc_Stilda = new MLABecLaplacian({geom}, {grids}, {dmap}, info);
      pnp_pc_diff = new MLABecCecLaplacian({geom}, {grids}, {dmap}, info);

      pnp_pc_Stilda->setMaxOrder(ef_PoissonMaxOrder);
      pnp_pc_diff->setMaxOrder(ef_PoissonMaxOrder);
      pnp_pc_drift->setMaxOrder(ef_PoissonMaxOrder);

      PCLinOp_needUpdate = 0;
   }

   // Set diff/drift operator
   {
      pnp_pc_diff->setScalars(-nE_scale/FnE_scale, -dt_lcl*nE_scale/FnE_scale, dt_lcl*nE_scale/FnE_scale);
      Real omega = 1.0;
      pnp_pc_diff->setRelaxation(omega);
      pnp_pc_diff->setACoeffs(0, 1.0);
      std::array<const MultiFab*,AMREX_SPACEDIM> bcoeffs{AMREX_D_DECL(De_ec[0],De_ec[1],De_ec[2])};
      pnp_pc_diff->setBCoeffs(0, bcoeffs);
      std::array<const MultiFab*,AMREX_SPACEDIM> ccoeffs{AMREX_D_DECL(&elec_Ueff[0],&elec_Ueff[1],&elec_Ueff[2])};
      pnp_pc_diff->setCCoeffs(0, ccoeffs);
   }

   // Stilda and Drift LinOp
   {
      // Upwinded edge neKe values
      MultiFab nEKe(grids,dmap,1,1);
      MultiFab nE_a(a_nl_state,amrex::make_alias,0,1);  // State is not scale at this point
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(nEKe, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx = mfi.growntilebox();
         auto const& neke   = nEKe.array(mfi);
         auto const& ne_arr = nE_a.const_array(mfi);
         amrex::ParallelFor(gbx, [neke,ne_arr]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getKappaE(i,j,k,neke);
            neke(i,j,k) *= ne_arr(i,j,k);
         });
      }

      FluxBoxes edge_fb(this, 1, 1);
      MultiFab** neKe_ec = edge_fb.get();
      auto math_bc = fetchBCArray(State_Type,nE,1);
      const Box& domain = geom.Domain();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(nEKe, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
         {
            const Box ebx = mfi.nodaltilebox(dir);
            const Box& edomain = amrex::surroundingNodes(domain,dir);
            const auto& neke_c  = nEKe.array(mfi);
            const auto& neke_ed = neKe_ec[dir]->array(mfi);
            const auto& ueff_ed = elec_Ueff[dir].array(mfi);
            const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
            const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
            amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, neke_c, neke_ed, ueff_ed, edomain]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
               bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
               cen2edg_upwind_cpp( i, j, k, dir, 1, on_lo, on_hi, ueff_ed, neke_c, neke_ed);
            });
         }
      }

      pnp_pc_drift->setScalars(0.0,0.5*phiV_scale/FnE_scale*dt_lcl);
      {
         std::array<const MultiFab*,AMREX_SPACEDIM> bcoeffs{AMREX_D_DECL(neKe_ec[0],neKe_ec[1],neKe_ec[2])};
         pnp_pc_drift->setBCoeffs(0, bcoeffs);
      }
      if ( ef_debug ) {
         VisMF::Write(*neKe_ec[0],"PC_Drift_nEKe_edgeX_lvl"+std::to_string(level));
         VisMF::Write(*neKe_ec[1],"PC_Drift_nEKe_edgeY_lvl"+std::to_string(level));
      }

      // Get the diagonal of pnp_pc_diff

      //pnp_pc_Stilda->setScalars(0.0,-phiV_scale/FphiV_scale);
      pnp_pc_Stilda->setScalars(0.0,-1.0);
      Real scalLap = EFConst::eps0 * EFConst::epsr / EFConst::elemCharge;
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
         neKe_ec[dir]->mult(0.5*dt_lcl,0,1);
         neKe_ec[dir]->plus(scalLap,0,1);
      }
      if ( ef_debug ) {
         VisMF::Write(*neKe_ec[0],"PC_Stilda_nEKepLap_edgeX_lvl"+std::to_string(level));
         VisMF::Write(*neKe_ec[1],"PC_Stilda_nEKepLap_edgeY_lvl"+std::to_string(level));
      }
      {
         std::array<const MultiFab*,AMREX_SPACEDIM> bcoeffs{AMREX_D_DECL(neKe_ec[0],neKe_ec[1],neKe_ec[2])};
         pnp_pc_Stilda->setBCoeffs(0, bcoeffs);
      }
   }

   // Set up the domainBCs
   std::array<LinOpBCType,AMREX_SPACEDIM> ne_lobc, ne_hibc;
   std::array<LinOpBCType,AMREX_SPACEDIM> phiV_lobc, phiV_hibc;
   ef_set_PoissonBC(phiV_lobc, phiV_hibc);
   ef_set_neBC(ne_lobc,ne_hibc);
   pnp_pc_Stilda->setDomainBC(phiV_lobc, phiV_hibc);
   pnp_pc_diff->setDomainBC(ne_lobc, ne_hibc);
   pnp_pc_drift->setDomainBC(phiV_lobc, phiV_hibc);

   // Trigger update of the MLMGs
   PCMLMG_needUpdate = 1;

}

void PeleLM::ef_applyPrecond (const MultiFab  &v,
                                    MultiFab  &Pv) {
   BL_PROFILE("PLM_EF::ef_applyPrecond()");

   //MultiFab::Copy(Pv, v, 0, 0, 2, 0);
   //return;

   // Set up some aliases to make things easier
   MultiFab a_ne(v,amrex::make_alias,0,1);
   MultiFab a_phiV(v,amrex::make_alias,1,1);
   MultiFab a_Pne(Pv,amrex::make_alias,0,1);
   MultiFab a_PphiV(Pv,amrex::make_alias,1,1);

   // TODO: I need to initialize the result to zero otherwise MLMG goes nuts
   // or do I ?
   a_Pne.setVal(0.0,0,1,1);
   a_PphiV.setVal(0.0,0,1,1);

   // Set up the linear solvers BCs
   pnp_pc_diff->setLevelBC(0, &a_Pne);
   pnp_pc_drift->setLevelBC(0, &a_PphiV);
   pnp_pc_Stilda->setLevelBC(0, &a_PphiV);

   // Set Coarse/Fine BCs
   // Assumes it should be at zero.
   if ( level > 0 ) {
      pnp_pc_diff->setCoarseFineBC(nullptr, crse_ratio[0]);
      pnp_pc_drift->setCoarseFineBC(nullptr, crse_ratio[0]);
      pnp_pc_Stilda->setCoarseFineBC(nullptr, crse_ratio[0]);
   }

   // Create MLMGs
   if ( PCMLMG_needUpdate ) {
      if ( mg_diff != nullptr ) {
         delete mg_diff;
         delete mg_drift;
         delete mg_Stilda;
      }
      mg_diff = new MLMG(*pnp_pc_diff);
      mg_drift = new MLMG(*pnp_pc_drift);
      mg_Stilda = new MLMG(*pnp_pc_Stilda);

      PCMLMG_needUpdate = 0;
   }

   mg_diff->setVerbose(0);
   mg_drift->setVerbose(0);
   mg_Stilda->setVerbose(0);
   if ( ef_PC_fixedIter > 0 ) {
      mg_diff->setFixedIter(ef_PC_fixedIter);
      mg_drift->setFixedIter(ef_PC_fixedIter);
      mg_Stilda->setFixedIter(ef_PC_fixedIter);
   }


   Real S_tol     = ef_PC_MG_Tol;
   Real S_tol_abs = a_ne.norm0() * ef_PC_MG_Tol;

   // Most inner mat
   // --                --
   // | [dtD-I]^-1     0 |
   // |                  |
   // |       0        I |
   // --                --
   mg_diff->solve({&a_Pne}, {&a_ne}, S_tol, S_tol_abs);
   MultiFab::Copy(a_PphiV,a_phiV,0,0,1,0);

   // Assembling mat
   // --       --
   // |  I    0 |
   // |         |
   // | -Ie   I |
   // --       --
   MultiFab::Saxpy(a_PphiV,nE_scale/FphiV_scale,a_Pne,0,0,1,0);

   // PhiV estimate mat
   // --         --
   // | I     0   |
   // |           |
   // | 0   S*^-1 |
   // --         --
   S_tol_abs = a_PphiV.norm0() * ef_PC_MG_Tol;
   MultiFab temp(grids,dmap,1,1);
   temp.setVal(0.0,0,1,0);
   // Scale the solve RHS
   a_PphiV.mult(FphiV_scale/phiV_scale);
   mg_Stilda->solve({&temp},{&a_PphiV}, S_tol, S_tol_abs);
   MultiFab::Copy(a_PphiV, temp, 0, 0, 1, 0);


   // Final mat
   // --                          --
   // | I       -[dtD - I]^-1 dtDr |
   // |                            |
   // | 0                I         |
   // --                          --
   mg_drift->apply({&temp},{&a_PphiV});
   S_tol_abs = temp.norm0() * ef_PC_MG_Tol;
   MultiFab temp2(grids,dmap,1,1);
   temp2.setVal(0.0,0,1,0);
   mg_diff->solve({&temp2},{&temp}, S_tol, S_tol_abs);
   temp2.mult(-1.0);
   MultiFab::Add(a_Pne,temp2,0,0,1,0);

}

void PeleLM::ef_calcGradPhiV(const Real&    time_lcl,
                                   MultiFab &a_phiv,
                                   MultiFab *grad_phiV[AMREX_SPACEDIM]) {

   // Set-up Poisson operator
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   info.setMaxCoarseningLevel(0);
   MLPoisson poisson({geom}, {grids}, {dmap}, info);
   poisson.setMaxOrder(ef_PoissonMaxOrder);

   // BCs
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
   ef_set_PoissonBC(mlmg_lobc, mlmg_hibc);
   poisson.setDomainBC(mlmg_lobc, mlmg_hibc);

   MultiFab phiV_crse;
   if (level > 0) {
      auto& crselev = getLevel(level-1);
      phiV_crse.define(crselev.grids, crselev.dmap, 1, 0, MFInfo(), crselev.Factory());
      FillPatch(crselev,phiV_crse,0,time_lcl, State_Type, PhiV, 1, 0);
      poisson.setCoarseFineBC(&phiV_crse, crse_ratio[0]);
   }
   poisson.setLevelBC(0, &a_phiv);

   // Linear solver
   MLMG mlmg(poisson);
   std::array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(grad_phiV[0],grad_phiV[1],grad_phiV[2])};
   mlmg.getFluxes({fp},{&a_phiv});

   if ( ef_debug ) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
         VisMF::Write(*grad_phiV[d],"GradPhiV_Dir"+std::to_string(d)+"_lvl"+std::to_string(level));
      }
   }

}

void PeleLM::ef_calcUDriftSpec(const Real &time,
                               int   use_instant_coeff){

   // Get CC t^n+1/2 species mobilities
   MultiFab KSpec_half(grids,dmap,NUM_SPECIES,1);
   if ( time > 0.0 && !use_instant_coeff ) {   // Not initial time step or using instantaneous diff coeff
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(KSpec_old,TilingIfNotGPU()); mfi.isValid();++mfi)
      {
         const Box& gbx = mfi.growntilebox();
         const auto& Ksp_o  = KSpec_old.array(mfi);
         const auto& Ksp_n  = KSpec_new.array(mfi);
         const auto& Ksp_h  = KSpec_half.array(mfi);
         amrex::ParallelFor(gbx, NUM_SPECIES, [Ksp_o,Ksp_n,Ksp_h]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            Ksp_h(i,j,k,n) = 0.5 * ( Ksp_o(i,j,k,n) + Ksp_n(i,j,k,n) );
         });
      }
   } else {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(KSpec_new,TilingIfNotGPU()); mfi.isValid();++mfi)
      {
         const Box& gbx = mfi.growntilebox();
         const auto& Ksp_n  = KSpec_new.array(mfi);
         const auto& Ksp_h  = KSpec_half.array(mfi);
         amrex::ParallelFor(gbx, NUM_SPECIES, [Ksp_n,Ksp_h]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            Ksp_h(i,j,k,n) = Ksp_n(i,j,k,n);
         });
      }
   }

   // CC -> EC
   FluxBoxes fb(this, NUM_SPECIES, 0);
   MultiFab  **KSpec_ec = fb.get();
   auto math_bc = fetchBCArray(State_Type,nE,1);
   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(KSpec_old,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         const Box ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         const auto& Ksp_h  = KSpec_half.array(mfi);
         const auto& Ksp_ec = KSpec_ec[dir]->array(mfi);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, use_harmonic_avg, Ksp_h, Ksp_ec, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, NUM_SPECIES, use_harmonic_avg, on_lo, on_hi, Ksp_h, Ksp_ec);
         });
      }
   }

   // Get old and new gPhiV
   MultiFab& S = get_new_data(State_Type);
   Real prev_time = state[State_Type].prevTime();
   FillPatchIterator fpi_old(*this,S,1,prev_time,State_Type,PhiV,1);
   MultiFab& phiV_old = fpi_old.get_mf();
   FillPatchIterator fpi_new(*this,S,1,time,State_Type,PhiV,1);
   MultiFab& phiV_new = fpi_new.get_mf();

   FluxBoxes fb_old(this, 1, 0);
   FluxBoxes fb_new(this, 1, 0);
   MultiFab** gphiV_o = fb_old.get();
   MultiFab** gphiV_n = fb_new.get();
   ef_calcGradPhiV(prev_time, phiV_old, gphiV_o);
   ef_calcGradPhiV(time, phiV_new, gphiV_n);

   // Build drift velo with time centered gPhiV and Ksp
   for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
   {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Udrift_spec[dir],TilingIfNotGPU()); mfi.isValid();++mfi)
      {
         const Box bx = mfi.tilebox();
         const auto& Ksp_ec = KSpec_ec[dir]->const_array(mfi);
         const auto& gp_o   = gphiV_o[dir]->const_array(mfi);
         const auto& gp_n   = gphiV_n[dir]->const_array(mfi);
         const auto& Ud_Sp  = Udrift_spec[dir].array(mfi);
         amrex::ParallelFor(bx, NUM_SPECIES, [Ksp_ec,gp_o,gp_n,Ud_Sp]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            Ud_Sp(i,j,k,n) = Ksp_ec(i,j,k,n) * 0.5 * ( gp_o(i,j,k) + gp_n(i,j,k) );
         });
      }
   }
   if ( ef_debug ) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
         VisMF::Write(Udrift_spec[d],"SpecUdrift_Dir"+std::to_string(d)+"_lvl"+std::to_string(level));
      }
   }
}

Real
PeleLM::ef_estTimeStep()
{
   Real ion_estdt = 1.0e20;
   const Real  cur_time = state[State_Type].curTime();
   MultiFab&   U_new = get_new_data(State_Type);
   ef_calc_transport(cur_time);
   ef_calcUDriftSpec(cur_time,1);      // Use instantaneous diffusion coefficients.

   // Get the cell centered max of effective velocity accross all directions
   // for each species
   MultiFab UeffMax_ccavg(grids, dmap, NUM_SPECIES, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(UeffMax_ccavg,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      const Box bx = mfi.tilebox();
      AMREX_D_TERM( const auto& udr  = Udrift_spec[0].const_array(mfi);,
                    const auto& vdr  = Udrift_spec[1].const_array(mfi);,
                    const auto& wdr  = Udrift_spec[2].const_array(mfi););
      AMREX_D_TERM( const auto& u    = U_new.const_array(mfi,Xvel);,
                    const auto& v    = U_new.const_array(mfi,Yvel);,
                    const auto& w    = U_new.const_array(mfi,Zvel););
      const auto& Ueff_arr = UeffMax_ccavg.array(mfi);
      amrex::ParallelFor(bx, NUM_SPECIES, [Ueff_arr,AMREX_D_DECL(udr,vdr,wdr),AMREX_D_DECL(u,v,w)]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
         Ueff_arr(i,j,k,n) = 0.0;
         Ueff_arr(i,j,k,n) = amrex::max(Ueff_arr(i,j,k,n),std::abs(u(i,j,k)+0.5*(udr(i,j,k,n)+udr(i+1,j,k,n))));
         Ueff_arr(i,j,k,n) = amrex::max(Ueff_arr(i,j,k,n),std::abs(v(i,j,k)+0.5*(vdr(i,j,k,n)+vdr(i,j+1,k,n))));
#if AMREX_SPACEDIM == 3
         Ueff_arr(i,j,k,n) = amrex::max(Ueff_arr(i,j,k,n),std::abs(w(i,j,k)+0.5*(wdr(i,j,k,n)+wdr(i,j,k+1,n))));
#endif
      });
   }

   const auto dx = geom.CellSizeArray();
   amrex::Real cfl_lcl = cfl;
   ion_estdt = amrex::ReduceMin(UeffMax_ccavg, 0, [dx,cfl_lcl]
      AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& ueffm) noexcept -> Real
      {
         using namespace amrex::literals;
         const auto lo = amrex::lbound(bx);
         const auto hi = amrex::ubound(bx);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
         amrex::Real velmax = std::numeric_limits<amrex::Real>::min();
#else
         amrex::Real velmax = -1.e37_rt;
#endif
         for       (int k = lo.z; k <= hi.z; ++k) {
            for    (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  // max accross species
                  amrex::Real velmaxsp = -1.e20_rt;
                  for ( int n = 0; n < NUM_SPECIES; n++ ) {
                     velmaxsp = amrex::max(velmaxsp,ueffm(i,j,k,n));
                  }
                  velmax = amrex::max(velmax,velmaxsp);
               }
            }
         }
         return dx[0]/velmax*cfl_lcl;
      });

   ParallelDescriptor::ReduceRealMin(ion_estdt);

   return ion_estdt;
}

void PeleLM::ef_bg_chrg(const Real      &dt_lcl,
                        const MultiFab  &Dn,
                        const MultiFab  &Dnp1,
                        const MultiFab  &Dhat) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(bg_charge,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhoYold  = get_old_data(State_Type).const_array(mfi,first_spec);
      auto const& adv_arr  = aofs->const_array(mfi,first_spec);
      auto const& dn_arr   = Dn.const_array(mfi);
      auto const& dnp1_arr = Dnp1.const_array(mfi);
      auto const& dhat_arr = Dhat.const_array(mfi);
      auto const& rhoYdot  = get_new_data(RhoYdot_Type).const_array(mfi);
      auto const& charge   = bg_charge.array(mfi);
      Real        factor = 1.0 / EFConst::elemCharge;
      amrex::ParallelFor(bx, [dt_lcl, rhoYold, adv_arr, dn_arr, dnp1_arr, dhat_arr, rhoYdot, charge, factor]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         charge(i,j,k) = 0.0;
         for (int n = 0; n < NUM_SPECIES; n++) {
            Real rhoYprov = rhoYold(i,j,k,n) + dt_lcl * ( adv_arr(i,j,k,n) +
                                                          0.5 * ( dn_arr(i,j,k,n) - dnp1_arr(i,j,k,n) ) +
                                                          dhat_arr(i,j,k,n) +
                                                          rhoYdot(i,j,k,n) );
            rhoYprov = amrex::max(rhoYprov,0.0);
            charge(i,j,k) += zk[n] * rhoYprov;
         }
         charge(i,j,k) *= factor;
      });
   }
}

void PeleLM::ef_normMF(const MultiFab &a_vec,
                             Real &norm){
   norm = 0.0;
   for ( int comp = 0; comp < a_vec.nComp(); comp++ ) {
      norm += MultiFab::Dot(a_vec,comp,a_vec,comp,1,0);
   }
   norm = std::sqrt(norm);
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

// Setup BC conditions for linear Poisson solve on PhiV. Directly copied from the diffusion one ...
void PeleLM::ef_set_neBC(std::array<LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                         std::array<LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc) {

    const BCRec& bc = get_desc_lst()[State_Type].getBC(nE);

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
