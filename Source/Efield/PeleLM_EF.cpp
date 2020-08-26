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
   pp.query("debug",ef_debug);
   pp.query("Poisson_tol",ef_PoissonTol);
   pp.query("Poisson_maxiter",ef_PoissonMaxIter);
   pp.query("Poisson_verbose",ef_PoissonVerbose);
   pp.query("Poisson_maxorder",ef_PoissonMaxOrder);
   pp.query("JFNK_lambda",ef_lambda_jfnk);
   pp.query("JFNK_difftype",ef_diffT_jfnk);
   pp.query("JFNK_maxNewton",ef_maxNewtonIter);
   pp.query("JFNK_newtonTol",ef_newtonTol);
   pp.query("GMRES_restart_size",ef_GMRES_size);
   pp.query("GMRES_rel_tol",ef_GMRES_reltol);
   pp.query("GMRES_max_restart",ef_GMRES_maxRst);
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

   ef_state_old.define(grids,dmap,2,1); 
   ef_state_refGhostCell.define(grids,dmap,2,Godunov::hypgrow()); 
   bg_charge.define(grids,dmap,1,1); 
   nl_state.define(grids,dmap,2,Godunov::hypgrow());
   nl_resid.define(grids,dmap,2,Godunov::hypgrow());
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
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, use_harmonic_avg, diff_c, diff_ed, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, 1, use_harmonic_avg, on_lo, on_hi, diff_c, diff_ed);
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
      MultiFab S_pert(grids,dmap,2,Godunov::hypgrow());
      MultiFab::Copy(S_pert, nl_state, 0, 0, 2, Godunov::hypgrow());
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
      MultiFab S_pertm(grids,dmap,2,Godunov::hypgrow());
      MultiFab S_pertp(grids,dmap,2,Godunov::hypgrow());
      MultiFab::Copy(S_pertp, nl_state, 0, 0, 2, Godunov::hypgrow());
      MultiFab::Copy(S_pertm, nl_state, 0, 0, 2, Godunov::hypgrow());
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
   FillPatchIterator nEfpi(*this,ef_state_refGhostCell,Godunov::hypgrow(),prev_time,State_Type,nE,2);
   MultiFab& nEfpi_mf = nEfpi.get_mf();
   MultiFab::Copy(ef_state_old,nEfpi_mf,0,0,2,1);                            // Copy into a storage to hold the old state
   MultiFab::Copy(ef_state_refGhostCell,nEfpi_mf,0,0,2,Godunov::hypgrow());  // Copy into a storage to hold the reference ghost cell value

   // Non-linear state & residual from FPI
   MultiFab::Copy(nl_state, ef_state_refGhostCell, 0, 0, 2, Godunov::hypgrow());
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
      //PrecondFunc prec = &PeleLM::applyPrec;
      //gmres.setPrecond(prec);
      gmres.setVerbose(1);
      gmres.setMaxRestart(ef_GMRES_maxRst);
   }

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
         newtonDir.setVal(0.0);
         if ( !ef_use_PETSC_direct ) {
            const Real S_tol     = ef_GMRES_reltol;
            const Real S_tol_abs = max_nlres * ef_GMRES_reltol;
            GMRES_tot_count += gmres.solve(newtonDir,nl_resid,S_tol_abs,S_tol);
            if ( ef_debug ) VisMF::Write(newtonDir,"NLDir_NewtIte"+std::to_string(NK_ite)+"_Lvl"+std::to_string(level));
         } else {
         }

         // Linesearch & update state: TODO
         newtonDir.mult(-1.0);
         nl_state.plus(newtonDir,0,2,Godunov::hypgrow());
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
         auto const& I_R_nE   = get_new_data(RhoYdot_Type).const_array(mfi);
         auto const& force    = ForcingnE.array(mfi);
         Real scaling         = nE_scale;
         Real dtinv           = 1.0 / dtsub;
         amrex::ParallelFor(bx, [old_nE, new_nE, I_R_nE, force, dtinv, scaling]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            force(i,j,k) += (new_nE(i,j,k) * scaling - old_nE(i,j,k)) * dtinv; // - I_R_nE(i,j,k);
         });
      }

      // If not last substep, unscale nl_state and update 'old' state
      if ( sstep != ef_substep - 1 ) {
         nl_state.mult(nE_scale,0,1);
         nl_state.mult(phiV_scale,1,1);
         MultiFab::Copy(ef_state_old,nl_state,0,0,2,1);
      }
      
   }

   // Update the state
   MultiFab&  S = get_new_data(State_Type);
   MultiFab::Copy(S,nl_state,0,nE,2,0);
   if ( ef_debug ) VisMF::Write(nl_state,"PostPNPState_"+std::to_string(level));

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

   if ( level == 2 ) Abort("Stop right there");

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

   if ( newtonIter >= ef_maxNewtonIter ) {
      exit = 1;
      amrex::Print() << " Max Newton iteration reached without convergence !!! \n";
   }

   return exit;
}

void PeleLM::ef_nlResidual(const Real      &dt_lcl,
                           const MultiFab  &a_nl_state,
                                 MultiFab  &a_nl_resid,
                                 int       update_res_scaling,
                                 int       update_precond){
   BL_PROFILE("PLM_EF::ef_nlResidual()");

   // Get the unscaled non-linear state
   MultiFab nl_state_usc(grids,dmap,2,Godunov::hypgrow());
   MultiFab::Copy(nl_state_usc, a_nl_state, 0, 0, 2, Godunov::hypgrow());
   nl_state_usc.mult(nE_scale,0,1);
   nl_state_usc.mult(phiV_scale,1,1);

   // Get aliases to make it easier
   MultiFab nE_a(nl_state_usc,amrex::make_alias,0,1);
   MultiFab phi_a(nl_state_usc,amrex::make_alias,1,1);

   // Lap(PhiV) and grad(PhiV)
   MultiFab laplacian_term(grids, dmap, 1, 0);
   compPhiVLap(phi_a,laplacian_term);
   if ( ef_debug ) VisMF::Write(laplacian_term,"NLRes_phiVLap_"+std::to_string(level));

   // Diffusion term nE
   MultiFab diffnE(grids, dmap, 1, 0);
   compElecDiffusion(nE_a,diffnE);
   if ( ef_debug ) VisMF::Write(diffnE,"NLRes_ElecDiff_"+std::to_string(level));

   // Advective term nE

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
      auto const& lapPhiV  = laplacian_term.const_array(mfi);
      auto const& ne_diff  = diffnE.const_array(mfi);
      auto const& ne_curr  = nE_a.const_array(mfi);
      auto const& ne_old   = ef_state_old.const_array(mfi,0);
      auto const& charge   = bg_charge.const_array(mfi);
      auto const& res_nE   = a_nl_resid.array(mfi,0);
      auto const& res_phiV = a_nl_resid.array(mfi,1);
      Real scalLap         = EFConst::eps0 * EFConst::epsr / EFConst::elemCharge;;
      amrex::ParallelFor(bx, [ne_curr,ne_old,lapPhiV,ne_diff,charge,res_nE,res_phiV,dt_lcl,scalLap]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         res_nE(i,j,k) = ne_old(i,j,k) - ne_curr(i,j,k) + dt_lcl * ( ne_diff(i,j,k) );
         res_phiV(i,j,k) = 0.0;//lapPhiV(i,j,k) * scalLap - ne_curr(i,j,k) + charge(i,j,k);
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
//      ef_setUpPrecond(dt, pnp_U, diffElec_ec);
   }
}

void PeleLM::compPhiVLap(MultiFab& phi,
                         MultiFab& phiLap){

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
//      if ( ef_debug ) VisMF::Write(phiV_crse,"NLRes_phiVCrse_Lap_"+std::to_string(level));
   }
   phiV_poisson.setLevelBC(0, &phi);

// LinearSolver to get divergence
   MLMG solver(phiV_poisson);
   solver.apply({&phiLap},{&phi});

// TODO: will need the flux (grad(phi))
   FluxBoxes fluxb(this, 1, 0);
   MultiFab** gradPhi = fluxb.get();
   Array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(gradPhi[0],gradPhi[1],gradPhi[2])};
   solver.getFluxes({fp},{&phi});
//   if ( ef_debug ) VisMF::Write(*gradPhi[0],"NLRes_gradPhiX_"+std::to_string(level));
//   if ( ef_debug ) VisMF::Write(*gradPhi[1],"NLRes_gradPhiY_"+std::to_string(level));
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
//      if ( ef_debug ) VisMF::Write(nE_crse,"NLRes_nECrse_nELap_"+std::to_string(level));
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
            charge(i,j,k) += zk[n] * rhoYprov * factor;
         }
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
