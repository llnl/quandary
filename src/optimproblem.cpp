#include "optimproblem.hpp"

OptimProblem::OptimProblem(const Config& config, OptimTarget* optim_target_, TimeStepper* timestepper_, MasterEq* mastereq_, MPI_Comm comm_init_, MPI_Comm comm_optim_, Output* output_, bool quietmode_){

  optim_target = optim_target_;
  timestepper = timestepper_;
  mastereq = mastereq_;
  ninit = config.getNInitialConditions();
  output = output_;
  quietmode = quietmode_;
  output_optimization_stride = config.getOutputOptimizationStride();
  optim_solver_type = config.getOptimSolverType();

  /* Reset */
  objective = 0.0;
  ksp_iters_last = 0;
  nonlinear_forward_valid = false;

  /* Store communicators */
  comm_init = comm_init_;
  comm_optim = comm_optim_;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize_world);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  MPI_Comm_rank(comm_init, &mpirank_init);
  MPI_Comm_size(comm_init, &mpisize_init);
  MPI_Comm_rank(comm_optim, &mpirank_optim);
  MPI_Comm_size(comm_optim, &mpisize_optim);

  /* Store number of initial conditions per init-processor group */
  ninit_local = ninit / mpisize_init;

  /* Store number of design parameters */
  int n = 0;
  for (size_t ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      n += mastereq->getOscillator(ioscil)->getNParams(); 
  }
  ndesign = n;
  if (mpirank_world == 0 && !quietmode) std::cout<< "Number of control parameters: " << ndesign << std::endl;

  /* Allocate adjoint terminal state */
  VecDuplicate(optim_target->getInitialState(), &rho_t0_bar);
  VecZeroEntries(rho_t0_bar);
  VecAssemblyBegin(rho_t0_bar); VecAssemblyEnd(rho_t0_bar);

  /* Get weights for the objective function (weighting the different initial conditions */
  obj_weights = config.getOptimWeights();

  /* Store other optimization parameters */
  gamma_tikhonov = config.getOptimTikhonovCoeff();
  tikhonov_use_x0 = config.getOptimTikhonovUseX0();

  // Get tolerance settings
  tol_grad_abs = config.getOptimTolGradAbs();
  tol_final_cost = config.getOptimTolFinalCost();
  tol_infidelity = config.getOptimTolInfidelity();
  tol_grad_rel = config.getOptimTolGradRel();
  maxiter = config.getOptimMaxiter();

  // Get penalty settings
  gamma_penalty_leakage = config.getOptimPenaltyLeakage();
  gamma_penalty_weightedcost = config.getOptimPenaltyWeightedCost();
  gamma_penalty_energy = config.getOptimPenaltyEnergy();
  gamma_penalty_dpdm = config.getOptimPenaltyDpdm();
  gamma_penalty_variation = config.getOptimPenaltyVariation();

  if (gamma_penalty_dpdm > 1e-13 && mastereq->decoherence_type != DecoherenceType::NONE){
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Disabling DpDm penalty term because it is not implemented for the Lindblad solver.\n");
    }
    gamma_penalty_dpdm = 0.0;
  }

  /* Store optimization bounds */
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xlower);
  VecSetFromOptions(xlower);
  VecDuplicate(xlower, &xupper);
  int col = 0;
  for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++){
    // Drive bounds (existing p/q controls)
    double drive_bound = config.getControlAmplitudeBound(iosc);
    drive_bound = drive_bound / (sqrt(2) * mastereq->getOscillator(iosc)->getNCarrierfrequencies());
    drive_bound = drive_bound * 2.0 * M_PI;
    for (size_t i = 0; i < mastereq->getOscillator(iosc)->getNDriveParams(); i++) {
      VecSetValue(xupper, col + i, drive_bound, INSERT_VALUES);
      VecSetValue(xlower, col + i, -1.0 * drive_bound, INSERT_VALUES);
    }
    col += mastereq->getOscillator(iosc)->getNDriveParams();

    // Flux bounds (independent f controls)
    double flux_bound = config.getControlFluxAmplitudeBound(iosc) * 2.0 * M_PI;
    for (size_t i = 0; i < mastereq->getOscillator(iosc)->getNFluxParams(); i++) {
      VecSetValue(xupper, col + i, flux_bound, INSERT_VALUES);
      VecSetValue(xlower, col + i, -1.0 * flux_bound, INSERT_VALUES);
    }
    col += mastereq->getOscillator(iosc)->getNFluxParams();
  }
  VecAssemblyBegin(xlower); VecAssemblyEnd(xlower);
  VecAssemblyBegin(xupper); VecAssemblyEnd(xupper);

  /* Create TAO optimization solver */
  TaoCreate(PETSC_COMM_SELF, &tao);
  /* Set optimization type and parameters */
  TaoSetType(tao,TAOBQNLS);         // Optim type: taoblmvm vs BQNLS ??
  TaoSetMaximumIterations(tao, maxiter);
  TaoSetTolerances(tao, tol_grad_abs, PETSC_DEFAULT, tol_grad_rel);
  TaoMonitorSet(tao, TaoMonitor, (void*)this, NULL);
  TaoSetVariableBounds(tao, xlower, xupper);
  TaoSetFromOptions(tao);
  /* Set user-defined objective and gradient evaluation routines */
  TaoSetObjective(tao, TaoEvalObjective, (void *)this);
  TaoSetGradient(tao, NULL, TaoEvalGradient,(void *)this);
  TaoSetObjectiveAndGradient(tao, NULL, TaoEvalObjectiveAndGradient, (void*) this);

  /* Allocate auxiliary vector */
  mygrad = new double[ndesign];

  /* Allocat xinit, xtmp */
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xinit);
  VecSetFromOptions(xinit);
  VecZeroEntries(xinit);
  VecCreateSeq(PETSC_COMM_SELF, ndesign, &xtmp);
  VecSetFromOptions(xtmp);
  VecZeroEntries(xtmp);
  VecDuplicate(xinit, &x_GN);
  VecZeroEntries(x_GN);
  VecDuplicate(xinit, &xprev);
  VecZeroEntries(xprev);
  VecDuplicate(xinit, &prev_solution_GN);
  VecZeroEntries(prev_solution_GN);
  VecAssemblyBegin(prev_solution_GN); VecAssemblyEnd(prev_solution_GN);

  /* Create MatShell for Gauss-Newton A=L^*L */
  MatCreateShell(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ndesign, ndesign, this, &GaussNewtonMatShell);
  MatShellSetOperation(GaussNewtonMatShell, MATOP_MULT, (void(*) (void)) applyGaussNewtonMatShell);
  VecDuplicate(xinit, &xeval_GN);
  VecZeroEntries(xeval_GN);
  VecAssemblyBegin(xeval_GN); VecAssemblyEnd(xeval_GN);

  // Include hessian of generalized J_inf wrt U in GaussNewton Approximation (for Jtrace only)
  if (optim_solver_type == OptimSolverType::GAUSS_NEWTON){
    if (optim_target->getObjectiveType() == ObjectiveType::JTRACE) {
      includeHessUJ = true;
    }
  }

  /* Create linear solver for solving Gauss-Newton Ax=b */
  KSPCreate(PETSC_COMM_SELF, &ksp_GN);
  KSPSetOperators(ksp_GN, GaussNewtonMatShell, GaussNewtonMatShell);
  KSPSetOptionsPrefix(ksp_GN, "gn_");  // Prefix for Gauss-Newton specific options

  // Set KSP type from config (default "cg", can be overridden by -gn_ksp_type)
  std::string ksp_type = config.getOptimGnKspType();
  KSPSetType(ksp_GN, ksp_type.c_str());

  // Set warm-start option from config (default false)
  bool warmstart = config.getOptimGnKspWarmstart();
  KSPSetInitialGuessNonzero(ksp_GN, warmstart ? PETSC_TRUE : PETSC_FALSE);
  KSPSetTolerances(ksp_GN, config.getOptimGnKspRtol(), PETSC_DEFAULT, PETSC_DEFAULT, config.getOptimGnKspMaxiter());

  // Configure preconditioner from config (default "none", can be overridden by -gn_pc_type)
  PC  pc;
  KSPGetPC(ksp_GN, &pc);
  std::string pc_type_str = config.getOptimGnPcType();
  PCSetType(pc, pc_type_str.c_str());

  // Allow command-line options to override defaults (must be called AFTER setting defaults)
  KSPSetFromOptions(ksp_GN);

  // If user requested shell PC (either from config or command-line), register the Jacobi functions
  PCType pc_type;
  PCGetType(pc, &pc_type);
  if (pc_type && strcmp(pc_type, PCSHELL) == 0) {
    PCShellSetContext(pc, this);
    PCShellSetSetUp(pc, jacobiShellSetup);
    PCShellSetApply(pc, jacobiShellApply);
  }

  /* Create eigenvalues solver for Gauss-Newton Ax=b */
#ifdef WITH_SLEPC
  EPSCreate(PETSC_COMM_SELF, &eps_GN);
  EPSSetOperators(eps_GN, GaussNewtonMatShell, NULL);
  EPSSetProblemType(eps_GN, EPS_HEP); // Hermitian
  EPSSetWhichEigenpairs(eps_GN, EPS_LARGEST_REAL); // largest eigenvalues
  neigvals = mastereq->getDim()*mastereq->getDim() - 1;
  ncv = neigvals + 2; // Max Krylov dimension. How to set??
  EPSSetDimensions(eps_GN, neigvals, ncv, PETSC_DEFAULT);
  EPSSetTolerances(eps_GN, eps_tol, eps_maxiter);
  EPSSetFromOptions(eps_GN);
#endif
}


OptimProblem::~OptimProblem() {
  delete [] mygrad;
  VecDestroy(&rho_t0_bar);

  VecDestroy(&xlower);
  VecDestroy(&xupper);
  VecDestroy(&xinit);
  VecDestroy(&xtmp);
  VecDestroy(&x_GN);
  VecDestroy(&xprev);

  MatDestroy(&GaussNewtonMatShell);
  VecDestroy(&xeval_GN);
  if (diag_GN != NULL) VecDestroy(&diag_GN);
  if (prev_solution_GN != NULL) VecDestroy(&prev_solution_GN);
  KSPDestroy(&ksp_GN);
#ifdef WITH_SLEPC
  EPSDestroy(&eps_GN);
#endif
  TaoDestroy(&tao);
}



double OptimProblem::evalF(const Vec x, bool writeTrajectoryDataFiles) {
  if (mpirank_world == 0 && !quietmode) printf("EVAL F... \n");

  /* Pass design vector x to oscillators */
  mastereq->setControlAmplitudes(x); 

  // Reset storage of final states in the target
  optim_target->resetFinalStates();

  /*  Iterate over initial condition */
  obj_cost  = 0.0;
  obj_regul = 0.0;
  obj_penal_leakage = 0.0;
  obj_penal_weightedcost = 0.0;
  obj_penal_dpdm = 0.0;
  obj_penal_energy = 0.0;
  obj_penal_variation = 0.0;
  fidelity = 0.0;
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  double fidelity_re = 0.0;
  double fidelity_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
      
    /* Prepare the initial condition in [rank * ninit_local, ... , (rank+1) * ninit_local - 1] */
    int initid = optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);

    /* Run forward with initial condition initid */
    if (mpirank_optim == 0 && !quietmode) printf("%d: Initial condition id=%d ...\n", mpirank_init, initid);
    Vec finalstate = timestepper->solveODE(initid, iinit, optim_target->getInitialState(), writeTrajectoryDataFiles, false);

    /* Store the final state for Geodesic Riemannian objective function */
    optim_target->storeFinalUnitaryColumn(iinit_global, finalstate);

    /* Add to leakage penalty term */
    obj_penal_leakage += obj_weights[iinit_global] * gamma_penalty_leakage * timestepper->getLeakageIntegral();

    /* Add to running cost penalty term */
    obj_penal_weightedcost += obj_weights[iinit_global] * gamma_penalty_weightedcost * timestepper->getWeightedCostIntegral();

    /* Add to second derivative penalty term */
    obj_penal_dpdm += obj_weights[iinit_global] * gamma_penalty_dpdm * timestepper->getDPDMIntegral();
    
    /* Add to energy integral penalty term */
    obj_penal_energy += obj_weights[iinit_global] * gamma_penalty_energy* timestepper->getEnergyIntegral();

    /* Evaluate J(finalstate) and add to final-time cost */
    double obj_iinit_re = 0.0;
    double obj_iinit_im = 0.0;
    optim_target->evalJ(finalstate,  &obj_iinit_re, &obj_iinit_im);
    obj_cost_re += obj_weights[iinit_global] * obj_iinit_re;
    obj_cost_im += obj_weights[iinit_global] * obj_iinit_im;

    /* Add to final-time fidelity */
    double fidelity_iinit_re = 0.0;
    double fidelity_iinit_im = 0.0;
    optim_target->HilbertSchmidtOverlap(finalstate, false, &fidelity_iinit_re, &fidelity_iinit_im);
    fidelity_re += 1./ ninit * fidelity_iinit_re;
    fidelity_im += 1./ ninit * fidelity_iinit_im;

    // printf("%d, %d: iinit obj_iinit: %f * (%1.14e + i %1.14e, Overlap=%1.14e + i %1.14e\n", mpirank_world, mpirank_init, obj_weights[iinit_global], obj_iinit_re, obj_iinit_im, fidelity_iinit_re, fidelity_iinit_im);
  }

  /* Sum up from initial conditions processors */
  double mypen_leak = obj_penal_leakage;
  double mypen_wcost = obj_penal_weightedcost;
  double mypen_dpdm = obj_penal_dpdm;
  double mypenen = obj_penal_energy;
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  double myfidelity_re = fidelity_re;
  double myfidelity_im = fidelity_im;
  MPI_Allreduce(&mypen_leak, &obj_penal_leakage, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_wcost, &obj_penal_weightedcost, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_dpdm, &obj_penal_dpdm, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypenen, &obj_penal_energy, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_re, &fidelity_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_im, &fidelity_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  /* Set the fidelity: If Schroedinger, need to compute the absolute value: Fid= |\sum_i \phi^\dagger \phi_target|^2 */
  if (mastereq->decoherence_type == DecoherenceType::NONE) {
    fidelity = pow(fidelity_re, 2.0) + pow(fidelity_im, 2.0);
  } else {
    fidelity = fidelity_re; 
  }
 
  /* Finalize the objective function */
  obj_cost = optim_target->finalizeJ(obj_cost_re, obj_cost_im, comm_init);

  /* Evaluate Tikhonov regularization term: gamma/2 * ||x-x0||^2*/
  double xnorm;
  if (!tikhonov_use_x0){  // ||x||^2
    VecNorm(x, NORM_2, &xnorm);
  } else {
    VecCopy(x, xtmp);
    VecAXPY(xtmp, -1.0, xinit);    // xtmp =  x - x_0
    VecNorm(xtmp, NORM_2, &xnorm);
  }
  obj_regul = gamma_tikhonov / 2. * pow(xnorm,2.0);

  /* Evaluate penality term for control variation */
  double var_reg = 0.0;
  for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++){
    var_reg += mastereq->getOscillator(iosc)->evalControlVariation(); // uses Oscillator::params instead of 'x'
  }
  obj_penal_variation = 0.5*gamma_penalty_variation*var_reg; 

  /* Sum, store and return objective value */
  objective = obj_cost + obj_regul + obj_penal_leakage + obj_penal_dpdm + obj_penal_energy + obj_penal_variation + obj_penal_weightedcost;

  /* Output */
  if (mpirank_world == 0 && !quietmode) {
    std::cout<< "Objective = " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal_leakage << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation << " + " << obj_penal_weightedcost << std::endl;
    std::cout<< "Fidelity = " << fidelity  << std::endl;
  }

  return objective;
}



void OptimProblem::evalGradF(const Vec x, Vec G, bool writeTrajectoryDataFiles){
  if (mpirank_world == 0 && !quietmode) std::cout<< "EVAL GRAD F... " << std::endl;

  /* Pass design vector x to oscillators */
  mastereq->setControlAmplitudes(x); 

  // DEBUG
  // output->writeControl(x, mastereq, timestepper->ntime, timestepper->dt);
  output->writeControlParams(x);

  /* Reset Gradient */
  VecZeroEntries(G);

  // Reset U_final
  optim_target->resetFinalStates();

  /* Derivative of regulatization terms (ADD ON ONE PROC ONLY!) */
  // if (mpirank_init == 0 && mpirank_optim == 0) { // TODO: Which one?? 
  if (mpirank_init == 0 ) {

    // Derivative of Tikhonov 0.5 * gamma * ||x||^2 
    VecAXPY(G, gamma_tikhonov, x);   // + gamma * x
    if (tikhonov_use_x0){
      VecAXPY(G, -1.0*gamma_tikhonov, xinit); // -gamma * xinit
    }

    // Derivative of penalization of control variation 
    double var_reg_bar = 0.5*gamma_penalty_variation;
    int skip_to_oscillator = 0;
    for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++){
      Oscillator* osc = mastereq->getOscillator(iosc);
      osc->evalControlVariationDiff(G, var_reg_bar, skip_to_oscillator);
      skip_to_oscillator += osc->getNParams();
    }
  }

  /*  Iterate over initial condition */
  obj_cost = 0.0;
  obj_regul = 0.0;
  obj_penal_leakage = 0.0;
  obj_penal_weightedcost = 0.0;
  obj_penal_dpdm = 0.0;
  obj_penal_energy = 0.0;
  obj_penal_variation = 0.0;
  fidelity = 0.0;
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  double fidelity_re = 0.0;
  double fidelity_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    // printf("%d: Initial condition id=%d ...\n", mpirank_init, iinit_global);

    /* Prepare the initial and target state */
    int initid = optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);

    /* --- Solve primal --- */
    // if (mpirank_optim == 0) printf("%d: %d FWD. ", mpirank_init, initid);

    /* Run forward with initial condition */
    Vec finalstate = timestepper->solveODE(initid, iinit, optim_target->getInitialState(), writeTrajectoryDataFiles, true);

    /* Store the final state for Riemannian objective function */
    optim_target->storeFinalUnitaryColumn(iinit_global, finalstate);

    /* Add to leakage penalty term */
    obj_penal_leakage += obj_weights[iinit_global] * gamma_penalty_leakage * timestepper->getLeakageIntegral();

    /* Add to running cost penalty term */
    obj_penal_weightedcost += obj_weights[iinit_global] * gamma_penalty_weightedcost * timestepper->getWeightedCostIntegral();

    /* Add to second derivative dpdm integral penalty term */
    obj_penal_dpdm += obj_weights[iinit_global] * gamma_penalty_dpdm * timestepper->getDPDMIntegral();
    /* Add to energy integral penalty term */
    obj_penal_energy += obj_weights[iinit_global] * gamma_penalty_energy * timestepper->getEnergyIntegral();

    /* Evaluate J(finalstate) and add to final-time cost */
    double obj_iinit_re = 0.0;
    double obj_iinit_im = 0.0;
    optim_target->evalJ(finalstate,  &obj_iinit_re, &obj_iinit_im);
    obj_cost_re += obj_weights[iinit_global] * obj_iinit_re;
    obj_cost_im += obj_weights[iinit_global] * obj_iinit_im;

    /* Add to final-time fidelity */
    double fidelity_iinit_re = 0.0;
    double fidelity_iinit_im = 0.0;
    optim_target->HilbertSchmidtOverlap(finalstate, false, &fidelity_iinit_re, &fidelity_iinit_im);
    fidelity_re += 1./ ninit * fidelity_iinit_re;
    fidelity_im += 1./ ninit * fidelity_iinit_im;
  }

  /* Sum up from initial conditions processors */
  double mypen_leak = obj_penal_leakage;
  double mypen_wcost = obj_penal_weightedcost;
  double mypen_dpdm = obj_penal_dpdm;
  double mypenen = obj_penal_energy;
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  double myfidelity_re = fidelity_re;
  double myfidelity_im = fidelity_im;
  MPI_Allreduce(&mypen_leak, &obj_penal_leakage, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_wcost, &obj_penal_weightedcost, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypen_dpdm, &obj_penal_dpdm, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mypenen, &obj_penal_energy, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_re, &fidelity_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&myfidelity_im, &fidelity_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  /* Set the fidelity: If Schroedinger, need to compute the absolute value: Fid= |\sum_i \phi^\dagger \phi_target|^2 */
  if (mastereq->decoherence_type == DecoherenceType::NONE) {
    fidelity = pow(fidelity_re, 2.0) + pow(fidelity_im, 2.0);
  } else {
    fidelity = fidelity_re; 
  }
 
  /* Finalize the objective function Jtrace to get the infidelity. 
     If Schroedingers solver, need to take the absolute value */
  obj_cost = optim_target->finalizeJ(obj_cost_re, obj_cost_im, comm_init);

  /* Evaluate Tikhonov regularization term += gamma/2 * ||x||^2*/
  double xnorm;
  if (!tikhonov_use_x0){  // ||x||^2
    VecNorm(x, NORM_2, &xnorm);
  } else {
    VecCopy(x, xtmp);
    VecAXPY(xtmp, -1.0, xinit);    // xtmp =  x_k - x_0
    VecNorm(xtmp, NORM_2, &xnorm);
  }
  obj_regul = gamma_tikhonov / 2. * pow(xnorm,2.0);

  /* Evaluate penalty term for control parameter variation */
  double var_reg = 0.0;
  for (size_t iosc = 0; iosc < mastereq->getNOscillators(); iosc++){
    var_reg += mastereq->getOscillator(iosc)->evalControlVariation(); // uses Oscillator::params instead of 'x'
  }
  obj_penal_variation = 0.5*gamma_penalty_variation*var_reg; 

  /* Sum, store and return objective value */
  objective = obj_cost + obj_regul + obj_penal_leakage + obj_penal_dpdm + obj_penal_energy + obj_penal_variation + obj_penal_weightedcost;

  double obj_cost_re_bar, obj_cost_im_bar;
  optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);

  /* Solve adjoint equations for all initial conditions . */
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;

    /* Recompute the initial state and target */
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
   
    /* Reset adjoint */
    VecZeroEntries(rho_t0_bar);

    /* Terminal condition for adjoint variable: Derivative of final time objective J */
    optim_target->evalJ_diff(timestepper->getFinalState(iinit), rho_t0_bar, obj_weights[iinit_global]*obj_cost_re_bar, obj_weights[iinit_global]*obj_cost_im_bar);

    // Derivative of storing final unitary (adds a column of U_final_bar into rho_t0_bar)
    optim_target->storeFinalUnitaryColumn_diff(iinit_global, rho_t0_bar);

    /* Derivative of time-stepping */
    timestepper->solveAdjointODE(iinit, rho_t0_bar, obj_weights[iinit_global] * gamma_penalty_leakage, obj_weights[iinit_global]*gamma_penalty_weightedcost, obj_weights[iinit_global]*gamma_penalty_dpdm, obj_weights[iinit_global]*gamma_penalty_energy);

    /* Add to optimizers's gradient */
    VecAXPY(G, 1.0, timestepper->getReducedGradient());
  } // end of initial condition loop 

  /* Sum up the gradient from all initial condition processors */
  PetscScalar* grad; 
  VecGetArray(G, &grad);
  for (int i=0; i<ndesign; i++) {
    mygrad[i] = grad[i];
  }
  MPI_Allreduce(mygrad, grad, ndesign, MPI_DOUBLE, MPI_SUM, comm_init);
  VecRestoreArray(G, &grad);

  /* Compute and store gradient norm */
  VecNorm(G, NORM_2, &(gnorm));

  /* Output */
  // if (mpirank_world == 0 && !quietmode) {
  //   std::cout<< "Objective = " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal_leakage << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation << " + " << obj_penal_weightedcost <<  std::endl;
  //   std::cout<< "Fidelity = " << fidelity << std::endl;
  // }
}



void OptimProblem::evalLinearizedForward(const Vec x, const Vec v){
  // if (mpirank_world == 0 && !quietmode) std::cout<< "EVAL LINEARIZED FWD ... " << std::endl;

  /* Pass design vector x to oscillators */
  mastereq->setControlAmplitudes(x); 
 
  /* Solve ODE and linearized ODE forward in time */
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    // printf("Solving ODE for initial condition %d (global index %d), total ninit_local = %d\n", iinit, iinit_global, ninit_local);

    int initid = optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);

    // Nonlinear forward at x is identical for every MatVec within the same KSP/EPS solve, so only
    // (re-)solve and store it once per xeval_GN; reuse the stored trajectory_states otherwise.
    if (!nonlinear_forward_valid) {
      bool writeTrajectoryDataFiles = false;
      bool storeStates = true;
      timestepper->solveODE(initid, iinit, optim_target->getInitialState(), writeTrajectoryDataFiles, storeStates);
    }

    // Solve linearized forward ODE in direction v while storing linearized states
    bool storeLinearizedStates = true;
    timestepper->solveLinearizedODE(iinit, v, storeLinearizedStates); 
  }
  nonlinear_forward_valid = true;
}

void OptimProblem::applyGaussNewtonMatShell(Mat A, const Vec v, Vec Av){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);
  // if (self->mpirank_world == 0) printf("APPLYING GAUSS-NEWTON...\n");
  self->GN_MatVec_counter++;

  // Grab the point of evaluation from the shell
  Vec x = self->xeval_GN;

  //  Reset output 
  VecZeroEntries(Av);
  
  // Apply linearized forward to get dU/dalpha x v
  // assumes that the timestepper's trajectory_states are already populated and computes all lin_trajectory_states.
  self->evalLinearizedForward(x, v);

  // Optionally, include \nabla^2_U J(U) in the terminal adjoint condition. 
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  if (self->includeHessUJ){
    // Only available for J_Inf (Generalized)
    assert(self->optim_target->getObjectiveType() == ObjectiveType::JTRACE);

    // Compute J_inf(W(T))
    for (int iinit=0; iinit<self->ninit_local; iinit++) {
      // Need to recompute the target (and initial) state
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      // eval J_inf(w(T))
      Vec wT = self->timestepper->getLinearizedFinalState(iinit);
      double obj_lin_iinit_re = 0.0;
      double obj_lin_iinit_im = 0.0;
      self->optim_target->evalJ(wT,  &obj_lin_iinit_re, &obj_lin_iinit_im);
      obj_cost_re += self->obj_weights[iinit] * obj_lin_iinit_re;
      obj_cost_im += self->obj_weights[iinit] * obj_lin_iinit_im;
    }
    // Gather result J_inf(W(T))
    double mycost_re = obj_cost_re;
    double mycost_im = obj_cost_im;
    MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
    MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
    
  }

  // Solve the adjoint ODE for initial condition 
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Set terminal condition for adjoint
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);

    // Set the terminal condition 
    VecCopy(lin_final_state, self->rho_t0_bar);  // w_i(T)

    if (self->includeHessUJ){
      // scale w_i(t) by 2/n (contribution from 1/n||U||^2)
      VecScale(self->rho_t0_bar, 2.0 / self->mastereq->getDim());
      // Add derivative of Jinf -2/n^2 Re tr(...)
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      double obj_cost_re_bar, obj_cost_im_bar;
      self->optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
      self->optim_target->evalJ_diff(lin_final_state, self->rho_t0_bar, self->obj_weights[iinit]*obj_cost_re_bar, self->obj_weights[iinit]*obj_cost_im_bar); 
    }

    // Solve adjoint backward ODE
    self->timestepper->solveAdjointODE(iinit, self->rho_t0_bar, 0.0, 0.0, 0.0, 0.0);

    // Add gradient to output
    VecAXPY(Av, 1.0, self->timestepper->getReducedGradient());
  }

  /* Sum up the gradient from all initial condition processors */
  PetscScalar* Av_data;
  VecGetArray(Av, &Av_data);
  MPI_Allreduce(MPI_IN_PLACE, Av_data, self->ndesign, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(Av, &Av_data);

  // Add damping term: Av += damping * v (for Levenberg-Marquardt regularization)
  if (self->ksp_damping > 0.0) {
    VecAXPY(Av, self->ksp_damping, v);
  }

}


PetscErrorCode OptimProblem::jacobiShellSetup(PC pc){
  OptimProblem *self;
  PCShellGetContext(pc, (void**)&self);

  // If diagonal already computed, skip
  if (self->diag_GN != NULL) return 0;

  // Create diagonal vector
  VecDuplicate(self->xeval_GN, &self->diag_GN);

  // Get the size
  PetscInt n;
  VecGetSize(self->diag_GN, &n);

  // Create temporary vectors
  Vec e_i, Ae_i;
  VecDuplicate(self->diag_GN, &e_i);
  VecDuplicate(self->diag_GN, &Ae_i);

  // For each parameter i, compute diag(A)_i = (A*e_i)_i
  for (PetscInt i = 0; i < n; i++) {
    // Create unit vector e_i
    VecZeroEntries(e_i);
    VecSetValue(e_i, i, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e_i);
    VecAssemblyEnd(e_i);

    // Compute A * e_i (includes damping)
    applyGaussNewtonMatShell(self->GaussNewtonMatShell, e_i, Ae_i);

    // Extract diagonal element (A*e_i)_i
    PetscScalar diag_i;
    VecGetValues(Ae_i, 1, &i, &diag_i);
    VecSetValue(self->diag_GN, i, diag_i, INSERT_VALUES);
  }

  VecAssemblyBegin(self->diag_GN);
  VecAssemblyEnd(self->diag_GN);

  // Cleanup
  VecDestroy(&e_i);
  VecDestroy(&Ae_i);

  return 0;
}


PetscErrorCode OptimProblem::jacobiShellApply(PC pc, Vec x, Vec y){
  OptimProblem *self;
  PCShellGetContext(pc, (void**)&self);

  // Apply Jacobi preconditioner: y = x ./ diag
  VecPointwiseDivide(y, x, self->diag_GN);

  return 0;
}


void OptimProblem::solveGaussNewtonKSP(Vec xinit, const Vec b, Vec Ainv_b, int outer_iter){

  // Store the point of evaluation for the Gauss-Newton matrix shell A(xinit)
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Set the matrix again, just in case, for reset.
  KSPSetOperators(ksp_GN, GaussNewtonMatShell, GaussNewtonMatShell);

  // Set initial guess: either warm-start with previous solution or use zero
  PetscBool warmstart;
  KSPGetInitialGuessNonzero(ksp_GN, &warmstart);
  if (warmstart) {
    VecCopy(prev_solution_GN, Ainv_b);
    VecAssemblyBegin(Ainv_b); VecAssemblyEnd(Ainv_b);
  } else {
    VecZeroEntries(Ainv_b);
  }

  // Monitor residual and solution norm at every iteration
  KSPMonitorCancel(ksp_GN);
  KSPMonitorSet(ksp_GN, KSPMonitorResidualAndSolution, (void*)this, NULL);

  // Open file for KSP convergence history
  char ksp_filename[256] = "";
  if (mpirank_world == 0) {
    // Get KSP type and tolerances
    KSPType ksp_type;
    PetscReal rtol, atol, dtol;
    PetscInt maxits;
    KSPGetType(ksp_GN, &ksp_type);
    KSPGetTolerances(ksp_GN, &rtol, &atol, &dtol, &maxits);

    // Create filename with KSP type and iteration number
    if (outer_iter >= 0) {
      snprintf(ksp_filename, sizeof(ksp_filename), "ksp_convergence_%s_iter%04d.dat", ksp_type, outer_iter);
    } else {
      snprintf(ksp_filename, sizeof(ksp_filename), "ksp_convergence_%s.dat", ksp_type);
    }

    ksp_history_file = fopen(ksp_filename, "w");
    if (ksp_history_file != NULL) {
      fprintf(ksp_history_file, "# KSP convergence history\n");
      fprintf(ksp_history_file, "# KSP_TYPE: %s\n", ksp_type);
      fprintf(ksp_history_file, "# RTOL: %e\n", rtol);
      fprintf(ksp_history_file, "# MAXITER: %d\n", maxits);
      fprintf(ksp_history_file, "# iteration  residual_norm  solution_norm\n");
    }
  }

  // Solve the linear system (L^*L + damping*I) x = b
  // Note: damping is applied inside applyGaussNewtonMatShell and getGaussNewtonDiagonal
  GN_MatVec_counter = 0;
  KSPSolve(ksp_GN, b, Ainv_b);

  // Report convergence and write to file
  KSPConvergedReason reason;
  int iters;
  double rnorm;
  KSPGetConvergedReason(ksp_GN, &reason);

  // Write convergence reason to file
  if (mpirank_world == 0 && ksp_history_file != NULL) {
    const char* reason_str = KSPConvergedReasons[reason];
    fprintf(ksp_history_file, "# CONVERGED_REASON: %s\n", reason_str);
  }

  // Close KSP history file and report
  if (mpirank_world == 0 && ksp_history_file != NULL) {
    fclose(ksp_history_file);
    ksp_history_file = NULL;
    if (!quietmode) {
      printf("KSP convergence history written to: %s\n", ksp_filename);
    }
  }

  KSPGetIterationNumber(ksp_GN, &iters);
  KSPGetResidualNorm(ksp_GN, &rnorm);
  ksp_iters_last = iters;
  if (mpirank_world == 0 && !quietmode) {
    KSPType ksp_type_used;
    KSPGetType(ksp_GN, &ksp_type_used);
    printf("Gauss-Newton %s stats: iterations = %d, MatVec counter = %d, residual norm = %1.14e\n", ksp_type_used, iters, GN_MatVec_counter, rnorm);
  }

  // Save solution for warm-starting next iteration
  VecCopy(Ainv_b, prev_solution_GN);
  VecAssemblyBegin(prev_solution_GN); VecAssemblyEnd(prev_solution_GN);
}


#ifdef WITH_SLEPC
void OptimProblem::solveGaussNewtonEPS(Vec xinit, const Vec b, Vec Ainv_b){

  // Get eigenvalues and eigenvectors of the Gauss-Newton matrix
  Mat evecs;
  std::vector<double> evals = computeGaussNewtonEvals(xinit, &evecs);

  // Print the eigenvalues
  // for (int i=0; i<evals.size(); i++) {
  //   if (mpirank_world == 0) printf("Eigenvalue %d: %1.14e\n", i, evals[i]);
  // }

  // Project the rhs onto evals: b_proj = evals^Tb
  Vec tmp;
  MatCreateVecs(evecs, &tmp, NULL);
  MatMultTranspose(evecs, b, tmp);

  // Apply inverse eigenvalue scaling: tmp = evals^-1 * tmp
  PetscScalar* tmp_data;
  VecGetArray(tmp, &tmp_data);
  for (int i=0; i<evals.size(); i++) {
    if (evals[i] > evals_cutoff) {
      tmp_data[i] /= evals[i];
    } else {
      tmp_data[i] = 0.0;
    }
  }
  VecRestoreArray(tmp, &tmp_data);

  // Transform back to the original space: Ainv_b = evecs * tmp
  MatMult(evecs, tmp, Ainv_b);

  MatDestroy(&evecs);
  VecDestroy(&tmp);
}
#endif

PetscErrorCode KSPMonitorResidualAndSolution(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx){
  OptimProblem* self = (OptimProblem*) ctx;

  Vec x;
  KSPBuildSolution(ksp, NULL, &x); // Current iterate, valid at this point in the KSP solve
  PetscReal xnorm;
  VecNorm(x, NORM_2, &xnorm);

  if (self->getMPIrank_world() == 0 && !self->getQuietmode()) {
    printf("KSP it %d: residual norm = %1.14e, solution norm = %1.14e\n", (int)it, (double)rnorm, (double)xnorm);
  }

  // Write to file if file pointer is set
  if (self->getMPIrank_world() == 0 && self->ksp_history_file != NULL) {
    fprintf(self->ksp_history_file, "%d %1.14e %1.14e\n", (int)it, (double)rnorm, (double)xnorm);
  }

  return 0;
}

#ifdef WITH_SLEPC
std::vector<double> OptimProblem::computeGaussNewtonEvals(Vec xinit, Mat* evecs_out){

  // Store xinit so the MatShell can use it as point of evaluation.
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // CAREFUL: EPS solver might need a reset if called multiple times!?

  // Solve the eigenvalue problem for the Gauss-Newton matrix
  GN_MatVec_counter = 0;
  EPSSolve(eps_GN);
  PetscInt numConv;
  PetscInt iters_taken;
  EPSGetConverged(eps_GN, &numConv);
  EPSGetIterationNumber(eps_GN,&iters_taken);
  if (mpirank_world == 0 && !quietmode) printf("Gauss-Newton EPS converged %d eigenvalues in %d iterations. MatVec counter = %d\n", numConv, iters_taken, GN_MatVec_counter);
  if (numConv < neigvals) {
      if (mpirank_world==0 && !quietmode) printf("WARNING: Only %d eigenvalues out of %d eigenvalues converged.\n", numConv, neigvals);
  }

  // Set up storage for eigenvalues and eigenvectors (should be real!)
  std::vector<double> evals_re(neigvals);
  std::vector<Vec> evec_re(neigvals);
  for (int ix = 0; ix < neigvals; ix++) {
    MatCreateVecs(GaussNewtonMatShell, &evec_re[ix], NULL);
  }

  // Retrieve eigenpairs of M and compute error. 
  for (PetscInt i = 0; i < numConv && i < neigvals; i++) {

    // Retrieve the eigenvalue (is real) and eigenvector
    EPSGetEigenpair(eps_GN, i, &evals_re[i], NULL, evec_re[i], NULL);
    // EPSGetEigenvalue(eps, i, &evals_re[i], NULL);

    // // Estimate the errror (needs one more application of A)
    // double error = 0.0;
    // EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
    // if (error > 1e-12) {
    //     if (mpirank_world==0) printf("WARNING: Relative error of eigenpair %d is large (error=%1.4e)\n", i, error);
    // }
  }

  // Resize to the number of converged eigenvalues. 
  PetscInt nconv = std::min(numConv, neigvals);
  evals_re.resize(nconv);

  // Assemble eigenvectors into a dense matrix, one eigenvector per column.
  MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ndesign, nconv, NULL, evecs_out);
  MatSetUp(*evecs_out);
  for (PetscInt col = 0; col < nconv; col++) {
    const PetscScalar *varr;
    VecGetArrayRead(evec_re[col], &varr);
    for (int row = 0; row < ndesign; row++) {
      MatSetValue(*evecs_out, row, col, varr[row], INSERT_VALUES);
    }
    VecRestoreArrayRead(evec_re[col], &varr);
  }
  MatAssemblyBegin(*evecs_out, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*evecs_out, MAT_FINAL_ASSEMBLY);

  // Cleanup
  for (int ix = 0; ix < neigvals; ix++) {
    VecDestroy(&evec_re[ix]);
  }

  return evals_re;
}
#endif


double OptimProblem::armijoLineSearch(Vec x, double f, Vec grad, Vec dir, Vec xnew, Vec step) {
  // Linesearch iteration with initial step size alpha = 1.0
  double alpha = 1.0;
  for (int ls = 0; ls < max_ls_iter; ls++) {

    // Trial point xnew = x - alpha*dir
    VecWAXPY(xnew, -alpha, dir, x);  

    // Projected onto the bound constraints
    VecPointwiseMax(xnew, xnew, xlower);
    VecPointwiseMin(xnew, xnew, xupper);

    // Actual (possibly clipped) step and its descent condition g^T step < 0
    VecWAXPY(step, -1.0, x, xnew); // step = xnew - x
    double gts;
    VecDot(grad, step, &gts);

    // Fall back to steepest descent if projection destroyed the descent property
    if (gts > 0.0) {
      // xnew = x - alpha * grad
      VecWAXPY(xnew, -alpha, grad, x); 
      // Projected onto the bound constraints
      VecPointwiseMax(xnew, xnew, xlower);
      VecPointwiseMin(xnew, xnew, xupper);
      // Recompute the step and its directional derivative after projection
      VecWAXPY(step, -1.0, x, xnew);
      VecDot(grad, step, &gts);
    }

    double fnew = evalF(xnew);
    if (fnew <= f + c1 * gts || ls == max_ls_iter - 1) {
      if (fnew > f + c1 * gts && mpirank_world == 0) {
        printf("Warning: Armijo line search did not find sufficient decrease within %d backtracks, accepting smallest step.\n", max_ls_iter);
      }
      return alpha;
    }
    alpha *= rho_backtrack;
  }

  return alpha; // unreachable
}

void OptimProblem::solve(Vec xinit) {

  switch (optim_solver_type) {
    case OptimSolverType::TAO_LBFGS:
      TaoSetSolution(tao, xinit);
      TaoSolve(tao);
      break;

    case OptimSolverType::GAUSS_NEWTON: {
      // Set the initial guess for the Gauss-Newton solver
      VecCopy(xinit, x_GN);

      // Work vectors for gradient, preconditioned search direction, and line search
      Vec G, Gprec, xnew, step;
      VecDuplicate(x_GN, &G);
      VecDuplicate(x_GN, &Gprec);
      VecDuplicate(x_GN, &xnew);
      VecDuplicate(x_GN, &step);

      bool stop = false;
      for (int iter = 0; !stop; iter++) {
        // Compute the gradient (and objective) at the current iterate
        evalGradF(x_GN, G);
        double f = objective;
        VecNorm(G, NORM_2, &gnorm);

        // Precondition the gradient: solve the Gauss-Newton system A(x)*Gprec = G via KSP
        solveGaussNewtonKSP(x_GN, G, Gprec, iter);
        // solveGaussNewtonEPS(x_GN, G, Gprec);
        // VecCopy(G, Gprec); // Steepest descent, no preconditioner

        // Backtracking Armijo line search along -Gprec, projected onto the bound constraints
        double alpha = armijoLineSearch(x_GN, f, G, Gprec, xnew, step);

        // Accept the step
        VecCopy(x_GN, xprev);
        VecCopy(xnew, x_GN);

        // Monitor progress and check stopping criteria
        stop = monitor(iter, f, gnorm, alpha);
      }

      VecDestroy(&G);
      VecDestroy(&Gprec);
      VecDestroy(&xnew);
      VecDestroy(&step);
      break;
    }

    default:
      printf("Unsupported optimization solver type.\n");
      break;
  }
}

void OptimProblem::getStartingPoint(Vec xinit){

  // Grab parameters from oscillators
  PetscScalar* xptr;
  VecGetArray(xinit, &xptr);
  int shift = 0;
  for (size_t ioscil = 0; ioscil<mastereq->getNOscillators(); ioscil++){
    mastereq->getOscillator(ioscil)->getControlParams(xptr + shift);
    shift += mastereq->getOscillator(ioscil)->getNParams();
  }
  VecRestoreArray(xinit, &xptr);
  
  /* Assemble initial guess */
  VecAssemblyBegin(xinit);
  VecAssemblyEnd(xinit);

  /* Pass to oscillator */
  mastereq->setControlAmplitudes(xinit);

  // Store it in the optimProblem
  VecCopy(xinit, this->xinit);
}


void OptimProblem::getSolution(Vec* xopt){
  
  /* Get ref to optimized parameters */
  if (optim_solver_type == OptimSolverType::TAO_LBFGS) {
    Vec params;
    TaoGetSolution(tao, &params);
    *xopt = params;
  } else if (optim_solver_type == OptimSolverType::GAUSS_NEWTON) {
    *xopt = x_GN;
  } else {
    printf("Unsupported optimization solver type.\n");
  }
}

bool OptimProblem::monitor(int iter, double f, double gnorm, double deltax){

  double F_avg = getFidelity();

  // Switch objective functions
  // if (1.0 - F_avg < 0.73) {
  // if (iter > 50) {
  //   printf("Switching to infidelity measure.\n");
  //   setRiemannianDistance(false);
  // }

  // // Freeze theta_avg if fidelity is sufficiently high
  // if (F_avg > 0.80) {
  //   getOptimTarget()->freeze_theta_avg = true;
  // }

  /* Additional Stopping criteria */
  std::string finalReason_str = "";
  if (1.0 - F_avg <= getTolInfidelity()) {
    finalReason_str = "Optimization converged with small infidelity.";
  // } else if (obj_cost <= getTolFinalCost()) {
  //   finalReason_str = "Optimization converged with small final time cost.";
  } else if (iter == getMaxIter()) {
    finalReason_str = "Optimization stopped at maximum number of iterations.";
  } else if (gnorm < getTolGradAbs()) {
    finalReason_str = "OPtimization converged with small gradient norm.";
  }
  bool lastIter = (finalReason_str.length() > 0);

  /* First iteration: Header for screen output of optimization history */
  if (iter == 0 && getMPIrank_world() == 0) {
    std::cout<<  "    Objective             Tikhonov               Penalty-Leakage        Penalty-StateVar       Penalty-TotalEnergy    Penalty-CtrlVar        Penalty-WeightedCost" << std::endl;
  }

  /* Every <output_optimization_stride> iterations: Output of optimization history */
  if (iter % getOutputOptimizationStride() == 0 || lastIter) {
    // Add to optimization history file 
    getOutput()->writeOptimFile(iter, f, gnorm, deltax, F_avg, obj_cost, obj_regul, obj_penal_leakage, obj_penal_dpdm, obj_penal_energy, obj_penal_variation, obj_penal_weightedcost, ksp_iters_last);
    // Screen output 
    if (getMPIrank_world() == 0) {
      std::cout<< iter <<  "  " << std::scientific<<std::setprecision(14) << obj_cost << " + " << obj_regul << " + " << obj_penal_leakage << " + " << obj_penal_dpdm << " + " << obj_penal_energy << " + " << obj_penal_variation << " + " << obj_penal_weightedcost;
      std::cout<< "  Fidelity = " << F_avg;
      std::cout<< "  ||Grad|| = " << gnorm;
      std::cout<< std::endl;
    }
  }

  /* Print last iteration stopping reason */
  if (lastIter && getMPIrank_world() == 0) {
    std::cout<< finalReason_str << std::endl;
  }

  return lastIter;
}

PetscErrorCode TaoMonitor(Tao tao,void*ptr){
  OptimProblem* ctx = (OptimProblem*) ptr;

  /* Get information from Tao optimization */
  PetscInt iter;
  PetscScalar deltax;
  TaoConvergedReason reason;
  PetscScalar f, gnorm;
  TaoGetSolutionStatus(tao, &iter, &f, &gnorm, NULL, &deltax, &reason);

  bool lastIter = ctx->monitor(iter, f, gnorm, deltax);
  if (lastIter) {
    TaoSetConvergedReason(tao, TAO_CONVERGED_USER);
  }

  return 0;
}


PetscErrorCode TaoEvalObjectiveAndGradient(Tao tao, Vec x, PetscReal *f, Vec G, void*ptr){

  TaoEvalGradient(tao, x, G, ptr);
  OptimProblem* ctx = (OptimProblem*) ptr;
  *f = ctx->getObjective();
  // *f = 1.0 - ctx->getFidelity();

  return 0;
}

PetscErrorCode TaoEvalObjective(Tao /*tao*/, Vec x, PetscReal *f, void*ptr){

  OptimProblem* ctx = (OptimProblem*) ptr;
  *f = ctx->evalF(x, false);
  
  return 0;
}


PetscErrorCode TaoEvalGradient(Tao /*tao*/, Vec x, Vec G, void*ptr){

  OptimProblem* ctx = (OptimProblem*) ptr;
  ctx->evalGradF(x, G, false);
  
  return 0;
}
