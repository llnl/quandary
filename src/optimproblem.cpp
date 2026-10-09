#include "optimproblem.hpp"
#include "gellmann.hpp"
#include <algorithm>

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

  /* Create MatShell for Gauss-Newton least-squares problem */
  // Create MatShell for GNLeastSquares solve
  PetscInt M = mastereq->getDim()*mastereq->getDim() - 1; 
  MatCreateShell(PETSC_COMM_SELF, M, ndesign, M, ndesign, this, &GNLeastSquaresShell);
  MatShellSetOperation(GNLeastSquaresShell, MATOP_MULT, (void(*)(void))GNLeastSquaresShell_MatMult);
  MatShellSetOperation(GNLeastSquaresShell, MATOP_MULT_TRANSPOSE, (void(*)(void))GNLeastSquaresShell_MatMultTranspose);
  MatShellSetOperation(GNLeastSquaresShell, MATOP_CREATE_VECS,     (void(*)(void))GNLeastSquaresShell_MatCreateVecs);


  /* Create MatShell for Gauss-Newton A=L^*L */
  MatCreateShell(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ndesign, ndesign, this, &GaussNewtonMatShell);
  MatShellSetOperation(GaussNewtonMatShell, MATOP_MULT, (void(*) (void)) applyGaussNewtonMatShell);
  VecDuplicate(xinit, &xeval_GN);
  VecZeroEntries(xeval_GN);
  VecAssemblyBegin(xeval_GN); VecAssemblyEnd(xeval_GN);

  /* Create dense matrix for Gauss-Newton */
  MatCreateDense(PETSC_COMM_SELF, ndesign, ndesign, ndesign, ndesign, NULL, &GaussNewtonMatDense);
  MatZeroEntries(GaussNewtonMatDense);

  // Include hessian of generalized J_inf wrt U in GaussNewton Approximation (for Jtrace only)
  includeHessUJ = false;
  if (optim_solver_type == OptimSolverType::GAUSS_NEWTON){
    if (optim_target->getObjectiveType() == ObjectiveType::JTRACE) {
      includeHessUJ = true;
    }
  }
  // DISABLE includeHessUJ: Using L^*L instead of L^*\nabla_U^2JL!
  includeHessUJ = false;

  // Create the KSP solver for the Gauss-Newton least-squares problem
  KSPCreate(PETSC_COMM_SELF, &ksp_LeastSquares);
  KSPSetOperators(ksp_LeastSquares, GNLeastSquaresShell, GNLeastSquaresShell);
  // KSPSetNormType(ksp_LeastSquares, KSP_NORM_UNPRECONDITIONED); // Unconditioned ressidual norm
  KSPSetType(ksp_LeastSquares, KSPLSQR);
  KSPSetTolerances(ksp_LeastSquares,config.getOptimKSPRtol(),PETSC_DEFAULT,PETSC_DEFAULT,config.getOptimKSPMaxiter());
  // KSPSetComputeSingularValues(ksp_LeastSquares, PETSC_TRUE);
  KSPSetFromOptions(ksp_LeastSquares);

  // Initialize Tao BRGN solver for least-squares problem if configured
  ls_solver = config.getLeastSquaresSolver();
  if (ls_solver == "BRGN") {
    // Create Tao solver on PETSC_COMM_SELF (same as ksp_LeastSquares)
    TaoCreate(PETSC_COMM_SELF, &tao_brgn);
    TaoSetType(tao_brgn, TAOBRGN);

    // Create storage for RHS vector (used in residual callback)
    GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, nullptr, &brgn_rhs);
    VecDuplicate(brgn_rhs, &brgn_residual);

    // Set convergence tolerances (reuse KSP settings)
    TaoSetTolerances(tao_brgn, config.getOptimKSPRtol(), PETSC_DEFAULT, PETSC_DEFAULT);
    TaoSetMaximumIterations(tao_brgn, config.getOptimKSPMaxiter());

    // Set damping parameter for regularization term λ||v||²
    brgn_damping = config.getBrgnDamping();
    TaoBRGNSetRegularizerWeight(tao_brgn, brgn_damping);

    // Set Residual routine F(v) = Lv - b 
    TaoSetResidualRoutine(tao_brgn, brgn_residual, TaoBRGN_EvalResidual, (void*)this);

    // Set Jacobian (the matrix L = GNLeastSquaresShell, constant for this problem)
    TaoSetJacobianResidualRoutine(tao_brgn, GNLeastSquaresShell, GNLeastSquaresShell,
                                   TaoBRGN_EvalJacobianResidual, (void*)this);

    // Set monitor function to track progress
    TaoMonitorSet(tao_brgn, TaoBRGN_Monitor, (void*)this, NULL);

    // Allow command-line overrides
    TaoSetFromOptions(tao_brgn);
  } else {
    tao_brgn = nullptr;
    brgn_rhs = nullptr;
    brgn_residual = nullptr;
  }

  // Cache state dimension for efficiency
  VecGetSize(rho_t0_bar, &state_dim_cached);

  // Pre-allocate workspace vectors for GNLeastSquaresShell operations
  wsub_workspace = new Vec[ninit_local];
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    VecDuplicate(rho_t0_bar, &wsub_workspace[iinit]);
  }

  /* Create linear solver for solving Gauss-Newton Ax=b */
  ksp_damping = config.getGnKspDamping();
  KSPCreate(PETSC_COMM_SELF, &ksp_GN);
  KSPSetOperators(ksp_GN, GaussNewtonMatShell, GaussNewtonMatShell);

  // Set KSP type from configuration (case-insensitive)
  std::string ksp_type_str = config.getGnKspType();
  std::string ksp_type_lower = ksp_type_str;
  std::transform(ksp_type_lower.begin(), ksp_type_lower.end(), ksp_type_lower.begin(), ::tolower);

  if (ksp_type_lower == "cg") {
    KSPSetType(ksp_GN, KSPCG);
  } else if (ksp_type_lower == "minres") {
    KSPSetType(ksp_GN, KSPMINRES);
  } else if (ksp_type_lower == "gmres") {
    KSPSetType(ksp_GN, KSPGMRES);
  } else {
    // Default or try to use the string directly for other types
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Unknown gn_ksp_type '%s', using it directly with PETSc\n", ksp_type_str.c_str());
    }
    KSPSetType(ksp_GN, ksp_type_str.c_str());
  }

  KSPSetInitialGuessNonzero(ksp_GN, PETSC_FALSE);

  // Enable MINRES-QLP if requested
  if (ksp_type_lower == "minres" && config.getGnMinresQlp()) {
    PetscBool set;
    PetscOptionsHasName(NULL, NULL, "-ksp_minres_qlp", &set);
    if (!set) {
      PetscOptionsSetValue(NULL, "-ksp_minres_qlp", NULL);
    }
  }

  KSPSetFromOptions(ksp_GN);
  PC  pc;
  KSPGetPC(ksp_GN, &pc);
  PCSetType(pc, PCNONE); // Disable preconditioner
  KSPSetNormType(ksp_GN, KSP_NORM_UNPRECONDITIONED); // Unconditioned ressidual norm
  KSPSetTolerances(ksp_GN,config.getOptimKSPRtol(),PETSC_DEFAULT,PETSC_DEFAULT,config.getOptimKSPMaxiter());

  /* Create eigenvalues solver for Gauss-Newton Ax=b */
  EPSCreate(PETSC_COMM_SELF, &eps_GN);
  EPSSetOperators(eps_GN, GaussNewtonMatShell, NULL);
  EPSSetProblemType(eps_GN, EPS_HEP); // Hermitian
  EPSSetWhichEigenpairs(eps_GN, EPS_LARGEST_REAL); // largest eigenvalues
  double eps_thresh = eps_evals_cutoff;
  EPSSetThreshold(eps_GN, eps_thresh, PETSC_FALSE);  // absolute threshold
  neigvals = mastereq->getDim()*mastereq->getDim() - 1;
  // ncv = neigvals + 2; // Max Krylov dimension. How to set??
  ncv = 2*neigvals ; // Max Krylov dimension. How to set??
  EPSSetDimensions(eps_GN, neigvals, ncv, PETSC_DEFAULT);
  EPSSetTolerances(eps_GN, eps_tol, config.getOptimKSPMaxiter());
  EPSSetFromOptions(eps_GN);

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

  MatDestroy(&GaussNewtonMatDense);
  MatDestroy(&GaussNewtonMatShell);
  MatDestroy(&GNLeastSquaresShell);
  VecDestroy(&xeval_GN);
  KSPDestroy(&ksp_GN);
  EPSDestroy(&eps_GN);
  KSPDestroy(&ksp_LeastSquares);
  if (tao_brgn) {
    TaoDestroy(&tao_brgn);
    VecDestroy(&brgn_rhs);
    VecDestroy(&brgn_residual);
  }
  TaoDestroy(&tao);

  // Clean up workspace vectors
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    VecDestroy(&wsub_workspace[iinit]);
  }
  delete[] wsub_workspace;
}


void OptimProblem::setGaussNewtonMaxiter(int maxiter) {
  if (mpirank_world == 0 && !quietmode) {
    printf("Setting max iterations for all Gauss-Newton solvers to: %d\n", maxiter);
  }

  // 1. Set for KSP solver (used in solveGaussNewtonKSP)
  PetscReal rtol, abstol, dtol;
  PetscInt current_maxits;
  KSPGetTolerances(ksp_GN, &rtol, &abstol, &dtol, &current_maxits);
  KSPSetTolerances(ksp_GN, rtol, abstol, dtol, maxiter);

  // 2. Set for Least Squares solver (used in solveGaussNewtonLeastSquares)
  if (ls_solver == "BRGN") {
    // For BRGN, this applies to the sub-solver, since BRGN does Newton steps (so only one anyways)
    Tao tao_sub;
    TaoBRGNGetSubsolver(tao_brgn, &tao_sub);
    TaoSetMaximumIterations(tao_sub, maxiter);
  } else {
    KSPGetTolerances(ksp_LeastSquares, &rtol, &abstol, &dtol, &current_maxits);
    KSPSetTolerances(ksp_LeastSquares, rtol, abstol, dtol, maxiter);
  }

  // 3. Set for EPS eigenvalue solver (used in solveGaussNewtonEPS)
  PetscReal eps_tol;
  EPSGetTolerances(eps_GN, &eps_tol, &current_maxits);
  EPSSetTolerances(eps_GN, eps_tol, maxiter);
}


int OptimProblem::getGaussNewtonMaxiter() {
  // Return the current max iterations from ksp_GN (they should all be the same)
  PetscReal rtol, abstol, dtol;
  PetscInt maxits;
  KSPGetTolerances(ksp_GN, &rtol, &abstol, &dtol, &maxits);
  return (int)maxits;
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

}


void OptimProblem::solveGaussNewtonKSP(Vec xinit, const Vec initial_guess, const Vec b, Vec Ainv_b){

  // Store the point of evaluation for the Gauss-Newton matrix shell A(xinit)
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Set the matrix again, just in case, for reset.
  KSPSetOperators(ksp_GN, GaussNewtonMatShell, GaussNewtonMatShell);

  // Set initial guess from initial_guess parameter
  VecCopy(initial_guess, Ainv_b);
  KSPSetInitialGuessNonzero(ksp_GN, PETSC_TRUE);

  // Monitor residual and solution norm at every iteration
  KSPMonitorCancel(ksp_GN);
  KSPMonitorSet(ksp_GN, KSPMonitorResidualAndSolution, (void*)this, NULL);

  // Optional Levenberg-Marquardt damping: A += mu I
  if (ksp_damping > 0.0) MatShift(GaussNewtonMatShell, ksp_damping);

  // Solve the linear system L^*L x = b
  GN_MatVec_counter = 0;
  KSPSolve(ksp_GN, b, Ainv_b);

  // Revert the optional scaling
  if (ksp_damping > 0.0) MatShift(GaussNewtonMatShell, -ksp_damping);

  // Report convergence
  KSPConvergedReason reason;
  int iters;
  double rnorm;
  KSPGetConvergedReason(ksp_GN, &reason);
  KSPGetIterationNumber(ksp_GN, &iters);
  KSPGetResidualNorm(ksp_GN, &rnorm);
  ksp_iters_last = iters;
  if (mpirank_world == 0 && !quietmode) {
    KSPType ksp_type;
    KSPGetType(ksp_GN, &ksp_type);
    printf("Gauss-Newton KSP (%s) stats: iterations = %d, MatVec counter = %d, residual norm = %1.14e\n", ksp_type, iters, GN_MatVec_counter, rnorm);
  }
}


void OptimProblem::solveGaussNewtonEPS(Vec xinit, const Vec b, Vec Ainv_b){

  // Get eigenvalues and eigenvectors of the Gauss-Newton matrix
  Mat evecs;
  std::vector<double> evals = computeGaussNewtonEvals(xinit, &evecs);

  // Print the eigenvalues
  // for (int i=0; i<evals.size(); i++) {
  //   if (mpirank_world == 0) printf("Eigenvalue %d: %1.14e\n", i, evals[i]);
  // }

  // Project the rhs onto evals: b_proj = evecs^Tb
  Vec tmp;
  MatCreateVecs(evecs, &tmp, NULL);
  MatMultTranspose(evecs, b, tmp);

  // Apply inverse eigenvalue scaling: tmp = evals^-1 * tmp
  PetscScalar* tmp_data;
  VecGetArray(tmp, &tmp_data);
  for (int i=0; i<evals.size(); i++) {
    if (evals[i] > eps_evals_cutoff) {
      tmp_data[i] /= (evals[i] + eps_damping);
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

PetscErrorCode KSPMonitorResidualAndSolution(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx){
  OptimProblem* self = (OptimProblem*) ctx;

  Vec x;
  KSPBuildSolution(ksp, NULL, &x); // Current iterate, valid at this point in the KSP solve
  PetscReal xnorm;
  VecNorm(x, NORM_2, &xnorm);

  if (self->getMPIrank_world() == 0 && !self->getQuietmode()) {
    printf("KSP it %d: residual norm = %1.14e, solution norm = %1.14e\n", (int)it, (double)rnorm, (double)xnorm);
  }

  return 0;
}

std::vector<double> OptimProblem::computeGaussNewtonEvals(Vec xinit, Mat* evecs_out){

  // Store xinit so the MatShell can use it as point of evaluation.
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Update EPS for a fresh solve on the new xeval_GN. 
  if (!GN_densemat) { 
    // use MatShell 
    EPSSetOperators(eps_GN, GaussNewtonMatShell, NULL);
  } else {
    // use dense Matrix representation of the Gauss-Newton matrix
    updateGaussNewtonMatDense(); // does N2-1 applications of the MatShell
    EPSSetOperators(eps_GN, GaussNewtonMatDense, NULL);
  }

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

    // Estimate the errror (needs one more application of A)
    double error = 0.0;
    EPSComputeError(eps_GN,i,EPS_ERROR_RELATIVE,&error);
    if (error > eps_tol) {
      if (mpirank_world==0) printf("ERROR: Relative error of eigenpair %d is large (error=%1.4e)\n", i, error);
    }
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
      Vec G, Gprec, xnew, step, v_zero;
      VecDuplicate(x_GN, &G);
      VecDuplicate(x_GN, &Gprec);
      VecDuplicate(x_GN, &xnew);
      VecDuplicate(x_GN, &step);
      VecDuplicate(x_GN, &v_zero);
      VecZeroEntries(v_zero);

      bool stop = false;
      for (int iter = 0; !stop; iter++) {
        // Compute the gradient (and objective) at the current iterate
        evalGradF(x_GN, G);
        double f = objective;
        VecNorm(G, NORM_2, &gnorm);

        // Precondition the gradient: solve the Gauss-Newton system A(x)*Gprec = -G via KSP
        // solveGaussNewtonKSP(x_GN, v_zero, G, Gprec);
        // solveGaussNewtonEPS(x_GN, G, Gprec);
        solveGaussNewtonLeastSquares(x_GN, v_zero, Gprec);
        VecScale(Gprec, -1.0); // For some reason, for the LeastSquares solver, the direction is -Gprec. TODO: Check. 

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
      VecDestroy(&v_zero);
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


void OptimProblem::updateGaussNewtonMatDense(){
  MatZeroEntries(GaussNewtonMatDense);

  Vec e;
  VecDuplicate(xeval_GN, &e);
  Vec Av;
  VecDuplicate(e, &Av);

  // Parallelize of comm_optim threads
  int ncols_per_rank = ndesign / mpisize_optim;
  int ncols_local = ncols_per_rank;
  if (mpirank_optim == mpisize_optim - 1) {
    ncols_local = ndesign - mpirank_optim * ncols_per_rank;
  }
  // printf("%d: Number of local columns = %d\n", mpirank_optim, ncols_local);

  // iterate over local columns
  for (int ix_local = 0; ix_local < ncols_local; ++ix_local) {
    int ix = ix_local + mpirank_optim * ncols_per_rank;
    // if (mpirank_init == 0 && !quietmode) printf("%d: Eval A*e_%d / %d \n", mpirank_optim, ix, ndesign);

    VecSet(e, 0.0);
    VecSetValue(e, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);
    // applyGaussNewtonMatShell(GaussNewtonMatDense, e, Av);
    MatMult(this->getGaussNewtonMatShell(), e, Av);
    for (int jx = 0; jx < ndesign; ++jx) {
      PetscScalar val;
      VecGetValues(Av, 1, &jx, &val);
      MatSetValue(GaussNewtonMatDense, jx, ix, val, INSERT_VALUES);
    }
  }
  VecDestroy(&e);
  VecDestroy(&Av);

  MatAssemblyBegin(GaussNewtonMatDense, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(GaussNewtonMatDense, MAT_FINAL_ASSEMBLY);

  // Allreduce to combine contributions from all MPI ranks
  PetscScalar *data;
  MatDenseGetArray(GaussNewtonMatDense, &data);
  int size = ndesign * ndesign;
  MPI_Allreduce(MPI_IN_PLACE, data, size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(GaussNewtonMatDense, &data);

}


void OptimProblem::solveGaussNewtonLeastSquares(const Vec xinit, const Vec initial_guess, Vec v_LeastSquares){

  // Store the point of evaluation (design parameters)
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Fill the RHS: b = -\nabla_U J projected to tangent space
  // \nabla_J = -2/n U - 2/n^2Tr(V^dU)V = 2/n(I-P)U

  Vec b; // RHS 
  GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, nullptr, &b);

  PetscInt tangent_dim;;
  VecGetSize(b, &tangent_dim);

  // Get U(T) from optim_target
  Mat U_final_re = optim_target->getFinalUnitaryRe();
  Mat U_final_im = optim_target->getFinalUnitaryIm();

  // first compute J(U(T)) needed for \nabla_UJ
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
    double obj_lin_iinit_re = 0.0;
    double obj_lin_iinit_im = 0.0;
    Vec ui = timestepper->getFinalState(iinit);
    optim_target->evalJ(ui,  &obj_lin_iinit_re, &obj_lin_iinit_im);
    obj_cost_re += obj_weights[iinit] * obj_lin_iinit_re;
    obj_cost_im += obj_weights[iinit] * obj_lin_iinit_im;
  }
  // Gather result J_inf(U(T))
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  PetscScalar *b_data;
  VecGetArray(b, &b_data);
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    Vec ui = timestepper->getFinalState(iinit);

    // Use pre-allocated workspace vector to compute -nabla_U J
    Vec b_iinit = wsub_workspace[iinit];

    // Copy ui to b_iinit
    VecCopy(ui, b_iinit);

    double obj_cost_re_bar, obj_cost_im_bar;
    VecScale(b_iinit, 2.0 / mastereq->getDim());
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
    optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
    double scale = 1.0 / mastereq->getDim();
    optim_target->evalJ_diff(ui, b_iinit, scale*obj_cost_re_bar, scale*obj_cost_im_bar);

    GellMann::projectVecToTangentSpace(b_iinit, U_final_re, U_final_im, mastereq->getDim(), iinit_global, b_data);
  }

  // Allreduce b across all comm_init processors
  MPI_Allreduce(MPI_IN_PLACE, b_data, tangent_dim, MPIU_SCALAR, MPI_SUM, comm_init);
  VecRestoreArray(b, &b_data);

  // Scale 
  VecScale(b, -1.0);

  // /*------ TEST setup of \nabla_U J and tangent space projection ----- */
  // if (mpisize_world > 1) {
  //   printf("Need mpisize_world == 1! \n");
  //   exit(1);
  // }

  // Mat V_re = optim_target->getTargetUnitaryRe();
  // Mat V_im = optim_target->getTargetUnitaryIm();

  // // Compute pre-factor -2/n^2 * <V,U> 
  // double factor_re, factor_im;
  // Mat VdU_re, VdU_im;
  // ComputeAdagB(V_re, V_im, U_final_re, U_final_im, &VdU_re, &VdU_im);
  // ComputeTrace(VdU_re, VdU_im, &factor_re, &factor_im);
  // factor_re = -2.0 / (mastereq->getDim() * mastereq->getDim()) * factor_re;
  // factor_im = -2.0 / (mastereq->getDim() * mastereq->getDim()) * factor_im;

  // // Compute nabla_U J = 2/n U + factor * V
  // // nablaUJ_re = 2/n U_re + factor_re*V_re - factor_im*V_im
  // // nablaUJ_im = 2/n U_im + factor_re*V_im + factor_im*V_re
  // Mat nablaUJ_re, nablaUJ_im;
  // MatDuplicate(U_final_re, MAT_COPY_VALUES, &nablaUJ_re);
  // MatDuplicate(U_final_im, MAT_COPY_VALUES, &nablaUJ_im);
  // MatScale(nablaUJ_re, 2.0 / mastereq->getDim());
  // MatScale(nablaUJ_im, 2.0 / mastereq->getDim());
  // MatAXPY(nablaUJ_re, factor_re, V_re, DIFFERENT_NONZERO_PATTERN);
  // MatAXPY(nablaUJ_re, -factor_im, V_im, DIFFERENT_NONZERO_PATTERN);
  // MatAXPY(nablaUJ_im, factor_re, V_im, DIFFERENT_NONZERO_PATTERN);
  // MatAXPY(nablaUJ_im, factor_im, V_re, DIFFERENT_NONZERO_PATTERN);

 
  // // Project onto tangent space
  // GellMann::projectMatToTangentSpace(nablaUJ_re, nablaUJ_im, U_final_re, U_final_im, mastereq->getDim(), b);

  // printf("MAT Projected RHS:\n");
  // VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // MatDestroy(&nablaUJ_re);
  // MatDestroy(&nablaUJ_im);

  // ----- END TANGENT SPACE PROJECTION SECTION -----


  // Note: includeHessUJ is disabled (set to false in constructor)
  if (includeHessUJ) {
    // TODO: Implement W^1/2 in tangent space if needed
    if (mpirank_world == 0) {
      printf("WARNING: includeHessUJ is not implemented for tangent space projection.\n");
    }
  }

  // // Adjoint test: <w, A v> should equal <A^T w, v>  (all real)
  // printf("GN Least-Squares: Adjoint test\n");
  // Vec v_test, w_test, Av, ATw;
  // GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, &v_test, nullptr);
  // GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, nullptr, &w_test);
  // VecSetRandom(v_test, nullptr);
  // VecSetRandom(w_test, nullptr);
  // GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, nullptr, &Av);
  // GNLeastSquaresShell_MatMult(GNLeastSquaresShell, v_test, Av);
  // GNLeastSquaresShell_MatCreateVecs(GNLeastSquaresShell, &ATw, nullptr);
  // GNLeastSquaresShell_MatMultTranspose(GNLeastSquaresShell, w_test, ATw);

  // PetscScalar lhs, rhs;
  // VecDot(w_test, Av, &lhs);
  // VecDot(ATw, v_test, &rhs);
  // PetscScalar max_val = PetscMax(PetscAbsScalar(lhs), PetscAbsScalar(rhs));
  // PetscScalar diff = lhs - rhs;
  // PetscScalar diff_rel = (lhs - rhs) / max_val;
  // PetscPrintf(PETSC_COMM_SELF, "adjoint check: %g vs %g, diff: %g, diff_rel: %g\n", (double)PetscRealPart(lhs), (double)PetscRealPart(rhs), (double)PetscRealPart(diff), (double)PetscRealPart(diff_rel));
  // exit(1);


  // Solve the Least Squares system with either TAOBRGN or KSPLSQR
  if (ls_solver == "BRGN") {

    // Cache the RHS for residual evaluations (store -b for F(v) = L*v - b)
    VecCopy(b, brgn_rhs);

    // Set initial guess from initial_guess parameter
    VecCopy(initial_guess, v_LeastSquares);
    TaoSetSolution(tao_brgn, v_LeastSquares);
    TaoSetUp(tao_brgn);
    // Forward initial guess to the inner solver  
    Tao tao_sub;
    TaoBRGNGetSubsolver(tao_brgn, &tao_sub);
    TaoSetSolution(tao_sub, v_LeastSquares);

    // Solve the regularized least-squares problem: min_v ||L*v - b||² + λ||v||²
    TaoSolve(tao_brgn);

    // Report convergence. Tao iterations are outer Newton steps (one suffices for a linear problem), the actual cost is in the inner linear solver iterations.
    TaoConvergedReason reason;
    PetscInt newton_iters, lin_iters;
    TaoGetConvergedReason(tao_brgn, &reason);
    TaoGetIterationNumber(tao_brgn, &newton_iters);
    TaoGetLinearSolveIterations(tao_sub, &lin_iters);
    ksp_iters_last = lin_iters;

    // TaoGetResidualNorm is not updated by TAOBRGN, so compute ||L*v - b|| directly
    PetscReal rnorm;
    MatMult(GNLeastSquaresShell, v_LeastSquares, brgn_residual);
    VecAXPY(brgn_residual, -1.0, brgn_rhs);
    VecNorm(brgn_residual, NORM_2, &rnorm);

    if (mpirank_world == 0 && !quietmode) {
      printf("Tao BRGN Least-Squares stats: newton iterations = %d, linear iterations = %d, residual norm = %1.14e, reason = %d\n", (int)newton_iters, (int)lin_iters, rnorm, reason);
    }

  } else if (ls_solver == "LSQR"){
    // Set the matrix again, just in case, for reset.
    KSPSetOperators(ksp_LeastSquares, GNLeastSquaresShell, GNLeastSquaresShell);

    // Set initial guess from initial_guess parameter
    VecCopy(initial_guess, v_LeastSquares);
    KSPSetInitialGuessNonzero(ksp_LeastSquares, PETSC_TRUE);

    // Monitor residual and solution norm at every iteration
    KSPMonitorCancel(ksp_LeastSquares);
    KSPMonitorSet(ksp_LeastSquares, KSPMonitorResidualAndSolution, (void*)this, NULL);

    // Solve the Gauss-Newton least-squares problem
    KSPSolve(ksp_LeastSquares, b, v_LeastSquares);

    // Report convergence
    KSPConvergedReason reason;
    int iters;
    double rnorm;
    KSPGetConvergedReason(ksp_LeastSquares, &reason);
    KSPGetIterationNumber(ksp_LeastSquares, &iters);
    KSPGetResidualNorm(ksp_LeastSquares, &rnorm);
    ksp_iters_last = iters;
    if (mpirank_world == 0 && !quietmode) {
      printf("Gauss-Newton Least-Squares stats: iterations = %d, residual norm = %1.14e\n", iters, rnorm);
    }
  }

  // Clean up
  VecDestroy(&b);

}

// Callback: Evaluate residual F(v) = L*v - b
PetscErrorCode OptimProblem::TaoBRGN_EvalResidual(Tao tao, Vec v, Vec F, void *ctx) {
  PetscFunctionBeginUser;

  OptimProblem* self = (OptimProblem*)ctx;

  // Compute residual: F(v) = L*v - b
  MatMult(self->GNLeastSquaresShell, v, F);  // F = L*v (use MatMult, not direct call)
  VecAXPY(F, -1.0, self->brgn_rhs);          // F = L*v - b

  PetscFunctionReturn(0);
}

// Callback: Evaluate Jacobian dF/dv = L (constant matrix)
PetscErrorCode OptimProblem::TaoBRGN_EvalJacobianResidual(Tao tao, Vec v, Mat J, Mat Jpre, void *ctx) {
  PetscFunctionBeginUser;

  // The Jacobian is L = GNLeastSquaresShell, which is constant
  // Already set via TaoSetJacobianResidualRoutine, nothing to compute

  PetscFunctionReturn(0);
}

// Monitor function for TaoBRGN iterations
PetscErrorCode OptimProblem::TaoBRGN_Monitor(Tao tao, void *ctx) {
  PetscFunctionBeginUser;

  OptimProblem* self = (OptimProblem*)ctx;

  // Get current status (iteration, objective, gradient norm, etc.)
  PetscInt iter;
  PetscReal f, gnorm, cnorm, xdiff;
  TaoConvergedReason reason;
  TaoGetSolutionStatus(tao, &iter, &f, &gnorm, &cnorm, &xdiff, &reason);

  // Get current solution norm
  Vec x;
  TaoGetSolution(tao, &x);
  PetscReal xnorm;
  VecNorm(x, NORM_2, &xnorm);

  // Get current residual norm
  PetscReal rnorm;
  TaoGetResidualNorm(tao, &rnorm);

  // Print progress (objective is (1/2)||F||^2 + (lambda/2)||x||^2)
  if (self->mpirank_world == 0 && !self->quietmode) {
    printf("  TaoBRGN it %d: objective = %1.8e, residual norm = %1.8e, solution norm = %1.8e\n",
           (int)iter, f, rnorm, xnorm);
  }

  PetscFunctionReturn(0);
}


void OptimProblem::GNLeastSquaresShell_MatMult(Mat A, Vec v, Vec y)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(y);

  PetscInt tangent_dim;
  VecGetSize(y, &tangent_dim);

  // Apply linearized forward to get dU/dalpha x v
  // assumes that the timestepper's trajectory_states are already populated and computes all lin_trajectory_states.
  self->evalLinearizedForward(self->xeval_GN, v);

  // Get U_final matrices. This needs that the matrices are synchronized across comm_init (happens in finalizeJ)
  Mat U_final_re = self->optim_target->getFinalUnitaryRe();
  Mat U_final_im = self->optim_target->getFinalUnitaryIm();

  // Set up Omega_j = U(T)^d *(Lv)_j and project onto tangent space
  PetscScalar *y_data;
  VecGetArray(y, &y_data);
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);

    GellMann::projectVecToTangentSpace(lin_final_state, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, y_data);
  }

  // Allreduce the data in y across all comm_init processors
  MPI_Allreduce(MPI_IN_PLACE, y_data, tangent_dim, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(y, &y_data);

  //   /* ----- TEST: Project Lv onto tangent space  ----*/
  // assert(self->mpisize_world == 1);
  // int dim = self->mastereq->getDim();

  // // First set up W = Lv as a Petsc matrix
  // Mat W_re, W_im;
  // MatDuplicate(U_final_re, MAT_DO_NOT_COPY_VALUES, &W_re);
  // MatDuplicate(U_final_im, MAT_DO_NOT_COPY_VALUES, &W_im);
  // MatZeroEntries(W_re);
  // MatZeroEntries(W_im);
  // for (int iinit = 0; iinit < self->ninit_local; iinit++) {
  //   int iinit_global = self->mpirank_init * self->ninit_local + iinit;
  //   Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);
  //   const PetscScalar *lin_final_state_ptr;
  //   VecGetArrayRead(lin_final_state, &lin_final_state_ptr);
  //   for (size_t row = 0; row < dim; row++) {
  //     int row_re = row;
  //     int row_im = row + dim;
  //     int col = iinit_global;
  //     MatSetValue(W_re, row, col, lin_final_state_ptr[row_re], INSERT_VALUES);
  //     MatSetValue(W_im, row, col, lin_final_state_ptr[row_im], INSERT_VALUES);
  //   }
  //   VecRestoreArrayRead(lin_final_state, &lin_final_state_ptr);
  // }
  // MatAssemblyBegin(W_re, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(W_re, MAT_FINAL_ASSEMBLY);
  // MatAssemblyBegin(W_im, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(W_im, MAT_FINAL_ASSEMBLY);

  // /* --- Set the tangent space projection from Omega ----*/
  // GellMann::projectMatToTangentSpace(W_re, W_im, U_final_re, U_final_im, dim, y);

  // printf("MAT Projected Lv:\n");
  // VecView(y, PETSC_VIEWER_STDOUT_WORLD);

  // MatDestroy(&W_re);
  // MatDestroy(&W_im);

  // exit(1);
  // /* ----- END TEST: Project Lv onto tangent space  ----*/

  // Apply W^1/2 to local buffer if needed
  // NOTE: includeHessUJ is disabled (set to false in constructor), so this branch is unlikely
  if (self->includeHessUJ){
    // Only available for J_Inf (Generalized)
    assert(self->optim_target->getObjectiveType() == ObjectiveType::JTRACE);

    // TODO: Need to implement W^1/2 in tangent space representation
    // For now, this is disabled as includeHessUJ = false in the constructor
    if (self->mpirank_world == 0) {
      printf("WARNING: includeHessUJ is not implemented for tangent space projection. Skipping W^1/2 application.\n");
    }
  }
}

void OptimProblem::GNLeastSquaresShell_MatMultTranspose(Mat A, Vec w, Vec vout)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(vout);

  // Get read access to w (in reduced tangent space) Note: w has size tangent_dim, 
  const PetscScalar *w_data;
  VecGetArrayRead(w, &w_data);

  // Optional, apply W^1/2 to w before reconstruction
  // NOTE: includeHessUJ is disabled (set to false in constructor)
  if (self->includeHessUJ){
    // Only available for J_Inf (Generalized)
    assert(self->optim_target->getObjectiveType() == ObjectiveType::JTRACE);

    // TODO: Need to implement W^1/2 in tangent space representation
    // For now, this is disabled as includeHessUJ = false in the constructor
    if (self->mpirank_world == 0) {
      printf("WARNING: includeHessUJ is not implemented for tangent space projection. Skipping W^1/2 application.\n");
    }
  }

  // Reconstruct full state vectors from tangent space coefficients for this processor's initial conditions
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    // Get U(T) from optim_target
    Mat U_final_re = self->optim_target->getFinalUnitaryRe();
    Mat U_final_im = self->optim_target->getFinalUnitaryIm();

    // Reconstruct full state vector from tangent space coefficients. This is the adjoint of projectVecToTangentSpace
    GellMann::reconstructVecFromTangentSpace(w_data, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, self->wsub_workspace[iinit]);
  }

  VecRestoreArrayRead(w, &w_data);

  // Solve adjoint backward ODE for each terminal condition in workspace
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Solve backwards
    self->timestepper->solveAdjointODE(iinit, self->wsub_workspace[iinit], 0.0, 0.0, 0.0, 0.0);

    // Add gradient to output
    VecAXPY(vout, 1.0, self->timestepper->getReducedGradient());
  }

  /* Sum up the gradient from all initial condition processors */
  PetscScalar* vout_data;
  VecGetArray(vout, &vout_data);
  MPI_Allreduce(MPI_IN_PLACE, vout_data, self->ndesign, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(vout, &vout_data);
}

void OptimProblem::GNLeastSquaresShell_MatCreateVecs(Mat A, Vec *right, Vec *left)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  if (right) VecCreateSeq(PETSC_COMM_SELF, self->ndesign, right);
  if (left) {
    // Create a reduced tangent space
    PetscInt N = self->mastereq->getDim();
    PetscInt tangent_dim = N * N - 1;
    VecCreateSeq(PETSC_COMM_SELF, tangent_dim, left);
  }
}
