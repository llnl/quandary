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
  TaoSetType(tao, TAOBQNLS);
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

  // Include hessian of generalized J_inf wrt U in GaussNewton Approximation (for Jtrace only)
  includeHessUJ = false;
  if (optim_solver_type == OptimSolverType::GAUSS_NEWTON ||
      optim_solver_type == OptimSolverType::LEAST_SQUARES){
    if (optim_target->getObjectiveType() == ObjectiveType::JTRACE) {
      includeHessUJ = true;
    }
  }
  includeHessUJ = true;
  // includeHessUJ = false;

  /* Create MatShell for Gauss-Newton least-squares problem */
  leastsquares_damping = config.getLeastSquaresDamping();
  int N = ndesign;
  PetscInt M = mastereq->getDim()*mastereq->getDim() - 1;
  // enlarge the output space if damping is enabled (only for least-squares mode)
  if (optim_solver_type == OptimSolverType::LEAST_SQUARES && leastsquares_damping > 0.0) {
    M += ndesign;
    if (mpirank_world == 0 && !quietmode) {
      printf("Least-squares damping enabled: λ = %e, augmented matrix size = %d × %d\n",
             leastsquares_damping, (int)M, N);
    }
  }
  MatCreateShell(PETSC_COMM_SELF, M, N, M, N, this, &GN_LeastSquares_MatShell);
  MatShellSetOperation(GN_LeastSquares_MatShell, MATOP_MULT, (void(*)(void))GN_LeastSquares_MatMult);
  MatShellSetOperation(GN_LeastSquares_MatShell, MATOP_MULT_TRANSPOSE, (void(*)(void))GN_LeastSquares_MatMultTranspose);
  MatShellSetOperation(GN_LeastSquares_MatShell, MATOP_CREATE_VECS,     (void(*)(void))GN_LeastSquares_MatCreateVecs);

  // Create the KSP solver for the Gauss-Newton least-squares problem
  KSPCreate(PETSC_COMM_SELF, &ksp_LeastSquares);
  KSPSetOperators(ksp_LeastSquares, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell);
  // KSPSetNormType(ksp_LeastSquares, KSP_NORM_UNPRECONDITIONED); // Unconditioned ressidual norm
  KSPSetType(ksp_LeastSquares, KSPLSQR);
  KSPSetTolerances(ksp_LeastSquares,config.getGnRtol(),PETSC_DEFAULT,PETSC_DEFAULT,config.getGnMaxiter());
  // KSPSetComputeSingularValues(ksp_LeastSquares, PETSC_TRUE);
  KSPSetFromOptions(ksp_LeastSquares);

  // Initialize Tao BRGN solver for least-squares problem if configured
  leastsquares_solver = config.getLeastSquaresSolver();
  if (leastsquares_solver == "brgn") {
    // Create Tao solver on PETSC_COMM_SELF (same as ksp_LeastSquares)
    TaoCreate(PETSC_COMM_SELF, &tao_brgn);
    TaoSetType(tao_brgn, TAOBRGN);

    // Create storage for RHS vector (used in residual callback)
    GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &brgn_rhs);
    VecDuplicate(brgn_rhs, &brgn_residual);

    // Set BRGN convergence tolerances. TODO: THIS NEEDS FIX! BRGN is a newton solver, so it will only need one outer newton step. The inner linear solver is the actual workhorse. Need to pass appropriate tolerances to the subsolver, not the BRGN solver. However, somehow, the below crashes when running on multiple processes. TODO: Check parallel implementation.
    // Tao tao_sub;
    // TaoBRGNGetSubsolver(tao_brgn, &tao_sub);
    // TaoSetTolerances(tao_sub, config.getGnRtol(), PETSC_DEFAULT, PETSC_DEFAULT);
    // TaoSetMaximumIterations(tao_sub, config.getGnMaxiter());
    TaoSetTolerances(tao_brgn, config.getGnRtol(), PETSC_DEFAULT, PETSC_DEFAULT);
    TaoSetMaximumIterations(tao_brgn, config.getGnMaxiter());

    // Set damping parameter for regularization term λ||v||²
    TaoBRGNSetRegularizerWeight(tao_brgn, leastsquares_damping);

    // Set Residual routine F(v) = Lv - b 
    TaoSetResidualRoutine(tao_brgn, brgn_residual, TaoBRGN_EvalResidual, (void*)this);

    // Set Jacobian (the matrix L = GN_LeastSquares_MatShell, constant for this problem)
    TaoSetJacobianResidualRoutine(tao_brgn, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell, TaoBRGN_EvalJacobianResidual, (void*)this);

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

  // Pre-allocate workspace vectors for GN_LeastSquares_MatShell operations
  wsub_workspace = new Vec[ninit_local];
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    VecDuplicate(rho_t0_bar, &wsub_workspace[iinit]);
  }

  /* Create MatShell for Gauss-Newton normal equation A=L^*L */
  MatCreateShell(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ndesign, ndesign, this, &GN_NormalEq_MatShell);
  MatShellSetOperation(GN_NormalEq_MatShell, MATOP_MULT, (void(*) (void)) GN_NormalEq_MatMult);
  VecDuplicate(xinit, &xeval_GN);
  VecZeroEntries(xeval_GN);
  VecAssemblyBegin(xeval_GN); VecAssemblyEnd(xeval_GN);

  /* Create dense matrix for Gauss-Newton normal equation */
  MatCreateDense(PETSC_COMM_SELF, ndesign, ndesign, ndesign, ndesign, NULL, &GN_NormalEq_MatDense);
  MatZeroEntries(GN_NormalEq_MatDense);

  gn_normaleq_damping = config.getGnNormaleqDamping();
  gn_rtol = config.getGnRtol();
  gn_solver_mode = config.getGnSolverMode();

  /* Create GN Normal equation linear solver: L^*Lv=b */
  KSPCreate(PETSC_COMM_SELF, &ksp_GN_NormalEq);
  KSPSetOperators(ksp_GN_NormalEq, GN_NormalEq_MatShell, GN_NormalEq_MatShell);
  std::string ksp_normaleq_solver_str = config.getGnNormaleqSolver();
  gn_normaleq_solver = ksp_normaleq_solver_str;  // Store as member variable
  std::string ksp_normaleq_solver_lower = ksp_normaleq_solver_str;
  std::transform(ksp_normaleq_solver_lower.begin(), ksp_normaleq_solver_lower.end(), ksp_normaleq_solver_lower.begin(), ::tolower);

  // Skip KSP configuration for projected solver (only applies to dual mode)
  if (ksp_normaleq_solver_lower != "proj") {
    if (ksp_normaleq_solver_lower == "cg") {
      KSPSetType(ksp_GN_NormalEq, KSPCG);
    } else if (ksp_normaleq_solver_lower == "minres") {
      KSPSetType(ksp_GN_NormalEq, KSPMINRES);
    } else if (ksp_normaleq_solver_lower == "gmres") {
      KSPSetType(ksp_GN_NormalEq, KSPGMRES);
    } else {
      // Default or try to use the string directly for other types
      if (mpirank_world == 0 && !quietmode) {
        printf("Warning: Unknown gn_normaleq_solver '%s', using it directly with PETSc\n", ksp_normaleq_solver_str.c_str());
      }
      KSPSetType(ksp_GN_NormalEq, ksp_normaleq_solver_str.c_str());
    }
    // Allow for non-zero initial guess
    KSPSetInitialGuessNonzero(ksp_GN_NormalEq, PETSC_TRUE);

    // Enable MINRES-QLP if requested
    if (ksp_normaleq_solver_lower == "minres" && config.getGnNormaleqMinresQlp()) {
      PetscBool set;
      PetscOptionsHasName(NULL, NULL, "-ksp_minres_qlp", &set);
      if (!set) {
        PetscOptionsSetValue(NULL, "-ksp_minres_qlp", NULL);
      }
    }
  } else {
    // For projected solver, set a default type (won't be used, but needs to be valid)
    KSPSetType(ksp_GN_NormalEq, KSPMINRES);
  }
  KSPSetFromOptions(ksp_GN_NormalEq);
  PC  pc;
  KSPGetPC(ksp_GN_NormalEq, &pc);
  PCSetType(pc, PCNONE); // Disable preconditioner
  KSPSetNormType(ksp_GN_NormalEq, KSP_NORM_UNPRECONDITIONED); // Unconditioned ressidual norm
  KSPSetTolerances(ksp_GN_NormalEq,config.getGnRtol(),PETSC_DEFAULT,PETSC_DEFAULT,config.getGnMaxiter());

  /* Dual mode uses pure tangent space (not augmented) */
  PetscInt M_dual = mastereq->getDim() * mastereq->getDim() - 1;

  /* Create MatShell for dual Gauss-Newton problem (L L^T + λI) */
  MatCreateShell(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, this, &GN_Dual_MatShell);
  MatShellSetOperation(GN_Dual_MatShell, MATOP_MULT, (void(*)(void))GN_Dual_MatMult);
  MatShellSetOperation(GN_Dual_MatShell, MATOP_GET_DIAGONAL, (void(*)(void))GN_Dual_GetDiagonal);
  MatShellSetOperation(GN_Dual_MatShell, MATOP_CREATE_VECS, (void(*)(void))GN_Dual_MatCreateVecs);

  /* Create dense matrices for dual Gauss-Newton problem */
  MatCreateDense(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, NULL, &GN_Dual_MatDense);
  MatZeroEntries(GN_Dual_MatDense);
  MatCreateDense(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, NULL, &GN_Dual_Precond_MatDense);
  MatZeroEntries(GN_Dual_Precond_MatDense);

  /* Create workspace vectors for dual solver */
  VecCreateSeq(PETSC_COMM_SELF, M_dual, &y_dual_workspace);  // tangent space
  VecDuplicate(xinit, &v_dual_workspace);  // design space
  VecCreateSeq(PETSC_COMM_SELF, M_dual, &gn_dual_precond_diag);  // D^{-1/2} for symmetric preconditioning
  VecCreateSeq(PETSC_COMM_SELF, M_dual, &gn_dual_temp);  // temp vector for preconditioning ops

  /* Create MatShell for symmetrically preconditioned operator D^{-1/2} A D^{-1/2} */
  MatCreateShell(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, this, &GN_Dual_Precond_MatShell);
  MatShellSetOperation(GN_Dual_Precond_MatShell, MATOP_MULT, (void(*)(void))GN_Dual_Precond_MatMult);

  /* Store preconditioner settings */
  gn_prec_update_interval = config.getGnPrecUpdateInterval();
  gn_dual_precond_type = config.getGnDualPrecondType();
  gn_dual_precond_rank = config.getGnDualPrecondRank();
  gn_maxiter = config.getGnMaxiter();
  gn_dual_ksp_baseline = 0;  // Will be set after first preconditioner formation
  gn_dual_ksp_excess_accum = 0;  // Accumulated excess KSP iterations

  /* Create low-rank preconditioner data structures */
  MatCreateDense(PETSC_COMM_SELF, M_dual, gn_dual_precond_rank, M_dual, gn_dual_precond_rank, NULL, &gn_dual_lowrank_U);
  VecCreateSeq(PETSC_COMM_SELF, gn_dual_precond_rank, &gn_dual_lowrank_sigma_invsqrt);
  VecCreateSeq(PETSC_COMM_SELF, gn_dual_precond_rank, &gn_dual_lowrank_temp);

  /* Create MatShell for low-rank preconditioned operator M^{-1/2} A M^{-1/2} */
  MatCreateShell(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, this, &GN_Dual_LowRank_MatShell);
  MatShellSetOperation(GN_Dual_LowRank_MatShell, MATOP_MULT, (void(*)(void))GN_Dual_LowRank_MatMult);

  /* Create KSP solver for dual Gauss-Newton problem */
  KSPCreate(PETSC_COMM_SELF, &ksp_GN_Dual);
  // Operators will be set at solve time depending on whether preconditioning is enabled
  KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);

  // Use same solver type configuration as normal equation solver (skip if using projected solver)
  if (ksp_normaleq_solver_lower != "proj") {
    if (ksp_normaleq_solver_lower == "cg") {
      KSPSetType(ksp_GN_Dual, KSPCG);
    } else if (ksp_normaleq_solver_lower == "minres") {
      KSPSetType(ksp_GN_Dual, KSPMINRES);
    } else if (ksp_normaleq_solver_lower == "gmres") {
      KSPSetType(ksp_GN_Dual, KSPGMRES);
    } else {
      if (mpirank_world == 0 && !quietmode) {
        printf("Warning: Unknown gn_normaleq_solver '%s' for dual solver, using it directly with PETSc\n", ksp_normaleq_solver_str.c_str());
      }
      KSPSetType(ksp_GN_Dual, ksp_normaleq_solver_str.c_str());
    }
    KSPSetInitialGuessNonzero(ksp_GN_Dual, PETSC_TRUE);

    // Apply MINRES-QLP if configured
    if (ksp_normaleq_solver_lower == "minres" && config.getGnNormaleqMinresQlp()) {
      PetscBool set;
      PetscOptionsHasName(NULL, NULL, "-ksp_minres_qlp", &set);
      if (!set) {
        PetscOptionsSetValue(NULL, "-ksp_minres_qlp", NULL);
      }
    }
  } else {
    // For projected solver, set a default type (won't be used, but needs to be valid)
    KSPSetType(ksp_GN_Dual, KSPMINRES);
  }

  KSPSetFromOptions(ksp_GN_Dual);
  PC pc_dual;
  KSPGetPC(ksp_GN_Dual, &pc_dual);
  PCSetType(pc_dual, PCNONE); // Disable preconditioner
  KSPSetNormType(ksp_GN_Dual, KSP_NORM_UNPRECONDITIONED);
  KSPSetTolerances(ksp_GN_Dual, config.getGnRtol(), PETSC_DEFAULT, PETSC_DEFAULT, config.getGnMaxiter());

  /* Create eigenvalues solver for Gauss-Newton normal equation Ax=b */
  EPSCreate(PETSC_COMM_SELF, &eps_GN_NormalEq);
  EPSSetOperators(eps_GN_NormalEq, GN_NormalEq_MatShell, NULL);
  EPSSetProblemType(eps_GN_NormalEq, EPS_HEP); // Hermitian
  EPSSetWhichEigenpairs(eps_GN_NormalEq, EPS_LARGEST_REAL); // largest eigenvalues
  double eps_thresh = eps_evals_cutoff;
  EPSSetThreshold(eps_GN_NormalEq, eps_thresh, PETSC_FALSE);  // absolute threshold
  neigvals = mastereq->getDim()*mastereq->getDim() - 1;
  // ncv = neigvals + 2; // Max Krylov dimension. How to set??
  ncv = 2*neigvals ; // Max Krylov dimension. How to set??
  EPSSetDimensions(eps_GN_NormalEq, neigvals, ncv, PETSC_DEFAULT);
  EPSSetTolerances(eps_GN_NormalEq, eps_tol, config.getGnMaxiter());
  EPSSetFromOptions(eps_GN_NormalEq);

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

  MatDestroy(&GN_NormalEq_MatDense);
  MatDestroy(&GN_NormalEq_MatShell);
  MatDestroy(&GN_NormalEq_MatDense);
  MatDestroy(&GN_LeastSquares_MatShell);
  MatDestroy(&GN_Dual_MatShell);
  MatDestroy(&GN_Dual_Precond_MatShell);
  MatDestroy(&GN_Dual_LowRank_MatShell);
  MatDestroy(&GN_Dual_MatDense);
  MatDestroy(&GN_Dual_Precond_MatDense);
  VecDestroy(&xeval_GN);
  VecDestroy(&y_dual_workspace);
  VecDestroy(&v_dual_workspace);
  VecDestroy(&gn_dual_precond_diag);
  VecDestroy(&gn_dual_temp);
  MatDestroy(&gn_dual_lowrank_U);
  VecDestroy(&gn_dual_lowrank_sigma_invsqrt);
  VecDestroy(&gn_dual_lowrank_temp);
  KSPDestroy(&ksp_GN_NormalEq);
  KSPDestroy(&ksp_GN_Dual);
  EPSDestroy(&eps_GN_NormalEq);
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


void OptimProblem::setGaussNewtonEvalPoint(Vec xinit) {
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false;
}

void OptimProblem::setGaussNewtonMaxiter(int maxiter) {
  if (mpirank_world == 0 && !quietmode) {
    printf("Setting max iterations for all Gauss-Newton solvers to: %d\n", maxiter);
  }

  // 1. Set for for solveGaussNewtonNormalEqKSP KSP solver
  PetscReal rtol, abstol, dtol;
  PetscInt current_maxits;
  KSPGetTolerances(ksp_GN_NormalEq, &rtol, &abstol, &dtol, &current_maxits);
  KSPSetTolerances(ksp_GN_NormalEq, rtol, abstol, dtol, maxiter);

  // 2. Set for solveGaussNewtonLeastSquares solvers
  if (leastsquares_solver == "brgn") {
    // For BRGN, this should apply to the sub-solver, since BRGN does Newton steps (so only one anyways). This crashes though. TODO: Fix. 
    // Tao tao_sub;
    // TaoBRGNGetSubsolver(tao_brgn, &tao_sub);
    // TaoSetMaximumIterations(tao_sub, maxiter);
    TaoSetMaximumIterations(tao_brgn, maxiter); // no effect
  } else {
    KSPGetTolerances(ksp_LeastSquares, &rtol, &abstol, &dtol, &current_maxits);
    KSPSetTolerances(ksp_LeastSquares, rtol, abstol, dtol, maxiter);
  }

  // 3. Set for dual Gauss-Newton KSP solver
  KSPGetTolerances(ksp_GN_Dual, &rtol, &abstol, &dtol, &current_maxits);
  KSPSetTolerances(ksp_GN_Dual, rtol, abstol, dtol, maxiter);

  // 4. Set for EPS eigenvalue solver (used in solveGaussNewtonNormalEqEPS)
  PetscReal eps_tol;
  EPSGetTolerances(eps_GN_NormalEq, &eps_tol, &current_maxits);
  EPSSetTolerances(eps_GN_NormalEq, eps_tol, maxiter);
}


int OptimProblem::getGaussNewtonMaxiter() {
  // Return the current max iterations from ksp_GN_NormalEq (they should all be the same)
  PetscReal rtol, abstol, dtol;
  PetscInt maxits;
  KSPGetTolerances(ksp_GN_NormalEq, &rtol, &abstol, &dtol, &maxits);
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

PetscErrorCode OptimProblem::GN_NormalEq_MatMult(Mat A, const Vec v, Vec Av){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);
  // if (self->mpirank_world == 0) printf("APPLYING GAUSS-NEWTON...\n");

  // Grab the point of evaluation from the shell
  Vec x = self->xeval_GN;

  //  Reset output 
  VecZeroEntries(Av);
  
  // Apply linearized forward to get dU/dalpha x v
  // assumes that the timestepper's trajectory_states are already populated and computes all lin_trajectory_states.
  self->evalLinearizedForward(x, v);

  // If include W^1/2, precompute z = 1/n<V,LinStates> 
  double z_re = 0.0;
  double z_im = 0.0;
  if (self->includeHessUJ){
    for (int iinit=0; iinit<self->ninit_local; iinit++) {
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      double tmp_re = 0.0;
      double tmp_im = 0.0;
      Vec wT = self->timestepper->getLinearizedFinalState(iinit);
      self->optim_target->HilbertSchmidtOverlap(wT, false, &tmp_re, &tmp_im);
      z_re +=  tmp_re / self->mastereq->getDim();
      z_im +=  tmp_im / self->mastereq->getDim();
    }
    // Gather result J_inf(W(T))
    double myz_re = z_re;
    double myz_im = z_im;
    MPI_Allreduce(&myz_re, &z_re, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
    MPI_Allreduce(&myz_im, &z_im, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
  }

  // Solve the adjoint ODE for initial condition 
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Set terminal condition for adjoint
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);

    // Set the terminal condition 
    VecCopy(lin_final_state, self->rho_t0_bar);  // w_i(T)

    // Include W = 2/n(I-<V,cdot>V), if requested
    if (self->includeHessUJ){
      double scale = 2.0 / self->mastereq->getDim();
      VecScale(self->rho_t0_bar, scale);
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      self->optim_target->HilbertSchmidtOverlap_diff(self->rho_t0_bar, false, -scale*z_re, -scale*z_im);
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

  return 0;
}

PetscErrorCode OptimProblem::GN_Dual_MatMult(Mat A, const Vec y, Vec Ay){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(Ay);

  // Compute Ay = L L^T y (using pure Jacobian operations for dual mode)
  // Note: We need to pass a matrix for context, but the functions use VecGetSize to get dimensions
  // Step 1: Compute v_temp = L^T y (tangent space → design space)
  GN_Jacobian_MatMultTranspose(A, y, self->v_dual_workspace);

  // Step 2: Compute Ay = L v_temp (design space → tangent space)
  GN_Jacobian_MatMult(A, self->v_dual_workspace, Ay);

  // Note: Damping (λI) is applied separately via MatShift on the shell

  return 0;
}

PetscErrorCode OptimProblem::GN_Dual_GetDiagonal(Mat A, Vec diag){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Return the stored diagonal from gn_dual_precond_diag
  VecCopy(self->gn_dual_precond_diag, diag);

  return 0;
}

PetscErrorCode OptimProblem::GN_Dual_Precond_MatMult(Mat A, const Vec z, Vec output){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Apply D^{-1/2} A D^{-1/2} z
  // Step 1: temp = D^{-1/2} * z
  VecPointwiseMult(self->gn_dual_temp, self->gn_dual_precond_diag, z);

  // Step 2: Apply A (the unpreconditioned operator GN_Dual_MatShell)
  // We need to apply (L L^T + λI), but MatShell already has the shift applied
  MatMult(self->GN_Dual_MatShell, self->gn_dual_temp, output);

  // Step 3: output = D^{-1/2} * output
  VecPointwiseMult(output, self->gn_dual_precond_diag, output);

  return 0;
}


void OptimProblem::solveGaussNewtonNormalEqKSP(Vec xinit, const Vec initial_guess, const Vec b, Vec Ainv_b){

  // Store the point of evaluation for the Gauss-Newton matrix shell A(xinit)
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Set the matrix again, just in case, for reset.
  KSPSetOperators(ksp_GN_NormalEq, GN_NormalEq_MatShell, GN_NormalEq_MatShell);

  // Set initial guess from initial_guess parameter
  VecCopy(initial_guess, Ainv_b);

  // Monitor residual and solution norm at every iteration
  KSPMonitorCancel(ksp_GN_NormalEq);
  KSPMonitorSet(ksp_GN_NormalEq, KSPMonitorResidualAndSolution, (void*)this, NULL);

  // Optional Levenberg-Marquardt damping: A += mu I
  if (gn_normaleq_damping > 0.0) MatShift(GN_NormalEq_MatShell, gn_normaleq_damping);

  // Solve the linear system L^*L x = b
  KSPSolve(ksp_GN_NormalEq, b, Ainv_b);

  // Revert the optional scaling
  if (gn_normaleq_damping > 0.0) MatShift(GN_NormalEq_MatShell, -gn_normaleq_damping);

  // Report convergence
  KSPConvergedReason reason;
  int iters;
  double rnorm;
  KSPGetConvergedReason(ksp_GN_NormalEq, &reason);
  KSPGetIterationNumber(ksp_GN_NormalEq, &iters);
  KSPGetResidualNorm(ksp_GN_NormalEq, &rnorm);
  ksp_iters_last = iters;
  if (mpirank_world == 0 && !quietmode) {
    KSPType ksp_type;
    KSPGetType(ksp_GN_NormalEq, &ksp_type);
    printf("Gauss-Newton KSP (%s) stats: iterations = %d, residual norm = %1.14e\n", ksp_type, iters, rnorm);
  }
}

void OptimProblem::solveGaussNewtonDual(Vec xinit, const Vec initial_guess, Vec v_out){

  // Check if using projected subspace solver instead of iterative KSP
  std::string solver = gn_normaleq_solver;
  std::string solver_lower;
  for (char c : solver) solver_lower += std::tolower(c);

  if (solver_lower == "proj") {
    // Use projected subspace solver (no KSP iterations)
    solveGaussNewtonDualProjected(xinit, v_out);
    return;
  }

  // Otherwise proceed with standard KSP-based solver

  // Store the point of evaluation for the Gauss-Newton matrix shell
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Compute RHS: b = -∇_U J projected to tangent space
  Vec b;
  GN_Dual_MatCreateVecs(GN_Dual_MatShell, nullptr, &b);

  PetscInt tangent_dim;
  VecGetSize(b, &tangent_dim);

  // Get U(T) from optim_target
  Mat U_final_re = optim_target->getFinalUnitaryRe();
  Mat U_final_im = optim_target->getFinalUnitaryIm();

  // Compute J(U(T)) needed for ∇_U J
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
    double obj_lin_iinit_re = 0.0;
    double obj_lin_iinit_im = 0.0;
    Vec ui = timestepper->getFinalState(iinit);
    optim_target->evalJ(ui, &obj_lin_iinit_re, &obj_lin_iinit_im);
    obj_cost_re += obj_weights[iinit] * obj_lin_iinit_re;
    obj_cost_im += obj_weights[iinit] * obj_lin_iinit_im;
  }
  // Gather result J_inf(U(T))
  double mycost_re = obj_cost_re;
  double mycost_im = obj_cost_im;
  MPI_Allreduce(&mycost_re, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(&mycost_im, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  // Compute -∇_U J and project to tangent space
  PetscScalar *b_data;
  VecGetArray(b, &b_data);
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    Vec ui = timestepper->getFinalState(iinit);

    // Use pre-allocated workspace vector to compute -∇_U J
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

  // Scale to get b = -∇_U J
  VecScale(b, -1.0);

  // Set the matrix operators
  KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);

  // Set initial guess (in tangent space)
  // Check if low-rank preconditioner will be used (defined later but needed here)
  bool use_lowrank_check = (gn_prec_update_interval >= 0 || gn_prec_update_interval == -1) && (gn_dual_precond_type == "low_rank");

  if (use_lowrank_check && gn_dual_lowrank_U != nullptr) {
    // Convert stored eigendata to format needed by computeEigenspaceApproximation
    // Eigenvalues: need Σ_r (not Σ_r^{-1/2}) for the function
    int rank = gn_dual_precond_rank;
    PetscScalar *eigenvalues = new PetscScalar[rank];

    // Extract eigenvalues from stored Σ_r^{-1/2}: Σ_r = (Σ_r^{-1/2})^{-2}
    PetscScalar *sigma_invsqrt_data;
    VecGetArray(gn_dual_lowrank_sigma_invsqrt, &sigma_invsqrt_data);
    for (int i = 0; i < rank; i++) {
        eigenvalues[i] = 1.0 / (sigma_invsqrt_data[i] * sigma_invsqrt_data[i]);
    }
    VecRestoreArray(gn_dual_lowrank_sigma_invsqrt, &sigma_invsqrt_data);

    // Check if complement is available (μ > 0)
    bool have_complement = (gn_dual_precond_mu > 0.0);

    // Compute eigenspace approximation from RHS b and stored eigendecomposition
    computeEigenspaceApproximation(b, gn_dual_lowrank_U, eigenvalues, rank,
                                   have_complement, gn_dual_precond_mu, y_dual_workspace);

    delete[] eigenvalues;

    if (mpirank_world == 0 && !quietmode) {
        printf("  Using precomputed eigenspace approximation as initial guess (rank %d)\n", rank);
    }
  } else {
    // Use provided initial guess (typically zero)
    VecCopy(initial_guess, y_dual_workspace);
  }
  KSPSetInitialGuessNonzero(ksp_GN_Dual, PETSC_TRUE);

  // Monitor residual and solution norm at every iteration
  KSPMonitorCancel(ksp_GN_Dual);
  KSPMonitorSet(ksp_GN_Dual, KSPMonitorResidualAndSolution, (void*)this, NULL);

  // Optional Levenberg-Marquardt damping: (L L^T + λI)
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, gn_normaleq_damping);
    // Also shift the preconditioned operator
    if (gn_prec_update_interval >= 0 || gn_prec_update_interval == -1) {
      MatShift(GN_Dual_Precond_MatShell, gn_normaleq_damping);
    }
  }

  // Apply preconditioning if enabled (interval >= 0 or adaptive mode = -1)
  Vec b_scaled, z_solution;
  bool use_precond = (gn_prec_update_interval >= 0 || gn_prec_update_interval == -1);
  bool use_symmetric_precond = use_precond && (gn_dual_precond_type == "symmetric");
  bool use_left_precond = use_precond && (gn_dual_precond_type == "left");
  bool use_lowrank_precond = use_precond && (gn_dual_precond_type == "low_rank");

  if (use_lowrank_precond) {
    // Low-rank preconditioning: M^{-1/2} (L L^T + λI) M^{-1/2} z = M^{-1/2} b
    // where M ≈ U Σ U^T is rank-r approximation

    // Scale RHS: b_scaled = M^{-1/2} * b
    VecDuplicate(b, &b_scaled);
    GN_Dual_LowRank_ApplyPrecond(GN_Dual_LowRank_MatShell, b, b_scaled);

    // Transform initial guess from y-space to z-space: z0 = M^{1/2} * y0
    // Mathematical insight: If y0 ≈ M^{-1} * b, then z0 = M^{1/2} * M^{-1} * b = M^{-1/2} * b
    // Therefore z0 should be b_scaled (the preconditioned RHS)
    VecDuplicate(y_dual_workspace, &z_solution);
    VecCopy(b_scaled, z_solution);  // Use preconditioned RHS as initial guess in z-space

    // Use the low-rank preconditioned operator
    KSPSetOperators(ksp_GN_Dual, GN_Dual_LowRank_MatShell, GN_Dual_LowRank_MatShell);

    // Solve: M^{-1/2} A M^{-1/2} z = M^{-1/2} b
    KSPSolve(ksp_GN_Dual, b_scaled, z_solution);

    // Recover y = M^{-1/2} * z
    GN_Dual_LowRank_ApplyPrecond(GN_Dual_LowRank_MatShell, z_solution, y_dual_workspace);

    VecDestroy(&b_scaled);
    VecDestroy(&z_solution);

    // Reset operators to unpreconditioned for next time
    KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);
  } else if (use_symmetric_precond) {
    // Symmetric preconditioning: D^{-1/2} (L L^T + λI) D^{-1/2} z = D^{-1/2} b
    // where D = diag(L L^T + λI), and we stored D^{-1/2} in gn_dual_precond_diag

    // Scale RHS: b_scaled = D^{-1/2} * b
    VecDuplicate(b, &b_scaled);
    VecPointwiseMult(b_scaled, gn_dual_precond_diag, b);

    // Scale initial guess: z0 = D^{1/2} * y0 = y0 / D^{-1/2}
    VecDuplicate(y_dual_workspace, &z_solution);
    VecPointwiseDivide(z_solution, y_dual_workspace, gn_dual_precond_diag);

    // Use the preconditioned operator
    KSPSetOperators(ksp_GN_Dual, GN_Dual_Precond_MatShell, GN_Dual_Precond_MatShell);

    // Solve: D^{-1/2} A D^{-1/2} z = D^{-1/2} b
    KSPSolve(ksp_GN_Dual, b_scaled, z_solution);

    // Recover y = D^{-1/2} * z
    VecPointwiseMult(y_dual_workspace, gn_dual_precond_diag, z_solution);

    VecDestroy(&b_scaled);
    VecDestroy(&z_solution);

    // Reset operators to unpreconditioned for next time
    KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);
  } else if (use_left_precond) {
    // Left preconditioning: D^{-1} (L L^T + λI) y = D^{-1} b
    // where D = diag(L L^T + λI), and we stored D^{-1} in gn_dual_precond_diag

    // Scale RHS: b_scaled = D^{-1} * b
    VecDuplicate(b, &b_scaled);
    VecPointwiseMult(b_scaled, gn_dual_precond_diag, b);

    // Use the preconditioned operator (applies D^{-1} on left)
    KSPSetOperators(ksp_GN_Dual, GN_Dual_Precond_MatShell, GN_Dual_Precond_MatShell);

    // Solve: D^{-1} A y = D^{-1} b
    KSPSolve(ksp_GN_Dual, b_scaled, y_dual_workspace);

    VecDestroy(&b_scaled);

    // Reset operators to unpreconditioned for next time
    KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);
  } else {
    // No preconditioning - solve directly
    KSPSolve(ksp_GN_Dual, b, y_dual_workspace);
  }

  // Revert the damping
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, -gn_normaleq_damping);
    if (gn_prec_update_interval >= 0 || gn_prec_update_interval == -1) {
      MatShift(GN_Dual_Precond_MatShell, -gn_normaleq_damping);
    }
  }

  // Cleanup
  VecDestroy(&b);

  // Recover primal solution: v = L^T y (use pure Jacobian transpose for dual mode)
  GN_Jacobian_MatMultTranspose(GN_Dual_MatShell, y_dual_workspace, v_out);

  // Report convergence
  KSPConvergedReason reason;
  int iters;
  double rnorm;
  KSPGetConvergedReason(ksp_GN_Dual, &reason);
  KSPGetIterationNumber(ksp_GN_Dual, &iters);
  KSPGetResidualNorm(ksp_GN_Dual, &rnorm);
  ksp_iters_last = iters;
  if (mpirank_world == 0 && !quietmode) {
    KSPType ksp_type;
    KSPGetType(ksp_GN_Dual, &ksp_type);
    printf("Gauss-Newton Dual KSP (%s) stats: iterations = %d, residual norm = %1.14e\n", ksp_type, iters, rnorm);
  }
}


void OptimProblem::solveGaussNewtonDualProjected(Vec xinit, Vec v_out){
  // Projected subspace solver for dual Gauss-Newton system
  // Algorithm:
  //   1. Compute eigendecomposition: (L L^T + λI) V_r = Σ_r V_r (r largest eigenvalues)
  //   2. Form RHS: b = -∇_U J in tangent space
  //   3. Project: b̃ = V_r^T b
  //   4. Solve in subspace: ỹ = Σ_r^{-1} b̃ (direct diagonal solve)
  //   5. Lift to full space: y = V_r ỹ
  //   6. Recover primal: v = L^T y

  // Store the point of evaluation
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false;

  // Compute RHS: b = -∇_U J
  Vec b;
  GN_Dual_MatCreateVecs(GN_Dual_MatShell, nullptr, &b);

  PetscInt tangent_dim;
  VecGetSize(b, &tangent_dim);

  // Get U(T) from optim_target
  Mat U_final_re = optim_target->getFinalUnitaryRe();
  Mat U_final_im = optim_target->getFinalUnitaryIm();

  // Compute J(U(T)) needed for ∇_U J
  double obj_cost_re = 0.0;
  double obj_cost_im = 0.0;
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
    double obj_lin_iinit_re = 0.0;
    double obj_lin_iinit_im = 0.0;
    Vec ui = timestepper->getFinalState(iinit);
    optim_target->evalJ(ui, &obj_lin_iinit_re, &obj_lin_iinit_im);
    obj_cost_re += obj_weights[iinit] * obj_lin_iinit_re;
    obj_cost_im += obj_weights[iinit] * obj_lin_iinit_im;
  }
  MPI_Allreduce(MPI_IN_PLACE, &obj_cost_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  MPI_Allreduce(MPI_IN_PLACE, &obj_cost_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  // Compute -∇_U J and project to tangent space
  PetscScalar *b_data;
  VecGetArray(b, &b_data);
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    int iinit_global = mpirank_init * ninit_local + iinit;
    Vec ui = timestepper->getFinalState(iinit);
    Vec b_iinit = wsub_workspace[iinit];
    VecCopy(ui, b_iinit);
    double obj_cost_re_bar, obj_cost_im_bar;
    VecScale(b_iinit, 2.0 / mastereq->getDim());
    optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
    optim_target->finalizeJ_diff(obj_cost_re, obj_cost_im, &obj_cost_re_bar, &obj_cost_im_bar);
    double scale = 1.0 / mastereq->getDim();
    optim_target->evalJ_diff(ui, b_iinit, scale*obj_cost_re_bar, scale*obj_cost_im_bar);
    GellMann::projectVecToTangentSpace(b_iinit, U_final_re, U_final_im, mastereq->getDim(), iinit_global, b_data);
  }
  MPI_Allreduce(MPI_IN_PLACE, b_data, tangent_dim, MPIU_SCALAR, MPI_SUM, comm_init);
  VecRestoreArray(b, &b_data);
  VecScale(b, -1.0);

  // Apply damping to the operator
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, gn_normaleq_damping);
  }

  // Compute eigendecomposition: (L L^T + λI) V_{r+1} = Σ_{r+1} V_{r+1}
  // We compute r+1 eigenvalues to get μ = λ_{r+1} for the complement space
  int rank = gn_dual_precond_rank;
  int num_eigs = (rank < tangent_dim) ? rank + 1 : rank;  // Compute r+1 unless at full rank

  if (mpirank_world == 0 && !quietmode) {
    printf("Computing %d largest eigenvalues/eigenvectors of dual GN matrix (n=%d)...\n", num_eigs, (int)tangent_dim);
  }

  EPS eps;
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, GN_Dual_MatShell, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
  EPSSetDimensions(eps, num_eigs, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSSetFromOptions(eps);
  EPSSolve(eps);

  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  if (mpirank_world == 0 && !quietmode) {
    printf("EPS converged %d/%d eigenvalues\n", (int)nconv, num_eigs);
  }

  // Check if we got enough eigenvalues
  bool have_complement = false;
  PetscScalar mu = 0.0;  // μ = λ_{r+1} for complement space

  if (nconv < rank) {
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Only %d eigenvalues converged (requested rank %d). Using rank %d.\n", (int)nconv, rank, (int)nconv);
    }
    rank = nconv;
    have_complement = false;
  } else if (rank < tangent_dim && nconv >= rank + 1) {
    // We have r+1 eigenvalues, use the (r+1)-th for complement
    have_complement = true;
  }

  // Extract eigenvectors into V_r (only first r)
  Mat V_r;
  MatCreateSeqDense(PETSC_COMM_SELF, tangent_dim, rank, NULL, &V_r);

  Vec eigvec_real, eigvec_imag;
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &eigvec_real);
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &eigvec_imag);

  // Allocate eigenvalues array for actual number converged (may be > num_eigs requested)
  PetscScalar *eigenvalues = new PetscScalar[nconv];

  for (int i = 0; i < nconv; i++) {
    PetscScalar eigval_r, eigval_i;
    EPSGetEigenpair(eps, i, &eigval_r, &eigval_i, eigvec_real, eigvec_imag);
    eigenvalues[i] = PetscRealPart(eigval_r);

    // Store only first r eigenvectors in V_r
    if (i < rank) {
      Vec v_col;
      MatDenseGetColumnVecWrite(V_r, i, &v_col);
      VecCopy(eigvec_real, v_col);
      MatDenseRestoreColumnVecWrite(V_r, i, &v_col);
    } else if (i == rank && have_complement) {
      // This is λ_{r+1}, use it for μ
      mu = eigenvalues[i];
    }
  }

  MatAssemblyBegin(V_r, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(V_r, MAT_FINAL_ASSEMBLY);

  VecDestroy(&eigvec_real);
  VecDestroy(&eigvec_imag);
  EPSDestroy(&eps);

  if (mpirank_world == 0 && !quietmode) {
    printf("Eigenvalue range: λ_1 = %e (largest), λ_%d = %e (smallest)\n",
           eigenvalues[0], rank, eigenvalues[rank-1]);
    if (have_complement) {
      printf("Using μ = λ_%d = %e for complement space approximation\n", rank+1, mu);
    }
  }

  // Compute eigenspace approximation: y = M^{-1} * b
  computeEigenspaceApproximation(b, V_r, eigenvalues, rank, have_complement, mu, y_dual_workspace);

  // Optional iterative refinement: refine the projected solution with KSP iterations
  // Refinement is enabled when gn_maxiter > 0
  // Reuse the eigendecomposition (V_r, Σ_r, μ) to construct low-rank preconditioner at zero cost
  int refine_iters = 0;
  if (gn_maxiter > 0) {
    // Normalize the RHS to avoid small ||b|| issues near optimum
    // Solve: (L L^T + λI) y = b_normalized, then scale: y_final = y * bnorm
    double bnorm;
    VecNorm(b, NORM_2, &bnorm);

    // Only refine if RHS is non-negligible
    if (bnorm > 1e-14) {
      // Normalize RHS: b_normalized = b / ||b||
      Vec b_normalized;
      VecDuplicate(b, &b_normalized);
      VecCopy(b, b_normalized);
      VecScale(b_normalized, 1.0 / bnorm);

      // Store projected solution as initial guess (also needs normalization)
      Vec y_proj;
      VecDuplicate(y_dual_workspace, &y_proj);
      VecCopy(y_dual_workspace, y_proj);
      VecScale(y_proj, 1.0 / bnorm);  // Normalize to match normalized RHS

      // Construct low-rank preconditioner M = V_r Σ_r V_r^T + μ(I - V_r V_r^T)
      // using the eigendecomposition already computed above (zero additional cost!)

      // Store eigenvalues^{-1/2} for preconditioning
      Vec sigma_invsqrt_temp;
      VecCreateSeq(PETSC_COMM_SELF, rank, &sigma_invsqrt_temp);
      PetscScalar *sigma_data;
      VecGetArray(sigma_invsqrt_temp, &sigma_data);
      for (int i = 0; i < rank; i++) {
        sigma_data[i] = 1.0 / sqrt(eigenvalues[i]);
      }
      VecRestoreArray(sigma_invsqrt_temp, &sigma_data);

      // Temporarily store the eigendecomposition in the preconditioner data structures
      // Save current state
      Mat V_r_saved = gn_dual_lowrank_U;
      Vec sigma_saved = gn_dual_lowrank_sigma_invsqrt;
      PetscScalar mu_saved = gn_dual_precond_mu;

      // Set preconditioner data to current eigendecomposition
      gn_dual_lowrank_U = V_r;  // V_r already computed above
      gn_dual_lowrank_sigma_invsqrt = sigma_invsqrt_temp;
      gn_dual_precond_mu = have_complement ? mu : eigenvalues[rank-1];  // Use μ if available

      // Apply symmetric low-rank preconditioning: M^{-1/2} A M^{-1/2} z = M^{-1/2} b_normalized
      Vec b_scaled, z_solution;
      VecDuplicate(b_normalized, &b_scaled);
      VecDuplicate(y_proj, &z_solution);

      // Scale normalized RHS: b_scaled = M^{-1/2} * b_normalized
      GN_Dual_LowRank_ApplyPrecond(GN_Dual_LowRank_MatShell, b_normalized, b_scaled);

      // Scale initial guess
      VecCopy(y_proj, z_solution);

      // Configure KSP for refinement with preconditioned operator
      KSPSetOperators(ksp_GN_Dual, GN_Dual_LowRank_MatShell, GN_Dual_LowRank_MatShell);
      KSPSetInitialGuessNonzero(ksp_GN_Dual, PETSC_TRUE);

      // Use pure relative tolerance - no absolute tolerance needed since ||b_normalized|| = 1
      KSPSetTolerances(ksp_GN_Dual, gn_rtol, PETSC_DEFAULT, PETSC_DEFAULT, gn_maxiter);

      // Solve preconditioned system: M^{-1/2} A M^{-1/2} z = M^{-1/2} b_normalized
      KSPSolve(ksp_GN_Dual, b_scaled, z_solution);

      // Recover y = M^{-1/2} * z
      GN_Dual_LowRank_ApplyPrecond(GN_Dual_LowRank_MatShell, z_solution, y_dual_workspace);

      // Scale back by ||b|| to get the actual solution
      VecScale(y_dual_workspace, bnorm);

    // Get actual number of iterations performed
    KSPGetIterationNumber(ksp_GN_Dual, &refine_iters);

      // Get actual number of iterations performed
      KSPGetIterationNumber(ksp_GN_Dual, &refine_iters);

      // Restore preconditioner data structures
      gn_dual_lowrank_U = V_r_saved;
      gn_dual_lowrank_sigma_invsqrt = sigma_saved;
      gn_dual_precond_mu = mu_saved;

      // Cleanup
      VecDestroy(&b_scaled);
      VecDestroy(&z_solution);
      VecDestroy(&sigma_invsqrt_temp);
      VecDestroy(&y_proj);
      VecDestroy(&b_normalized);

      // Reset KSP operators to unpreconditioned for next time
      KSPSetOperators(ksp_GN_Dual, GN_Dual_MatShell, GN_Dual_MatShell);
    }
  }

  // Revert damping
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, -gn_normaleq_damping);
  }

  // Cleanup
  VecDestroy(&b);
  MatDestroy(&V_r);
  delete[] eigenvalues;

  // Recover primal solution: v = L^T y (use pure Jacobian transpose for dual mode)
  GN_Jacobian_MatMultTranspose(GN_Dual_MatShell, y_dual_workspace, v_out);

  // Store KSP iterations for reporting in optimization history
  ksp_iters_last = refine_iters;

  if (mpirank_world == 0 && !quietmode) {
    if (refine_iters > 0) {
      printf("Gauss-Newton Dual Projected: solved using rank-%d subspace + %d refinement iterations\n", rank, refine_iters);
    } else if (have_complement) {
      printf("Gauss-Newton Dual Projected: solved using rank-%d subspace + complement approximation (no KSP iterations)\n", rank);
    } else {
      printf("Gauss-Newton Dual Projected: solved in rank-%d subspace (no KSP iterations)\n", rank);
    }
  }
}


void OptimProblem::solveGaussNewtonNormalEqEPS(Vec xinit, const Vec b, Vec Ainv_b){

  // Get eigenvalues and eigenvectors of the Gauss-Newton matrix
  Mat evecs;
  std::vector<double> evals = computeGaussNewtonNormalEqEvals(xinit, &evecs);

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
      tmp_data[i] /= (evals[i] + eps_normaleq_damping);
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

std::vector<double> OptimProblem::computeGaussNewtonNormalEqEvals(Vec xinit, Mat* evecs_out){

  // Store xinit so the MatShell can use it as point of evaluation.
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Update EPS for a fresh solve on the new xeval_GN. 
  if (!GN_NormalEq_usedensemat) { 
    // use MatShell 
    EPSSetOperators(eps_GN_NormalEq, GN_NormalEq_MatShell, NULL);
  } else {
    // use dense Matrix representation of the Gauss-Newton matrix
    updateGaussNewtonMatDense(); // does N2-1 applications of the MatShell
    EPSSetOperators(eps_GN_NormalEq, GN_NormalEq_MatDense, NULL);
  }

  // Solve the eigenvalue problem for the Gauss-Newton matrix
  EPSSolve(eps_GN_NormalEq);
  PetscInt numConv;
  PetscInt iters_taken;
  EPSGetConverged(eps_GN_NormalEq, &numConv);
  EPSGetIterationNumber(eps_GN_NormalEq,&iters_taken);
  if (mpirank_world == 0 && !quietmode) printf("Gauss-Newton EPS converged %d eigenvalues in %d iterations.\n", numConv, iters_taken);
  if (numConv < neigvals) {
      if (mpirank_world==0 && !quietmode) printf("WARNING: Only %d eigenvalues out of %d eigenvalues converged.\n", numConv, neigvals);
  }

  // Set up storage for eigenvalues and eigenvectors (should be real!)
  std::vector<double> evals_re(neigvals);
  std::vector<Vec> evec_re(neigvals);
  for (int ix = 0; ix < neigvals; ix++) {
    MatCreateVecs(GN_NormalEq_MatShell, &evec_re[ix], NULL);
  }

  // Retrieve eigenpairs of M and compute error. 
  for (PetscInt i = 0; i < numConv && i < neigvals; i++) {

    // Retrieve the eigenvalue (is real) and eigenvector
    EPSGetEigenpair(eps_GN_NormalEq, i, &evals_re[i], NULL, evec_re[i], NULL);

    // Estimate the errror (needs one more application of A)
    double error = 0.0;
    EPSComputeError(eps_GN_NormalEq,i,EPS_ERROR_RELATIVE,&error);
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
      Vec G, Gprec, xnew, step, v_zero, y_zero_dual;
      VecDuplicate(x_GN, &G);
      VecDuplicate(x_GN, &Gprec);
      VecDuplicate(x_GN, &xnew);
      VecDuplicate(x_GN, &step);

      // Note: Preconditioner initialization moved to first iteration (after gradient eval)
      // because low-rank preconditioner needs valid Jacobian matrix

      bool stop = false;
      for (int iter = 0; !stop; iter++) {
        // Initialize dual preconditioner on first iteration (after gradient is computed)
        // Skip preconditioner initialization when using projected solver
        bool enable_precond = (gn_prec_update_interval >= 0 || gn_prec_update_interval == -1) && (gn_normaleq_solver != "proj");
        if (iter == 0 && gn_solver_mode == "dual" && enable_precond && gn_dual_precond_type != "low_rank") {
          // Diagonal preconditioner can be initialized before gradient eval
          if (mpirank_world == 0 && !quietmode) {
            printf("Initializing dual Gauss-Newton Jacobi preconditioner...\n");
          }
          setGaussNewtonEvalPoint(x_GN);
          updateGaussNewtonDualPreconditioner();
        }
        // Update dual preconditioner: fixed interval or adaptive
        bool should_update = false;
        if (gn_solver_mode == "dual" && iter > 0) {
          if (gn_prec_update_interval > 0 && iter % gn_prec_update_interval == 0) {
            // Fixed interval update
            should_update = true;
            if (mpirank_world == 0 && !quietmode) {
              printf("Iteration %d: Updating preconditioner (fixed interval = %d)\n", iter, gn_prec_update_interval);
            }
          } else if (gn_prec_update_interval == -1 && gn_dual_ksp_baseline > 0) {
            // Adaptive update: accumulate excess KSP iterations beyond n_0
            int excess = ksp_iters_last - gn_dual_ksp_baseline;
            if (excess > 0) {
              gn_dual_ksp_excess_accum += excess;
              if (mpirank_world == 0 && !quietmode) {
                printf("Iteration %d: KSP excess = %d (accumulated = %d, threshold = %d)\n",
                       iter, excess, gn_dual_ksp_excess_accum, gn_dual_precond_rank);
              }
            }

            // Update when accumulated excess exceeds r
            if (gn_dual_ksp_excess_accum >= gn_dual_precond_rank) {
              should_update = true;
              if (mpirank_world == 0 && !quietmode) {
                printf("Iteration %d: Updating preconditioner (accumulated excess %d >= r = %d)\n",
                       iter, gn_dual_ksp_excess_accum, gn_dual_precond_rank);
              }
            }
          }

          if (should_update) {
            setGaussNewtonEvalPoint(x_GN);
            if (gn_dual_precond_type == "low_rank") {
              updateGaussNewtonDualPreconditionerLowRank(gn_dual_precond_rank);
            } else {
              updateGaussNewtonDualPreconditioner();
            }
            // Reset accumulator after update
            if (gn_prec_update_interval == -1) {
              gn_dual_ksp_excess_accum = 0;
            }
          }
        }

        // Compute the gradient (and objective) at the current iterate
        evalGradF(x_GN, G);
        double f = objective;
        VecNorm(G, NORM_2, &gnorm);

        // Initialize low-rank preconditioner on first iteration (after gradient is available)
        if (iter == 0 && gn_solver_mode == "dual" && enable_precond && gn_dual_precond_type == "low_rank") {
          if (mpirank_world == 0 && !quietmode) {
            printf("Initializing dual Gauss-Newton low-rank preconditioner (rank=%d)...\n", gn_dual_precond_rank);
          }
          setGaussNewtonEvalPoint(x_GN);
          updateGaussNewtonDualPreconditionerLowRank(gn_dual_precond_rank);
        }

        // Precondition the gradient: solve the Gauss-Newton system A(x)*Gprec = -G
        // Choose solver based on gn_solver_mode configuration ("primal" or "dual")
        if (gn_solver_mode == "dual") {
          // Use dual formulation: (L L^T + λI) y = -∇_U J, then v = L^T y
          MatCreateVecs(GN_Dual_MatShell, nullptr, &y_zero_dual);
          VecZeroEntries(y_zero_dual);
          solveGaussNewtonDual(x_GN, y_zero_dual, Gprec);
          VecDestroy(&y_zero_dual);
          VecScale(Gprec, -1.0);

          // Update baseline KSP iterations after preconditioner formation
          if ((iter == 0 && enable_precond) || should_update) {
            gn_dual_ksp_baseline = ksp_iters_last;
            if (mpirank_world == 0 && !quietmode) {
              printf("  Baseline KSP iterations set to n_0 = %d\n", gn_dual_ksp_baseline);
            }
          }
        } else {
          // Use primal formulation via normal equations KSP solver: (L^T L + λI) v = -G
          VecDuplicate(x_GN, &v_zero);
          Vec neg_G;
          VecDuplicate(G, &neg_G);
          VecCopy(G, neg_G);
          VecScale(neg_G, -1.0);
          VecZeroEntries(v_zero);
          solveGaussNewtonNormalEqKSP(x_GN, v_zero, neg_G, Gprec);
          VecDestroy(&v_zero);
          VecDestroy(&neg_G);
          VecScale(Gprec, -1.0);
        } 

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

    case OptimSolverType::LEAST_SQUARES: {
      // Least-squares Gauss-Newton solver: min ||L*v - b||²
      // Uses either LSQR or TAO BRGN depending on leastsquares_solver config

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

        // Solve the least-squares Gauss-Newton system: min ||L*v - b||²
        solveGaussNewtonLeastSquares(x_GN, v_zero, Gprec);
        VecScale(Gprec, -1.0);

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
  } else if (optim_solver_type == OptimSolverType::GAUSS_NEWTON ||
             optim_solver_type == OptimSolverType::LEAST_SQUARES) {
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
  MatZeroEntries(GN_NormalEq_MatDense);

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
    MatMult(this->getGN_NormalEq_MatShell(), e, Av);
    for (int jx = 0; jx < ndesign; ++jx) {
      PetscScalar val;
      VecGetValues(Av, 1, &jx, &val);
      MatSetValue(GN_NormalEq_MatDense, jx, ix, val, INSERT_VALUES);
    }
  }
  VecDestroy(&e);
  VecDestroy(&Av);

  MatAssemblyBegin(GN_NormalEq_MatDense, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(GN_NormalEq_MatDense, MAT_FINAL_ASSEMBLY);

  // Allreduce to combine contributions from all MPI ranks
  PetscScalar *data;
  MatDenseGetArray(GN_NormalEq_MatDense, &data);
  int size = ndesign * ndesign;
  MPI_Allreduce(MPI_IN_PLACE, data, size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(GN_NormalEq_MatDense, &data);

}


void OptimProblem::updateGaussNewtonDualMatDense(){
  MatZeroEntries(GN_Dual_MatDense);

  // Get tangent space dimension (dim²-1)
  PetscInt tangent_dim;
  MatGetSize(GN_Dual_MatShell, &tangent_dim, NULL);

  Vec e;
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &e);
  Vec Ay;
  VecDuplicate(e, &Ay);

  // Parallelize over comm_optim threads
  int ncols_per_rank = tangent_dim / mpisize_optim;
  int ncols_local = ncols_per_rank;
  if (mpirank_optim == mpisize_optim - 1) {
    ncols_local = tangent_dim - mpirank_optim * ncols_per_rank;
  }

  // Iterate over local columns: apply GN_Dual_MatShell to each unit vector
  for (int ix_local = 0; ix_local < ncols_local; ++ix_local) {
    int ix = ix_local + mpirank_optim * ncols_per_rank;

    VecSet(e, 0.0);
    VecSetValue(e, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);
    MatMult(GN_Dual_MatShell, e, Ay);
    for (int jy = 0; jy < tangent_dim; ++jy) {
      PetscScalar val;
      VecGetValues(Ay, 1, &jy, &val);
      MatSetValue(GN_Dual_MatDense, jy, ix, val, INSERT_VALUES);
    }
  }
  VecDestroy(&e);
  VecDestroy(&Ay);

  MatAssemblyBegin(GN_Dual_MatDense, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(GN_Dual_MatDense, MAT_FINAL_ASSEMBLY);

  // Allreduce to combine contributions from all MPI ranks
  PetscScalar *dual_data;
  MatDenseGetArray(GN_Dual_MatDense, &dual_data);
  int dual_size = tangent_dim * tangent_dim;
  MPI_Allreduce(MPI_IN_PLACE, dual_data, dual_size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(GN_Dual_MatDense, &dual_data);

}


void OptimProblem::updateGaussNewtonDualPrecondMatDense(){
  MatZeroEntries(GN_Dual_Precond_MatDense);

  // Get tangent space dimension (dim²-1)
  PetscInt tangent_dim;
  MatGetSize(GN_Dual_Precond_MatShell, &tangent_dim, NULL);

  Vec e;
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &e);
  Vec Ay;
  VecDuplicate(e, &Ay);

  // Parallelize over comm_optim threads
  int ncols_per_rank = tangent_dim / mpisize_optim;
  int ncols_local = ncols_per_rank;
  if (mpirank_optim == mpisize_optim - 1) {
    ncols_local = tangent_dim - mpirank_optim * ncols_per_rank;
  }

  // Iterate over local columns: apply preconditioned MatShell to each unit vector
  for (int ix_local = 0; ix_local < ncols_local; ++ix_local) {
    int ix = ix_local + mpirank_optim * ncols_per_rank;

    VecSet(e, 0.0);
    VecSetValue(e, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);
    MatMult(GN_Dual_Precond_MatShell, e, Ay);
    for (int jy = 0; jy < tangent_dim; ++jy) {
      PetscScalar val;
      VecGetValues(Ay, 1, &jy, &val);
      MatSetValue(GN_Dual_Precond_MatDense, jy, ix, val, INSERT_VALUES);
    }
  }
  VecDestroy(&e);
  VecDestroy(&Ay);

  MatAssemblyBegin(GN_Dual_Precond_MatDense, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(GN_Dual_Precond_MatDense, MAT_FINAL_ASSEMBLY);

  // Allreduce to combine contributions from all MPI ranks
  PetscScalar *precond_data;
  MatDenseGetArray(GN_Dual_Precond_MatDense, &precond_data);
  int precond_size = tangent_dim * tangent_dim;
  MPI_Allreduce(MPI_IN_PLACE, precond_data, precond_size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(GN_Dual_Precond_MatDense, &precond_data);

}


void OptimProblem::updateGaussNewtonDualPreconditioner(){
  // Extract diagonal of (L L^T + λI) by applying MatShell to unit vectors

  // Get tangent space dimension (dim²-1)
  PetscInt tangent_dim;
  VecGetSize(gn_dual_precond_diag, &tangent_dim);

  Vec e;
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &e);
  Vec Ay;
  VecDuplicate(e, &Ay);

  // Parallelize over comm_optim threads (though for diagonal we only need diagonal entries)
  int ncols_per_rank = tangent_dim / mpisize_optim;
  int ncols_local = ncols_per_rank;
  if (mpirank_optim == mpisize_optim - 1) {
    ncols_local = tangent_dim - mpirank_optim * ncols_per_rank;
  }

  // Extract diagonal: d_ii = (GN_Dual_MatShell * e_i)_i
  PetscScalar *diag_data;
  VecGetArray(gn_dual_precond_diag, &diag_data);

  for (int ix_local = 0; ix_local < ncols_local; ++ix_local) {
    int ix = ix_local + mpirank_optim * ncols_per_rank;

    // Create unit vector e_i
    VecSet(e, 0.0);
    VecSetValue(e, ix, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);

    // Apply (L L^T + λI) to e_i
    MatMult(GN_Dual_MatShell, e, Ay);

    // Extract diagonal entry
    PetscScalar val;
    VecGetValues(Ay, 1, &ix, &val);
    diag_data[ix] = val;
  }

  // Allreduce BEFORE restoring the array (while we still have the pointer)
  MPI_Allreduce(MPI_IN_PLACE, diag_data, tangent_dim, MPIU_SCALAR, MPI_SUM, comm_optim);

  // Compute preconditioner: D^{-1/2} (symmetric) or D^{-1} (left)
  bool use_symmetric = (gn_dual_precond_type == "symmetric");

  for (PetscInt i = 0; i < tangent_dim; i++) {
    if (PetscAbsScalar(diag_data[i]) > 1e-14) {
      if (use_symmetric) {
        diag_data[i] = 1.0 / PetscSqrtScalar(diag_data[i]);  // D^{-1/2}
      } else {
        diag_data[i] = 1.0 / diag_data[i];  // D^{-1}
      }
    } else {
      diag_data[i] = 1.0;  // Avoid division by zero
    }
  }

  VecRestoreArray(gn_dual_precond_diag, &diag_data);

  VecDestroy(&e);
  VecDestroy(&Ay);

  if (mpirank_world == 0 && !quietmode) {
    if (use_symmetric) {
      printf("Updated dual Gauss-Newton symmetric diagonal preconditioner (D^{-1/2})\n");
    } else {
      printf("Updated dual Gauss-Newton left diagonal preconditioner (D^{-1})\n");
    }
  }
}


void OptimProblem::updateGaussNewtonDualPreconditionerLowRank(int rank){
  // Construct full-rank preconditioner M = V_r Σ_r V_r^T + μ(I - V_r V_r^T)
  // Algorithm:
  //   1. Use EPS to compute r+1 largest eigenvalues/eigenvectors of A
  //   2. Store V_r (first r eigenvectors), Σ_r (first r eigenvalues), μ = σ_{r+1}
  //   3. Apply M^{-1/2} = V_r Σ_r^{-1/2} V_r^T + μ^{-1/2}(I - V_r V_r^T)

  PetscInt tangent_dim;
  MatGetSize(GN_Dual_MatShell, &tangent_dim, NULL);

  if (rank >= tangent_dim) {
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: rank %d >= tangent_dim %d, setting rank = tangent_dim - 1\n", rank, tangent_dim);
    }
    rank = tangent_dim - 1;
  }

  // 1. Use EPS to compute r+1 largest eigenvalues/eigenvectors of A
  int nev = rank + 1;  // Request r+1 eigenvalues

  if (mpirank_world == 0 && !quietmode) {
    printf("Computing %d largest eigenvalues/eigenvectors of dual GN matrix (n=%d)...\n", nev, (int)tangent_dim);
  }

  EPS eps;
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, GN_Dual_MatShell, NULL);
  EPSSetProblemType(eps, EPS_HEP);  // Hermitian eigenvalue problem
  EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);  // Largest eigenvalues
  EPSSetDimensions(eps, nev, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSSetFromOptions(eps);
  EPSSolve(eps);

  // Check convergence
  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  if (mpirank_world == 0 && !quietmode) {
    printf("EPS converged %d/%d eigenvalues\n", (int)nconv, nev);
  }

  if (nconv < nev) {
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: Only %d eigenvalues converged (requested %d). Adjusting rank.\n", (int)nconv, nev);
    }
    if (nconv < 2) {
      printf("Error: Need at least 2 eigenvalues for low-rank preconditioner. Disabling.\n");
      EPSDestroy(&eps);
      return;
    }
    rank = nconv - 1;  // Use nconv-1 for preconditioner, last one for μ
  }

  // 2. Extract first r eigenvectors into V_r (stored in gn_dual_lowrank_U)
  MatZeroEntries(gn_dual_lowrank_U);

  Vec eigvec_real, eigvec_imag;
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &eigvec_real);
  VecCreateSeq(PETSC_COMM_SELF, tangent_dim, &eigvec_imag);

  // Store eigenvalues temporarily
  PetscScalar *eigenvalues = new PetscScalar[rank + 1];

  for (int i = 0; i < rank + 1; i++) {
    PetscScalar eigval_r, eigval_i;
    EPSGetEigenpair(eps, i, &eigval_r, &eigval_i, eigvec_real, eigvec_imag);
    eigenvalues[i] = PetscRealPart(eigval_r);

    // Store first r eigenvectors in gn_dual_lowrank_U
    if (i < rank) {
      Vec v_col;
      MatDenseGetColumnVecWrite(gn_dual_lowrank_U, i, &v_col);
      VecCopy(eigvec_real, v_col);
      MatDenseRestoreColumnVecWrite(gn_dual_lowrank_U, i, &v_col);
    }
  }

  MatAssemblyBegin(gn_dual_lowrank_U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(gn_dual_lowrank_U, MAT_FINAL_ASSEMBLY);

  VecDestroy(&eigvec_real);
  VecDestroy(&eigvec_imag);
  EPSDestroy(&eps);

  // 3. Store Σ_r^{-1/2} and μ^{-1/2}
  PetscScalar *sigma_invsqrt_data;
  VecGetArray(gn_dual_lowrank_sigma_invsqrt, &sigma_invsqrt_data);

  // First r eigenvalues → Σ_r^{-1/2}
  for (int i = 0; i < rank; i++) {
    if (eigenvalues[i] > 1e-14) {
      sigma_invsqrt_data[i] = 1.0 / PetscSqrtScalar(eigenvalues[i]);
    } else {
      printf("Warning: Eigenvalue %d = %e is too small, regularizing\n", i, eigenvalues[i]);
      sigma_invsqrt_data[i] = 1.0 / PetscSqrtScalar(1e-14);
    }
  }

  VecRestoreArray(gn_dual_lowrank_sigma_invsqrt, &sigma_invsqrt_data);

  // Store μ = σ_{r+1} in gn_dual_precond_mu (need to add this member variable)
  gn_dual_precond_mu = eigenvalues[rank];

  if (mpirank_world == 0 && !quietmode) {
    printf("Low-rank preconditioner updated:\n");
    printf("  Rank r = %d\n", rank);
    printf("  λ_1 = %e (largest)\n", eigenvalues[0]);
    printf("  λ_r = %e\n", eigenvalues[rank-1]);
    printf("  μ = λ_{r+1} = %e\n", gn_dual_precond_mu);
    printf("  μ^{-1/2} = %e\n", 1.0/PetscSqrtScalar(gn_dual_precond_mu));
  }

  delete[] eigenvalues;
}


// Apply M^{-1/2} where M ≈ U Σ U^T
// M^{-1/2} x = U * (Σ^{-1/2} .* (U^T * x))
void OptimProblem::GN_Dual_LowRank_ApplyPrecond(Mat A, const Vec x, Vec y){
  // Apply M^{-1/2} = V_r Σ_r^{-1/2} V_r^T + μ^{-1/2}(I - V_r V_r^T)
  // y = V_r (Σ_r^{-1/2} .* (V_r^T x)) + μ^{-1/2} (x - V_r (V_r^T x))

  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Step 1: temp = V_r^T * x (project onto top-r subspace)
  MatMultTranspose(self->gn_dual_lowrank_U, x, self->gn_dual_lowrank_temp);

  // Step 2: temp = Σ_r^{-1/2} .* temp (scale by inverse sqrt eigenvalues)
  VecPointwiseMult(self->gn_dual_lowrank_temp, self->gn_dual_lowrank_sigma_invsqrt, self->gn_dual_lowrank_temp);

  // Step 3: y = V_r * temp (part 1: V_r Σ_r^{-1/2} V_r^T x)
  MatMult(self->gn_dual_lowrank_U, self->gn_dual_lowrank_temp, y);

  // Step 4: Compute complement part: μ^{-1/2} (x - V_r V_r^T x)
  // First compute V_r^T x again (already in lowrank_temp before scaling)
  Vec vr_proj;
  VecDuplicate(x, &vr_proj);
  MatMultTranspose(self->gn_dual_lowrank_U, x, self->gn_dual_lowrank_temp);  // Recompute projection
  MatMult(self->gn_dual_lowrank_U, self->gn_dual_lowrank_temp, vr_proj);  // V_r V_r^T x

  // complement = x - V_r V_r^T x
  Vec complement;
  VecDuplicate(x, &complement);
  VecWAXPY(complement, -1.0, vr_proj, x);  // complement = x - vr_proj

  // Add μ^{-1/2} * complement to y
  PetscScalar mu_invsqrt = 1.0 / PetscSqrtScalar(self->gn_dual_precond_mu);
  VecAXPY(y, mu_invsqrt, complement);

  VecDestroy(&vr_proj);
  VecDestroy(&complement);
}


// Compute eigenspace approximation: y_approx = M^{-1} * b
// where M = V_r Σ_r V_r^T + μ(I - V_r V_r^T)
void OptimProblem::computeEigenspaceApproximation(
    const Vec b,                      // RHS vector (tangent space)
    const Mat V_r,                    // Eigenvectors (tangent_dim × rank)
    const PetscScalar* eigenvalues,   // Eigenvalues Σ_r (rank), NOT inverse
    int rank,                         // Number of eigenpairs
    bool have_complement,             // Whether we have μ = λ_{r+1}
    PetscScalar mu,                   // Complement eigenvalue
    Vec y_approx)                     // Output: eigenspace approximation
{
    // Step 1: Project RHS: b̃ = V_r^T b
    Vec b_tilde;
    VecCreateSeq(PETSC_COMM_SELF, rank, &b_tilde);
    MatMultTranspose(V_r, b, b_tilde);

    // Step 2: Solve in subspace: ỹ = Σ_r^{-1} b̃ (direct diagonal solve)
    Vec y_tilde;
    VecDuplicate(b_tilde, &y_tilde);

    PetscScalar *b_tilde_data, *y_tilde_data;
    VecGetArray(b_tilde, &b_tilde_data);
    VecGetArray(y_tilde, &y_tilde_data);

    for (int i = 0; i < rank; i++) {
        if (eigenvalues[i] > 1e-14) {
            y_tilde_data[i] = b_tilde_data[i] / eigenvalues[i];
        } else {
            y_tilde_data[i] = 0.0;
            if (mpirank_world == 0 && !quietmode) {
                printf("Warning: Eigenvalue %d = %e is too small, setting component to zero\n",
                       i, eigenvalues[i]);
            }
        }
    }

    VecRestoreArray(b_tilde, &b_tilde_data);
    VecRestoreArray(y_tilde, &y_tilde_data);

    // Step 3: Lift to full space: y = V_r ỹ (top-r subspace part)
    MatMult(V_r, y_tilde, y_approx);

    // Step 4: Add complement space contribution: y += μ^{-1} (I - V_r V_r^T) b
    if (have_complement && mu > 1e-14) {
        Vec temp;
        VecDuplicate(b, &temp);

        // Compute V_r (V_r^T b) to get projection onto top-r subspace
        Vec b_proj;
        VecDuplicate(b, &b_proj);
        MatMult(V_r, b_tilde, b_proj);  // Reuse b_tilde = V_r^T b from Step 1

        // Compute complement: (I - V_r V_r^T) b = b - V_r (V_r^T b)
        VecWAXPY(temp, -1.0, b_proj, b);  // temp = b - b_proj

        // Add scaled complement to solution: y += μ^{-1} * temp
        VecAXPY(y_approx, 1.0/mu, temp);

        VecDestroy(&temp);
        VecDestroy(&b_proj);
    }

    // Cleanup
    VecDestroy(&b_tilde);
    VecDestroy(&y_tilde);
}


// Apply M^{-1/2} A M^{-1/2} z
void OptimProblem::GN_Dual_LowRank_MatMult(Mat A, const Vec z, Vec output){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  Vec temp;
  VecDuplicate(z, &temp);

  // Step 1: temp = M^{-1/2} * z
  GN_Dual_LowRank_ApplyPrecond(A, z, temp);

  // Step 2: temp = A * temp
  MatMult(self->GN_Dual_MatShell, temp, output);

  // Step 3: output = M^{-1/2} * output
  VecCopy(output, temp);  // Save intermediate result
  GN_Dual_LowRank_ApplyPrecond(A, temp, output);

  VecDestroy(&temp);
}



void OptimProblem::solveGaussNewtonLeastSquares(const Vec xinit, const Vec initial_guess, Vec v_LeastSquares){

  // Store the point of evaluation (design parameters)
  VecCopy(xinit, xeval_GN);
  nonlinear_forward_valid = false; // Force a fresh nonlinear forward solve for the new xeval_GN

  // Fill the RHS: b = -\nabla_U J projected to tangent space
  // \nabla_J = -2/n U - 2/n^2Tr(V^dU)V = 2/n(I-P)U

  Vec b; // RHS
  GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &b);
  VecZeroEntries(b);  // Ensure bottom part (for damping rows) is initialized to zero

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
  PetscInt tangent_dim = mastereq->getDim()*mastereq->getDim()-1;
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

  if (includeHessUJ) {
    // TODO: Implement W^1/2 in tangent space if needed
    if (mpirank_world == 0) {
      printf("WARNING: W^1/2 is not implemented for LeastSquares solver tangent space projection.\n");
    }
  }

  // // Adjoint test: <w, A v> should equal <A^T w, v>  (all real)
  // printf("GN Least-Squares: Adjoint test\n");
  // Vec v_test, w_test, Av, ATw;
  // GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, &v_test, nullptr);
  // GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &w_test);
  // VecSetRandom(v_test, nullptr);
  // VecSetRandom(w_test, nullptr);
  // GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &Av);
  // GN_LeastSquares_MatMult(GN_LeastSquares_MatShell, v_test, Av);
  // GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, &ATw, nullptr);
  // GN_LeastSquares_MatMultTranspose(GN_LeastSquares_MatShell, w_test, ATw);

  // PetscScalar lhs, rhs;
  // VecDot(w_test, Av, &lhs);
  // VecDot(ATw, v_test, &rhs);
  // PetscScalar max_val = PetscMax(PetscAbsScalar(lhs), PetscAbsScalar(rhs));
  // PetscScalar diff = lhs - rhs;
  // PetscScalar diff_rel = (lhs - rhs) / max_val;
  // PetscPrintf(PETSC_COMM_SELF, "adjoint check: %g vs %g, diff: %g, diff_rel: %g\n", (double)PetscRealPart(lhs), (double)PetscRealPart(rhs), (double)PetscRealPart(diff), (double)PetscRealPart(diff_rel));
  // exit(1);


  // Solve the Least Squares system with either TAOBRGN or KSPLSQR
  if (leastsquares_solver == "brgn") {

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
    MatMult(GN_LeastSquares_MatShell, v_LeastSquares, brgn_residual);
    VecAXPY(brgn_residual, -1.0, brgn_rhs);
    VecNorm(brgn_residual, NORM_2, &rnorm);

    if (mpirank_world == 0 && !quietmode) {
      printf("Tao BRGN Least-Squares stats: newton iterations = %d, linear iterations = %d, residual norm = %1.14e, reason = %d\n", (int)newton_iters, (int)lin_iters, rnorm, reason);
    }

  } else if (leastsquares_solver == "lsqr"){
    // Set the matrix again, just in case, for reset.
    KSPSetOperators(ksp_LeastSquares, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell);

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
  } else {
    if (mpirank_world == 0 && !quietmode) {
      printf("Unknown Gauss-Newton least-squares solver: %s\n", leastsquares_solver.c_str());
      exit(1);
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
  MatMult(self->GN_LeastSquares_MatShell, v, F);  // F = L*v (use MatMult, not direct call)
  VecAXPY(F, -1.0, self->brgn_rhs);          // F = L*v - b

  PetscFunctionReturn(0);
}

// Callback: Evaluate Jacobian dF/dv = L (constant matrix)
PetscErrorCode OptimProblem::TaoBRGN_EvalJacobianResidual(Tao tao, Vec v, Mat J, Mat Jpre, void *ctx) {
  PetscFunctionBeginUser;

  // The Jacobian is L = GN_LeastSquares_MatShell, which is constant
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


PetscErrorCode OptimProblem::GN_LeastSquares_MatMult(Mat A, Vec v, Vec y)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  PetscInt M, N;
  MatGetSize(A, &M, &N);

  // Reset output
  VecZeroEntries(y);

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
  MPI_Allreduce(MPI_IN_PLACE, y_data, M, MPIU_SCALAR, MPI_SUM, self->comm_init);

  // Fill the bottom of y with damping: lambda * v
  int tangentdim = self->mastereq->getDim()*self->mastereq->getDim() - 1;
  const PetscScalar *v_data;
  VecGetArrayRead(v, &v_data);
  for (int i=0; i<M - tangentdim; i++) {
    y_data[i + tangentdim] = self->leastsquares_damping * v_data[i];
  }
  VecRestoreArrayRead(v, &v_data);
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

  if (self->includeHessUJ){
    // TODO: Need to implement W^1/2 in tangent space representation
    if (self->mpirank_world == 0) {
      printf("WARNING: W^1/2 is not implemented for Least Squares solver tangent space projection. Skipping.\n");
    }
  }

  return 0;
}

PetscErrorCode OptimProblem::GN_LeastSquares_MatMultTranspose(Mat A, Vec w_long, Vec vout)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(vout);

  int M,N;
  MatGetSize(A, &M, &N);

  // Create a tangent space vector w and fill it will the first part of w_long
  Vec w;
  int tangentdim = self->mastereq->getDim()*self->mastereq->getDim() - 1;
  VecCreateSeq(PETSC_COMM_SELF, tangentdim, &w);
  VecSetFromOptions(w);
  VecAssemblyBegin(w); VecAssemblyEnd(w);
  PetscScalar *w_ptr;
  const PetscScalar *wl_data;
  VecGetArrayRead(w_long, &wl_data);
  VecGetArray(w, &w_ptr);
  for (int i=0; i<tangentdim; i++) {
    w_ptr[i] = wl_data[i];
  }
  VecRestoreArray(w, &w_ptr);
  VecRestoreArrayRead(w_long, &wl_data);

  // Get read access to w (in reduced tangent space) Note: w has size tangent_dim, 
  const PetscScalar *w_data;
  VecGetArrayRead(w, &w_data);

  // Optional, apply W^1/2 to w before reconstruction
  if (self->includeHessUJ){
    if (self->mpirank_world == 0) {
      printf("WARNING: W^1/2 is not implemented for Least Squarestangent space projection. Skipping.\n");
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

  // Add damping lambda * w_long[tangentdim:end] contribution to the output
  VecGetArrayRead(w_long, &wl_data);
  for (int i=0; i<self->ndesign; i++) {
    vout_data[i] += self->leastsquares_damping * wl_data[i + tangentdim];
  }
  VecRestoreArray(vout, &vout_data);
  VecRestoreArrayRead(w_long, &wl_data);
  VecDestroy(&w);
  return 0;
}

PetscErrorCode OptimProblem::GN_LeastSquares_MatCreateVecs(Mat A, Vec *right, Vec *left)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  PetscInt M,N;
  MatGetSize(A, &M, &N);

  if (right) {
    // Create a control pulse vector
    VecCreateSeq(PETSC_COMM_SELF, N, right);
  }
  if (left) {
    // Create a reduced tangent space (augmented with ndesign rows for least-squares damping)
    VecCreateSeq(PETSC_COMM_SELF, M, left);
  }
  return 0;
}

PetscErrorCode OptimProblem::GN_Dual_MatCreateVecs(Mat A, Vec *right, Vec *left)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Always use pure tangent space dimension (never augmented)
  PetscInt tangent_dim = self->mastereq->getDim() * self->mastereq->getDim() - 1;
  PetscInt N = self->ndesign;

  if (right) {
    // Create a control pulse vector
    VecCreateSeq(PETSC_COMM_SELF, N, right);
  }
  if (left) {
    // Create a pure tangent space vector (not augmented)
    VecCreateSeq(PETSC_COMM_SELF, tangent_dim, left);
  }
  return 0;
}

void OptimProblem::GN_Jacobian_MatMult(Mat A, Vec v, Vec y)
{
  // Jacobian forward: y = P_G W^1/2 Lv 

  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(y);

  PetscInt tangent_dim;
  VecGetSize(y, &tangent_dim);

  // Apply linearized forward to get dU/dalpha * v
  self->evalLinearizedForward(self->xeval_GN, v);

  // If include W^1/2, precompute z = 1/n<V,W> 
  double z_re = 0.0;
  double z_im = 0.0;
  if (self->includeHessUJ){
    for (int iinit=0; iinit<self->ninit_local; iinit++) {
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      double tmp_re = 0.0;
      double tmp_im = 0.0;
      Vec wT = self->timestepper->getLinearizedFinalState(iinit);
      self->optim_target->HilbertSchmidtOverlap(wT, false, &tmp_re, &tmp_im);
      z_re +=  tmp_re / self->mastereq->getDim();
      z_im +=  tmp_im / self->mastereq->getDim();
    }
    // Gather result J_inf(W(T))
    double myz_re = z_re;
    double myz_im = z_im;
    MPI_Allreduce(&myz_re, &z_re, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
    MPI_Allreduce(&myz_im, &z_im, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
  }

  // Get U_final matrices
  Mat U_final_re = self->optim_target->getFinalUnitaryRe();
  Mat U_final_im = self->optim_target->getFinalUnitaryIm();

  // For each initial condition, apply W^1/2 and Project onto tangent space
  PetscScalar *y_data;
  VecGetArray(y, &y_data);
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;

    // w_i(T). Use rho_t0_bar as temporary storage
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);
    VecCopy(lin_final_state, self->rho_t0_bar);  

    // Include W^1/2 = sqrt(2/n)(I-z*V), if requested
    if (self->includeHessUJ){
      double scale = sqrt(2.0 / self->mastereq->getDim());
      VecScale(self->rho_t0_bar, scale);
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      self->optim_target->HilbertSchmidtOverlap_diff(self->rho_t0_bar, false, -scale*z_re, -scale*z_im);
    }

    // Project onto tangent space -> fills the appropriate entries in y_data
    GellMann::projectVecToTangentSpace(self->rho_t0_bar, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, y_data);
  }

  // Allreduce across comm_init processors
  MPI_Allreduce(MPI_IN_PLACE, y_data, tangent_dim, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(y, &y_data);

  // No damping rows - pure Jacobian only
}

void OptimProblem::GN_Jacobian_MatMultTranspose(Mat A, Vec w, Vec vout)
{
  // Jacobian transpose: vout = L^T W^1/2 P_G^{-1} w 

  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  // Reset output
  VecZeroEntries(vout);

  // Get read access to w (tangent space)
  const PetscScalar *w_data;
  VecGetArrayRead(w, &w_data);

  // P_G^{-1}: For each initial condition, reconstruct full state vectors from tangent space coefficients
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    Mat U_final_re = self->optim_target->getFinalUnitaryRe();
    Mat U_final_im = self->optim_target->getFinalUnitaryIm();

    // Reconstruct full state vector (adjoint of projectVecToTangentSpace) -> fills wsub_workspace[iinit]
    GellMann::reconstructVecFromTangentSpace(w_data, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, self->wsub_workspace[iinit]);
  }
  VecRestoreArrayRead(w, &w_data);

  // If include W^1/2, precompute z = 1/n<V,wsub_workspace> 
  double z_re = 0.0;
  double z_im = 0.0;
  if (self->includeHessUJ){
    for (int iinit=0; iinit<self->ninit_local; iinit++) {
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      double tmp_re = 0.0;
      double tmp_im = 0.0;
      self->optim_target->HilbertSchmidtOverlap(self->wsub_workspace[iinit], false, &tmp_re, &tmp_im);
      z_re +=  tmp_re / self->mastereq->getDim();
      z_im +=  tmp_im / self->mastereq->getDim();
    }
    // Gather result J_inf(W(T))
    double myz_re = z_re;
    double myz_im = z_im;
    MPI_Allreduce(&myz_re, &z_re, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
    MPI_Allreduce(&myz_im, &z_im, 1, MPI_DOUBLE, MPI_SUM, self->comm_init);
  }

  // Solve adjoint backward ODE for each terminal condition
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Include W = sqrt(2/n)(I-<V,cdot>V), if requested
    if (self->includeHessUJ){
      double scale = sqrt(2.0 / self->mastereq->getDim());
      VecScale(self->wsub_workspace[iinit], scale);
      int iinit_global = self->mpirank_init * self->ninit_local + iinit;
      self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
      self->optim_target->HilbertSchmidtOverlap_diff(self->wsub_workspace[iinit], false, -scale*z_re, -scale*z_im);
    }

    self->timestepper->solveAdjointODE(iinit, self->wsub_workspace[iinit], 0.0, 0.0, 0.0, 0.0);
    VecAXPY(vout, 1.0, self->timestepper->getReducedGradient());
  }

  // Sum up the gradient from all initial condition processors
  PetscScalar* vout_data;
  VecGetArray(vout, &vout_data);
  MPI_Allreduce(MPI_IN_PLACE, vout_data, self->ndesign, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(vout, &vout_data);

  // No damping contribution - pure Jacobian transpose only
}
