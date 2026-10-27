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

  /* Create TAO optimization solver. */
  // Each processor runs an independent TAO solver on the same replicated control vector. The control vector x is synchronized via the TaoSynchronize callback which is called at the start of each evaluation.
  TaoCreate(PETSC_COMM_SELF, &tao_lbfgs);
  TaoSetType(tao_lbfgs, TAOBQNLS);
  TaoSetMaximumIterations(tao_lbfgs, maxiter);
  TaoSetTolerances(tao_lbfgs, tol_grad_abs, PETSC_DEFAULT, tol_grad_rel);
  TaoSetUpdate(tao_lbfgs, TaoSynchronize, (void*)this);  
  TaoMonitorSet(tao_lbfgs, TaoMonitor, (void*)this, NULL);
  TaoSetVariableBounds(tao_lbfgs, xlower, xupper);
  TaoSetFromOptions(tao_lbfgs);
  TaoSetObjectiveAndGradient(tao_lbfgs, NULL, TaoEvalObjectiveAndGradient, (void*) this);

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

  /* Initialize the Gauss-Newton solvers */
  VecDuplicate(xinit, &xeval_GN);
  VecZeroEntries(xeval_GN);
  VecAssemblyBegin(xeval_GN); VecAssemblyEnd(xeval_GN);

  // Get for the inner Gauss-Newton solvers from config 
  gn_maxiter = config.getGnMaxiter();
  gn_rtol = config.getGnRtol();
  gn_solver_mode = config.getGnSolverMode();
  gn_normaleq_damping = config.getGnNormaleqDamping();
  gn_normaleq_solver = config.getGnNormaleqSolver();  
  gn_leastsquares_damping = config.getGnLeastSquaresDamping();
  gn_leastsquares_solver = config.getGnLeastSquaresSolver();
  gn_prec_update_interval = config.getGnPrecUpdateInterval();
  gn_precond_rank = config.getGnPrecondRank();

  gn_precond_enabled = (gn_prec_update_interval >= -1);
  tangent_dim = mastereq->getDim()*mastereq->getDim()-1;

  // Create MatShell for least-squares problem
  int N = ndesign;
  PetscInt M = tangent_dim;
  if (gn_leastsquares_damping > 0.0) {
    M += ndesign;
  }
  MatCreateShell(PETSC_COMM_SELF, M, N, M, N, this, &GN_LeastSquares_MatShell);
  MatShellSetOperation(GN_LeastSquares_MatShell, MATOP_MULT, (void(*)(void))GN_LeastSquares_MatMult);
  MatShellSetOperation(GN_LeastSquares_MatShell, MATOP_MULT_TRANSPOSE, (void(*)(void))GN_LeastSquares_MatMultTranspose);
  // Pre-allocate workspace vectors for GN_LeastSquares_MatShell operations
  wsub_workspace = new Vec[ninit_local];
  for (int iinit = 0; iinit < ninit_local; iinit++) {
    VecDuplicate(rho_t0_bar, &wsub_workspace[iinit]);
  }

  // Create the LSQR and BRGN solvers for least-squares problem
  KSPCreate(PETSC_COMM_SELF, &ksp_gn_leastsquares);
  KSPSetOperators(ksp_gn_leastsquares, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell);
  KSPSetType(ksp_gn_leastsquares, KSPLSQR);
  KSPSetTolerances(ksp_gn_leastsquares,gn_rtol,PETSC_DEFAULT,PETSC_DEFAULT,gn_maxiter);
  // KSPSetComputeSingularValues(ksp_gn_leastsquares, PETSC_TRUE);
  KSPSetFromOptions(ksp_gn_leastsquares);
  KSPMonitorCancel(ksp_gn_leastsquares);
  KSPMonitorSet(ksp_gn_leastsquares, KSPMonitorResidualAndSolution, (void*)this, NULL);
  TaoCreate(PETSC_COMM_SELF, &tao_gn_leastsquares_brgn);
  TaoSetType(tao_gn_leastsquares_brgn, TAOBRGN);
  MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &brgn_rhs);
  VecDuplicate(brgn_rhs, &brgn_residual);
  // Set BRGN convergence tolerances. TODO: THIS NEEDS FIX! BRGN is a newton solver, so it will only need one outer newton step. The inner linear solver is the actual workhorse. Need to pass appropriate tolerances to the subsolver, not the BRGN solver. However, somehow, the below crashes when running on multiple processes. TODO: Check parallel implementation.
  // Tao tao_sub;
  // TaoBRGNGetSubsolver(tao_gn_leastsquares_brgn, &tao_sub);
  // TaoSetTolerances(tao_sub, config.getGnRtol(), PETSC_DEFAULT, PETSC_DEFAULT);
  // TaoSetMaximumIterations(tao_sub, config.getGnMaxiter());
  TaoSetTolerances(tao_gn_leastsquares_brgn, gn_rtol, PETSC_DEFAULT, PETSC_DEFAULT);
  TaoSetMaximumIterations(tao_gn_leastsquares_brgn, gn_maxiter);
  TaoBRGNSetRegularizerWeight(tao_gn_leastsquares_brgn, gn_leastsquares_damping);
  TaoSetResidualRoutine(tao_gn_leastsquares_brgn, brgn_residual, TaoBRGN_EvalResidual, (void*)this);
  TaoSetJacobianResidualRoutine(tao_gn_leastsquares_brgn, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell, TaoBRGN_EvalJacobianResidual, (void*)this);
  TaoMonitorSet(tao_gn_leastsquares_brgn, TaoBRGN_Monitor, (void*)this, NULL);
  TaoSetFromOptions(tao_gn_leastsquares_brgn);

  // Create MatShell for Gauss-Newton primal equation
  MatCreateShell(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, ndesign, ndesign, this, &GN_Primal_MatShell);
  MatShellSetOperation(GN_Primal_MatShell, MATOP_MULT, (void(*) (void)) GN_Primal_MatMult);
  MatCreateDense(PETSC_COMM_SELF, ndesign, ndesign, ndesign, ndesign, NULL, &GN_Primal_MatDense); // only used for debug and output 
  MatZeroEntries(GN_Primal_MatDense);

  // Create MatShell for dual Gauss-Newton problem
  PetscInt M_dual = tangent_dim;
  MatCreateShell(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, this, &GN_Dual_MatShell);
  MatShellSetOperation(GN_Dual_MatShell, MATOP_MULT, (void(*)(void))GN_Dual_MatMult);
  MatCreateDense(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, NULL, &GN_Dual_MatDense); // only debug and output
  MatZeroEntries(GN_Dual_MatDense);

  // Create MatShell for dual and primal low-rank preconditioned 
  MatCreateShell(PETSC_COMM_SELF, M_dual, M_dual, M_dual, M_dual, this, &GN_Dual_LowRank_MatShell);
  MatShellSetOperation(GN_Dual_LowRank_MatShell, MATOP_MULT, (void(*)(void))GN_Dual_LowRank_MatMult);
  MatCreateShell(PETSC_COMM_SELF, ndesign, ndesign, ndesign, ndesign, this, &GN_Primal_LowRank_MatShell);
  MatShellSetOperation(GN_Primal_LowRank_MatShell, MATOP_MULT, (void(*)(void))GN_Primal_LowRank_MatMult);

  // Create Gauss-Newton normal equation linear solver (used for both the primal and the dual equations) 
  KSPCreate(PETSC_COMM_SELF, &ksp_gn_normaleq);
  KSPSetInitialGuessNonzero(ksp_gn_normaleq, PETSC_TRUE);
  if (gn_normaleq_solver == "cg") {
    KSPSetType(ksp_gn_normaleq, KSPCG);
  } else if (gn_normaleq_solver == "minres") {
    KSPSetType(ksp_gn_normaleq, KSPMINRES);
    // Apply MINRES-QLP if configured
    if (config.getGnNormaleqMinresQlp()) {
      PetscBool set;
      PetscOptionsHasName(NULL, NULL, "-ksp_minres_qlp", &set);
      if (!set) {
        PetscOptionsSetValue(NULL, "-ksp_minres_qlp", NULL);
      }
    }
  } else if (gn_normaleq_solver == "gmres") {
    KSPSetType(ksp_gn_normaleq, KSPGMRES);
  } 
  KSPSetFromOptions(ksp_gn_normaleq);
  // KSPSetNormType(ksp_gn_normaleq, KSP_NORM_UNPRECONDITIONED); // Uncomment if unpreconditioned norm is desired, but is slower. 
  KSPSetTolerances(ksp_gn_normaleq, gn_rtol, PETSC_DEFAULT, PETSC_DEFAULT, gn_maxiter);
  KSPMonitorCancel(ksp_gn_normaleq);
  KSPMonitorSet(ksp_gn_normaleq, KSPMonitorResidualAndSolution, (void*)this, NULL);

  // Pre-allocate workspace for dual solver and reset
  VecCreateSeq(PETSC_COMM_SELF, M_dual, &y_dual_workspace);  // tangent space
  VecDuplicate(xinit, &v_workspace);  // design space

  // preallocate preconditioner workspace
  int M_precond = 0;
  if (gn_solver_mode == "dual") {
    M_precond = M_dual;
  } else if (gn_solver_mode == "primal") {
    M_precond = ndesign;
  }
  MatCreateDense(PETSC_COMM_SELF, M_precond, gn_precond_rank, M_precond, gn_precond_rank, NULL, &gn_dual_lowrank_U);
  VecCreateSeq(PETSC_COMM_SELF, gn_precond_rank, &gn_dual_lowrank_sigma_invsqrt);
  VecCreateSeq(PETSC_COMM_SELF, gn_precond_rank, &gn_dual_lowrank_temp);


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

  MatDestroy(&GN_Primal_MatDense);
  MatDestroy(&GN_Primal_MatShell);
  MatDestroy(&GN_Primal_MatDense);
  MatDestroy(&GN_LeastSquares_MatShell);
  MatDestroy(&GN_Dual_MatShell);
  MatDestroy(&GN_Dual_LowRank_MatShell);
  MatDestroy(&GN_Primal_LowRank_MatShell);
  MatDestroy(&GN_Dual_MatDense);
  VecDestroy(&xeval_GN);
  VecDestroy(&y_dual_workspace);
  VecDestroy(&v_workspace);
  MatDestroy(&gn_dual_lowrank_U);
  VecDestroy(&gn_dual_lowrank_sigma_invsqrt);
  VecDestroy(&gn_dual_lowrank_temp);
  KSPDestroy(&ksp_gn_normaleq);
  KSPDestroy(&ksp_gn_leastsquares);
  if (tao_gn_leastsquares_brgn) {
    TaoDestroy(&tao_gn_leastsquares_brgn);
    VecDestroy(&brgn_rhs);
    VecDestroy(&brgn_residual);
  }
  TaoDestroy(&tao_lbfgs);

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

  // 1. Set for for solveGaussNewtonPrimalKSP KSP solver
  PetscReal rtol, abstol, dtol;
  PetscInt current_maxits;

  // For solveGaussNewtonLeastSquares solvers
  if (gn_leastsquares_solver == "brgn") {
    // For BRGN, this should apply to the sub-solver, since BRGN does Newton steps (so only one anyways). This crashes though. TODO: Fix. 
    // Tao tao_sub;
    // TaoBRGNGetSubsolver(tao_gn_leastsquares_brgn, &tao_sub);
    // TaoSetMaximumIterations(tao_sub, maxiter);
    TaoSetMaximumIterations(tao_gn_leastsquares_brgn, maxiter); // no effect
  } else {
    KSPGetTolerances(ksp_gn_leastsquares, &rtol, &abstol, &dtol, &current_maxits);
    KSPSetTolerances(ksp_gn_leastsquares, rtol, abstol, dtol, maxiter);
  }

  // For Gauss-Newton normal equation KSP solvers
  KSPGetTolerances(ksp_gn_normaleq, &rtol, &abstol, &dtol, &current_maxits);
  KSPSetTolerances(ksp_gn_normaleq, rtol, abstol, dtol, maxiter);
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

    // Solve the ODE only once per new xeval_GN; reuse the stored trajectory_states otherwise.
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

PetscErrorCode OptimProblem::GN_Primal_MatMult(Mat A, const Vec v, Vec Av){
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

  // Solve the adjoint ODE for initial condition 
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Set terminal condition for adjoint
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);

    // Set the terminal condition 
    VecCopy(lin_final_state, self->rho_t0_bar);  // w_i(T)

    // Include W = 2/n(I-<V,cdot>V)
    double scale = 2.0 / self->mastereq->getDim();
    VecScale(self->rho_t0_bar, scale);
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
    self->optim_target->HilbertSchmidtOverlap_diff(self->rho_t0_bar, false, -scale*z_re, -scale*z_im);

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

  // v_temp = L^T W^1/2 P_G y (tangent space -> design space)
  self->apply_GN_Jacobian_Tangentspace_Adjoint(y, self->v_workspace);
  // Ay = P_G W^1/2 L v_temp (design space -> tangent space)
  self->apply_GN_Jacobian_Tangentspace(self->v_workspace, Ay);

  return 0;
}

Vec OptimProblem::computeGaussNewtonRHS(const std::string& gn_solver_mode, const Vec G){
  assert(gn_solver_mode == "dual" || gn_solver_mode == "primal" || gn_solver_mode == "least_squares");

  Vec b = nullptr;

  if (gn_solver_mode == "primal") {
    // For primal solver, use the negative reduced gradient (\nabla_alpha J)

    MatCreateVecs(GN_Primal_MatShell, nullptr, &b);
    VecCopy(G, b);
    VecScale(b, -1.0);
    
  } else if (gn_solver_mode == "dual" || gn_solver_mode == "least_squares") {
    // For dual and least-squares solver: - W^1/2 U = -sqrt(n/2)\nabla_U J projected to tangent space 

    // Create the RHS vector of appropriate size 
    if (gn_solver_mode == "dual") {
      MatCreateVecs(GN_Dual_MatShell, nullptr, &b);
      VecZeroEntries(b);
    } else  { // least-squares
      MatCreateVecs(GN_LeastSquares_MatShell, nullptr, &b);
      VecZeroEntries(b);
    }

    // Precompute z = 1/n tr(V^dagger U), needed for ∇_U J
    double z_re = 0.0;
    double z_im = 0.0;
    for (int iinit = 0; iinit < ninit_local; iinit++) {
      int iinit_global = mpirank_init * ninit_local + iinit;
      optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
      double obj_lin_iinit_re = 0.0;
      double obj_lin_iinit_im = 0.0;
      Vec ui = timestepper->getFinalState(iinit);
      optim_target->evalJ(ui, &obj_lin_iinit_re, &obj_lin_iinit_im);
      z_re += obj_weights[iinit] * obj_lin_iinit_re;
      z_im += obj_weights[iinit] * obj_lin_iinit_im;
    }
    double myz_re = z_re;
    double myz_im = z_im;
    MPI_Allreduce(&myz_re, &z_re, 1, MPI_DOUBLE, MPI_SUM, comm_init);
    MPI_Allreduce(&myz_im, &z_im, 1, MPI_DOUBLE, MPI_SUM, comm_init);

    // Compute -sqrt(n/2)∇_U J and project to tangent space
    PetscScalar *b_data;
    VecGetArray(b, &b_data);
    for (int iinit = 0; iinit < ninit_local; iinit++) {
      int iinit_global = mpirank_init * ninit_local + iinit;
      Vec ui = timestepper->getFinalState(iinit);

      // Use pre-allocated workspace vector
      VecCopy(ui,  wsub_workspace[iinit]);

      // compute wsub = alpha * (ui - z * V), where z = 1/n tr(V^dagger U)
      double scale = 2.0 / mastereq->getDim(); // dual: alpha = 2/n, no W^1/2!
      if (gn_solver_mode == "least_squares") { // for least-squares: alpha = sqrt(2/n) includes W^1/2 factor
        scale = sqrt(2.0 / mastereq->getDim());
      }
      VecScale(wsub_workspace[iinit], scale);
      optim_target->prepareInitialAndTargetState(iinit_global, ninit, mastereq->nlevels, mastereq->nessential);
      optim_target->HilbertSchmidtOverlap_diff(wsub_workspace[iinit], false, -scale*z_re, -scale*z_im);

      // Project onto tangent space 
      Mat U_final_re = optim_target->getFinalUnitaryRe();
      Mat U_final_im = optim_target->getFinalUnitaryIm();
      GellMann::projectVecToTangentSpace(wsub_workspace[iinit], U_final_re, U_final_im, mastereq->getDim(), iinit_global, b_data);
    }

    // Allreduce b across all comm_init processors
    MPI_Allreduce(MPI_IN_PLACE, b_data, tangent_dim, MPIU_SCALAR, MPI_SUM, comm_init);
    VecRestoreArray(b, &b_data);

    VecScale(b, -1.0);

  } else {
    printf("Wrong GN solver mode in computeGaussNewtonRHS. This should never happen!\n");
    exit(1);
  }
  return b;
}

Vec OptimProblem::computeGaussNewtonInitialGuess(const std::string& gn_solver_mode, const Vec b_rhs){
  assert(gn_solver_mode == "dual" || gn_solver_mode == "primal" || gn_solver_mode == "least_squares");

  Vec y_initial = nullptr;

  // Create vector of appropriate size and set default zero initial guess
  if (gn_solver_mode == "least_squares") {
    MatCreateVecs(GN_LeastSquares_MatShell, &y_initial, nullptr);
    VecZeroEntries(y_initial);
  } else if (gn_solver_mode == "dual") {
    MatCreateVecs(GN_Dual_MatShell, &y_initial, nullptr);
    VecZeroEntries(y_initial);
  } else if (gn_solver_mode == "primal") {
    MatCreateVecs(GN_Primal_MatShell, &y_initial, nullptr);
    VecZeroEntries(y_initial);
  }

  // If preconditioner enabled, use M^{-1/2} b as initial guess.
  if (gn_precond_enabled) {
    apply_GN_LowRank_Precond(b_rhs, y_initial);
  }

  return y_initial;
}

bool OptimProblem::solveGaussNewtonPrecondKSP(const std::string& gn_solver_mode, Vec xinit, const Vec initial_guess, const Vec b_rhs, Vec v_out){

  // Store the point of evaluation for the Gauss-Newton matrix shell
  setGaussNewtonEvalPoint(xinit); 
  // Forces a fresh nonlinear forward solve during the next call to evalLinearizedForward() for the new xeval_GN
  // TODO: This can probably be dropped here, since it should be set before any solveGaussNewtonXYZ calls. 

  Mat SystemMatrix = (gn_solver_mode == "dual") ? GN_Dual_MatShell : GN_Primal_MatShell;
  Mat Preconditioner = (gn_solver_mode == "dual") ? GN_Dual_LowRank_MatShell : GN_Primal_LowRank_MatShell;


  // Optional Levenberg-Marquardt damping: (L L^T + λI)
  if (gn_normaleq_damping > 0.0) {
    MatShift(SystemMatrix, gn_normaleq_damping);
  }

  Vec vec_workspace = nullptr;
  if (gn_solver_mode == "dual") {
    vec_workspace = y_dual_workspace;
  } else {
    vec_workspace = v_workspace;
  }

  if (!gn_precond_enabled) {
    // No preconditioning - solve directly
    VecCopy(initial_guess, vec_workspace);
    KSPReset(ksp_gn_normaleq);
    KSPSetOperators(ksp_gn_normaleq, SystemMatrix, SystemMatrix);
    KSPSolve(ksp_gn_normaleq, b_rhs, vec_workspace);

  } else { 
    // Low-rank preconditioning: M^{-1/2} A M^{-1/2} z = M^{-1/2} b  where M ≈ U Σ U^T is rank-r approximation

    // Scale RHS: b_scaled = M^{-1/2} * b
    Vec b_scaled;
    VecDuplicate(b_rhs, &b_scaled);
    apply_GN_LowRank_Precond(b_rhs, b_scaled);

    // Copy initial guess into solution vector
    Vec z_solution;
    VecDuplicate(initial_guess, &z_solution);
    VecCopy(initial_guess, z_solution);  

    // Solve: M^{-1/2} A M^{-1/2} z = M^{-1/2} b
    KSPReset(ksp_gn_normaleq);
    KSPSetOperators(ksp_gn_normaleq, Preconditioner, Preconditioner);
    KSPSolve(ksp_gn_normaleq, b_scaled, z_solution);

    // Recover y = M^{-1/2} * z
    apply_GN_LowRank_Precond(z_solution, vec_workspace);

    VecDestroy(&b_scaled);
    VecDestroy(&z_solution);
  } 

  // Revert the damping
  if (gn_normaleq_damping > 0.0) {
    MatShift(SystemMatrix, -gn_normaleq_damping);
  }

  // Recover primal solution: v = L^* y 
  if (gn_solver_mode == "dual") {
    apply_GN_Jacobian_Tangentspace_Adjoint(y_dual_workspace, v_out);
  } else {
    VecCopy(vec_workspace, v_out);
  }

  // Report convergence
  KSPConvergedReason reason;
  int iters;
  double rnorm;
  KSPGetConvergedReason(ksp_gn_normaleq, &reason);
  KSPGetIterationNumber(ksp_gn_normaleq, &iters);
  KSPGetResidualNorm(ksp_gn_normaleq, &rnorm);
  ksp_iters_last = iters;
  if (mpirank_world == 0 && !quietmode) {
    KSPType ksp_type;
    KSPGetType(ksp_gn_normaleq, &ksp_type);
    printf("Gauss-Newton Dual KSP (%s) stats: iterations = %d, residual norm = %1.14e\n", ksp_type, iters, rnorm);

    // Print convergence reason for debugging
    const char* reason_str = KSPConvergedReasons[reason];
    printf("  KSP Convergence Reason: %s (code %d)\n", reason_str, (int)reason);
  }

  // Detect stale eigendecomposition and return flag to caller
  bool eigendecomp_is_stale = false;
  std::string stale_reason = "";

  // Check 1: KSP diverged
  if (reason < 0 && gn_precond_enabled) {
    eigendecomp_is_stale = true;
    stale_reason = "KSP diverged";
  }

  // Check 2: Large initial residual (only if we used low-rank preconditioner)
  if (!eigendecomp_is_stale && gn_precond_enabled) {
    // Check if initial guess led to huge residual
    // For well-conditioned eigenspace approximation, initial residual should be moderate
    // If r_0 is very large compared to ||b||, eigendecomposition is likely stale
   // Using Heuristic: if initial residual > 100 * ||b||, eigendecomposition is suspect
    double b_norm;
    VecNorm(b_rhs, NORM_2, &b_norm);
    if (iters == 0 && rnorm > 100.0 * b_norm) {
      eigendecomp_is_stale = true;
      stale_reason = "Initial residual too large";
    }
  }

  return eigendecomp_is_stale;
}


void OptimProblem::solveGaussNewtonDualProjected(Vec xinit, const Vec b_rhs, Vec v_out){
  // TODO: Understand this function and how it interacts with the low-rank preconditioner and the Gauss-Newton dual system.
  // TODO: Remove the KSP iterations. Make this a full-rank rank solver in the projected subspace.

  // Store the point of evaluation. This should have been done already! 
  setGaussNewtonEvalPoint(xinit); 
  // Forces a fresh nonlinear forward solve during the next call to evalLinearizedForward() for the new xeval_GN
  // TODO: This can probably be dropped here, since it should be set before any solveGaussNewtonXYZ calls. 

  // Use passed RHS
  const Vec b = b_rhs;

  // Apply damping to the operator
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, gn_normaleq_damping);
  }

  // Compute eigendecomposition: (L L^T + λI) V_{r+1} = Σ_{r+1} V_{r+1}
  // We compute r+1 eigenvalues to get μ = λ_{r+1} for the complement space
  int rank = gn_precond_rank;
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
  applyEigenspaceApproximation(b, V_r, eigenvalues, rank, have_complement, mu, y_dual_workspace);

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

      // Apply low-rank preconditioning: M^{-1/2} A M^{-1/2} z = M^{-1/2} b_normalized
      Vec b_scaled, z_solution;
      VecDuplicate(b_normalized, &b_scaled);
      VecDuplicate(y_proj, &z_solution);

      // Scale normalized RHS: b_scaled = M^{-1/2} * b_normalized
      apply_GN_LowRank_Precond(b_normalized, b_scaled);

      // Scale initial guess
      VecCopy(y_proj, z_solution);

      // Solve preconditioned system: M^{-1/2} A M^{-1/2} z = M^{-1/2} b_normalized
      KSPReset(ksp_gn_normaleq);
      KSPSetOperators(ksp_gn_normaleq, GN_Dual_LowRank_MatShell, GN_Dual_LowRank_MatShell);
      KSPSolve(ksp_gn_normaleq, b_scaled, z_solution);

      // Recover y = M^{-1/2} * z
      apply_GN_LowRank_Precond(z_solution, y_dual_workspace);

      // Scale back by ||b|| to get the actual solution
      VecScale(y_dual_workspace, bnorm);

      // Get actual number of iterations performed
      KSPGetIterationNumber(ksp_gn_normaleq, &refine_iters);

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
      KSPSetOperators(ksp_gn_normaleq, GN_Dual_MatShell, GN_Dual_MatShell);
    }
  }

  // Revert damping
  if (gn_normaleq_damping > 0.0) {
    MatShift(GN_Dual_MatShell, -gn_normaleq_damping);
  }

  // Cleanup (b_rhs is destroyed by caller)
  MatDestroy(&V_r);
  delete[] eigenvalues;

  // Recover primal solution: v = L^T y (use pure Jacobian transpose for dual mode)
  apply_GN_Jacobian_Tangentspace_Adjoint(y_dual_workspace, v_out);

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

bool OptimProblem::solveGaussNewtonSystem(const Vec x_GN, const Vec initial_guess, const Vec b, Vec Ainv_b){

  // Default return: Preconditioner is NOT state. Currently only "dual" has a precondioner that could change this return.
  bool precond_is_stale = false; 

  if (gn_solver_mode == "primal") {
    // Primal normal equations: (L^*WL + λI) v = - \nabla_alpha J
    if (gn_normaleq_solver == "proj") {
      solveGaussNewtonPrimalProjected(x_GN, b, Ainv_b);
    } else {
      precond_is_stale = solveGaussNewtonPrecondKSP("primal", x_GN, initial_guess, b, Ainv_b);
    }
  } 
  else if (gn_solver_mode == "dual") {
    // Dual normal equation: (L L^* + λI) y = -\nabla_U J and v = L^* y
    if (gn_normaleq_solver == "proj") {
      solveGaussNewtonDualProjected(x_GN, b, Ainv_b);
    } else {
      precond_is_stale = solveGaussNewtonPrecondKSP("dual", x_GN, initial_guess, b, Ainv_b);
    } 
  } 
  else if (gn_solver_mode == "least_squares") {
    // Least squares problem min ||L v - G||^2
    solveGaussNewtonLeastSquares(x_GN, initial_guess, b, Ainv_b);
  } 
  else {
    if (mpirank_world == 0 && !quietmode) {
      printf("Unknown Gauss-Newton solver mode: %s. This should never happen.\n", gn_solver_mode.c_str());
      exit(1);
    }
  }

  return precond_is_stale;
}

void OptimProblem::solveGaussNewtonPrimalProjected(Vec xinit, const Vec b, Vec Ainv_b){

  // Get eigenvalues and eigenvectors of the Gauss-Newton matrix
  Mat evecs;
  std::vector<double> evals = computeGaussNewtonPrimalEvals(xinit, &evecs);

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
    if (evals[i] > 1e-14) {
      tmp_data[i] /= (evals[i] + gn_normaleq_damping);
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

std::vector<double> OptimProblem::computeGaussNewtonPrimalEvals(Vec xinit, Mat* evecs_out){

  // Store xinit so the MatShell can use it as point of evaluation.
  setGaussNewtonEvalPoint(xinit); 
  // Forces a fresh nonlinear forward solve during the next call to evalLinearizedForward() for the new xeval_GN
  // TODO: This can probably be dropped here, since it should be set before any solveGaussNewtonXYZ calls. 


  /* Create eigenvalues solver for Gauss-Newton primal equation Ax=b */
  EPS eps_gn_primal; // EPS solver for the Gauss-Newton primal equation 
  EPSCreate(PETSC_COMM_SELF, &eps_gn_primal);
  EPSSetOperators(eps_gn_primal, GN_Primal_MatShell, NULL);
  EPSSetProblemType(eps_gn_primal, EPS_HEP); // Hermitian
  EPSSetWhichEigenpairs(eps_gn_primal, EPS_LARGEST_REAL); // largest eigenvalues
  int neigvals = tangent_dim;
  int ncv = 2*neigvals ; // Max Krylov dimension. How to set??
  EPSSetDimensions(eps_gn_primal, neigvals, ncv, PETSC_DEFAULT);
  EPSSetTolerances(eps_gn_primal, gn_rtol, gn_maxiter);
  EPSSetFromOptions(eps_gn_primal);

  // Update EPS for a fresh solve on the new xeval_GN. 
  if (!gn_primal_usedensemat) { 
    // use MatShell 
    EPSSetOperators(eps_gn_primal, GN_Primal_MatShell, NULL);
  } else {
    // use dense Matrix representation of the Gauss-Newton matrix
    updateGaussNewtonPrimalMatDense(); // does N2-1 applications of the MatShell
    EPSSetOperators(eps_gn_primal, GN_Primal_MatDense, NULL);
  }

  // Solve the eigenvalue problem for the Gauss-Newton matrix
  EPSSolve(eps_gn_primal);
  PetscInt numConv;
  PetscInt iters_taken;
  EPSGetConverged(eps_gn_primal, &numConv);
  EPSGetIterationNumber(eps_gn_primal,&iters_taken);
  if (mpirank_world == 0 && !quietmode) printf("Gauss-Newton EPS converged %d eigenvalues in %d iterations.\n", numConv, iters_taken);
  if (numConv < neigvals) {
      if (mpirank_world==0 && !quietmode) printf("WARNING: Only %d eigenvalues out of %d eigenvalues converged.\n", numConv, neigvals);
  }

  // Set up storage for eigenvalues and eigenvectors (should be real!)
  std::vector<double> evals_re(neigvals);
  std::vector<Vec> evec_re(neigvals);
  for (int ix = 0; ix < neigvals; ix++) {
    MatCreateVecs(GN_Primal_MatShell, &evec_re[ix], NULL);
  }

  // Retrieve eigenpairs of M and compute error. 
  for (PetscInt i = 0; i < numConv && i < neigvals; i++) {

    // Retrieve the eigenvalue (is real) and eigenvector
    EPSGetEigenpair(eps_gn_primal, i, &evals_re[i], NULL, evec_re[i], NULL);

    // Estimate the errror (needs one more application of A)
    double error = 0.0;
    EPSComputeError(eps_gn_primal,i,EPS_ERROR_RELATIVE,&error);
    if (error > gn_rtol) {
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
  EPSDestroy(&eps_gn_primal);

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
      TaoSetSolution(tao_lbfgs, xinit);
      TaoSolve(tao_lbfgs);
      break;

    case OptimSolverType::GAUSS_NEWTON: {
      // Work vectors for gradient, preconditioned search direction, and line search
      Vec G, Gprec, xnew, step;
      VecDuplicate(x_GN, &G);
      VecDuplicate(x_GN, &Gprec);
      VecDuplicate(x_GN, &xnew);
      VecDuplicate(x_GN, &step);

      // Reset KSP counters. TODO: Missing some? 
      gn_dual_ksp_baseline = 0;  
      gn_dual_ksp_excess_accum = 0;  
      ksp_iters_last = 0;
      nonlinear_forward_valid = false;

      // Set the initial iterate for Gauss-Newton optimization
      VecCopy(xinit, x_GN);

      bool stop = false;
      for (int iter = 0; !stop; iter++) {
        // Set the current evaluation point for the Gauss-Newton matrix
        setGaussNewtonEvalPoint(x_GN);
        
        // Compute gradient and objective at the current iterate
        evalGradF(x_GN, G);
        double f = objective;

        // Update the preconditioner. Currently only for the dual mode.
        bool did_update_PC = updateGaussNewtonPreconditioner(iter, false);

        // Set RHS and initial guess for the Gauss-Newton system
        Vec b_rhs = computeGaussNewtonRHS(gn_solver_mode, G); 
        Vec x_initial = computeGaussNewtonInitialGuess(gn_solver_mode, b_rhs);

        // Solve for the Gauss-Newton system for the search direction (Gprec) 
        bool precond_is_stale = solveGaussNewtonSystem(x_GN, x_initial, b_rhs, Gprec);
          
        // If preconditioner is stale, force immediate recomputation and retry
        if (precond_is_stale) {  
          if (mpirank_world == 0 && !quietmode) {
            printf("  Recomputing the preconditioner at the current iterate and forcing new GN solve...\n");
          }
          did_update_PC = updateGaussNewtonPreconditioner(iter, true); // Force the update
          VecDestroy(&x_initial);
          x_initial = computeGaussNewtonInitialGuess(gn_solver_mode, b_rhs); 
          precond_is_stale = solveGaussNewtonSystem(x_GN, x_initial, b_rhs, Gprec);
          if (precond_is_stale && mpirank_world == 0 && !quietmode) {
            printf("  WARNING: Solve still problematic after eigendecomposition update\n");
          }
        }
        VecDestroy(&x_initial);
        VecDestroy(&b_rhs);

        // Update baseline KSP iterations 
        if (iter == 0 || did_update_PC) {
          gn_dual_ksp_baseline = ksp_iters_last;
        }
        
        // Backtracking Armijo line search along -Gprec, projected onto the bound constraints
        VecScale(Gprec, -1.0);
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
    TaoGetSolution(tao_lbfgs, &params);
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
    // Compute norm of control vector, normalized by sqrt(D)
    double xnorm;
    if (optim_solver_type == OptimSolverType::TAO_LBFGS) {
      Vec params;
      TaoGetSolution(tao_lbfgs, &params);
      VecNorm(params, NORM_2, &xnorm);
    } else {
      VecNorm(x_GN, NORM_2, &xnorm);
    }
    // Normalize by sqrt(D) where D is the dimension
    xnorm /= std::sqrt(static_cast<double>(ndesign));

    // Add to optimization history file
    getOutput()->writeOptimFile(iter, f, gnorm, deltax, F_avg, obj_cost, obj_regul, obj_penal_leakage, obj_penal_dpdm, obj_penal_energy, obj_penal_variation, obj_penal_weightedcost, ksp_iters_last, xnorm);
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

PetscErrorCode TaoSynchronize(Tao tao, PetscInt iter, void*ptr){
  OptimProblem* ctx = (OptimProblem*) ptr;

  /* Get solution vector */
  Vec x;
  TaoGetSolution(tao, &x);

  /* Synchronize x across all processors: Average solution from each rank. */
  PetscScalar* xdata;
  VecGetArray(x, &xdata);
  MPI_Allreduce(MPI_IN_PLACE, xdata, ctx->getNdesign(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  VecRestoreArray(x, &xdata);
  VecScale(x, 1.0 / ctx->getMPIsize_world());
  return 0;
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
  OptimProblem* ctx = (OptimProblem*) ptr;

  ctx->evalGradF(x, G, false);
  *f = ctx->getObjective();

  return 0;
}

void OptimProblem::updateGaussNewtonPrimalMatDense(){
  MatZeroEntries(GN_Primal_MatDense);

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
    MatMult(this->getGN_Primal_MatShell(), e, Av);
    for (int jx = 0; jx < ndesign; ++jx) {
      PetscScalar val;
      VecGetValues(Av, 1, &jx, &val);
      MatSetValue(GN_Primal_MatDense, jx, ix, val, INSERT_VALUES);
    }
  }
  VecDestroy(&e);
  VecDestroy(&Av);

  MatAssemblyBegin(GN_Primal_MatDense, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(GN_Primal_MatDense, MAT_FINAL_ASSEMBLY);

  // Allreduce to combine contributions from all MPI ranks
  PetscScalar *data;
  MatDenseGetArray(GN_Primal_MatDense, &data);
  int size = ndesign * ndesign;
  MPI_Allreduce(MPI_IN_PLACE, data, size, MPIU_SCALAR, MPI_SUM, comm_optim);
  MatDenseRestoreArray(GN_Primal_MatDense, &data);

}


void OptimProblem::updateGaussNewtonDualMatDense(){
  MatZeroEntries(GN_Dual_MatDense);

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


bool OptimProblem::updateGaussNewtonPreconditioner(int iter, bool force_update){

  // preconditioner only available for dual or primal solver modes
  if (gn_solver_mode != "dual" && gn_solver_mode != "primal") {
    return false; 
  }

  // Don't update if preconditioner is disabled
  if (!gn_precond_enabled) {
    return false; 
  }

  // Don't update if this is a direct projected eigensolve 
  if (gn_normaleq_solver == "proj") {
    return false;
  }

  // Always update the preconditioner at the first iteration, or if forced
  if (iter == 0 || force_update){
    updateGaussNewtonPreconditionerLowRank(gn_precond_rank);
    // Reset accumulator after adaptive update
    if (gn_prec_update_interval == -1) {
      gn_dual_ksp_excess_accum = 0;
    }
    return true;
  }

  // At larger iterations, decide whether to update the preconditioner
  bool should_update = false;
  if (gn_prec_update_interval > 0 && iter % gn_prec_update_interval == 0) {
    // Fixed interval update
    should_update = true;
    if (mpirank_world == 0 && !quietmode) {
      printf("Iteration %d: Updating preconditioner (fixed interval = %d)\n", iter, gn_prec_update_interval);
    }
  } else if (gn_prec_update_interval == -1 && gn_dual_ksp_baseline > 0) {
    // Adaptive: Update if accumulate excess exceeds rank 
    int excess = ksp_iters_last - gn_dual_ksp_baseline;
    if (excess > 0) {
      gn_dual_ksp_excess_accum += excess;
    }
    if (gn_dual_ksp_excess_accum >= gn_precond_rank) {
      should_update = true;
      if (mpirank_world == 0 && !quietmode) {
        printf("Iteration %d: Updating preconditioner (accumulated excess %d >= r = %d)\n", iter, gn_dual_ksp_excess_accum, gn_precond_rank);
      }
    }
  }

  if (should_update) {
    updateGaussNewtonPreconditionerLowRank(gn_precond_rank);
    // Reset accumulator after adaptive update
    if (gn_prec_update_interval == -1) {
      gn_dual_ksp_excess_accum = 0;
    }
  }

  return should_update;
}


void OptimProblem::updateGaussNewtonPreconditionerLowRank(int rank){
  // Construct full-rank preconditioner M = V_r Σ_r V_r^T + μ(I - V_r V_r^T)
  // Algorithm:
  //   1. Use EPS to compute r+1 largest eigenvalues/eigenvectors of A
  //   2. Store V_r (first r eigenvectors), Σ_r (first r eigenvalues), μ = σ_{r+1}
  //   3. Apply M^{-1/2} = V_r Σ_r^{-1/2} V_r^T + μ^{-1/2}(I - V_r V_r^T)

  if (rank >= tangent_dim) {
    if (mpirank_world == 0 && !quietmode) {
      printf("Warning: rank %d >= tangent_dim %d, setting rank = tangent_dim - 1\n", rank, tangent_dim);
    }
    rank = tangent_dim - 1;
  }

  // 1. Use EPS to compute r+1 largest eigenvalues/eigenvectors of A
  int nev = rank + 1;  // Request r+1 eigenvalues

  if (mpirank_world == 0 && !quietmode) {
    printf("Computing %d largest eigenvalues/eigenvectors of the GN matrix (n=%d)...\n", nev, (int)tangent_dim);
  }

  EPS eps;
  EPSCreate(PETSC_COMM_SELF, &eps);
  if (gn_solver_mode == "dual") {
    EPSSetOperators(eps, GN_Dual_MatShell, NULL);
  } else if (gn_solver_mode == "primal") {
    EPSSetOperators(eps, GN_Primal_MatShell, NULL);
  }
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
  int vecdim;
  MatGetSize(gn_dual_lowrank_U, &vecdim, NULL);
  VecCreateSeq(PETSC_COMM_SELF, vecdim, &eigvec_real);
  VecCreateSeq(PETSC_COMM_SELF, vecdim, &eigvec_imag);

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

void OptimProblem::apply_GN_LowRank_Precond(const Vec x, Vec y){
  OptimProblem *self = this;

  // temp = V_r^T * x (project onto top-r subspace)
  MatMultTranspose(self->gn_dual_lowrank_U, x, self->gn_dual_lowrank_temp);

  // temp = Σ_r^{-1/2} .* temp (scale by inverse sqrt eigenvalues)
  VecPointwiseMult(self->gn_dual_lowrank_temp, self->gn_dual_lowrank_sigma_invsqrt, self->gn_dual_lowrank_temp);

  // y = V_r * temp (part 1: V_r Σ_r^{-1/2} V_r^T x)
  MatMult(self->gn_dual_lowrank_U, self->gn_dual_lowrank_temp, y);

  // Compute complement part: μ^{-1/2} (x - V_r V_r^T x)
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

void OptimProblem::applyEigenspaceApproximation(const Vec b, const Mat V_r, const PetscScalar* eigenvalues, int rank, bool have_complement, PetscScalar mu, Vec y_approx) {

  // Project RHS: b̃ = V_r^T b
  Vec b_tilde;
  VecCreateSeq(PETSC_COMM_SELF, rank, &b_tilde);
  MatMultTranspose(V_r, b, b_tilde);

  // Solve in subspace: ỹ = Σ_r^{-1} b̃ (direct diagonal solve)
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
          printf("Warning: Eigenvalue %d = %e is too small, setting component to zero\n", i, eigenvalues[i]);
      }
    }
  }
  VecRestoreArray(b_tilde, &b_tilde_data);
  VecRestoreArray(y_tilde, &y_tilde_data);

  // Lift to full space: y = V_r ỹ (top-r subspace part)
  MatMult(V_r, y_tilde, y_approx);

  // Add complement space contribution: y += μ^{-1} (I - V_r V_r^T) b
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

PetscErrorCode OptimProblem::GN_Dual_LowRank_MatMult(Mat A, const Vec z, Vec output){
  // Applies M^{-1/2} A M^{-1/2} z where A=GN_dual_MatShell = P_G W^1/2 L L^* W^1/2 P_G 

  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  Vec temp;
  VecDuplicate(z, &temp); // TODO: Should probably reuse a preallocated workspace!

  self->apply_GN_LowRank_Precond(z, temp);
  MatMult(self->GN_Dual_MatShell, temp, output);
  VecCopy(output, temp);  // Save intermediate result. TODO: Probably not needed. Check.
  self->apply_GN_LowRank_Precond(temp, output);

  VecDestroy(&temp);
  return 0;
}


PetscErrorCode OptimProblem::GN_Primal_LowRank_MatMult(Mat A, const Vec z, Vec output){
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  Vec temp;
  VecDuplicate(z, &temp); // TODO: Should probably reuse a preallocated workspace!

  self->apply_GN_LowRank_Precond(z, temp);
  MatMult(self->GN_Primal_MatShell, temp, output);
  VecCopy(output, temp);  // Save intermediate result. TODO: Probably not needed. Check.
  self->apply_GN_LowRank_Precond(temp, output);

  VecDestroy(&temp);
  return 0;
}


void OptimProblem::solveGaussNewtonLeastSquares(const Vec xinit, const Vec initial_guess, const Vec b, Vec v_LeastSquares){

  // Store the point of evaluation (design parameters)
  setGaussNewtonEvalPoint(xinit); 
  // Forces a fresh nonlinear forward solve during the next call to evalLinearizedForward() for the new xeval_GN
  // TODO: This can probably be dropped here, since it should be set before any solveGaussNewtonXYZ calls. 

  // Solve the Least Squares system with either TAOBRGN or KSPLSQR
  if (gn_leastsquares_solver == "brgn") {

    // Cache the RHS for residual evaluations (store -b for F(v) = L*v - b)
    VecCopy(b, brgn_rhs);

    // Set initial guess from initial_guess parameter
    VecCopy(initial_guess, v_LeastSquares);
    TaoSetSolution(tao_gn_leastsquares_brgn, v_LeastSquares);
    TaoSetUp(tao_gn_leastsquares_brgn);
    // Forward initial guess to the inner solver  
    Tao tao_sub;
    TaoBRGNGetSubsolver(tao_gn_leastsquares_brgn, &tao_sub);
    TaoSetSolution(tao_sub, v_LeastSquares);

    // Solve the regularized least-squares problem: min_v ||L*v - b||² + λ||v||²
    TaoSolve(tao_gn_leastsquares_brgn);

    // Report convergence. Tao iterations are outer Newton steps (one suffices for a linear problem), the actual cost is in the inner linear solver iterations.
    TaoConvergedReason reason;
    PetscInt newton_iters, lin_iters;
    TaoGetConvergedReason(tao_gn_leastsquares_brgn, &reason);
    TaoGetIterationNumber(tao_gn_leastsquares_brgn, &newton_iters);
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

  } else if (gn_leastsquares_solver == "lsqr"){

    // Set initial guess from initial_guess parameter
    VecCopy(initial_guess, v_LeastSquares);
    KSPSetInitialGuessNonzero(ksp_gn_leastsquares, PETSC_TRUE);

    // Solve the Gauss-Newton least-squares problem
    KSPReset(ksp_gn_leastsquares);
    KSPSetOperators(ksp_gn_leastsquares, GN_LeastSquares_MatShell, GN_LeastSquares_MatShell);
    KSPSolve(ksp_gn_leastsquares, b, v_LeastSquares);

    // Report convergence
    KSPConvergedReason reason;
    int iters;
    double rnorm;
    KSPGetConvergedReason(ksp_gn_leastsquares, &reason);
    KSPGetIterationNumber(ksp_gn_leastsquares, &iters);
    KSPGetResidualNorm(ksp_gn_leastsquares, &rnorm);
    ksp_iters_last = iters;
    if (mpirank_world == 0 && !quietmode) {
      printf("Gauss-Newton Least-Squares stats: iterations = %d, residual norm = %1.14e\n", iters, rnorm);
    }
  } else {
    if (mpirank_world == 0 && !quietmode) {
      printf("Unknown Gauss-Newton least-squares solver: %s\n", gn_leastsquares_solver.c_str());
      exit(1);
    }
  }
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

  // Fill the top part of y with Lv projected onto tangent space
  self->apply_GN_Jacobian_Tangentspace(v, y);

  // Fill the bottom of y with damping: lambda * v
  PetscInt M;
  MatGetSize(A, &M, nullptr);
  const PetscScalar *v_data;
  PetscScalar *y_data;
  VecGetArrayRead(v, &v_data);
  VecGetArray(y, &y_data);
  for (int i=0; i<M - self->tangent_dim; i++) {
    y_data[i + self->tangent_dim] = self->gn_leastsquares_damping * v_data[i];
  }
  VecRestoreArrayRead(v, &v_data);
  VecRestoreArray(y, &y_data);

  return 0;
}

PetscErrorCode OptimProblem::GN_LeastSquares_MatMultTranspose(Mat A, Vec w_long, Vec vout)
{
  OptimProblem *self;
  MatShellGetContext(A, (void**)&self);

  int M,N;
  MatGetSize(A, &M, &N);

  // Create a tangent space vector w and fill it will the first part of w_long
  Vec w;
  VecCreateSeq(PETSC_COMM_SELF, self->tangent_dim, &w);
  VecSetFromOptions(w);
  VecAssemblyBegin(w); VecAssemblyEnd(w);
  PetscScalar *w_ptr;
  const PetscScalar *wl_data;
  VecGetArrayRead(w_long, &wl_data);
  VecGetArray(w, &w_ptr);
  for (int i=0; i<self->tangent_dim; i++) {
    w_ptr[i] = wl_data[i];
  }
  VecRestoreArray(w, &w_ptr);
  VecRestoreArrayRead(w_long, &wl_data);

  // Apply adjoint of GN Jacobian in tangent space
  self->apply_GN_Jacobian_Tangentspace_Adjoint(w, vout);

  // Add damping lambda * w_long[tangent_dim:end] contribution to the output
  PetscScalar* vout_data;
  VecGetArray(vout, &vout_data);
  VecGetArrayRead(w_long, &wl_data);
  for (int i=0; i<self->ndesign; i++) {
    vout_data[i] += self->gn_leastsquares_damping * wl_data[i + self->tangent_dim];
  }
  VecRestoreArray(vout, &vout_data);
  VecRestoreArrayRead(w_long, &wl_data);
  VecDestroy(&w);
  return 0;
}


void OptimProblem::apply_GN_Jacobian_Tangentspace(Vec v, Vec y)
{
  // Jacobian forward: y = P_G W^1/2 Lv 

  OptimProblem *self = this;

  // Reset output
  VecZeroEntries(y);

  // Apply linearized forward to get dU/dalpha * v
  self->evalLinearizedForward(self->xeval_GN, v);

  // If include W^1/2, precompute z = 1/n<V,W> 
  double z_re = 0.0;
  double z_im = 0.0;
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

  // For each initial condition, apply W^1/2 and Project onto tangent space
  PetscScalar *y_data;
  VecGetArray(y, &y_data);
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;

    // w_i(T). Use rho_t0_bar as temporary storage
    Vec lin_final_state = self->timestepper->getLinearizedFinalState(iinit);
    VecCopy(lin_final_state, self->rho_t0_bar);  

    // Include W^1/2 = sqrt(2/n)(I-z*V), if requested
    double scale = sqrt(2.0 / self->mastereq->getDim());
    VecScale(self->rho_t0_bar, scale);
    self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
    self->optim_target->HilbertSchmidtOverlap_diff(self->rho_t0_bar, false, -scale*z_re, -scale*z_im);

    // Project onto tangent space -> fills the appropriate entries in y_data
    Mat U_final_re = self->optim_target->getFinalUnitaryRe();
    Mat U_final_im = self->optim_target->getFinalUnitaryIm();
    GellMann::projectVecToTangentSpace(self->rho_t0_bar, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, y_data);
  }

  // Allreduce across comm_init processors
  MPI_Allreduce(MPI_IN_PLACE, y_data, tangent_dim, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(y, &y_data);

  // No damping rows - pure Jacobian only
}

void OptimProblem::apply_GN_Jacobian_Tangentspace_Adjoint(Vec w, Vec vout)
{
  // Jacobian transpose: vout = L^T W^1/2 P_G^{-1} w 

  OptimProblem *self = this;

  // Reset output
  VecZeroEntries(vout);

  // Get read access to w (tangent space)
  const PetscScalar *w_data;
  VecGetArrayRead(w, &w_data);

  // P_G^{-1}: For each initial condition, reconstruct full state vectors from tangent space coefficients
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    // Reconstruct full state vector (adjoint of projectVecToTangentSpace) -> fills wsub_workspace[iinit]
    Mat U_final_re = self->optim_target->getFinalUnitaryRe();
    Mat U_final_im = self->optim_target->getFinalUnitaryIm();
    GellMann::reconstructVecFromTangentSpace(w_data, U_final_re, U_final_im, self->mastereq->getDim(), iinit_global, self->wsub_workspace[iinit]);
  }
  VecRestoreArrayRead(w, &w_data);

  // If include W^1/2, precompute z = 1/n<V,wsub_workspace> 
  double z_re = 0.0;
  double z_im = 0.0;
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
  

  // Solve adjoint backward ODE for each terminal condition
  for (int iinit = 0; iinit < self->ninit_local; iinit++) {

    // Include W = sqrt(2/n)(I-<V,cdot>V), if requested
    double scale = sqrt(2.0 / self->mastereq->getDim());
    VecScale(self->wsub_workspace[iinit], scale);
    int iinit_global = self->mpirank_init * self->ninit_local + iinit;
    self->optim_target->prepareInitialAndTargetState(iinit_global, self->ninit, self->mastereq->nlevels, self->mastereq->nessential);
    self->optim_target->HilbertSchmidtOverlap_diff(self->wsub_workspace[iinit], false, -scale*z_re, -scale*z_im);

    self->timestepper->solveAdjointODE(iinit, self->wsub_workspace[iinit], 0.0, 0.0, 0.0, 0.0);
    VecAXPY(vout, 1.0, self->timestepper->getReducedGradient());
  }

  // Sum up the gradient from all initial condition processors
  PetscScalar* vout_data;
  VecGetArray(vout, &vout_data);
  MPI_Allreduce(MPI_IN_PLACE, vout_data, self->ndesign, MPIU_SCALAR, MPI_SUM, self->comm_init);
  VecRestoreArray(vout, &vout_data);
}
