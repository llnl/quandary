#include "math.h"
#include <assert.h>
#include <petsctao.h>
#include "defs.hpp"
#include "timestepper.hpp"
#include <iostream>
#include <algorithm>
#include "optimtarget.hpp"
#pragma once

/* if Petsc version < 3.17: Change interface for Tao Optimizer */
#if PETSC_VERSION_MAJOR<4 && PETSC_VERSION_MINOR<17
#define TaoSetObjective TaoSetObjectiveRoutine
#define TaoSetGradient(tao, NULL, TaoEvalGradient,this)  TaoSetGradientRoutine(tao, TaoEvalGradient,this) 
#define TaoSetObjectiveAndGradient(tao, NULL, TaoEvalObjectiveAndGradient,this) TaoSetObjectiveAndGradientRoutine(tao, TaoEvalObjectiveAndGradient,this)
#define TaoSetSolution(tao, xinit) TaoSetInitialVector(tao, xinit)
#define TaoGetSolution(tao, params) TaoGetSolutionVector(tao, params) 
#endif
/* if Petsc version < 3.21: Change interface for Tao Monitor */
#if PETSC_VERSION_MAJOR<4 && PETSC_VERSION_MINOR<21
#define TaoMonitorSet TaoSetMonitor 
#endif

/**
 * @brief Optimization problem solver for quantum optimal control.
 *
 * This class manages the optimization of quantum control pulses using PETSc's TAO
 * optimization library. It handles objective function evaluation, by propagating initial states 
 * forward in time solving the dynamical equation and computing the final-time objective cost function 
 * and integral penalty terms, as well as the gradient computation by backpropagating the adjoint 
 * terminal states backwards in time solving the adjoint dynamical equation and collecting gradient
 * contributions. It further defines the interface functions for PETSc's TAO optimization 
 * via L-BFGS, including a callback function to monitor optimization progress. 
 * 
 * Main functionality:
 *    - @ref evalF evaluates the objective function by calling @ref TimeStepper::solveODE to evolve initial states 
 *      to final time T and summing up the objective function measure over each target state
 *    - @ref evalGradF evaluates the objective function and its gradient with respect to the optimization parameters
 *       by calling @ref TimeStepper::solveODE and @ref TimeStepper::solveAdjointODE to propagate initial states
 *       forward and backward through the time domain while accumulating objective and gradient information.
 * 
 * This class contains references to:
 *    - @ref TimeStepper for handling the forward and backward time stepping process
 *    - @ref OptimTarget for evaluating the final-time cost for each initial condition
 *    - @ref Output      for writing monitored optimization convergence to file
 */
class OptimProblem {
  protected:

  size_t ninit; ///< Number of initial conditions to be considered (N^2, N, or 1)
  int ninit_local; ///< Local number of initial conditions on this processor
  Vec rho_t0_bar; ///< Storage for adjoint initial condition of the adjoint ODE (aka the terminal condition)

  OptimTarget* optim_target; ///< Pointer to the optimization target (gate or state)

  MPI_Comm comm_init; ///< MPI communicator for initial condition parallelization
  MPI_Comm comm_optim; ///< MPI communicator for optimization parallelization, currently not used (size 1)
  int mpirank_optim, mpisize_optim; ///< MPI rank and size for optimization communicator
  int mpirank_petsc, mpisize_petsc; ///< MPI rank and size for spatial parallelization (PETSc)
  int mpirank_world, mpisize_world; ///< MPI rank and size for global communicator
  int mpirank_init, mpisize_init; ///< MPI rank and size for initial condition communicator

  bool quietmode; ///< Flag for quiet mode operation

  OptimSolverType optim_solver_type; ///< Type of optimization solver to use
  std::vector<double> obj_weights; ///< Weights for averaging objective over initial conditions
  int ndesign; ///< Number of global design (optimization) parameters
  double objective = 0.0; ///< Current objective function value (sum over final-time cost, regularization terms and penalty terms)
  double obj_cost = 0.0; ///< Final-time measure J(T) in objective
  double obj_regul = 0.0; ///< Regularization term in objective
  double obj_penal_leakage = 0.0; ///< Penalty term for leakage into guard levels
  double obj_penal_weightedcost = 0.0; ///< Penalty term for weighted running cost 
  double obj_penal_dpdm = 0.0; ///< Penalty term second-order state derivatives (penalizes variations of the state evolution)
  double obj_penal_variation = 0.0; ///< Penalty term for variation of control parameters
  double obj_penal_energy = 0.0; ///< Energy penalty term in objective
  double fidelity = 0.0; ///< Final-time fidelity: 1/ninit sum_i Tr(rho_target^dag rho(T)) for Lindblad, |1/ninit sum_i phi_target^dag phi|^2 for Schrodinger
  double gnorm = 0.0; ///< Current norm of gradient
  double gamma_tikhonov; ///< Parameter for Tikhonov regularization
  bool tikhonov_use_x0; ///< Switch to use ||x - x0||^2 for Tikhonov regularization instead of ||x||^2
  double gamma_penalty_leakage; ///< Parameter multiplying integral leakage term
  double gamma_penalty_weightedcost; ///< Parameter multiplying integral weighted cost function 
  double gamma_penalty_dpdm; ///< Parameter multiplying integral penalty term for 2nd derivative of state variation
  double gamma_penalty_energy; ///< Parameter multiplying energy penalty
  double gamma_penalty_variation; ///< Parameter multiplying finite-difference squared regularization term
  double tol_grad_abs; ///< Stopping criterion based on absolute gradient norm
  double tol_grad_rel; ///< Stopping criterion based on relative gradient norm
  double tol_final_cost; ///< Stopping criterion based on objective function value
  double tol_infidelity; ///< Stopping criterion based on infidelity
  int maxiter; ///< Stopping criterion based on maximum number of iterations
  Tao tao_lbfgs; ///< PETSc's TAO LBFGS optimization solver
  double* mygrad; ///< Auxiliary gradient storage
  Vec xtmp; ///< Temporary vector storage
  int output_optimization_stride; ///< Write output files every N optimization iterations

  TimeStepper* timestepper; ///< Pointer to time-stepping scheme
  Output* output; ///< Pointer to output handler
  MasterEq* mastereq; ///< Pointer to master equation solver

  /* Storage for Gauss-Newton optimization solvers */
  Vec x_GN; ///< Current iterate for Gauss-Newton optimization. Holds solution after finished. 
  Vec xeval_GN; ///< Point of evaluation for Gauss-Newton (evaluate L(x) at this x) 
  bool nonlinear_forward_valid; ///< True once the nonlinear forward solveODE has been done and stored for the current xeval_GN, reset whenever xeval_GN changes!

  // Options for the Armijo line search 
  const double c1 = 1e-4;    //< Sufficient decrease parameter
  const double rho_backtrack = 0.5;   ///< Backtracking factor
  const int max_ls_iter = 20; ///< Maximum number of backtracking steps

  // All matshells/densemats for Gauss-Newton solvers
  Mat GN_Primal_MatShell; ///< Applies L^*WL  
  Mat GN_Primal_MatDense; ///< Dense L^*WL (debug only)
  bool gn_primal_usedensemat = false; ///< Always false.
  Mat GN_LeastSquares_MatShell; ///< Applies P_gL and, if damped, augments with lambda I. Misses W^1/2 at this point. TODO: Should probably be calling GN_Jacobian_Tangentspace plus augmentation. 
  Mat GN_Dual_MatShell; ///< Applies P_G W^1/2 L L^* W^1/2 P_G
  Mat GN_Dual_MatDense; ///< Dense P_G W^1/2 L L^* W^1/2 P_G
  Mat GN_Dual_LowRank_MatShell; ///< Applies low-rank preconditioned dual: M^{-1/2} (L L^T + λI) M^{-1/2}
  Mat GN_Dual_Precond_MatDense; // Preconditoning matrix in dense format

  // Solvers for the inner Gauss-Newton equations
  KSP ksp_gn_normaleq;  ///< Linear solver for Gauss-Newton normal equation (primal or dual)
  KSP ksp_gn_leastsquares; ///< Linear solver for Gauss-Newton least-squares problem
  Tao tao_gn_leastsquares_brgn;        ///< Tao BRGN solver for least squares
  EPS eps_gn_primal; // EPS solver for the Gauss-Newton primal equation 

  // Options for Gauss-Newton solvers
  std::string gn_solver_mode; ///< Gauss-Newton solver mode: "primal" or "dual" or "least-squares"
  std::string gn_leastsquares_solver; ///< Least-squares solver name: "brgn" or "lsqr"
  double gn_leastsquares_damping; ///< Damping parameter λ for Least Squares regularization
  std::string gn_normaleq_solver; ///< Normal equation solver: "minres", "cg", "gmres", or "proj"
  double gn_normaleq_damping; ///< Damping parameter for Gauss-Newton matrix shift
  int gn_prec_update_interval; ///< Preconditioner update interval (0 = startup only, -1 = adaptive)
  int gn_precond_rank; ///< Rank for low-rank preconditioner
  int gn_maxiter; ///< Maximum iterations for Gauss-Newton KSP solver
  double gn_rtol; ///< Relative tolerance for Gauss-Newton KSP solver
  PetscReal eps_tol = 1e-4; ///< Tolerance for EPS eigenvalue solver in primal mode
  double eps_evals_cutoff = 1e-5; ///< Cutoff for eigenvalues of the Gauss-Newton matrix in primal mode
  double eps_primal_damping = 1e-3; ///< Damping for eigenvalues of the Gauss-Newton matrix

  // Workspace for Gauss-Newton solvers
  int tangent_dim; ///< Dimension of the tangent space (dim²-1)
  bool gn_precond_enabled; ///< Flag whether Gauss-Newton preconditioner is enabled
  int ksp_iters_last; ///< Number of KSP iterations used in the most recent Gauss-Newton normal or least squares solvers 
  Vec *wsub_workspace; ///< Pre-allocated workspace to compute -\nabla_U J
  Vec brgn_rhs; ///< Cached RHS vector for BRGN residual evaluation
  Vec brgn_residual; ///< Cached Residual vector for BRGN solver
  Vec y_dual_workspace; ///< Workspace vector for dual solver (tangent space, size dim²-1)
  Vec v_dual_workspace; ///< Workspace vector for intermediate L^T y computation (design space)
  int gn_dual_ksp_baseline; ///< Baseline KSP iterations right after preconditioner formation (for adaptive update)
  int gn_dual_ksp_excess_accum; ///< Accumulated excess KSP iterations beyond n_0 (for adaptive update)
  Mat gn_dual_lowrank_U; ///< V_r: first r eigenvectors (tangent_dim × rank, dense)
  Vec gn_dual_lowrank_sigma_invsqrt; ///< Σ_r^{-1/2} (rank eigenvalues^{-1/2})
  Vec gn_dual_lowrank_temp; ///< Workspace for low-rank application (size = rank)
  PetscScalar gn_dual_precond_mu; ///< μ = σ_{r+1} for complement space

  // TODO: Move these to local functions
  int neigvals; ///< Number of eigenvalues to compute (=N^2-1)
  int ncv; ///< Number of Lanczos vectors to use in EPS solver. HOW TO CHOOSE?? 

  public: 
    Vec xlower, xupper; ///< Lower and upper bounds for optimization variables
    Vec xprev; ///< Design vector at previous iteration
    Vec xinit; ///< Initial design vector


  /**
   * @brief Constructor for optimization problem.
   *
   * @param config Configuration parameters from input file
   * @param optim_target_ Pointer to optimization target
   * @param timestepper_ Pointer to time-stepping scheme
   * @param mastereq_ Pointer to master equation solver
   * @param comm_init_ MPI communicator for initial condition parallelization
   * @param comm_optim MPI communicator for optimization parallelization
   * @param output_ Pointer to output handler
   * @param quietmode Flag for quiet operation (default: false)
   */
  OptimProblem(const Config& config, OptimTarget* optim_target_, TimeStepper* timestepper_, MasterEq* mastereq_, MPI_Comm comm_init_, MPI_Comm comm_optim, Output* output_, bool quietmode=false);

  ~OptimProblem();

  int getNdesign(){ return ndesign; };
  double getObjective(){ return objective; };
  double getCostT()    { return obj_cost; };
  double getRegul()    { return obj_regul; };
  double getPenaltyLeakage()  { return obj_penal_leakage; };
  double getPenaltyWeightedCost()  { return obj_penal_weightedcost; };
  double getPenaltyDpDm()  { return obj_penal_dpdm; };
  double getPenaltyVariation()  { return obj_penal_variation; };
  double getPenaltyEnergy()  { return obj_penal_energy; };
  double getFidelity() { return fidelity; };
  double getTolFinalCost()    { return tol_final_cost; };
  double getTolGradAbs()    { return tol_grad_abs; };
  double getTolInfidelity()   { return tol_infidelity; };
  int getMPIrank_world() { return mpirank_world;};
  int getMPIsize_world() { return mpisize_world;};
  int getMaxIter()     { return maxiter; };
  OptimTarget* getOptimTarget() { return optim_target; };
  Mat getGN_Primal_MatShell() { return GN_Primal_MatShell; };
  Mat getGN_Primal_MatDense() { return GN_Primal_MatDense; };
  Mat getGN_Dual_MatShell() { return GN_Dual_MatShell; };
  Mat getGN_Dual_LowRank_MatShell() { return GN_Dual_LowRank_MatShell; };
  Mat getGN_Dual_LowRank_U() { return gn_dual_lowrank_U; };
  Vec getGN_Dual_LowRank_Sigma_InvSqrt() { return gn_dual_lowrank_sigma_invsqrt; };
  PetscScalar getGN_Dual_Precond_Mu() { return gn_dual_precond_mu; };
  Mat getGN_Dual_MatDense() { return GN_Dual_MatDense; };
  Mat getGN_Dual_Precond_MatDense() { return GN_Dual_Precond_MatDense; };
  Mat getGNLeastSquaresShell() { return GN_LeastSquares_MatShell; };
  Mat getGN_LeastSquares_MatShell() { return GN_LeastSquares_MatShell; };
  bool getQuietmode() { return quietmode; };

  int getOutputOptimizationStride() { return output_optimization_stride; };
  Output* getOutput() { return output; };
  TimeStepper* getTimeStepper() { return timestepper; };

  /**
   * @brief Override the maximum number of iterations for all Gauss-Newton solvers.
   */
  void setGaussNewtonMaxiter(int maxiter);

  /**
   * @brief Sets the evaluation point for Gauss-Newton matrix operators 
   *
   * Copies xinit to xeval_GN and invalidates the nonlinear forward solve cache.
   *
   * @param xinit Design parameter vector to use as evaluation point
   */
  void setGaussNewtonEvalPoint(Vec xinit);

  /**
   * @brief Evaluates the objective function F(x).
   * 
   * Performs forward simulations for each initial conditions and
   * evaluates the objective function. 
   *
   * @param x Design vector
   * @param writeTrajectoryDataFiles Flag to determine whether trajectory data should be written during forward simulations to files (default: false)
   * @return double Objective function value
   */
  double evalF(const Vec x, bool writeTrajectoryDataFiles=false);

  /**
   * @brief Evaluates the gradient of the objective function with respect to the control parameters
   *
   * @param x Design (optimization) vector
   * @param G Gradient vector to store result
   */
  void evalGradF(const Vec x, Vec G, bool writeTrajectoryDataFiles=false);

  /**
   * @brief Evaluate linearized forward operator: Lv = \sum_k dU/dv_k v_k. 
   * 
   * This does one linearized ODE solve, preceded by a nonlinear ODE solve only if the nonlinear
   * trajectory for x has not already been computed and stored (see nonlinear_forward_valid). After
   * this, the timesteppers trajectory_states and lin_trajectory_states will be set.
   * 
   * @param[in] x Point of evaluation
   * @param[in] v Direction vector 
   */
  void evalLinearizedForward(const Vec x, const Vec v);

  /**
   * @brief Solves the Gauss-Newton least-squares problem min_v ||W^1/2 L v + W^-1/2 \nabla_UJ||_F^2
   *
   * First fills the VecNest RHS with the appropriate values based on the current point of evaluation xinit, then solves the least-squares problem to obtain v_LeastSquares.
   *
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param initial_guess Initial guess for the solution
   * @param b Right-hand side vector
   * @param v_LeastSquares Solution vector to store the result
   */
  void solveGaussNewtonLeastSquares(const Vec xinit, const Vec initial_guess, const Vec b, Vec v_LeastSquares);


  /**
   * @brief Callback for TaoBRGN residual evaluation: F(v) = L*v - b
   */
  static PetscErrorCode TaoBRGN_EvalResidual(Tao tao, Vec v, Vec F, void *ctx);

  /**
   * @brief Callback for TaoBRGN Jacobian evaluation: dF/dv = L (constant)
   */
  static PetscErrorCode TaoBRGN_EvalJacobianResidual(Tao tao, Vec v, Mat J, Mat Jpre, void *ctx);

  /**
   * @brief Monitor function for TaoBRGN iterations
   */
  static PetscErrorCode TaoBRGN_Monitor(Tao tao, void *ctx);

  /**
   * @brief MatMult for GN_Primal_MatShell: Applies L^*WL
   * 
   * Solves the linearized forward equation, followed by an adjoint solve with terminal condition WL. 
   * 
   * @note: Requires xeval_GN to be set before calling. Performs a nonlinear solveODE at xeval_GN if !nonlinear_forward_valid, so to store U(t;xeval_GN) for all t. 
   */
  static PetscErrorCode GN_Primal_MatMult(Mat A, const Vec v, Vec Av);

  /**
   * @brief MatMult for GN_LeastSquares_MatShell: Applies P_G L, and, if damped, augments with lambda_I 
   * 
   * @note: Requires xeval_GN to be set before calling. Performs a nonlinear solveODE at xeval_GN if !nonlinear_forward_valid, so to store U(t;xeval_GN) for all t. 
   * 
   * @note: TODO: Implement W^1/2 in tangent space representation for the least-squares solver. Probably use the apply_GN_Jacobian_Tangentspace function.
   */
  static PetscErrorCode GN_LeastSquares_MatMult(Mat A, Vec v, Vec y); 
  static PetscErrorCode GN_LeastSquares_MatMultTranspose(Mat A, Vec w, Vec vout); 

  /**
   * @brief MatMult GN_Dual_MatShell: Applies P_G W^1/2 L L^* W^1/2 P_G 
   * 
   * Calls apply_GN_Jacobian_Tangentspace_Adjoint (L*^ W^1/2 P_G)followed by apply_GN_Jacobian_Tangentspace (P_G W^1/2 L)
   * 
   * @note: Requires xeval_GN to be set before calling. Performs a nonlinear solveODE at xeval_GN if !nonlinear_forward_valid, so to store U(t;xeval_GN) for all t. 
   */
  static PetscErrorCode GN_Dual_MatMult(Mat A, const Vec y, Vec Ay);

  /**
   * @brief MatMult for GN_Dual_LowRank_MatShell: Applies M^{-1/2} A M^{-1/2} z where A=GN_dual_MatShell  
   */
  static PetscErrorCode GN_Dual_LowRank_MatMult(Mat A, const Vec z, Vec output);


  /**
   * @brief Applies the Gauss-Newton Jacobian projected onto the tangent space: Y = P_G W^1/2 Lv
   * 
   * @note: Requires xeval_GN to be set before calling. Performs a nonlinear solveODE at xeval_GN if !nonlinear_forward_valid, so to store U(t;xeval_GN) for all t. 
   * 
   * @param v Input vector in the design space
   * @param y Output vector in the tangent space
   */
  void apply_GN_Jacobian_Tangentspace(Vec v, Vec y); 

  /**
   * @brief Applies the adjoint of the Gauss-Newton Jacobian projected onto the tangent space: vout = (P_G W^1/2 L)^T w
   * 
   * @note: Requires that U(t; xeval_GN) is valid (i.e., a nonlinear forward solveODE has been performed at xeval_GN)
   * 
   * @param w Input vector in the tangent space
   * @param vout Output vector in the design space
   */
  void apply_GN_Jacobian_Tangentspace_Adjoint(Vec w, Vec vout); 

  /**
   * @brief Updates the dense primal Gauss-Newton matrix based on the current point of evaluation xeval_GN by calling GN_Primal_MatShell on each unit vectoor.
   * 
   * @note This operation can be expensive as it involves N^2-1 applications of the MatShell.
   */
  void updateGaussNewtonPrimalMatDense();

  /**
   * @brief Assembles the dense dual Gauss-Newton matrix using GN_Dual_MatShell
   *
   * @note This operation can be expensive as it involves (dim²-1) applications of the dual MatShell.
   */
  void updateGaussNewtonDualMatDense();

  /**
   * @brief Updates the Gauss-Newton preconditioner if needed
   * 
   * Currently only available for the DUAL system
   * 
   * @returns True if the preconditioner was updated, false otherwise.
   */
  bool updateGaussNewtonPreconditioner(int iter, bool force_update);

  /**
   * @brief Constructs rank-r low-rank approximation M ≈ U Σ U^T for dual Gauss-Newton preconditioner
   *
   * Uses EPS solver to compute r+1 largest eigenvalues and eigenvectors of the dual GN_Dual_MatShell. Populates the low-rank preconditioner data structures: V_r (gn_dual_lowrank_U), Σ_r (gn_dual_lowrank_sigma_invsqrt), and μ (gn_dual_precond_mu).
   * 
   * @param rank Number of basis vectors to use (r)
   */
  void updateGaussNewtonDualPreconditionerLowRank(int rank);


  /**
   * @brief Construct dense verson of the rank-r low-rank preconditoner M^{-1/2} ≈ U Σ^{-1/2} U^T + mu^{-1/2} (I - U U^T) for the dual Gauss-Newton preconditioner
   *
   * Uses EPS solver to compute r+1 largest eigenvalues and eigenvectors of the dual GN_Dual_MatShell. Populates the low-rank preconditioner data structures: V_r (gn_dual_lowrank_U), Σ_r (gn_dual_lowrank_sigma_invsqrt), and μ (gn_dual_precond_mu).
   * 
   */
  void updateGaussNewtonDualPreconditionerMatDense();

  /**
   * @brief Applies the low-rank dual Gauss-Newton preconditioner: y = M^{-1/2} x, where M ≈ U Σ U^T is the precomputed low-rank approximation.
   * 
   * This applies the low-rank preconditioner using the formula:
   *   M^{-1/2} = V_r Σ_r^{-1/2} V_r^T + μ^{-1/2}(I - V_r V_r^T)
   *   y = V_r (Σ_r^{-1/2} .* (V_r^T x)) + μ^{-1/2} (x - V_r (V_r^T x))
   *
   * @param[in] x Input vector (tangent space, size dim²-1)
   * @param[out] y Output vector after applying M^{-1/2} (tangent space)
   */
  void apply_GN_Dual_LowRank_Precond(const Vec x, Vec y);

  /**
   * @brief Solves the primal Gauss-Newton equation iteratively using the ksp_gn_normaleq solver
   *
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param initial_guess Initial guess for iterative solver
   * @param b Right-hand side vector
   * @param Ainv_b Solution vector to store the result
   */
  void solveGaussNewtonPrimalKSP(Vec xinit, const Vec initial_guess, const Vec b, Vec Ainv_b);

  /**
   * @brief Solves the primal Gauss-Newton equation via eigenvalue decomposition
   *
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param b Right-hand side vector
   * @param Ainv_b Solution vector to store the result
   */
  void solveGaussNewtonPrimalProjected(Vec xinit, const Vec b, Vec Ainv_b);

  /**
   * @brief Solves the dual Gauss-Newton problem iteratively using low-rank preconditioned KSP. 
   *
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param initial_guess Initial guess for iterative solver (in tangent space, size dim²-1)
   * @param b_rhs Right-hand side vector (tangent space, size dim²-1)
   * @param v_out Solution vector (in design space, size ndesign)
   * @return true if eigendecomposition appears stale, false otherwise
   */
  bool solveGaussNewtonDualPrecondKSP(Vec xinit, const Vec initial_guess, const Vec b_rhs, Vec v_out);

  /**
   * @brief Solves the dual Gauss-Newton problem by projection onto dominant eigenspace.
   *
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param b_rhs Right-hand side vector (tangent space, size dim²-1)
   * @param v_out Solution vector (in design space, size ndesign)
   */
  void solveGaussNewtonDualProjected(Vec xinit, const Vec b_rhs, Vec v_out);

  /**
   * @brief Compute RHS for Gauss-Newton problem depending on the solver mode. 
   * 
   * primal: b = -G in design space
   * least-squares and dual: b = -∇_U J projected to tangent space 
   * 
   * @return Newly created RHS vector for the Gauss-Newton system
   */
  Vec computeGaussNewtonRHS(const std::string& gn_solver_mode, const Vec G);

  /** 
   * @brief Computes an initial guess for the Gauss-Newton system depending on the solver mode.
   * 
   * @return Newly created initial guess vector for the Gauss-Newton system
   */
  Vec computeGaussNewtonInitialGuess(const std::string& gn_solver_mode, const Vec b_rhs);

  /**
   * @brief Top-level function to solve the Gauss-Newton system. 
   * 
   * Solves either "primal", "dual", or "least_squares" Gauss-Newton systems depending on gn_solver_mode. Delegates to the appropriate solver function.
   * 
   * @return true if the preconditioner is stale (if it is available), false otherwise
   */
  bool solveGaussNewtonSystem(const Vec xinit, const Vec initial_guess, const Vec b, Vec Ainv_b);

  /**
   * @brief Computes a low-rank approximation of the dual Gauss-Newton matrix using provided eigenvalues and eigenvectors.
   * 
   * Compute eigenspace approximation: y_approx = M^{-1} * b
   * where M = V_r Σ_r V_r^T + μ(I - V_r V_r^T)
   *
   * @param b Right-hand side vector
   * @param V_r Matrix of eigenvectors (ndesign x rank)
   * @param eigenvalues Array of eigenvalues corresponding to V_r
   * @param rank Number of eigenvectors/eigenvalues
   * @param have_complement Whether a complementary subspace is available
   * @param mu Regularization parameter for the complementary subspace
   * @param y_approx Output vector storing the approximation result
   */
  void applyEigenspaceApproximation(const Vec b, const Mat V_r, const PetscScalar* eigenvalues, int rank, bool have_complement, PetscScalar mu, Vec y_approx);

  /**
   * @brief Backtracking Armijo line search along a descent direction, projected onto bound constraints.
   *
   * Finds a step length alpha such that f(P[x - alpha*dir]) <= f(x) + c1 * grad^T (P[x - alpha*dir] - x), where P[.] clips onto bounds [xlower, xupper]. Falls back to the steepest descent direction if dir is not a descent direction (e.g. after projection).
   *
   * @param[in] x Current iterate
   * @param[in] f Objective function value at x
   * @param[in] grad Gradient at x
   * @param[in] dir Search direction (e.g. preconditioned gradient)
   * @param[out] xnew Accepted trial point
   * @param[out] step Work vector to store the accepted step xnew-x
   * @return Accepted step length alpha
   */
  double armijoLineSearch(Vec x, double f, Vec grad, Vec dir, Vec xnew, Vec step);

  /**
   * @brief Compute evals and evecs of Gauss-Newton Normal Equation A=L^*L matrix
   * 
   * @param[in] xinit Point of evaluation
   * @param[out] evecs_out Newly created dense matrix (ndesign x number-of-converged-evals) holding one eigenvector per column
   * @return Eigenvalues of Gauss-Newton matrix
   */
  std::vector<double> computeGaussNewtonPrimalEvals(Vec xinit, Mat* evecs_out);

  /**
   * @brief Runs the optimization solver.
   *
   * @param xinit Initial guess for design variables
   */
  void solve(Vec xinit);

  /**
   * @brief Computes initial guess for optimization variables.
   *
   * @param x Vector to store the initial guess
   */
  void getStartingPoint(Vec x);

  /**
   * @brief Retrieves the optimization solution and prints summary information.
   *
   * This method should be called after TaoSolve() has finished.
   *
   * @param opt Pointer to vector to store the optimal solution
   */
  void getSolution(Vec* opt);

  /**
   * @brief Monitor function called in each optimization iteration.
   * 
   * @param iter Current iteration number
   * @param f Current objective function value
   * @param gnorm Current gradient norm
   * @param deltax Current step size
   * @return True if this stopping criterion is satisfied (if last iteration)
   */
  bool monitor(int iter, double f, double gnorm, double deltax);
};

/**
 * @brief Monitors optimization progress during TAO optimization iterations.
 *
 * This callback function is called at each iteration of TaoSolve() to
 * track convergence and output progress information.
 *
 * @param tao TAO solver object
 * @param ptr Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode TaoMonitor(Tao tao,void*ptr);


/**
 * @brief Synchronizes the solution vector across all MPI ranks at the start of each iteration. 
 *
 * @param tao TAO solver object
 * @param iter Current iteration number
 * @param ptr Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode TaoSynchronize(Tao tao, PetscInt iter, void*ptr);

/**
 * @brief PETSc TAO interface routine for combined objective and gradient evaluation.
 *
 * @param tao TAO solver object
 * @param x Design vector
 * @param f Pointer to store objective function value
 * @param G Gradient vector
 * @param ptr Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode TaoEvalObjectiveAndGradient(Tao tao, Vec x, PetscReal *f, Vec G, void*ptr);

/**
 * @brief Monitors the Gauss-Newton KSP solve, printing residual and solution norms.
 *
 * @param ksp KSP solver object
 * @param it Iteration number
 * @param rnorm (Unpreconditioned) residual norm at this iteration
 * @param ctx Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode KSPMonitorResidualAndSolution(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx);
