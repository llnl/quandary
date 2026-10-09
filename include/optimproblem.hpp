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
  Tao tao; ///< PETSc's TAO optimization solver
  double* mygrad; ///< Auxiliary gradient storage
  Vec xtmp; ///< Temporary vector storage
  int output_optimization_stride; ///< Write output files every N optimization iterations

  TimeStepper* timestepper; ///< Pointer to time-stepping scheme
  Output* output; ///< Pointer to output handler
  MasterEq* mastereq; ///< Pointer to master equation solver

  Mat GN_NormalEq_MatShell; ///< MatShell for applying Gauss-Newton normal equations matrix A(x)=L(x)^* L(x) to a vector v 
  Mat GN_NormalEq_MatDense; ///< Dense matrix representation of the Gauss-Newton normal equation matrix
  bool GN_NormalEq_usedensemat = false; ///< Flag indicating if the dense matrix representation of the Gauss-Newton normal equation matrix is used. Always false.

  Vec xeval_GN; ///< Point of evaluation for Gauss-Newton (evaluate L(x) at this x) 

  Vec x_GN; ///< Current iterate for GN optimization. Holds solution after finished. 
  int ksp_iters_last; ///< Number of KSP iterations used in the most recent Gauss-Newton normal or least squares solvers 
  bool nonlinear_forward_valid; ///< True once the nonlinear forward has been solved and stored for the current xeval_GN, reset whenever xeval_GN changes
  bool includeHessUJ; ///< Flag to include Hessian of J(U) in the terminal adjoint condition

  // Options for the Armijo line search 
  const double c1 = 1e-4;    //< Sufficient decrease parameter
  const double rho_backtrack = 0.5;   ///< Backtracking factor
  const int max_ls_iter = 20; ///< Maximum number of backtracking steps

  // Gauss-Newton least squares solver
  Mat GN_LeastSquares_MatShell; ///< MatShell for the Gauss-Newton least-squares problem
  KSP ksp_LeastSquares;
  Vec *wsub_workspace; ///< Pre-allocated workspace vectors for MatMultTranspose (size: ninit_local)
  PetscInt state_dim_cached; ///< Cached state dimension for efficiency

  // Tao BRGN solver for least-squares problem
  Tao tao_brgn;        ///< Tao BRGN solver (alternative to KSPLSQR)
  Vec brgn_rhs;        ///< Cached RHS vector for BRGN residual evaluation
  Vec brgn_residual;   ///< Cached Residual vector for BRGN solver
  double gn_leastsquares_brgn_damping; ///< Damping parameter λ for BRGN regularization
  std::string gn_leastsquares_solver; ///< Least-squares solver name: "brgn" or "lsqr"

  // KSP linear solver
  KSP ksp_GN_NormalEq;  ///< Linear solver for Gauss-Newton Normal Equation solver
  double gn_normaleq_damping; ///< Damping parameter for Gauss-Newton matrix shift 

  EPS eps_GN_NormalEq; // EPS solver for the Gauss-Newton Normal Equation matrix
  PetscReal eps_tol = 1e-4; ///< Tolerance for EPS eigenvalue solver
  double eps_evals_cutoff = 1e-5; ///< Cutoff for eigenvalues of the Gauss-Newton matrix
  double eps_normaleq_damping = 1e-3; ///< Damping for eigenvalues of the Gauss-Newton matrix
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
  int getMaxIter()     { return maxiter; };
  OptimTarget* getOptimTarget() { return optim_target; };
  Mat getGN_NormalEq_MatShell() { return GN_NormalEq_MatShell; };
  Mat getGN_NormalEq_MatDense() { return GN_NormalEq_MatDense; };
  Mat getGN_LeastSquares_MatShell() { return GN_LeastSquares_MatShell; };
  bool getQuietmode() { return quietmode; };

  int getOutputOptimizationStride() { return output_optimization_stride; };
  Output* getOutput() { return output; };
  TimeStepper* getTimeStepper() { return timestepper; };

  void setXevalGN(const Vec x){ VecCopy(x, xeval_GN); }

  /**
   * @brief Override the maximum number of iterations for all Gauss-Newton solvers.
   *
   * This sets the maximum iterations for:
   * - KSP solver (used in solveGaussNewtonNormalEqKSP)
   * - Least Squares solver (used in solveGaussNewtonLeastSquares)
   * - EPS eigenvalue solver (used in solveGaussNewtonNormalEqEPS)
   *
   * @param maxiter Maximum number of iterations
   */
  void setGaussNewtonMaxiter(int maxiter);

  /**
   * @brief Get the current maximum number of iterations for Gauss-Newton solvers.
   *
   * @return int Current maximum number of iterations (from KSP solver)
   */
  int getGaussNewtonMaxiter();

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
   * @param v_LeastSquares Solution vector to store the result
   */
  void solveGaussNewtonLeastSquares(const Vec xinit, const Vec initial_guess, Vec v_LeastSquares);


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

  // For the Gauss-Newton least-squares MatShell operations
  static void GN_LeastSquares_MatMult(Mat A, Vec v, Vec y); // Apply y = W^1/2Lv
  static void GN_LeastSquares_MatMultTranspose(Mat A, Vec w, Vec vout); // Apply vout = L* W^1/2 w
  static void GN_LeastSquares_MatCreateVecs(Mat A, Vec *right, Vec *left);

  /**
   * @brief MatMult operation for MatShell Gauss-Newton Normal Equation L^*Lv: Linearized forward + adjoint operator. 
   * 
   * The point of evaluation xeval_GN must be set correctly in the OptimProblem before calling this.
   * 
   * @param[in] v Direction vector
   * @param[out] Av Resulting vector after applying the linearized forward and adjoint operators
   */
  static void GN_NormalEq_MatMult(Mat A, const Vec v, Vec Av);

  /**
   * @brief Updates the dense full  Gauss-Newton matrix based on the current point of evaluation xeval_GN by calling GN_NormalEqShell_MatMult on each unit vectoor.
   * @note This operation can be expensive as it involves N^2-1 applications of the MatShell.
   */
  void updateGaussNewtonMatDense();

  /**
   * @brief Solves the Gauss-Newton normal equation L^* L(x) v = b for v using CG iterations
   * 
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param b Right-hand side vector
   * @param Ainv_b Solution vector to store the result
   */
  void solveGaussNewtonNormalEqKSP(Vec xinit, const Vec initial_guess, const Vec b, Vec Ainv_b);

  /**
   * @brief Solves the Gauss-Newton normal equation L^* L(x) v = b via eigenvalue decomposition
   * 
   * @param xinit Point of evaluation for the Gauss-Newton matrix
   * @param b Right-hand side vector
   * @param Ainv_b Solution vector to store the result
   */
  void solveGaussNewtonNormalEqEPS(Vec xinit, const Vec b, Vec Ainv_b);

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
  std::vector<double> computeGaussNewtonNormalEqEvals(Vec xinit, Mat* evecs_out);

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
 * @brief PETSc TAO interface routine for objective function evaluation.
 *
 * @param tao TAO solver object
 * @param x Design vector
 * @param f Pointer to store objective function value
 * @param ptr Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode TaoEvalObjective(Tao tao, Vec x, PetscReal *f, void*ptr);

/**
 * @brief PETSc TAO interface routine for gradient evaluation.
 *
 * @param tao TAO solver object
 * @param x Design vector
 * @param G Gradient vector
 * @param ptr Pointer to user context (OptimProblem instance)
 * @return PetscErrorCode Error code
 */
PetscErrorCode TaoEvalGradient(Tao tao, Vec x, Vec G, void*ptr);

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
