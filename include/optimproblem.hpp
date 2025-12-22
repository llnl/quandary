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
  Vec rho_t0; ///< Storage for initial condition of the ODE
  Vec rho_t0_bar; ///< Storage for adjoint initial condition of the adjoint ODE (aka the terminal condition)
  std::vector<Vec> store_finalstates; ///< Storage for final states for each initial condition

  OptimTarget* optim_target; ///< Pointer to the optimization target (gate or state)

  MPI_Comm comm_init; ///< MPI communicator for initial condition parallelization
  MPI_Comm comm_optim; ///< MPI communicator for optimization parallelization, currently not used (size 1)
  int mpirank_optim, mpisize_optim; ///< MPI rank and size for optimization communicator
  int mpirank_petsc, mpisize_petsc; ///< MPI rank and size for spatial parallelization (PETSc)
  int mpirank_world, mpisize_world; ///< MPI rank and size for global communicator
  int mpirank_init, mpisize_init; ///< MPI rank and size for initial condition communicator

  bool quietmode; ///< Flag for quiet mode operation

  std::vector<double> obj_weights; ///< Weights for averaging objective over initial conditions
  int ndesign; ///< Number of global design (optimization) parameters
  double objective; ///< Current objective function value (sum over final-time cost, regularization terms and penalty terms)
  double obj_cost; ///< Final-time measure J(T) in objective
  double obj_regul; ///< Regularization term in objective
  double obj_penal_leakage; ///< Penalty term for leakage into guard levels
  double obj_penal_weightedcost; ///< Penalty term for weighted running cost 
  double obj_penal_dpdm; ///< Penalty term second-order state derivatives (penalizes variations of the state evolution)
  double obj_penal_variation; ///< Penalty term for variation of control parameters
  double obj_penal_energy; ///< Energy penalty term in objective
  double fidelity; ///< Final-time fidelity: 1/ninit sum_i Tr(rho_target^dag rho(T)) for Lindblad, |1/ninit sum_i phi_target^dag phi|^2 for Schrodinger
  double gnorm; ///< Current norm of gradient
  double gamma_tikhonov; ///< Parameter for Tikhonov regularization
  bool tikhonov_use_x0; ///< Switch to use ||x - x0||^2 for Tikhonov regularization instead of ||x||^2
  double gamma_penalty_leakage; ///< Parameter multiplying integral leakage term
  double gamma_penalty_weightedcost; ///< Parameter multiplying integral weighted cost function 
  double gamma_penalty_dpdm; ///< Parameter multiplying integral penalty term for 2nd derivative of state variation
  double gamma_penalty_energy; ///< Parameter multiplying energy penalty
  double gamma_penalty_variation; ///< Parameter multiplying finite-difference squared regularization term
  double tol_grad_abs; ///< Stopping criterion based on absolute gradient norm
  double tol_grad_rel; ///< Stopping criterion based on relative gradient norm
  double tol_finalcost; ///< Stopping criterion based on objective function value
  double tol_infidelity; ///< Stopping criterion based on infidelity
  int maxiter; ///< Stopping criterion based on maximum number of iterations
  Tao tao; ///< PETSc's TAO optimization solver
  std::vector<double> initguess_fromfile; ///< Initial guess read from file
  double* mygrad; ///< Auxiliary gradient storage
    
  Vec xtmp; ///< Temporary vector storage
  
  public: 
    Output* output; ///< Pointer to output handler
    TimeStepper* timestepper; ///< Pointer to time-stepping scheme
    Vec xlower, xupper; ///< Lower and upper bounds for optimization variables
    Vec xprev; ///< Design vector at previous iteration
    Vec xinit; ///< Initial design vector

  /**
   * @brief Constructor for optimization problem.
   *
   * @param config Configuration parameters from input file
   * @param timestepper_ Pointer to time-stepping scheme
   * @param comm_init_ MPI communicator for initial condition parallelization
   * @param comm_optim MPI communicator for optimization parallelization
   * @param output_ Pointer to output handler
   * @param quietmode Flag for quiet operation (default: false)
   */
  OptimProblem(const Config& config, TimeStepper* timestepper_, MPI_Comm comm_init_, MPI_Comm comm_optim, Output* output_, bool quietmode=false);

  ~OptimProblem();

  /**
   * @brief Retrieves the number of design variables.
   *
   * @return int Number of optimization parameters
   */
  int getNdesign(){ return ndesign; };

  /**
   * @brief Retrieves the current objective function value.
   *
   * @return double Total objective function value
   */
  double getObjective(){ return objective; };

  /**
   * @brief Retrieves the final-time cost term.
   *
   * @return double Final-time cost J(T)
   */
  double getCostT()    { return obj_cost; };

  /**
   * @brief Retrieves the Tikhonov regularization term.
   *
   * @return double Tikhonov regularization contribution
   */
  double getRegul()    { return obj_regul; };

  /**
   * @brief Retrieves the integral term for leakage.
   *
   * @return double Penalty term contribution
   */
  double getPenaltyLeakage()  { return obj_penal_leakage; };

  /**
   * @brief Retrieves the integral penalty term for weighted running cost.
   *
   * @return double Penalty term contribution
   */
  double getPenaltyWeightedCost()  { return obj_penal_weightedcost; };



  /**
   * @brief Retrieves the second-order state derivative penalty term.
   *
   * @return double Second-order penalty contribution
   */
  double getPenaltyDpDm()  { return obj_penal_dpdm; };

  /**
   * @brief Retrieves the control variation penalty term.
   *
   * @return double Control variation penalty contribution
   */
  double getPenaltyVariation()  { return obj_penal_variation; };

  /**
   * @brief Retrieves the control energy penalty term.
   *
   * @return double Energy penalty contribution
   */
  double getPenaltyEnergy()  { return obj_penal_energy; };

  /**
   * @brief Retrieves the current fidelity.
   *
   * @return double Quantum fidelity measure
   */
  double getFidelity() { return fidelity; };

  /**
   * @brief Retrieves the objective function tolerance.
   *
   * @return double Absolute tolerance for objective function convergence
   */
  double getTolFinalCost()    { return tol_finalcost; };

  /**
   * @brief Retrieves the gradient tolerance.
   *
   * @return double Absolute tolerance for gradient norm convergence
   */
  double getTolGradAbs()    { return tol_grad_abs; };

  /**
   * @brief Retrieves the infidelity tolerance.
   *
   * @return double Tolerance for infidelity convergence
   */
  double getTolInfidelity()   { return tol_infidelity; };

  /**
   * @brief Retrieves the MPI rank in the world communicator.
   *
   * @return int MPI rank
   */
  int getMPIrank_world() { return mpirank_world;};

  /**
   * @brief Retrieves the maximum number of iterations.
   *
   * @return int Maximum iteration limit
   */
  int getMaxIter()     { return maxiter; };

  /**
   * @brief Evaluates the objective function F(x).
   * 
   * Performs forward simulations for each initial conditions and
   * evaluates the objective function. 
   *
   * @param x Design vector
   * @return double Objective function value
   */
  double evalF(const Vec x);

  /**
   * @brief Evaluates the gradient of the objective function with respect to the control parameters
   *
   * @param x Design (optimization) vector
   * @param G Gradient vector to store result
   */
  void evalGradF(const Vec x, Vec G);

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
