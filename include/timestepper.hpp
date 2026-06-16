#include <petsc/private/tsimpl.h>        /*I "petscts.h"  I*/
#include <petscts.h>
#include <petscksp.h>
#include "mastereq.hpp"
#include <assert.h> 
#include <iostream> 
#include "defs.hpp"
#include "output.hpp"
#include "optimtarget.hpp"
#include <deque>
#pragma once

/**
 * @brief Base class for time integration schemes to evolve the quantum dynamics.
 *
 * This abstract class provides the interface for time-stepping methods used to
 * integrate the quantum evolution equations (Lindblad master equation or Schroedinger
 * equation). It supports both forward and adjoint time integration, handles output of evolution data, 
 * and evaluates penalty integral terms.
 * 
 * Main functionality:
 *    - @ref solveODE propagates an initial state at time t0 to the final time T, while writing evolution data to
 *      files and adding to integral penalty terms, if needed. 
 *    - @ref solveAdjointODE propagates a terminal (adjoint) condition through the time domain backwards in time 
 *      from T to 0, while updating the reduced gradient along the way. 
 * 
 * This class contains references to:
 *    - @ref MasterEq for evaluating and applying the right-hand-side system matrix of the real-valued, vectorized
 *      differential equation to a state vector at each time-step
 *    - @ref OptimTarget for evaluating integral penalty terms at each time step
 *    - @ref Output for writing evolution data to output files at each time step
 */
class TimeStepper{
  protected:
    Vec x; ///< Auxiliary vector for forward time stepping
    Vec xadj; ///< Auxiliary vector needed for adjoint (backward) time stepping
    Vec xprimal; ///< Auxiliary vector for backward time stepping
    std::vector<std::vector<Vec>> trajectory_states; ///< Storage for primal states during forward evolution, one trajectory for each local initial condition..
    std::vector<Vec> final_states; ///< Storage for final states for each local initial condition. Always filled after solveODE.
    std::vector<Vec> dpdm_states; ///< Storage for states needed for second-order derivative penalty
    int mpirank_world; ///< MPI rank in global communicator
    int mpisize_petsc; ///< MPI size in Petsc communicator
    int mpirank_petsc; ///< MPI rank in Petsc communicator
    PetscInt localsize_u; ///< Size of local sub vector u or v in state x=[u,v]
    PetscInt ilow; ///< First index of the local sub vector u,v
    PetscInt iupp; ///< Last index (+1) of the local sub vector u,v

    bool eval_leakage; ///< Flag to compute leakage integral term
    bool eval_energy; ///< Flag to compute energy integral term
    bool eval_dpdm; ///< Flag to compute second-order derivative integral term
    bool eval_weightedcost; ///< Flag to compute weighted cost integral term
    double leakage_integral=0.0; ///< Sums the integral over leakage 
    double energy_integral=0.0; ///< Sums the energy term
    double dpdm_integral=0.0; ///< Sums second-order derivative variation value
    double weightedcost_integral=0.0; ///< Sums the integral over weighted cost function
    double weightedcost_width; ///< Width parameter for weighted cost function

    int ntime; ///< Number of time steps
    double total_time; ///< Final evolution time
    double dt; ///< Time step size
    bool writeTrajectoryDataFiles;  ///< Flag to determine whether or not trajectory data will be written to files during forward simulation */

    Vec redgrad; ///< Reduced gradient vector for optimization
    OptimTarget* optim_target; ///< Pointer to optimization target specification
    Output* output; ///< Pointer to output handler
    MasterEq* mastereq; ///< Pointer to master equation solver

  public: 
    TimeStepper(); 

    /**
     * @brief Constructor for time stepper.
     *
     * @param config Configuration parameters from input file
     * @param mastereq_ Pointer to master equation solver
     * @param output_ Pointer to output handler
     * @param ninit_local Number of initial conditions local to this processor
     */
    TimeStepper(const Config& config, MasterEq* mastereq_, Output* output_, int ninit_local); 

    virtual ~TimeStepper(); 

    double getLeakageIntegral(){ return leakage_integral; };
    double getWeightedCostIntegral(){ return weightedcost_integral; };
    double getEnergyIntegral(){ return energy_integral; };
    double getDPDMIntegral(){ return dpdm_integral; };
    Vec getReducedGradient(){ return redgrad; };
    void setWriteTrajectoryDataFiles(bool write){ writeTrajectoryDataFiles = write; };

    void setOptimTarget(OptimTarget* optim_target_){ optim_target = optim_target_; };

    /**
     * @brief Retrieves the final state for a specific local initial condition.
     *
     * @param iinit_local Local index of the initial condition
     * @return Vec Final state vector for the specified initial condition
     */
    Vec getFinalState(size_t iinit_local){ return final_states[iinit_local]; }

    /**
     * @brief Get the smallest timestep size used during timestepping.
     *
     * For adaptive timestepping methods, returns the minimum timestep size chosen.
     * For fixed timestep methods, returns the regular timestep size (dt).
     *
     * @return double Smallest timestep size
     */
    virtual double getMinTimestepSize() const { return dt; }

    /**
     * @brief Solves the ODE forward in time.
     * 
     * This performs the time-stepping to propagate an initial condition to the final time.
     *
     * @param initid Initial condition identifier
     * @param iinit_local Local index of initial condition for this processor
     * @param rho_t0 Initial state vector
     * @return Vec Final state vector at time T
     */
    virtual Vec solveODE(int initid, int iinit_local, Vec rho_t0);

    /**
     * @brief Solves the adjoint ODE backward in time.
     * 
     * This performs backward time-stepping to backpropagate an adjoint initial condition at 
     * final time (aka a terminal condtion) to time t=0, while accumulating the reduced gradient. 
     *
     * @param iinit_local Local index of initial condition for this processor
     * @param rho_t0_bar Terminal condition for adjoint state
     * @param Jbar_leakage Adjoint of leakage integral term
     * @param Jbar_weightedcost Adjoint of weighted cost integral term
     * @param Jbar_dpdm Adjoint of second-order derivative variation
     * @param Jbar_energy Adjoint of energy integral term
     */
    virtual void solveAdjointODE(int iinit_local, Vec rho_t0_bar, double Jbar_leakage, double Jbar_weightedcost, double Jbar_dpdm, double Jbar_energy);

    /**
     * @brief Evaluates leakage into guard levels 
     *
     * @param x Current state vector
     * @return double Leakage term value / ntime
     */
    double evalLeakage(const Vec x);

    /**
     * @brief Computes derivative of leakage term.
     *
     * @param x Current state vector
     * @param xbar Adjoint state vector to update
     * @param Jbar Adjoint of leakage integral term
     */
    void evalLeakage_diff(const Vec x, Vec xbar, double Jbar);

    /**
     * @brief Evaluates the weighted cost function term.
     *
     * @param time Current time
     * @param x Current state vector
     * @return double Weighted cost integral term value
     */
    double evalWeightedCost(double time, const Vec x);

    /**
     * @brief Computes derivative of weighted cost function term.
     *
     * @param time Current time
     * @param x Current state vector
     * @param xbar Adjoint state vector to update
     * @param Jbar Adjoint of weighted cost term
     */
    void evalWeightedCost_diff(double time, const Vec x, Vec xbar, double Jbar);

    /**
     * @brief Evaluates second-order derivative variation for the state.
     *
     * @param x Current state vector
     * @param xm1 State vector at previous time step
     * @param xm2 State vector at two time steps ago
     * @return double Second-order variation value
     */
    double evalDpDm(Vec x, Vec xm1, Vec xm2);

    /**
     * @brief Computes derivative of second-order derivative variation term.
     *
     * @param n Time step index
     * @param xbar Adjoint state vector to update
     * @param Jbar Adjoint of dpdm term
     */
    void evalDpDm_diff(int n, Vec xbar, double Jbar);
    
    /**
     * @brief Evaluates energy integral term.
     *
     * @param time Current time
     * @return double Energy value
     */
    double evalEnergy(double time);

    /**
     * @brief Computes derivative of energy integral.
     *
     * @param time Current time
     * @param Jbar Adjoint of energy 
     * @param redgrad Reduced gradient vector to update
     */
    void evalEnergy_diff(double time, double Jbar, Vec redgrad);

    /**
     * @brief Evolves state forward by one time-step from tstart to tstop.
     *
     * Pure virtual function to be implemented by the derived time-stepping classes.
     *
     * @param tstart Start time
     * @param tstop Stop time
     * @param x State vector to evolve
     */
    virtual void evolveFWD(const double tstart, const double tstop, Vec x) = 0;

    /**
     * @brief Evolves adjoint state backward by one time-step and updates reduced gradient.
     * 
     * Abstract base-class implementation is empty. Derived classes that need backward time-stepping should implement this function.
     *
     * @param tstart Start time (backward evolution)
     * @param tstop Stop time (backward evolution)
     * @param x_stop State at stop time
     * @param x_adj Adjoint state vector
     * @param grad Gradient vector to update
     * @param compute_gradient Flag to compute gradient
     */
    virtual void evolveBWD(const double tstart, const double tstop, const Vec x_stop, Vec x_adj, Vec grad, bool compute_gradient);
};

/**
 * @brief Explicit Euler time integration scheme.
 *
 * First-order explicit time stepping method. Simple, requires small time steps
 * for stability. Mainly used for testing and comparison purposes.
 */
class ExplEuler : public TimeStepper {
  protected:
  Vec stage; ///< Intermediate vector
  public:
    /**
     * @brief Constructor for explicit Euler scheme.
     *
     * @param config Configuration parameters from input file
     * @param mastereq_ Pointer to master equation solver
     * @param output_ Pointer to output handler
     * @param ninit_local Number of initial conditions local to this processor
     */
    ExplEuler(const Config& config, MasterEq* mastereq_, Output* output_, int ninit_local);

    ~ExplEuler();

    /**
     * @brief Evolves state forward using explicit Euler method.
     *
     * @param tstart Start time
     * @param tstop Stop time
     * @param x State vector to evolve
     */
    void evolveFWD(const double tstart, const double tstop, Vec x);

    /**
     * @brief Evolves adjoint backward using explicit Euler method.
     *
     * @param tstart Start time (backward evolution)
     * @param tstop Stop time (backward evolution)
     * @param x_stop State at stop time
     * @param x_adj Adjoint state vector
     * @param grad Gradient vector to update
     * @param compute_gradient Flag to compute gradient
     */
    void evolveBWD(const double tstart, const double tstop, const Vec x_stop, Vec x_adj, Vec grad, bool compute_gradient);
};

/**
 * @brief Implicit midpoint rule time integration scheme.
 *
 * Second-order implicit method with symplectic properties. This is the default and recommended
 * time integration, due to its stability and energy conservation properties.
 *
 * Runge-Kutta tableau:
 * @code
 * 1/2 | 1/2
 * ----------
 *     |  1
 * @endcode
 */
class ImplMidpoint : public TimeStepper {
  protected:
  Vec stage, stage_adj; ///< Intermediate stage vectors for forward and adjoint
  Vec rhs, rhs_adj; ///< Right-hand side vectors for forward and adjoint
  KSP ksp; ///< PETSc's linear solver context for GMRES
  PC  preconditioner; ///< Preconditioner for linear solver
  LinearSolverType linsolve_type; ///< Linear solver type (GMRES or NEUMANN)
  int linsolve_maxiter; ///< Maximum number of linear solver iterations
  double linsolve_abstol; ///< Absolute tolerance for linear solver
  double linsolve_reltol; ///< Relative tolerance for linear solver
  int linsolve_iterstaken_avg; ///< Average number of linear solver iterations
  double linsolve_error_avg; ///< Average error of linear solver
  int linsolve_counter; ///< Counter for linear solve calls
  Vec tmp, err; ///< Auxiliary vectors for Neumann iterations

  public:
    /**
     * @brief Constructor for implicit midpoint scheme.
     *
     * @param config Configuration parameters from input file
     * @param mastereq_ Pointer to master equation solver
     * @param output_ Pointer to output handler
     * @param ninit_local Number of initial conditions local to this processor
     */
    ImplMidpoint(const Config& config, MasterEq* mastereq_, Output* output_, int ninit_local);

    ~ImplMidpoint();

    /**
     * @brief Evolves state forward using implicit midpoint rule.
     *
     * @param tstart Start time
     * @param tstop Stop time
     * @param x State vector to evolve
     */
    virtual void evolveFWD(const double tstart, const double tstop, Vec x);

    /**
     * @brief Evolves adjoint backward using implicit midpoint rule and adds to reduced gradient.
     *
     * @param tstart Start time (backward evolution)
     * @param tstop Stop time (backward evolution)
     * @param x_stop State at stop time
     * @param x_adj Adjoint state vector
     * @param grad Gradient vector to update
     * @param compute_gradient Flag to compute gradient
     */
    virtual void evolveBWD(const double tstart, const double tstop, const Vec x_stop, Vec x_adj, Vec grad, bool compute_gradient);

    /**
     * @brief Solves (I - alpha*A) * x = b using Neumann iterations.
     *
     * @param A Matrix A
     * @param b Right-hand side vector
     * @param x Solution vector
     * @param alpha Scaling parameter
     * @param transpose Flag to solve transposed system (I - alpha*A^T)*x = b
     * @return int Number of iterations taken
     */
    int NeumannSolve(Mat A, Vec b, Vec x, double alpha, bool transpose);
};

/**
 * @brief Compositional implicit midpoint rule for higher-order accuracy.
 *
 * Extends the implicit midpoint rule to higher order by using multiple
 * substeps with specific coefficients. Maintains symplectic properties
 * while achieving better accuracy for larger time steps.
 */
class CompositionalImplMidpoint : public ImplMidpoint {
  protected:

  std::vector<double> gamma; ///< Coefficients for compositional step sizes
  std::vector<Vec> x_stage; ///< Storage for primal states at intermediate stages
  Vec aux; ///< Auxiliary vector
  int order; ///< Order of the compositional method

  public:
    /**
     * @brief Constructor for compositional implicit midpoint scheme.
     *
     * @param config Configuration parameters from input file
     * @param mastereq_ Pointer to master equation solver
     * @param output_ Pointer to output handler
     * @param ninit_local Number of initial conditions local to this processor
     * @param order_ Order of the compositional method
     */
    CompositionalImplMidpoint(const Config& config, MasterEq* mastereq_, Output* output_, int ninit_local, int order_);

    ~CompositionalImplMidpoint();

    /**
     * @brief Evolves state forward using compositional implicit midpoint rule.
     *
     * @param tstart Start time
     * @param tstop Stop time
     * @param x State vector to evolve
     */
    void evolveFWD(const double tstart, const double tstop, Vec x);

    /**
     * @brief Evolves adjoint backward using compositional implicit midpoint rule and accumulates gradient.
     *
     * @param tstart Start time (backward evolution)
     * @param tstop Stop time (backward evolution)
     * @param x_stop State at stop time
     * @param x_adj Adjoint state vector
     * @param grad Gradient vector to update
     * @param compute_gradient Flag to compute gradient
     */
    void evolveBWD(const double tstart, const double tstop, const Vec x_stop, Vec x_adj, Vec grad, bool compute_gradient);
};


/**
 * @brief Petsc's adaptive timestepper
 */
class PetscTS : public TimeStepper {
  protected:
    TS ts;      ///< Backward-compatible alias to ts_pool[0].
    std::vector<TS> ts_pool; ///< One PETSc TS per initial condition.
    std::vector<Vec> q_pool; ///< One quadrature state vector per TS/initial condition.
    Vec redgrad_ts; ///< TS-internal gradient vector with PETSc communicator-compatible layout.

    Mat dRHSdp; ///< MatShell for applying derivative of the RHS to control parameters.
    Mat dIntegralCostdY; ///< Dense matrix for derivative of integral costs wrt state.
    Mat dIntegralCostdP; ///< Dense matrix for derivative of integral costs wrt parameters.
    PetscReal dRHSdp_time; ///< Cached time used by the RHSJacobianP MatShell callbacks.
    double adj_scale_leakage; ///< Per-solve scaling for leakage integral adjoint contribution.
    double adj_scale_weightedcost; ///< Per-solve scaling for weighted-cost integral adjoint contribution.
    double adj_scale_energy; ///< Per-solve scaling for energy integral adjoint contribution.
    double min_timestep_size; ///< Smallest timestep size chosen during adaptive timestepping

  public:
    /**
     * @brief Constructor for PetscTS.
     *
     * @param config Configuration parameters from input file
     * @param mastereq_ Pointer to master equation solver
     * @param output_ Pointer to output handler
     * @param ninit_local Number of initial conditions local to this processor
     */
    PetscTS(const Config& config, MasterEq* mastereq_, Output* output_, int ninit_local);
    ~PetscTS();

    /**
     * @brief Overwrites the default time-stepping by calling PETSc's TSSolve.
     * 
     * @param initid Initial condition identifier
     * @param iinit_local Local index of initial condition for this processor
     * @param rho_t0 Initial state vector
     * @return Vec Final state vector at time T
     */
    Vec solveODE(int initid, int iinit_local, Vec rho_t0) override;

    /** 
     * @brief Overwrites the default adjoint time-stepping by calling PETSc's TSSolve on the adjoint system.
     * 
     * @param iinit_local Local index of initial condition for this processor
     * @param rho_t0_bar Terminal condition for adjoint state
     * @param Jbar_leakage Adjoint seed of leakage integral term
     * @param Jbar_weightedcost Adjoint seed of weighted cost integral term
     * @param Jbar_dpdm Adjoint seed of second-order derivative variation term
     * @param Jbar_energy Adjoint seed of energy integral term
     */
    void solveAdjointODE(int iinit_local, Vec rho_t0_bar, double Jbar_leakage, double Jbar_weightedcost, double Jbar_dpdm, double Jbar_energy) override;

    /** 
     * @brief PETSc callback to update the RHS matrix for the forward system at time t.
      *
      * This is called by PETSc's TS at each time step to get the current system matrix. We use it to update the MasterEq's RHS MatShell with the current time and assemble it, which is needed for time-dependent Hamiltonian terms.
      * @param ts PETSc TS context
      * @param t Current time
      * @param x Current state vector (not used here, but required by PETSc's interface)
      * @param A Matrix to store the RHS (Jacobian) for the forward system
      * @param B Matrix to store the mass matrix (not used here, but required by PETSc's interface)
      * @param ptr Pointer to PetscTS object (self)
      * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode RHSMatrixUpdate(TS ts, PetscReal t, Vec x, Mat A, Mat B, void *ptr);

    /**
     * @brief PETSc callback to update the matrix for the derivative of the RHS with respect to control parameters at time t.
     * 
     * @param ts PETSc TS context
     * @param t Current time
     * @param x Current state vector (not used here, but required by PETSc's interface)
     * @param A Matrix to store the derivative of the RHS with respect to parameters
     * @param ptr Pointer to PetscTS object (self)
     * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode dRHSdpMatrixUpdate(TS ts, PetscReal t, Vec x, Mat A, void *ptr);

    /** 
     * @brief PETSc callback to compute the action of the derivative of the RHS with respect to parameters on a vector, used during adjoint solves.
     *  
     * Computes y = (dRHS/dp)^T x via MasterEq::compute_dRHS_dParams.
     * 
     * @param A Matrix representing the derivative of the RHS with respect to parameters
     * @param x Input vector
     * @param y Output vector
     * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode computedRHSdp(Mat A, Vec x, Vec y);

    /** 
     * @brief PETSc callback to evaluate integral cost functions at time t for the current state x.
     * 
     * @param ts PETSc TS context
     * @param t Current time
     * @param x Current state vector
     * @param F Vector to store the evaluated integral cost function values
     * @param ctx Pointer to PetscTS object (self)
     * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode IntegralCosts(TS ts, PetscReal t, Vec x, Vec F, void *ctx);

    /** 
     * @brief PETSc callback to update the matrix for the derivative of integral cost functions with respect to the state at time t.
     * 
     * @param ts PETSc TS context
     * @param t Current time
     * @param x Current state vector
     * @param A Matrix to store the derivative of integral cost functions with respect to the state
     * @param B Matrix to store additional information (not used here, but required by PETSc's interface)
     * @param ptr Pointer to PetscTS object (self)
     * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode dIntegralCostdYUpdate(TS ts, PetscReal t, Vec x, Mat A, Mat B, void *ptr);

    /** 
     * @brief PETSc callback to update the matrix for the derivative of integral cost functions with respect to the parameters at time t.
     * 
     * @param ts PETSc TS context
     * @param t Current time
     * @param x Current state vector
     * @param A Matrix to store the derivative of integral cost functions with respect to the parameters
     * @param ptr Pointer to PetscTS object (self)
     * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode dIntegralCostdPUpdate(TS ts, PetscReal t, Vec x, Mat A, void *ptr);

    // Callback function during TSSolve to evaluate trajectory data at each accepted time step 
    /**
     * @brief PETSc callback to monitor the trajectory during time-stepping and write data to output files.
      *
      * This is called by PETSc's TS at each accepted time step. We use it to write trajectory data to output files via the Output handler, and also to track the minimum timestep size chosen by the adaptive time-stepping.
      * @param ts PETSc TS context
      * @param step Current time step index
      * @param time Current time
      * @param state Current state vector
      * @param ctx Pointer to PetscTS object (self) 
      * @return PetscErrorCode PETSc error code
     */
    static PetscErrorCode monitorTrajectory(TS ts, PetscInt step, PetscReal time, Vec state, void *ctx);

    // THESE ARE NOT USED. Instead, solveODE and solveAdjointODE overwrites the default time-stepping by calling TSSolve. 
    void evolveFWD(const double, const double, Vec) override {};
    void evolveBWD(const double, const double, const Vec, Vec, Vec, bool) override {};

    /**
     * @brief Get the smallest timestep size chosen during adaptive timestepping.
     *
     * @return double Smallest timestep size from the last forward solve
     */
    double getMinTimestepSize() const override { return min_timestep_size; }
};