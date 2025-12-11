#include "defs.hpp"
#include "gate.hpp"
#include "mastereq.hpp"
#pragma once

/**
 * @brief Optimization target specification for quantum control.
 *
 * This class manages the target specification for quantum optimal control problems,
 * including gate optimization and pure state preparation. It handles target and initial state
 * preparation and the evaluation of the final-time objective function measure.
 * 
 * Main functionality: 
 *    - @ref prepareInitialState prepares and returns the states at time t=0, depending on what the optimization 
 *       target is
 *    - @ref prepareTargetState prepares and stores the corresponding target state for this initial state
 *    - @ref evalJ for computing the final-time objective function measure
 * 
 * This class contains references to:
 *    - @ref Gate for evaluating the target state for quantum gate optimization
 */
class OptimTarget{
    protected:

    PetscInt dim; ///< State dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
    PetscInt dim_rho; ///< Dimension of Hilbert space = N
    PetscInt dim_ess; ///< Dimension of essential level system = N_e
    int noscillators; ///< Number of oscillators in the system
 
    TargetType target_type; ///< Type of optimization target (pure state preparation or gate optimization)
    ObjectiveType objective_type; ///< Type of objective function measure (Frobenius, trace, pure-state measure)
    Gate *targetgate; ///< Pointer to target gate (if gate optimization)
    double purity_rho0; ///< Purity of initial state Tr(rho(0)^2)
    PetscInt purestateID; ///< For pure state preparation: integer m for preparing the target state \f$ e_m e_m^{\dagger}\f$
    std::string target_filename; ///< Filename if target state is read from file
    Vec targetstate; ///< Storage for the target state vector (NULL for pure states, \f$V\rho V^\dagger\f$ for gates, density matrix from file)
    InitialConditionType initcond_type; ///< Type of initial conditions
    std::vector<size_t> initcond_IDs; ///< Integer list for pure-state initialization
    LindbladType lindbladtype; ///< Type of Lindblad decoherence operators, or NONE for Schroedinger solver
    int mpisize_petsc; ///< Size of PETSc communicator
    int mpirank_petsc; ///< Rank of PETSc communicator
    PetscInt localsize_u; ///< Size of local sub vector u or v in state x=[u,v]
    PetscInt ilow; ///< First index of the local sub vector u,v
    PetscInt iupp; ///< Last index (+1) of the local sub vector u,v

    Vec aux; ///< Auxiliary vector for gate optimization objective computation
    bool quietmode; ///< Flag for quiet mode operation

  public:
    OptimTarget();

    /**
     * @brief Constructor with full target specification.
     *
     * @param target_str Vector of strings specifying the target
     * @param objective_str String specifying the objective function type
     * @param initcond_str Vector of strings specifying initial conditions
     * @param mastereq Pointer to master equation solver
     * @param total_time Total evolution time
     * @param read_gate_rot Gate rotation parameters
     * @param rho_t0 Initial state vector
     * @param quietmode_ Flag for quiet operation
     */
    OptimTarget(std::vector<std::string> target_str, const std::string& objective_str, std::vector<std::string> initcond_str, MasterEq* mastereq, double total_time, std::vector<double> read_gate_rot, Vec rho_t0, bool quietmode_);

    ~OptimTarget();

    /**
     * @brief Retrieves the target type.
     *
     * @return TargetType Type of optimization target
     */
    TargetType getTargetType(){ return target_type; };

    /**
     * @brief Retrieves the objective function type.
     *
     * @return ObjectiveType Type of objective function
     */
    ObjectiveType getObjectiveType(){ return objective_type; };

    /**
     * @brief Prepares the initial condition state.
     *
     * @param iinit Index in processor range [rank * ninit_local .. (rank+1) * ninit_local - 1]
     * @param ninit Total number of initial conditions
     * @param nlevels Number of levels per oscillator
     * @param nessential Number of essential levels per oscillator
     * @param rho0 Vector to store the initial condition
     * @return int Identifier for this initial condition (element number in matrix vectorization)
     */
    int prepareInitialState(const int iinit, const int ninit, const std::vector<size_t>& nlevels, const std::vector<size_t>& nessential,  Vec rho0);

    /**
     * @brief Prepares the target state for gate optimization.
     *
     * For gate optimization, computes the rotated target state \f$V \rho V^{\dagger}\f$
     * for a given initial state \f$\rho\f$ and stores it locally as a class member. 
     * Also stores the purity of \f$rho\f$ needed for scaling the Hilbert-Schmidt 
     * overlap in the trace objective function.
     *
     * @param rho Initial state vector
     */
    void prepareTargetState(const Vec rho);

    /**
     * @brief Evaluates the final-time objective function measure \f$J(\rho(T))\f$.
     *
     * The target state must be prepared and stored before calling this function.
     * Returns both real and imaginary parts of the final-time measure. The imaginary part
     * is generally zero except for Schroedinger solver with the trace objective 
     * function measure.
     *
     * @param[in] state Current state vector
     * @param[out] J_re_ptr Pointer to store real part of objective
     * @param[out] J_im_ptr Pointer to store imaginary part of objective
     */
    void evalJ(const Vec state, double* J_re_ptr, double* J_im_ptr);

    /**
     * @brief Computes derivative of the final-time objective function measure.
     *
     * Updates the adjoint state vector for gradient computation.
     *
     * @param state Final-time state vector 
     * @param statebar Adjoint state vector to update
     * @param J_re_bar Adjoint of real part of objective
     * @param J_im_bar Adjoint of imaginary part of objective
     */
    void evalJ_diff(const Vec state, Vec statebar, const double J_re_bar, const double J_im_bar);

    /**
     * @brief Finalizes the objective function computation.
     * 
     * Compute the infidelity (1-fidelity).
     * 
     * @param obj_cost_re Real part of objective cost
     * @param obj_cost_im Imaginary part of objective cost
     * @return double Final objective function value
     */
    double finalizeJ(const double obj_cost_re, const double obj_cost_im); 

    /**
     * @brief Computes the Riemannian distance between target and current unitary: J(U) = 1/2 || log(U^\dagger V) ||^2_F
     */
    double RiemannianDistance(const Mat U_final_re, const Mat U_final_im);

    /**
     * @brief Derivative of objective function finalization.
     *
     * @param[in] obj_cost_re Real part of objective cost
     * @param[in] obj_cost_im Imaginary part of objective cost
     * @param[out] obj_cost_re_bar Pointer to store adjoint of real part
     * @param[out] obj_cost_im_bar Pointer to store adjoint of imaginary part
     */
    void finalizeJ_diff(const double obj_cost_re, const double obj_cost_im, double* obj_cost_re_bar, double* obj_cost_im_bar); 

    /**
     * @brief Computes Frobenius distance between target and current state.
     *
     * Calculates \f$F = 1/2 || \rho_{target} - \rho||^2_F\f$
     *
     * @param state Current state vector
     * @return double Frobenius distance
     */
    double FrobeniusDistance(const Vec state);

    /**
     * @brief Derivative of Frobenius distance computation.
     *
     * @param state Current state vector
     * @param statebar Adjoint state vector to update
     * @param Jbar Adjoint seed 
     */
    void FrobeniusDistance_diff(const Vec state, Vec statebar, const double Jbar);

    /**
     * @brief Computes Hilbert-Schmidt overlap between state and target.
     *
     * Calculates \f$ Tr(\rho^\dagger \rho_{target})\f$, optionally scaled by the purity of the target state.
     *
     * @param state Current state vector
     * @param scalebypurity Flag to scale by purity of target state
     * @param HS_re_ptr Pointer to store real part of overlap
     * @param Hs_im_ptr Pointer to store imaginary part of overlap
     */
    void HilbertSchmidtOverlap(const Vec state, const bool scalebypurity, double* HS_re_ptr, double* Hs_im_ptr );

    /**
     * @brief Derivative of Hilbert-Schmidt overlap computation.
     *
     * @param statebar Adjoint state vector to update
     * @param scalebypurity Flag to scale by purity of target state
     * @param HS_re_bar Adjoint of real part of overlap
     * @param HS_im_bar Adjoint of imaginary part of overlap
     */
    void HilbertSchmidtOverlap_diff(Vec statebar, bool scalebypurity, const double HS_re_bar, const double HS_im_bar);
};

