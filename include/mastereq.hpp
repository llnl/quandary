#include "defs.hpp"
#include "oscillator.hpp"
#include "util.hpp"
#include <petscts.h>
#include <vector>
#include <assert.h>
#include <iostream> 
#include "gate.hpp"
#include "hamiltonianfilereader.hpp"

#pragma once


/**
 * @brief Matrix shell context containing data needed for applying the right-hand-side (RHS) system matrix to a vector.
 *
 * This structure holds all the necessary data for applying the real-valued 
 * and vectorized system matrix to a state vector.
 */
typedef struct {
  PetscInt dim; ///< Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
  std::vector<size_t> nlevels; ///< Number of levels per oscillator
  IS *isu, *isv; ///< Vector strides for accessing real and imaginary parts
  Oscillator** oscil_vec; ///< Array of pointers to the oscillators
  std::vector<double> crosskerr; ///< Cross-Kerr coupling coefficients
  std::vector<double> Jkl; ///< Dipole-dipole coupling strength
  std::vector<double> eta; ///< Frequency differences of the rotating frames
  LindbladType lindbladtype; ///< Type of Lindblad operators to include
  bool addT1, addT2; ///< Flags for T1 decay and T2 dephasing
  std::vector<double> control_Re;  ///< Real parts of control pulse \f$p(t)\f$
  std::vector<double> control_Im;  ///< Imaginary parts of control pulse \f$q(t)\f$
  std::vector<Mat> Ac_vec; ///< Vector of real parts of control matrices per oscillator
  std::vector<Mat> Bc_vec; ///< Vector of imaginary parts of control matrices per oscillator
  Mat *Ad; ///< Real parts of time-independent system matrix 
  Mat *Bd; ///< Imaginary parts of time-independent system matrix 
  std::vector<Mat> Ad_vec; ///< Vector of real parts of dipole-dipole coupling system matrices
  std::vector<Mat> Bd_vec; ///< Vector of imaginary parts of dipole-dipole coupling system matrices
  std::vector<double> Bd_coeffs;  //< Time-dependent coefficients for dipole-dipole coupling matrices: cos(eta_k*t)
  std::vector<double> Ad_coeffs;  //< Time-dependent coefficients for dipole-dipole coupling matrices: sin(eta_k*t)
  Vec *aux; ///< Auxiliary vector for computations
  double time; ///< Current time
} MatShellCtx;


/**
 * @brief Matrix-vector products to apply the RHS system matrix to a state vector.
 * 
 * Each function returns y = RHS*x. Matrix-free versions as well as sparse-matrix versions are implemented. Those 
 * functions will be passed to PETSc's MatMult operations to realize a Matrix-Vector Multiplication of the RHS with 
 * a state vector, such that one can call PETSc's MatMult(RHS, x, y) for the RHS. 
 * Note: The RHS matrix must be assembled before usage, see @ref MasterEq::assemble_RHS.
 */
int applyRHS_matfree_1Osc(Mat RHS, Vec x, Vec y); ///< Matrix-free MatMult for 1 oscillator
int applyRHS_matfree_transpose_1Osc(Mat RHS, Vec x, Vec y); ///< Transpose matrix-free MatMult for 1 oscillator
int applyRHS_matfree_2Osc(Mat RHS, Vec x, Vec y); ///< Matrix-free MatMult for 2 oscillators
int applyRHS_matfree_transpose_2Osc(Mat RHS, Vec x, Vec y); ///< Transpose matrix-free MatMult for 2 oscillators
int applyRHS_matfree_3Osc(Mat RHS, Vec x, Vec y); ///< Matrix-free MatMult for 3 oscillators
int applyRHS_matfree_transpose_3Osc(Mat RHS, Vec x, Vec y); ///< Transpose matrix-free MatMult for 3 oscillators
int applyRHS_matfree_4Osc(Mat RHS, Vec x, Vec y); ///< Matrix-free MatMult for 4 oscillators
int applyRHS_matfree_transpose_4Osc(Mat RHS, Vec x, Vec y); ///< Transpose matrix-free MatMult for 4 oscillators
int applyRHS_matfree_5Osc(Mat RHS, Vec x, Vec y); ///< Matrix-free MatMult for 5 oscillators
int applyRHS_matfree_transpose_5Osc(Mat RHS, Vec x, Vec y); ///< Transpose matrix-free MatMult for 5 oscillators
int applyRHS_sparsemat(Mat RHS, Vec x, Vec y); ///< Sparse matrix MatMult
int applyRHS_sparsemat_transpose(Mat RHS, Vec x, Vec y); ///< Transpose sparse matrix MatMult


/**
 * @brief Implementation of the real-valued right-hand-side (RHS) system matrix of the quantum dynamical equations.
 *
 * This class provides functionality for evaluating and applying the real-valued right-hand-side (RHS) 
 * system matrix of the vectorized Lindblad master equation (open quantum systems) or Schroedinger equation 
 * (closed systems). It supports a matrix-free and a sparse-matrix version for applying the RHS to a state vector. 
 * The system matrix (RHS) is stored as a Matrix Shell, that applies the real and imaginary parts of 
 * the system matrix to the corresponding parts of the state vector separately. The RHS needs to be assembled 
 * (evaluated at each time step t) before usage.  
 * 
 * Main functionality:
 *    - @ref assemble_RHS  for preparing the RHS system matrix shell with current control pulse values at given 
 *      time t. 
 *    - @ref compute_dRHS_dParams for updating the gradient with += x^T*RHS(t)^T*x_bar, where x & x_bar are the primal and 
 *      adjoint state at time t
 *    - @ref setControlAmplitudes for passing the global optimization vector to each oscillator
 *    - @ref expectedEnergy and @ref population for evaluating expected energy level and level occupations of the 
 *      full composite system.
 *    - definition of MatMult interface functions to multiply the RHS matrix shell to a state vector (sparse matrix
 *      multiplication or matrix-free multiplication)
 * 
 * This class contains references to:
 *    - Array of @ref Oscillator for defining the Hamiltonian matrices for each subsystem (oscillator) and for evaluating their time-dependent control pulses at each time step
 */
class MasterEq{

  protected:
    PetscInt dim; ///< Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
    PetscInt dim_rho; ///< Dimension of Hilbert space = N
    PetscInt dim_ess; ///< Dimension of essential levels = N_e
    int noscillators; ///< Number of oscillators in the system
    Oscillator** oscil_vec; ///< Array of pointers to oscillator objects

    Mat RHS; ///< MatShell for real-valued, vectorized system matrix (size 2N^2 x 2N^2)
    MatShellCtx RHSctx; ///< Context containing data for applying the RHS to a state

    std::vector<Mat> Ac_vec;  // Vector of constant mats for time-varying control term (real). One for each oscillators. 
    std::vector<Mat> Bc_vec;  // Vector of constant mats for time-varying control term (imag). One for each oscillators. 
    Mat  Ad, Bd;  // Real and imaginary part of constant system matrix
    std::vector<Mat> Ad_vec;  // Vector of constant mats for Dipole-Dipole coupling term in drift Hamiltonian (real)
    std::vector<Mat> Bd_vec;  // Vector of constant mats for Dipole-Dipole coupling term in drift Hamiltonian (imag)

    std::vector<double> crosskerr; ///< Cross-Kerr coefficients (rad/time) \f$\xi_{kl}\f$ for ZZ-coupling \f$a_k^\dagger a_k a_l^\dagger a_l\f$
    std::vector<double> Jkl; ///< Dipole-dipole coupling coefficients (rad/time), multiplies \f$a_k^\dagger a_l + a_k a_l^\dagger\f$
    std::vector<double> eta; ///< Frequency differences in rotating frame (rad/time) for dipole-dipole coupling
    bool addT1, addT2; ///< Flags for including T1 decay and T2 dephasing Lindblad operators

    int mpirank_world; ///< Rank of global MPI communicator
    int mpirank_petsc; ///< Rank of PETSc's communicator
    int mpisize_petsc; ///< Size of PETSc's communicator
    PetscInt localsize_u; ///< Size of local sub vector u or v in state x=[u,v]
    PetscInt ilow; ///< First index of the local sub vector u,v
    PetscInt iupp; ///< Last index (+1) of the local sub vector u,v

    IS isu, isv; ///< Vector strides for accessing real and imaginary parts u=Re(x), v=Im(x)
    Vec aux; ///< Auxiliary vector for computations
    bool quietmode; ///< Flag for quiet mode operation
    std::string hamiltonian_file_Hsys; ///< Filename if a custom system Hamiltonian is read from file ('none' if standard Hamiltonian is used)
    std::string hamiltonian_file_Hc; ///< Filename if a custom control Hamiltonians are read from file ('none' if standard Hamiltonian is used)

    bool transmon_resonator_system; ///< Flag for using transmon-resonator Hamiltonian
    double charge_offset; ///< Charge offset for transmon qubits in charge basis
    double Ec;  ///< Charging energy for transmon qubits
    double Ej;  ///< Josephson energy for transmon qubits
    double resonator_freq;  ///< Resonator frequency
    double resonator_rot_freq;  ///< Resonator rotating frame frequency


  public:
    std::vector<size_t> nlevels; ///< Number of levels per oscillator
    std::vector<size_t> nessential; ///< Number of essential levels per oscillator
    bool usematfree; ///< Flag for using matrix-free solver
    LindbladType lindbladtype; ///< Type of Lindblad operators to include (NONE means Schroedinger equation)

  public:
    MasterEq();

    /**
     * @brief Constructor with full system specification.
     *
     * @param config Configuration parameters
     * @param nlevels Number of levels per oscillator
     * @param nessential Number of essential levels per oscillator
     * @param oscil_vec_ Array of pointers to oscillator objects
     * @param crosskerr_ Cross-Kerr coupling coefficients
     * @param Jkl_ Dipole-dipole coupling coefficients
     * @param eta_ Frequency differences for rotating frame
     * @param lindbladtype_ Type of Lindblad operators to include
     * @param usematfree_ Flag to use matrix-free solver
     * @param hamiltonian_file_Hsys Filename for system Hamiltonian data
     * @param hamiltonian_file_Hc Filename for control Hamiltonian data
     * @param quietmode Flag for quiet operation (default: false)
     */
    MasterEq(Config config, const std::vector<size_t>& nlevels, const std::vector<size_t>& nessential, Oscillator** oscil_vec_, const std::vector<double>& crosskerr_, const std::vector<double>& Jkl_, const std::vector<double>& eta_, LindbladType lindbladtype_, bool usematfree_, const std::string& hamiltonian_file_Hsys, const std::string& hamiltonian_file_Hc, bool quietmode=false);

    ~MasterEq();

    /** 
     * @brief Pass MatMult operations for applying the RHS to PETSc.
     * 
     * The RHS is stored as a Matrix Shell. This routine passes the (transpose)
     * Matrix-Multiplication operation for applying the RHS to a state vector to PETSc. 
     */
    void set_RHS_MatMult_operation();

    /**
     * @brief Allocates and initializes matrices needed for the sparse matrix solver.
     */
    void initSparseMatSolver();

    void initStdHamiltonianMats();
    void initStdControlMats(int osc_id, int nlevels_osc, int nlevels_postosc);
    void initStdDecoherenceMats();


    /*
    * @brief Initializes sparse matrices for transmon-resonator system. 
    */
    void initTransmonResonatorSparseMats();

    /**
     * @brief Retrieves the i-th oscillator.
     *
     * @param i Index of the oscillator
     * @return Oscillator* Pointer to the oscillator object
     */
    Oscillator* getOscillator(const size_t i) { return oscil_vec[i]; }

    /**
     * @brief Retrieves the number of oscillators in the system.
     *
     * @return size_t Number of oscillators
     */
    size_t getNOscillators() { return noscillators; }

    /**
     * @brief Retrieves the dimension of the vectorized system.
     *
     * @return PetscInt \f$N^2\f$ for Lindblad solver, \f$N\f$ for Schroedinger solver
     */
    PetscInt getDim(){ return dim; }

    /**
     * @brief Retrieves the dimension of the essential level system.
     *
     * @return PetscInt Dimension N_e of essential levels
     */
    PetscInt getDimEss(){ return dim_ess; }

    /**
     * @brief Retrieves the dimension of the density matrix.
     *
     * @return PetscInt Dimension N of the Hilbert space
     */
    PetscInt getDimRho(){ return dim_rho; }

    /**
     * @brief Assembles the real-valued system matrix (RHS) at time t.
     *
     * Updates the time-dependent parameters in the RHS MatShell context. This must be 
     * called before applying the RHS to a state vector.
     *
     * @param t Current time
     * @return int Error code
     */
    int assemble_RHS(const double t);

    /**
     * @brief Retrieves the right-hand-side system matrix.
     *
     * @return Mat PETSc matrix shell object representing the RHS.
     */
    Mat getRHS();

    /**
     * @brief Computes gradient of RHS with respect to control parameters.
     *
     * Updates grad += alpha * x^T * (d RHS / d params)^T * x_bar. Supports both 
     * the sparse-matrix and the matrix-free version of the RHS. 
     *
     * @param t Current time
     * @param x State vector
     * @param x_bar Adjoint state vector
     * @param alpha Scaling factor
     * @param grad Gradient vector to update
     */
    void compute_dRHS_dParams(const double t,const Vec x,const Vec x_bar, const double alpha, Vec grad);

    /**
     * @brief Pass control parameters from global design vector to each oscillator.
     *
     * @param x Global design vector containing control parameters
     */
    void setControlAmplitudes(const Vec x);

    /**
     * @brief Computes expected energy of the full composite system.
     *
     * @param x State vector
     * @return double Expected energy value
     */
    double expectedEnergy(const Vec x);

    /**
     * @brief Computes population of the full composite system.
     *
     * @param x State vector
     * @param population_com Reference to vector to store population values
     */
    void population(const Vec x, std::vector<double> &population_com);

    // /* Compute reduced density operator for a sub-system defined by IDs in the oscilIDs vector */
    // void createReducedDensity(const Vec rho, Vec *reduced, const std::vector<int>& oscilIDs);
    // /* Derivative of reduced density computation */
    // void createReducedDensity_diff(Vec rhobar, const Vec reducedbar, const std::vector<int>& oscilIDs);
};


/**
 * @brief: Sparse-matrix version to compute gradient of RHS with respect to parameters 
 *
 * Updates grad += alpha * x^T * (d RHS / d params)^T * x_bar with sparse-matrix
 * version of RHS. Compare @ref compute_dRHS_dParams_matfree for the matrix-free version of this routine.
 *
 * @param[in] t Current time
 * @param[in] x State vector
 * @param[in] x_bar Adjoint state vector
 * @param[in] alpha Scaling factor
 * @param[out] grad Gradient vector to update
 * @param[in] nlevels Number of energy levels per subsystem 
 * @param[in] isu Index stride to access real parts of a state vector
 * @param[in] isv Index stride to access imaginar parts of a state vector
 * @param[in] Ac_vec Vector of real parts of control matrices per oscillator
 * @param[in] Bc_vec Vector of imaginary parts of control matrices per oscillator
 * @param[in] aux Auxiliary vector for computations 
 * @param[in] oscil_vec Vector of quantum oscilators
 */
void compute_dRHS_dParams_sparsemat(const double t,const Vec x,const Vec x_bar, const double alpha, Vec grad, std::vector<size_t>& nlevels, IS isu, IS isv, std::vector<Mat>& Ac_vec, std::vector<Mat>& Bc_vec, Vec aux, Oscillator** oscil_vec);

/**
 * @brief: Matrix free version to compute gradient of RHS with respect to parameters 
 *
 * Updates grad += alpha * x^T * (d RHS / d params)^T * x_bar in a matrix-free 
 * manner. See @ref compute_dRHS_dParams_sparsemat for the sparse-matrix version of this routine.
 *
 * @param[in] dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param[in] t Current time
 * @param[in] x State vector
 * @param[in] x_bar Adjoint state vector
 * @param[in] alpha Scaling factor
 * @param[out] grad Gradient vector to update
 * @param[in] nlevels Number of energy levels per subsystem
 * @param[in] lindbladtype Type of Lindblad decoherence operators, or NONE
 * @param[in] oscil_vec Vector of quantum oscillators 
 */
void compute_dRHS_dParams_matfree(const PetscInt dim, const double t,const Vec x,const Vec x_bar, const double alpha, Vec grad, std::vector<size_t>& nlevels, LindbladType lindbladtype, Oscillator** oscil_vec);



// Inline functions for the Matrix-free RHS application

/**
 * @brief Inline for Matrix-free RHS to Compute detuning for 1 oscillator.
 *
 * @param detuning0 Detuning frequency 
 * @param a Occupation number
 * @return double Detuning energy contribution
 */
inline double H_detune(const double detuning0, const int a) {
  return detuning0*a;
};

/**
 * @brief Inline for Matrix-free RHS to Compute self-Kerr coefficient for 1 oscillator.
 *
 * @param xi0 Self-Kerr coefficient for oscillator 0
 * @param a Occupation number
 * @return double Self-Kerr energy contribution
 */
inline double H_selfkerr(const double xi0, const int a) {
  return - xi0 / 2.0 * a * (a-1);
};

/**
 * @brief Inline for Matrix-free RHS to compute L2 dephasing coefficient for 1 oscillator.
 *
 * @param dephase0 Dephasing rate for oscillator 0
 * @param i0 Occupation number (bra)
 * @param i0p Occupation number (ket)
 * @return double L2 operator contribution
 */
inline double L2(const double dephase0, const int i0, const int i0p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) );
};

/**
 * @brief Inline for Matrix-free RHS to compute L1 decay coefficient for 1 oscillator.
 *
 * @param decay0 Decay rate for oscillator 0
 * @param i0 Occupation number (bra)
 * @param i0p Occupation number (ket)
 * @return double L1 diagonal contribution
 */
inline double L1diag(const double decay0, const int i0, const int i0p){
  return - decay0 / 2.0 * ( i0 + i0p ) ;
};

/**
 * @brief Inline for Matrix-free RHS to compute tensor product index for 1 oscillator system.
 *
 * @param nlevels0 Number of levels for oscillator 0
 * @param i0 Occupation number (bra)
 * @param i0p Occupation number (ket)
 * @return int Linear index in tensor product space
 */
inline int TensorGetIndex(const int nlevels0, const  int i0, const int i0p){
  return i0 + nlevels0 * i0p;
};

/**
 * @brief Inline for Matrix-free RHS, Computes detuning Hamiltonian term for 2 oscillators.
 *
 * @param detuning0 Detuning frequency for oscillator 0
 * @param detuning1 Detuning frequency for oscillator 1
 * @param a Occupation number for oscillator 0
 * @param b Occupation number for oscillator 1
 * @return double Total detuning energy contribution
 */
inline double H_detune(const double detuning0, const double detuning1, const int a, const int b) {
  return detuning0*a + detuning1*b;
};

/**
 * @brief Inline for Matrix-free RHS, Computes self-Kerr nonlinearity terms for 2 oscillators.
 *
 * @param xi0 Self-Kerr coefficient for oscillator 0
 * @param xi1 Self-Kerr coefficient for oscillator 1
 * @param a Occupation number for oscillator 0
 * @param b Occupation number for oscillator 1
 * @return double Total self-Kerr energy contribution
 */
inline double H_selfkerr(const double xi0, const double xi1, const int a, const int b) {
  return - xi0 / 2.0 * a * (a-1) - xi1 / 2.0 * b * (b-1);
};

/**
 * @brief Inline for Matrix-free RHS, Computes cross-Kerr coupling between 2 oscillators.
 *
 * @param xi01 Cross-Kerr coefficient between oscillators 0 and 1
 * @param a Occupation number for oscillator 0
 * @param b Occupation number for oscillator 1
 * @return double Cross-Kerr energy contribution
 */
inline double H_crosskerr(const double xi01, const int a, const int b) {
  return - xi01 * a * b;
};

/**
 * @brief Inline for Matrix-free RHS, Computes L2 dephasing Lindblad operator for 2 oscillators.
 *
 * @param dephase0 Dephasing rate for oscillator 0
 * @param dephase1 Dephasing rate for oscillator 1
 * @param i0 Occupation number for oscillator 0 (bra)
 * @param i1 Occupation number for oscillator 1 (bra)
 * @param i0p Occupation number for oscillator 0 (ket)
 * @param i1p Occupation number for oscillator 1 (ket)
 * @return double Total L2 operator contribution
 */
inline double L2(const double dephase0, const double dephase1, const int i0, const int i1, const int i0p, const int i1p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) ) + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) );
};

/**
 * @brief Inline for Matrix-free RHS, Computes L1 decay Lindblad operator diagonal term for 2 oscillators.
 *
 * @param decay0 Decay rate for oscillator 0
 * @param decay1 Decay rate for oscillator 1
 * @param i0 Occupation number for oscillator 0 (bra)
 * @param i1 Occupation number for oscillator 1 (bra)
 * @param i0p Occupation number for oscillator 0 (ket)
 * @param i1p Occupation number for oscillator 1 (ket)
 * @return double Total L1 diagonal contribution
 */
inline double L1diag(const double decay0, const double decay1, const int i0, const int i1, const int i0p, const int i1p){
  return - decay0 / 2.0 * ( i0 + i0p ) - decay1 / 2.0 * ( i1 + i1p );
};

/**
 * @brief Inline for Matrix-free RHS, Computes tensor product index for 2 oscillator system.
 *
 * @param nlevels0 Number of levels for oscillator 0
 * @param nlevels1 Number of levels for oscillator 1
 * @param i0 Occupation number for oscillator 0 (bra)
 * @param i1 Occupation number for oscillator 1 (bra)
 * @param i0p Occupation number for oscillator 0 (ket)
 * @param i1p Occupation number for oscillator 1 (ket)
 * @return int Linear index in tensor product space
 */
inline int TensorGetIndex(const int nlevels0, const int nlevels1,const  int i0, const int i1, const int i0p, const int i1p){
  return i0*nlevels1 + i1 + (nlevels0 * nlevels1) * ( i0p * nlevels1 + i1p);
};


// Inlines for Matrix-free RHS for 3 oscillators
// See documentation for 1-2 oscillator functions above for parameter descriptions
inline double H_detune(const double detuning0, const double detuning1, const double detuning2, const int i0, const int i1, const int i2) {
  return detuning0*i0 + detuning1*i1 + detuning2*i2;
};
inline double H_selfkerr(const double xi0, const double xi1, const double xi2, const int i0, const int i1, const int i2) {
  return - xi0 / 2.0 * i0 * (i0-1) - xi1 / 2.0 * i1 * (i1-1) - xi2 / 2.0 * i2 * (i2-1);
};
inline double H_crosskerr(const double xi01, const double xi02, const double xi12, const int i0, const int i1, const int i2) {
  return - xi01 * i0 * i1 - xi02 * i0 * i2 - xi12 * i1 * i2;
};
inline double L2(const double dephase0, const double dephase1, const double dephase2, const int i0, const int i1, const int i2, const int i0p, const int i1p, const int i2p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) ) 
       + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) )
       + dephase2 * ( i2*i2p - 1./2. * (i2*i2 + i2p*i2p) );
};
inline double L1diag(const double decay0, const double decay1, const double decay2, const int i0, const int i1, const int i2, const int i0p, const int i1p, const int i2p){
  return - decay0 / 2.0 * ( i0 + i0p ) 
         - decay1 / 2.0 * ( i1 + i1p )
         - decay2 / 2.0 * ( i2 + i2p );
};
inline int TensorGetIndex(const int nlevels0, const int nlevels1, const int nlevels2, const  int i0, const int i1, const int i2, const int i0p, const int i1p, const int i2p){
  return i0*nlevels1*nlevels2 + i1*nlevels2 + i2 + (nlevels0 * nlevels1 * nlevels2) * ( i0p * nlevels1*nlevels2 + i1p*nlevels2 + i2p);
};


// Matrix-free solver inlines for 4 oscillators
// See documentation for 1-2 oscillator functions above for parameter descriptions
inline double H_detune(const double detuning0, const double detuning1, const double detuning2, const double detuning3, const int i0, const int i1, const int i2, const int i3) {
  return detuning0*i0 + detuning1*i1 + detuning2*i2 + detuning3*i3;
};
inline double H_selfkerr(const double xi0, const double xi1, const double xi2, const double xi3, const int i0, const int i1, const int i2, const int i3) {
  return - xi0 / 2.0 * i0 * (i0-1) - xi1 / 2.0 * i1 * (i1-1) - xi2 / 2.0 * i2 * (i2-1) - xi3/2.0 * i3 * (i3-1);
};
inline double H_crosskerr(const double xi01, const double xi02, const double xi03, const double xi12, const double xi13, const double xi23, const int i0, const int i1, const int i2, const int i3) {
  return - xi01 * i0 * i1 - xi02 * i0 * i2  - xi03*i0*i3 - xi12 * i1 * i2 - xi13*i1*i3 - xi23*i2*i3;
};
inline double L2(const double dephase0, const double dephase1, const double dephase2, const double dephase3, const int i0, const int i1, const int i2, const int i3, const int i0p, const int i1p, const int i2p, const int i3p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) ) 
       + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) )
       + dephase2 * ( i2*i2p - 1./2. * (i2*i2 + i2p*i2p) )
       + dephase3 * ( i3*i3p - 1./2. * (i3*i3 + i3p*i3p) );
};
inline double L1diag(const double decay0, const double decay1, const double decay2, const double decay3, const int i0, const int i1, const int i2, const int i3, const int i0p, const int i1p, const int i2p, const int i3p){
  return - decay0 / 2.0 * ( i0 + i0p ) 
         - decay1 / 2.0 * ( i1 + i1p )
         - decay2 / 2.0 * ( i2 + i2p )
         - decay3 / 2.0 * ( i3 + i3p );
};
inline int TensorGetIndex(const int nlevels0, const int nlevels1, const int nlevels2, const int nlevels3, const  int i0, const int i1, const int i2, const int i3, const int i0p, const int i1p, const int i2p, const int i3p){
  return i0*nlevels1*nlevels2*nlevels3 + i1*nlevels2*nlevels3 + i2*nlevels3 + i3 + (nlevels0 * nlevels1 * nlevels2 * nlevels3) * ( i0p * nlevels1*nlevels2*nlevels3 + i1p*nlevels2*nlevels3 + i2p*nlevels3 + i3p);
}

// Matrix-free solver inlines for 5 oscillators
// See documentation for 1-2 oscillator functions above for parameter descriptions
inline double H_detune(const double detuning0, const double detuning1, const double detuning2, const double detuning3, const double detuning4, const int i0, const int i1, const int i2, const int i3, const int i4) {
  return detuning0*i0 + detuning1*i1 + detuning2*i2 + detuning3*i3 + detuning4*i4;
};
inline double H_selfkerr(const double xi0, const double xi1, const double xi2, const double xi3, const double xi4, const int i0, const int i1, const int i2, const int i3, const int i4) {
  return - xi0 / 2.0 * i0 * (i0-1) - xi1 / 2.0 * i1 * (i1-1) - xi2 / 2.0 * i2 * (i2-1) - xi3/2.0 * i3 * (i3-1) - xi4/2.0 * i4 * (i4-1);
};
inline double H_crosskerr(const double xi01, const double xi02, const double xi03, const double xi04, const double xi12, const double xi13, const double xi14, const double xi23, const double xi24, const double xi34, const int i0, const int i1, const int i2, const int i3, const int i4) {
  return - xi01 * i0 * i1 - xi02 * i0 * i2  - xi03*i0*i3 - xi04*i0*i4 - xi12 * i1 * i2 - xi13*i1*i3 - xi14*i1*i4 - xi23*i2*i3 - xi24*i2*i4 - xi34*i3*i4;
};
inline double L2(const double dephase0, const double dephase1, const double dephase2, const double dephase3, const double dephase4, const int i0, const int i1, const int i2, const int i3, const int i4, const int i0p, const int i1p, const int i2p, const int i3p, const int i4p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) ) 
       + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) )
       + dephase2 * ( i2*i2p - 1./2. * (i2*i2 + i2p*i2p) )
       + dephase3 * ( i3*i3p - 1./2. * (i3*i3 + i3p*i3p) )
       + dephase4 * ( i4*i4p - 1./2. * (i4*i4 + i4p*i4p) );
};
inline double L1diag(const double decay0, const double decay1, const double decay2, const double decay3, const double decay4, const int i0, const int i1, const int i2, const int i3, const int i4, const int i0p, const int i1p, const int i2p, const int i3p, const int i4p){
  return - decay0 / 2.0 * ( i0 + i0p ) 
         - decay1 / 2.0 * ( i1 + i1p )
         - decay2 / 2.0 * ( i2 + i2p )
         - decay3 / 2.0 * ( i3 + i3p )
         - decay4 / 2.0 * ( i4 + i4p );
};
inline int TensorGetIndex(const int nlevels0, const int nlevels1, const int nlevels2, const int nlevels3, const int nlevels4, const  int i0, const int i1, const int i2, const int i3, const int i4, const int i0p, const int i1p, const int i2p, const int i3p, const int i4p){
  return i0*nlevels1*nlevels2*nlevels3*nlevels4 + i1*nlevels2*nlevels3*nlevels4 + i2*nlevels3*nlevels4 + i3*nlevels4 + i4 + (nlevels0 * nlevels1 * nlevels2 * nlevels3*nlevels4) * ( i0p * nlevels1*nlevels2*nlevels3*nlevels4 + i1p*nlevels2*nlevels3*nlevels4 + i2p*nlevels3*nlevels4+ i3p*nlevels4 + i4p);
}


/**
 * @brief Inline for Matrix-free RHS for gradient updates.
 *
 * Computes coefficients for gradient computation with respect to control parameters.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param n Number of levels for current oscillator
 * @param np Number of levels for conjugate oscillator
 * @param i Current occupation number
 * @param ip Conjugate occupation number
 * @param stridei Stride for current oscillator
 * @param strideip Stride for conjugate oscillator
 * @param xptr Pointer to state vector data
 * @param res_p_re Pointer to store real part of p result
 * @param res_p_im Pointer to store imaginary part of p result
 * @param res_q_re Pointer to store real part of q result
 * @param res_q_im Pointer to store imaginary part of q result
 */
inline void dRHSdp_getcoeffs(const PetscInt dim, const int it, const int n, const int np, const int i, const int ip, const int stridei, const int strideip, const double* xptr, double* res_p_re, double* res_p_im, double* res_q_re, double* res_q_im) {

  *res_p_re = 0.0;
  *res_p_im = 0.0;
  *res_q_re = 0.0;
  *res_q_im = 0.0;

  /* ik+1..,ik'.. term */
  if (i < n-1) {
    int itx = it + stridei;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(i + 1);
    *res_p_re +=   sq * xim;
    *res_p_im += - sq * xre;
    *res_q_re +=   sq * xre;
    *res_q_im +=   sq * xim;
  }
  /* \rho(ik..,ik'+1..) */
  if (ip < np-1) {
    int itx = it + strideip;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(ip + 1);
    *res_p_re += - sq * xim;
    *res_p_im += + sq * xre;
    *res_q_re +=   sq * xre;
    *res_q_im +=   sq * xim;
  }
  /* \rho(ik-1..,ik'..) */
  if (i > 0) {
    int itx = it - stridei;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(i);
    *res_p_re += + sq * xim;
    *res_p_im += - sq * xre;
    *res_q_re += - sq * xre;
    *res_q_im += - sq * xim;
  }
  /* \rho(ik..,ik'-1..) */
  if (ip > 0) {
    int itx = it - strideip;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(ip);
    *res_p_re += - sq * xim;
    *res_p_im += + sq * xre;
    *res_q_re += - sq * xre;
    *res_q_im += - sq * xim;
  }
}

/**
 * @brief Matrix-free solver inline for dipole-dipole coupling between oscillators.
 *
 * Implements J_kl coupling terms in the rotating frame between oscillator i and j.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param ni Number of levels for oscillator i
 * @param nj Number of levels for oscillator j
 * @param nip Number of levels for oscillator i (conjugate)
 * @param njp Number of levels for oscillator j (conjugate)
 * @param i Occupation number for oscillator i
 * @param ip Occupation number for oscillator i (conjugate)
 * @param j Occupation number for oscillator j
 * @param jp Occupation number for oscillator j (conjugate)
 * @param stridei Stride for oscillator i
 * @param strideip Stride for oscillator i (conjugate)
 * @param stridej Stride for oscillator j
 * @param stridejp Stride for oscillator j (conjugate)
 * @param xptr Pointer to state vector data
 * @param Jij Coupling coefficient
 * @param cosij Cosine of frequency difference
 * @param sinij Sine of frequency difference
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void Jkl_coupling(const PetscInt dim, const int it, const int ni, const int nj, const int nip, const int njp, const int i, const int ip, const int j, const int jp, const int stridei, const int strideip, const int stridej, const int stridejp, const double* xptr, const double Jij, const double cosij, const double sinij, double* yre, double* yim) {
  if (fabs(Jij)>1e-10) {
    //  1) J_kl (-icos + sin) * ρ_{E−k+l i, i′}
    if (i > 0 && j < nj-1) {
      int itx = it - stridei + stridej;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(i * (j + 1));
      // sin u + cos v + i ( -cos u + sin v)
      *yre += Jij * sq * (   cosij * xim + sinij * xre);
      *yim += Jij * sq * ( - cosij * xre + sinij * xim);
    }
    // 2) J_kl (−icos − sin)sqrt(il*(ik +1)) ρ_{E+k−li,i′}
    if (i < ni-1 && j > 0) {
      int itx = it + stridei - stridej;  // E+k-l i, i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(j * (i + 1)); // sqrt( il*(ik+1))
      // -sin u + cos v + i (-cos u - sin v)
      *yre += Jij * sq * (   cosij * xim - sinij * xre);
      *yim += Jij * sq * ( - cosij * xre - sinij * xim);
    }
    // 3) J_kl ( icos + sin)sqrt(ik'*(il' +1)) ρ_{i,E-k+li'}
    if (ip > 0 && jp < njp-1) {
      int itx = it - strideip + stridejp;  // i, E-k+l i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(ip * (jp + 1)); // sqrt( ik'*(il'+1))
      //  sin u - cos v + i ( cos u + sin v)
      *yre += Jij * sq * ( - cosij * xim + sinij * xre);
      *yim += Jij * sq * (   cosij * xre + sinij * xim);
    }
    // 4) J_kl ( icos - sin)sqrt(il'*(ik' +1)) ρ_{i,E+k-li'}
    if (ip < nip-1 && jp > 0) {
      int itx = it + strideip - stridejp;  // i, E+k-l i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(jp * (ip + 1)); // sqrt( il'*(ik'+1))
      // - sin u - cos v + i ( cos u - sin v)
      *yre += Jij * sq * ( - cosij * xim - sinij * xre);
      *yim += Jij * sq * (   cosij * xre - sinij * xim);
    }
  }
}

/**
 * @brief Transpose of dipole-dipole coupling for adjoint computations.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param ni Number of levels for oscillator i
 * @param nj Number of levels for oscillator j
 * @param nip Number of levels for oscillator i (conjugate)
 * @param njp Number of levels for oscillator j (conjugate)
 * @param i Occupation number for oscillator i
 * @param ip Occupation number for oscillator i (conjugate)
 * @param j Occupation number for oscillator j
 * @param jp Occupation number for oscillator j (conjugate)
 * @param stridei Stride for oscillator i
 * @param strideip Stride for oscillator i (conjugate)
 * @param stridej Stride for oscillator j
 * @param stridejp Stride for oscillator j (conjugate)
 * @param xptr Pointer to state vector data
 * @param Jij Coupling coefficient
 * @param cosij Cosine of frequency difference
 * @param sinij Sine of frequency difference
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void Jkl_coupling_T(const PetscInt dim, const int it, const int ni, const int nj, const int nip, const int njp, const int i, const int ip, const int j, const int jp, const int stridei, const int strideip, const int stridej, const int stridejp, const double* xptr, const double Jij, const double cosij, const double sinij, double* yre, double* yim) {
  if (fabs(Jij)>1e-10) {
    //  1) [...] * \bar y_{E+k-l i, i′}
    if (i < ni-1 && j > 0) {
      int itx = it + stridei - stridej;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(j * (i + 1));
      *yre += Jij * sq * ( - cosij * xim + sinij * xre);
      *yim += Jij * sq * ( + cosij * xre + sinij * xim);
    }
    // 2) J_kl (−icos − sin)sqrt(ik*(il +1)) \bar y_{E-k+li,i′}
    if (i > 0 && j < nj-1) {
      int itx = it - stridei + stridej;  // E-k+l i, i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(i * (j + 1)); // sqrt( ik*(il+1))
      *yre += Jij * sq * ( - cosij * xim - sinij * xre);
      *yim += Jij * sq * ( + cosij * xre - sinij * xim);
    }
    // 3) J_kl ( icos + sin)sqrt(il'*(ik' +1)) \bar y_{i,E+k-li'}
    if (ip < nip-1 && jp > 0) {
      int itx = it + strideip - stridejp;  // i, E+k-l i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(jp * (ip + 1)); // sqrt( il'*(ik'+1))
      *yre += Jij * sq * (   cosij * xim + sinij * xre);
      *yim += Jij * sq * ( - cosij * xre + sinij * xim);
    }
    // 4) J_kl ( icos - sin)sqrt(ik'*(il' +1)) \bar y_{i,E-k+li'}
    if (ip > 0 && jp < njp-1) {
      int itx = it - strideip + stridejp;  // i, E-k+l i'
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(ip * (jp + 1)); // sqrt( ik'*(il'+1))
      *yre += Jij * sq * (   cosij * xim - sinij * xre);
      *yim += Jij * sq * ( - cosij * xre - sinij * xim);
    }
  }

}

/**
 * @brief Matrix-free solver inline for off-diagonal L1 decay term.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param n Number of levels
 * @param i Occupation number (bra)
 * @param ip Occupation number (ket)
 * @param stridei Stride for bra index
 * @param strideip Stride for ket index
 * @param xptr Pointer to state vector data
 * @param decayi Decay rate
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void L1decay(const PetscInt dim, const int it, const int n, const int i, const int ip, const int stridei, const int strideip, const double* xptr, const double decayi, double* yre, double* yim){
  if  (fabs(decayi) > 1e-12) {
    if (i < n-1 && ip < n-1) {
      double l1off = decayi * sqrt((i+1)*(ip+1));
      int itx = it + stridei + strideip;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      *yre += l1off * xre;
      *yim += l1off * xim;
    }
  }
}


/**
 * @brief Matrix-free inline Transpose of off-diagonal L1 decay for adjoint computations.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param i Occupation number (bra)
 * @param ip Occupation number (ket)
 * @param stridei Stride for bra index
 * @param strideip Stride for ket index
 * @param xptr Pointer to state vector data
 * @param decayi Decay rate
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void L1decay_T(const PetscInt dim, const int it, const int i, const int ip, const int stridei, const int strideip, const double* xptr, const double decayi, double* yre, double* yim){
  if (fabs(decayi) > 1e-12) {
      if (i > 0 && ip > 0) {
        double l1off = decayi * sqrt(i*ip);
        int itx = it - stridei - strideip;
        double xre = xptr[itx];
        double xim = xptr[itx+dim];
        *yre += l1off * xre;
        *yim += l1off * xim;
      }
    }
}

/**
 * @brief Matrix-free solver inline for control terms.
 *
 * Applies control Hamiltonian terms (ladder operators) to the state.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param n Number of levels
 * @param i Occupation number (bra)
 * @param np Number of levels (conjugate)
 * @param ip Occupation number (ket)
 * @param stridei Stride for bra index
 * @param strideip Stride for ket index
 * @param xptr Pointer to state vector data
 * @param pt Real part of control amplitude
 * @param qt Imaginary part of control amplitude
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void control(const PetscInt dim, const int it, const int n, const int i, const int np, const int ip, const int stridei, const int strideip, const double* xptr, const double pt, const double qt, double* yre, double* yim){
  /* \rho(ik+1..,ik'..) term */
  if (i < n-1) {
      int itx = it + stridei;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(i + 1);
      *yre += sq * (   pt * xim + qt * xre);
      *yim += sq * ( - pt * xre + qt * xim);
    }
    /* \rho(ik..,ik'+1..) */
    if (ip < np-1) {
      int itx = it + strideip;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(ip + 1);
      *yre += sq * ( -pt * xim + qt * xre);
      *yim += sq * (  pt * xre + qt * xim);
    }
    /* \rho(ik-1..,ik'..) */
    if (i > 0) {
      int itx = it - stridei;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(i);
      *yre += sq * (  pt * xim - qt * xre);
      *yim += sq * (- pt * xre - qt * xim);
    }
    /* \rho(ik..,ik'-1..) */
    if (ip > 0) {
      int itx = it - strideip;
      double xre = xptr[itx];
      double xim = xptr[itx+dim];
      double sq = sqrt(ip);
      *yre += sq * (- pt * xim - qt * xre);
      *yim += sq * (  pt * xre - qt * xim);
    }
}


/**
 * @brief Matrix-free Transpose of control terms for adjoint computations.
 *
 * @param dim Dimension of full vectorized system: N^2 if Lindblad, N if Schroedinger
 * @param it Current tensor index
 * @param n Number of levels
 * @param i Occupation number (bra)
 * @param np Number of levels (conjugate)
 * @param ip Occupation number (ket)
 * @param stridei Stride for bra index
 * @param strideip Stride for ket index
 * @param xptr Pointer to state vector data
 * @param pt Real part of control amplitude
 * @param qt Imaginary part of control amplitude
 * @param yre Pointer to store real part of result
 * @param yim Pointer to store imaginary part of result
 */
inline void control_T(const PetscInt dim, const int it, const int n, const int i, const int np, const int ip, const int stridei, const int strideip, const double* xptr, const double pt, const double qt, double* yre, double* yim){
  /* \rho(ik+1..,ik'..) term */
  if (i > 0) {
    int itx = it - stridei;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(i);
    *yre += sq * ( - pt * xim + qt * xre);
    *yim += sq * (   pt * xre + qt * xim);
  }
  /* \rho(ik..,ik'+1..) */
  if (ip > 0) {
    int itx = it - strideip;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(ip);
    *yre += sq * (  pt * xim + qt * xre);
    *yim += sq * ( -pt * xre + qt * xim);
  }
  /* \rho(ik-1..,ik'..) */
  if (i < n-1) {
    int itx = it + stridei;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(i+1);
    *yre += sq * (- pt * xim - qt * xre);
    *yim += sq * (  pt * xre - qt * xim);
  }
  /* \rho(ik..,ik'-1..) */
  if (ip < np-1) {
    int itx = it + strideip;
    double xre = xptr[itx];
    double xim = xptr[itx+dim];
    double sq = sqrt(ip+1);
    *yre += sq * (+ pt * xim - qt * xre);
    *yim += sq * (- pt * xre - qt * xim);
  }
}