#include <stdio.h>
#include <iostream> 
#include <math.h>
#include <assert.h>
#include <petscmat.h>
#include <vector>
#include "defs.hpp"
#include "util.hpp"
#pragma once

/**
 * @brief Base class for quantum gate operations.
 *
 * The Gate class provides the fundamental interface for quantum gate operations
 * in quantum optimal control. It handles the transformation \f$V \rho V^\dagger\f$ where \f$V\f$ is
 * the unitary gate matrix (target unitary transformation).
 */
class Gate {
  protected:
    Mat V_re, V_im; ///< Real and imaginary parts of V_target, non-vectorized, essential levels only.

    std::vector<size_t> nessential; ///< Number of essential levels per oscillator.
    std::vector<size_t> nlevels; ///< Total number of levels per oscillator.
    int mpirank_petsc; ///< MPI rank in PETSc communicator.
    int mpisize_petsc; ///< MPI size in PETSc communicator.
    int mpirank_world; ///< MPI rank in world communicator.

    bool quietmode; ///< Flag to suppress output messages.

    PetscInt dim_ess; ///< Dimension of target gate matrix (non-vectorized), essential levels only.
    PetscInt dim_rho; ///< Dimension of system matrix rho (non-vectorized), all levels, N.

    double final_time; ///< Final time T. Time of gate rotation.
    std::vector<double> gate_rot_freq; ///< Frequencies of gate rotation (rad/time). Default: 0.0.

    DecoherenceType decoherence_type; ///< Type of Lindblad operators for open system dynamics.


  private:
    Mat VxV_re, VxV_im;     ///< Real and imaginary parts of vectorized gate \f$\bar V \kron V\f$.
    Vec x;                  ///< Auxiliary state vector for computations.
    IS isu, isv;            ///< Vector strides for accessing real and imaginary parts of the state.

  public:
    Gate();

    /**
     * @brief Constructor for quantum gates.
     *
     * @param nlevels_ Number of levels per oscillator.
     * @param nessential_ Number of essential levels per oscillator.
     * @param time_ Final time for gate rotation.
     * @param gate_rot_freq Frequencies of gate rotation (rad/time).
     * @param decoherence_type_ Type of Lindblad operators for open system dynamics.
     * @param quietmode Flag to suppress output messages (default: false).
     */
    Gate(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time_, const std::vector<double>& gate_rot_freq, DecoherenceType decoherence_type_, bool quietmode=false);

    virtual ~Gate();

    /**
     * @brief Retrieves the dimension of the density matrix.
     *
     * @return PetscInt Dimension of the system density matrix (all levels).
     */
    PetscInt getDimRho() { return dim_rho; };

    /**
     * @brief Assembles the vectorized gate matrix.
     *
     * Computes VxV_re = \f$\text{Re}(\bar{V} \otimes V)\f$ and VxV_im = \f$\text{Im}(\bar{V} \otimes V)\f$
     * where V is the gate matrix.
     */
    void assembleGate();

    /**
     * @brief Applies the gate transformation to a quantum state.
     *
     * Computes \f$V\rho V^\dagger\f$ where \f$V\f$ is the gate matrix and \f$\rho\f$ is the input state.
     * The output vector VrhoV must be pre-allocated.
     *
     * @param state Input quantum state vector.
     * @param VrhoV Output vector storing the transformed state.
     */
    void applyGate(const Vec state, Vec VrhoV);
};

/**
 * @brief Pauli-X gate implementation.
 *
 * Implements the Pauli-X (NOT) gate for a single qubit:
 * @code
 * V = | 0  1 |
 *     | 1  0 |
 * @endcode
 */
class XGate : public Gate {
  public:
    XGate(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~XGate();
};

/**
 * @brief Pauli-Y gate implementation.
 *
 * Implements the Pauli-Y gate for a single qubit:
 * @code
 * V = | 0 -i |
 *     | i  0 |
 * @endcode
 */
class YGate : public Gate {
  public:
    YGate(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~YGate();
};

/**
 * @brief Pauli-Z gate implementation.
 *
 * Implements the Pauli-Z gate for a single qubit:
 * @code
 * V = | 1   0 |
 *     | 0  -1 |
 * @endcode
 */
class ZGate : public Gate {
  public:
    ZGate(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~ZGate();
};

/**
 * @brief Hadamard gate implementation.
 *
 * Implements the Hadamard gate for a single qubit:
 * @code
 * V = 1/sqrt(2) * | 1   1 |
 *                 | 1  -1 |
 * @endcode
 */
class HadamardGate : public Gate {
  public:
    HadamardGate(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~HadamardGate();
};

/**
 * @brief Controlled-NOT (CNOT) gate implementation.
 *
 * Implements the CNOT gate for two qubits:
 * @code
 * V = | 1  0  0  0 |
 *     | 0  1  0  0 |
 *     | 0  0  0  1 |
 *     | 0  0  1  0 |
 * @endcode
 */
class CNOT : public Gate {
    public:
    CNOT(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~CNOT();
};

/**
 * @brief SWAP gate implementation.
 *
 * Implements the SWAP gate for two qubits:
 * @code
 * V = | 1  0  0  0 |
 *     | 0  0  1  0 |
 *     | 0  1  0  0 |
 *     | 0  0  0  1 |
 * @endcode
 */
class SWAP: public Gate {
    public:
    SWAP(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~SWAP();
};


/**
 * @brief Multi-qubit SWAP gate implementation.
 *
 * Implements a SWAP gate for Q qubits, swapping qubit 0 with qubit Q-1
 * while leaving all other qubits unchanged.
 */
class SWAP_0Q: public Gate {
    public:
    SWAP_0Q(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~SWAP_0Q();
};

/**
 * @brief Multi-qubit CQNOT gate implementation.
 *
 * Implements a Q-qubit controlled-NOT gate where the NOT operation is applied
 * to qubit Q-1 and is controlled by all other qubits (qubits 0 through Q-2).
 */
class CQNOT: public Gate {
    public:
    CQNOT(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~CQNOT();
};

/**
 * @brief Quantum Fourier Transform gate implementation.
 *
 * Implements the discrete Quantum Fourier Transform gate for multiple qubits.
 */
class QFT: public Gate {
  public:
    QFT(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, bool quietmode=false);
    ~QFT();
};

/**
 * @brief Gate loaded from external file.
 *
 * Implements a custom gate by loading the gate matrix from an external file.
 */
class FromFile: public Gate {
  public:
    /**
     * @brief Constructor for gate loaded from file.
     *
     * @param nlevels_ Number of levels per oscillator.
     * @param nessential_ Number of essential levels per oscillator.
     * @param time Final time for gate rotation.
     * @param rotation_frequencies_ Frequencies of gate rotation.
     * @param decoherence_type_ Type of Lindblad operators.
     * @param filename Path to file containing gate matrix data.
     * @param quietmode Flag to suppress output messages (default: false).
     */
    FromFile(const std::vector<size_t>& nlevels_, const std::vector<size_t>& nessential_, double time, const std::vector<double>& rotation_frequencies_, DecoherenceType decoherence_type_, const std::string& filename, bool quietmode=false);
    ~FromFile();
};


/**
 * @brief Factory function to initialize a gate from configuration string.
 *
 * Creates and returns a pointer to the appropriate Gate subclass based on
 * the target specification string.
 *
 * @param target_gate target_gate Type of gate to initialize.
 * @param file Path to the file containing gate matrix data (if applicable).
 * @param nlevels Number of levels per oscillator.
 * @param nessential Number of essential levels per oscillator.
 * @param total_time Total time for gate operation.
 * @param decoherence_type Type of Lindblad operators for open system dynamics.
 * @param gate_rot_freq Frequencies of gate rotation.
 * @param quietmode Flag to suppress output messages.
 * @return Gate* Pointer to the initialized gate object.
 */
Gate* initTargetGate(GateType target_gate, const std::string& file, const std::vector<size_t>& nlevels, const std::vector<size_t>& nessential, double total_time, DecoherenceType decoherence_type, const std::vector<double>& gate_rot_freq, bool quietmode);
