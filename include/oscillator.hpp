#include <stdio.h>
#include "defs.hpp"
#include "controlbasis.hpp"
#include <fstream>
#include <iostream> 
#include <iomanip>
#include <petscmat.h>
#include <vector>
#include <assert.h>
#include "util.hpp"
#include "config.hpp"
#include <stdlib.h> 
#include<random>

#pragma once

/**
 * @brief Quantum oscillator (multi-level qubit) with control capabilities.
 *
 * This class represents a single quantum oscillator that is controlled by external 
 * control pulses. It stores the oscillator parameters, such as number of energy levels, 
 * frequency detuning, anharmonicity and Lindblad decay and dephasing times.  
 * It also manages this oscillator's control pulse parameterization and carrier wave frequencies.
 * 
 * Main functionality:
 *    - @ref evalControl computes the rotating-frame drive pulses p(t) & q(t) as well as the flux control f(t) at a given time t. 
 *      The drive pulses are products of fixed-frequency carrier waves multiplied with an outer envelope
 *      (e.g. spline) represented by @ref drive_basisfunctions_re and @ref drive_basisfunctions_im, for the coefficients alpha^1 and alpha^2, respectively.
 *      The flux control is a single scalar channel represented by @ref flux_basisfunctions.
 *    - @ref expectedEnergy and @ref population for computing this oscillators expected Energy and level occupations
 *      given a current state
 *    - @ref evalControlVariation for evaluating control-variation penalties used in optimization
 * 
 * This class stores @ref ControlBasis objects for evaluating control pulse envelopes at runtime.
 */
class Oscillator {
  protected:
    int myid; ///< Integer identifier for this oscillator
    size_t nlevels; ///< Number of energy levels for this oscillator
    double ground_freq; ///< Fundamental 0-1 transition frequency of this oscillator \f$\omega_k\f$
    double selfkerr; ///< Self-Kerr frequency \f$\xi_k\f$, multiplies \f$a_k^\dagger a_k^\dagger a_k a_k\f$

    double detuning_freq; ///< Detuning frequency, detuning = ground_freq - rotational_freq, multiplies \f$a_k^\dagger a_k\f$
    DecoherenceType decoherence_type; ///< Type of Lindblad decoherence operators to include
    double decay_time; ///< Characteristic time for T1 decay operations
    double dephase_time; ///< Characteristic time for T2 dephasing operations

    double total_time; ///< Final evolution time
    ControlBasis* drive_basisfunctions_re; ///< Real part of drive envelope basis (per carrier), alpha^1
    ControlBasis* drive_basisfunctions_im; ///< Imaginary part of drive envelope basis (per carrier), alpha^2
    ControlBasis* flux_basisfunctions; ///< Flux control parameterization
    std::vector<double> carrier_freq; ///< Frequencies of the carrier waves for this oscillator

    int mpirank_world; ///< Rank of MPI_COMM_WORLD
    int mpirank_petsc; ///< Rank of PETSc's communicator
    int mpisize_petsc; ///< Size of PETSc's communicator
    PetscInt localsize_u; ///< Size of local sub vector u or v in state x=[u,v]
    PetscInt ilow; ///< First index of the local sub vector u,v
    PetscInt iupp; ///< Last index (+1) of the local sub vector u,v

  public:
    PetscInt dim_preOsc; ///< Dimension of coupled subsystems preceding this oscillator
    PetscInt dim_postOsc; ///< Dimension of coupled subsystems following this oscillator

    Oscillator();

    /**
     * @brief Constructor with simplified configuration-based specification.
     *
     * @param config Configuration parameters containing all oscillator settings
     * @param id Oscillator identifier
     * @param rand_engine Random number generator engine
     * @param param_offset Offset for global control parameter indexing
     * @param quietmode Flag for quiet mode operation
     */
    Oscillator(const Config& config, size_t id, std::mt19937& rand_engine, int param_offset, bool quietmode);

    virtual ~Oscillator();

    /**
     * @brief Retrieves the number of control parameters.
     *
     * @return size_t Number of control parameters
     */
    size_t getNParams() { return getNDriveParams() + getNFluxParams(); };
    size_t getNFluxParams() { return flux_basisfunctions ? static_cast<size_t>(flux_basisfunctions->getNparams()) : 0; };
    size_t getNDriveParams() {
      size_t nparams = 0;
      if (drive_basisfunctions_re) nparams += static_cast<size_t>(drive_basisfunctions_re->getNparams());
      if (drive_basisfunctions_im) nparams += static_cast<size_t>(drive_basisfunctions_im->getNparams());
      return nparams;
    };

    /**
     * @brief Retrieves the number of energy levels.
     *
     * @return size_t Number of energy levels
     */
    size_t getNLevels() { return nlevels; };

    /**
     * @brief Retrieves the self-Kerr coefficient.
     *
     * @return double Self-Kerr frequency
     */
    double getSelfkerr() { return selfkerr; }; 

    /**
     * @brief Retrieves the detuning frequency.
     *
     * @return double Detuning frequency (rad/time)
     */
    double getDetuning() { return detuning_freq; }; 

    /**
     * @brief Retrieves the T1 decay time.
     *
     * @return double Decay time
     */
    double getDecayTime() {return decay_time; };

    /**
     * @brief Retrieves the T2 dephasing time.
     *
     * @return double Dephasing time
     */
    double getDephaseTime() {return dephase_time; };

    /**
     * @brief Retrieves the number of carrier frequencies.
     *
     * @return size_t Number of carrier frequencies
     */
    size_t getNCarrierfrequencies() {return carrier_freq.size(); };

    /**
     * @brief Retrieves the rotating frame frequency.
     *
     * @return double Rotating frame frequency (Hz)
     */
    double getRotFreq() {return (ground_freq - detuning_freq) / (2.0*M_PI); };

    /**
     * @brief Sets control parameters from a global storage.
     *
     * @param x Array of control parameter values
     */
    void setControlParams(const double* x);

    // Backward-compatible alias used by legacy call sites.
    void setParams(const double* x) { setControlParams(x); }

    /**
     * @brief Retrieves control parameters into a global storage.
     *
     * @param x Array to store control parameter values
     */
    void getControlParams(double* x);

    // Backward-compatible alias used by legacy call sites.
    void getParams(double* x) { getControlParams(x); }

    /**
     * @brief Evaluates drive-only controls p(t), q(t).
      *
      * @param[in] t Time at which to evaluate
      * @param[out] p_ptr Pointer to store p(t)
      * @param[out] q_ptr Pointer to store q(t)
      * @return int Error code
     */
    int evalDriveControl(const double t, double* p_ptr, double* q_ptr);

    /**
     * @brief Evaluates the rotating-frame drives p,q and flux control functions.
     *
     * Computes p(t), q(t), and f(t), where p(t) and q(t) are the real and imaginary 
     * parts of the drive control function and f(t) is the flux control function.
     *
     * @param[in] t Time at which to evaluate
     * @param[out] p_ptr Pointer to store real part, p(t)
     * @param[out] q_ptr Pointer to store imaginary part, q(t)
     * @param[out] flux_ptr Pointer to store flux control, f(t)
     * @return int Error code
     */
    int evalControl(const double t, double* p_ptr, double* q_ptr, double* flux_ptr);

    /**
     * @brief Computes derivatives of drive control functions p(t) and q(t) with respect to the parameters.
     *
     * @param[in] t Time at which to evaluate derivatives
     * @param[out] grad_for_this_oscillator Array to update the gradient
     * @param[in] pbar Adjoint scaling factor for the gradient of p (seed) 
     * @param[in] qbar Adjoint scaling factor for the gradient of q (seed) 
     * @return int Error code
     */
    int evalDriveControl_diff(const double t, double* grad_for_this_oscillator, const double pbar, const double qbar);

    /**
     * @brief Computes derivatives of drive and flux control functions p(t), q(t), and f(t) with respect to parameters.
     *
     * @param[in] t Time at which to evaluate derivatives
     * @param[out] grad_for_this_oscillator Array to update the gradient
     * @param[in] pbar Adjoint seed for p(t)
     * @param[in] qbar Adjoint seed for q(t)
     * @param[in] fbar Adjoint seed for f(t)
     * @return int Error code
     */
    int evalControl_diff(const double t, double* grad_for_this_oscillator, const double pbar, const double qbar, const double fbar);

    /**
     * @brief Computes expected energy for this oscillator.
     *
     * Returns the expected value of the number operator for this oscillator's subsystem.
     *
     * @param x State vector
     * @return double Expected energy value
     */
    double expectedEnergy(const Vec x);

    /**
     * @brief Computes derivative of expected energy computation.
     *
     * @param x State vector
     * @param x_bar Adjoint state vector to update
     * @param obj_bar Adjoint of expected energy
     */
    void expectedEnergy_diff(const Vec x, Vec x_bar, const double obj_bar);

    /**
     * @brief Computes population in each energy level of this oscillator.
     *
     * Extracts the diagonal elements of the reduced density matrix for this oscillator.
     *
     * @param x State vector
     * @param pop Reference to vector to store population values
     */
    void population(const Vec x, std::vector<double> &pop);

    /**
     * @brief Evaluates control variation penalty term.
     *
     * Computes finite-difference based regularization of control parameters.
     *
     * @return double Control variation penalty value
     */
    double evalControlVariation();

    /**
     * @brief Computes derivative of control variation penalty.
     *
     * @param G Gradient vector to update
     * @param var_reg_bar Adjoint of variation penalty
     * @param skip_to_oscillator Offset to this oscillator's parameters in global vector
     */
    void evalControlVariationDiff(Vec G, double var_reg_bar, int skip_to_oscillator);
};



