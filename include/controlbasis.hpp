#include "defs.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cassert>
#include "util.hpp"
#pragma once


/**
 * @brief Abstract base class for control parameterizations.
 *
 * The `ControlBasis` class defines the interface for various parameterizations of 
 * the control pulse envelopes, which are multiplied by carrier waves
 * in the Oscillator class. Derived classes implement specific control parameterizations, 
 * such as the most standard parameterization via 2nd order Bsplines. Specific control 
 * parameterizations are initialized in the constructor of the oscillator. 
 * 
 * Main functionality:
 *      - @ref evaluate        for computing the outer envelop shape of the pulses at given time t
 *      - @ref derivative      for updating the local gradient of the @ref evaluate function
 *      - @ref enforceBoundary for setting pulse shape envelopes to zero at initial and final time
 */
class ControlBasis {
    protected:
        int nparams; ///< Number of parameters that define the control pulse.
        double tstart; ///< Start time of the interval where the control basis is applied.
        double tstop; ///< Stop time of the interval where the control basis is applied.
        int skip; ///< Offset to the starting location for this basis inside the global control vector.
        ControlType controltype; ///< Type of control parameterization.
        bool enforceZeroBoundary; ///< Flag to enforce zero boundary conditions for control pulses.

    public: 
        ControlBasis();

        /**
         * @brief Constructor with parameters.
         *
         * @param nparams_ Number of parameters defining the controls.
         * @param tstart Start time of the interval.
         * @param tstop Stop time of the interval.
         * @param enforceZeroBoundary Flag to enforce zero boundary conditions. If true, ensures that the controls start and end at zero.
         */
        ControlBasis(int nparams_, double tstart, double tstop, bool enforceZeroBoundary);

        virtual ~ControlBasis();

        /**
         * @brief Retrieves the number of parameters defining the controls.
         *
         * @return int Number of parameters.
         */
        int getNparams() {return nparams; };

        /**
         * @brief Retrieves the start time of the interval.
         *
         * @return double Start time.
         */
        double getTstart() {return tstart; };

        /**
         * @brief Retrieves the stop time of the interval.
         *
         * @return double Stop time.
         */
        double getTstop() {return tstop; };

        /**
         * @brief Retrieves the type of control parameterization.
         *
         * @return ControlType Type of control.
         */
        ControlType getType() {return controltype;};

        /**
         * @brief Sets the offset for the starting location in the global control vector.
         *
         * @param skip_ Offset value.
         */
        void setSkip(int skip_) {skip = skip_;};

        /**
         * @brief Retrieves the offset for the starting location in the global control vector.
         *
         * @return int Offset value.
         */
        int getSkip(){return skip;};

        /**
         * @brief Retrieves the number of splines (default implementation returns 0).
         *
         * @return int Number of splines.
         */
        virtual int getNSplines() {return 0;};

        /**
         * @brief Computes the variation of control parameters (default implementation returns 0.0).
         *
         * Default implementation ignores all input parameters.
         *
         * @return double Variation value.
         */
        virtual double computeVariation(std::vector<double>& /*params*/, int /*carrierfreqID*/){return 0.0;};

        /**
         * @brief Computes the gradient of the variation (default implementation does nothing).
         *
         * Default implementation ignores all input parameters.
         */
        virtual void computeVariation_diff(double* /*grad*/, std::vector<double>& /*params*/, double /*var_bar*/, int /*carrierfreqID*/){};

        /**
         * @brief Enforces boundary conditions for controls (default implementation does nothing).
         *
         * For some control parameterizations, this can be used to enforce that the controls start and end at zero.
         * For example, the 2nd order Bspline parameterization  will overwrite the parameters of the first and last 
         * two splines by zero, ensuring that the resulting control pulses start and end at zero.
         *
         * Default implementation ignores all input parameters.
         */
        virtual void enforceBoundary(double* /*x*/, int /*carrier_id*/) {};

        /**
         * @brief Evaluates the control basis at a given time using the provided coefficients.
         *
         * @param[in] t Time at which to evaluate.
         * @param[in] coeff Vector of parameters (coefficients of the basis parameterization).
         * @param[in] carrier_freq_id ID of the carrier frequency, provided by the oscillator.
         * @param[out] Blt1 Pointer to store the real part of the control parameterization.
         * @param[out] Blt2 Pointer to store the imaginary part of the control parameterization.
         */
        virtual void evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Blt1, double*Blt2) = 0;

        /**
         * @brief Evaluates the derivative of the control basis at a given time.
         *
         * @param t Time at which to evaluate.
         * @param coeff Vector of coefficients.
         * @param coeff_diff Pointer to the derivative coefficients.
         * @param valbar1 Multiplier for the real derivative term.
         * @param valbar2 Multiplier for the imaginary derivative term.
         * @param carrier_freq_id ID of the carrier frequency.
         */
        virtual void derivative(const double t, const std::vector<double>& coeff, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id)= 0;
};

/**
 * @brief Control parameterization using quadratic (2nd order) Bspline basis functions.
 *
 * Bspline basis functions have local support with width = 3*dtknot,
 * where dtknot = T/(nsplines -2) is the time knot vector spacing.
 */
class BSpline2nd : public ControlBasis {
    protected:
        int nsplines; ///< Number of splines.
        double dtknot; ///< Spacing of time knot vector.
        double *tcenter; ///< Vector of basis function center positions.
        double width; ///< Support of each basis function (m*dtknot).

        /**
         * @brief Evaluate one basis function B_i(tau_i(t)).
         *
         * @param id Index of the basis function.
         * @param t Time at which to evaluate.
         * @return double Value of the basis function.
         */
        double basisfunction(int id, double t);

    public:
        /**
         * @brief Constructor for quadratic Bsplines.
         *
         * @param nsplines Number of splines.
         * @param tstart Start time of the interval.
         * @param tstop Stop time of the interval.
         * @param enforceZeroBoundary Flag to enforce zero boundary conditions.
         */
        BSpline2nd(int nsplines, double tstart, double tstop, bool enforceZeroBoundary);

        ~BSpline2nd();
        
        int getNSplines() {return nsplines;};

        /**
         * @brief Sets the first and last two spline coefficients in x to zero for this carrier wave, 
         * so that the controls start and end at zero.
         *
         * @param x Pointer to the control parameters.
         * @param carrier_id ID of the carrier wave.
         */
        void enforceBoundary(double* x, int carrier_id);

        void evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Blt1_ptr, double* Blt2_ptr);

        void derivative(const double t, const std::vector<double>& coeff, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id);
};

/**
 * @brief Control parameterization where only the pulse amplitude is parameterized 
 * by Bsplines, while the phase is time-independent.
 *
 * Discretization of the controls amplitudes using quadratic Bsplines.
 * Bspline basis functions have local support with width = 3 * dtknot, 
 * where dtknot = T / (nsplines - 2) is the time knot vector spacing.
 */
class BSpline2ndAmplitude : public ControlBasis {
    protected:
        int nsplines; ///< Number of splines.
        double dtknot; ///< Spacing of the time knot vector.
        double *tcenter; ///< Vector of basis function center positions.
        double width; ///< Support of each basis function (m * dtknot).
        double scaling; ///< Scaling for the phase.

        /**
         * @brief Evaluate one basis function B_i(tau_i(t)).
         *
         * @param id Index of the basis function.
         * @param t Time at which to evaluate.
         * @return double Value of the basis function.
         */
        double basisfunction(int id, double t);

    public:
        /**
         * @brief Constructor for quadratic Bsplines for amplitude parameterization.
         *
         * @param nsplines Number of splines.
         * @param scaling Scaling factor for the phase.
         * @param tstart Start time of the interval.
         * @param tstop Stop time of the interval.
         * @param enforceZeroBoundary Flag to enforce zero boundary conditions.
         */
        BSpline2ndAmplitude(int nsplines, double scaling, double tstart, double tstop, bool enforceZeroBoundary);

        ~BSpline2ndAmplitude();

        int getNSplines() {return nsplines;};

        /**
         * @brief Sets the first and last two spline coefficients in x to zero, so that the controls start and end at zero.
         *
         * @param x Pointer to the control parameters.
         * @param carrier_id ID of the carrier wave.
         */
        void enforceBoundary(double* x, int carrier_id);

        void evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Blt1_ptr, double* Blt2_ptr);

        void derivative(const double t, const std::vector<double>& coeff, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id);
};

/**
 * @brief Discretization of the Controls using piece-wise constant (0-th order) Bsplines.
 *
 * This class parameterizes the control pulse using 0-th order Bspline basis functions 
 * (hat functions), with local support of width = T/nsplines.
 */
class BSpline0 : public ControlBasis {
    protected:
        int nsplines; ///< Number of splines.
        double dtknot; ///< Spacing of time knot vector.
        double width; ///< Support of each basis function (m*dtknot).

    public:
        BSpline0(int nsplines, double tstart, double tstop, bool enforceZeroBoundary);
        ~BSpline0();

        int getNSplines() {return nsplines;};

        /**
         * @brief Sets the first and last parameter to zero for this carrier wave, 
         * so that the controls start and end at zero.
         *
         * @param x Pointer to the control parameters array
         * @param carrier_id ID of the carrier wave
         */
        void enforceBoundary(double* x, int carrier_id);

        /**
         * @brief Computes total variation of the control parameters.
         *
         * Computes \f$\frac{1}{n_{splines}} \sum_{splines} (\alpha_i - \alpha_{i-1})^2\f$.
         *
         * @param params Vector of control parameters
         * @param carrierfreqID ID of the carrier frequency
         * @return double Variation value
         */
        double computeVariation(std::vector<double>& params, int carrierfreqID);

        /**
         * @brief Computes derivative of control parameter variation.
         *
         * @param grad Pointer to gradient array to update
         * @param params Vector of control parameters
         * @param var_bar Adjoint of variation term
         * @param carrierfreqID ID of the carrier frequency
         */
        virtual void computeVariation_diff(double* grad, std::vector<double>&params, double var_bar, int carrierfreqID);

        void evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Blt1_ptr, double* Blt2_ptr);

        void derivative(const double t, const std::vector<double>& coeff, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id);
};

