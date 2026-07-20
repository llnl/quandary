#include "defs.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cassert>
#include "util.hpp"
#pragma once


/**
 * @brief Control parameterization interface and empty fallback implementation.
 *
 * The `ControlBasis` class defines the interface for various parameterizations of 
 * the control pulse envelopes, which are multiplied by carrier waves
 * in the Oscillator class. Derived classes implement specific control parameterizations, 
 * such as the most standard parameterization via 2nd order Bsplines. Specific control 
 * parameterizations are initialized in the constructor of the oscillator. 
 * The base class itself is also usable as a concrete empty basis that stores no
 * parameters and evaluates to zero.
 * 
 * Main functionality:
 *      - @ref evaluate        for computing the outer envelop shape of the pulses at given time t
 *      - @ref derivative      for updating the local gradient of the @ref evaluate function
 *      - @ref enforceBoundary for setting pulse shape envelopes to zero at initial and final time
 */
class ControlBasis {
    protected:
        std::vector<std::vector<double>> params; ///< Coefficients of the control parameterizations, one for each carrier wave
        double tstart; ///< Start time of the interval where the control basis is applied.
        double tstop; ///< Stop time of the interval where the control basis is applied.
        ControlType controltype; ///< Type of control parameterization.
        bool enforceZeroBoundary; ///< Flag to enforce zero boundary conditions for control pulses.

    public: 
        ControlBasis();

        /**
         * @brief Constructor with parameters.
         *
         * @param nbasisfunctions Number of basis functions for the control parameterization.
         * @param ncarrier Number of carrier waves for the control parameterization.
         * @param tstart_ Start time of the interval.
         * @param tstop_ Stop time of the interval.
         * @param control_zero_boundary_condition_ Flag to enforce zero boundary conditions. If true,
         * ensures that the controls start and end at zero.
         */
        ControlBasis(int nbasisfunctions, int ncarrier, double tstart_, double tstop_, bool control_zero_boundary_condition_);

        virtual ~ControlBasis();

        /**
         * @brief Retrieves the total number of parameters defining the controls.
         *
         * @return int Number of parameters.
         */
        int getNparams() { 
            int nparams = 0;
            for (const auto& p : params) {
                nparams += p.size();
            }
            return nparams;
        };

        /**
         * @brief Retrieves the number of parameters for a specific carrier wave.
         *
         * @param f Index of the carrier wave.
         * @return int Number of parameters for the specified carrier wave.
         */
        int getNparams(size_t f) { return params.size() > f ? params[f].size() : 0; };

        /**
         * @brief Set parameter coefficients
         * @param x Pointer to the array of parameter coefficients.
         * @param carrier_freq_id ID of the carrier frequency for which to set the parameters.
         */
        void setParams(const double* x, int carrier_freq_id);

        /**
         * @brief Get parameter coefficients
         * @param x Pointer to the array to store parameter coefficients.
         * @param carrier_freq_id ID of the carrier frequency for which to get the parameters.
         */
        void getParams(double* x, int carrier_freq_id);

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
         * @brief Computes the variation of control parameters (default implementation returns 0.0).
         *
         * Default implementation ignores all input parameters.
         *
         * @return double Variation value.
         */
        virtual double computeVariation(){return 0.0;};

        /**
         * @brief Computes the gradient of the variation (default implementation does nothing).
         *
         * Default implementation ignores all input parameters.
         */
        virtual void computeVariation_diff(double* /*grad*/,  double /*var_bar*/){};

        /**
         * @brief Enforces boundary conditions for controls (default implementation does nothing).
         *
         * For some control parameterizations, this can be used to enforce that the controls start and end at zero.
         * For example, the 2nd order Bspline parameterization  will overwrite the parameters of the first and last 
         * two splines by zero, ensuring that the resulting control pulses start and end at zero.
         *
         * Default implementation ignores all input parameters.
         */
        virtual void enforceBoundary() {};

        /**
         * @brief Evaluates the control basis at a given time using the provided coefficients.
         *
         * @param[in] t Time at which to evaluate.
         * @param[in] carrier_freq_id ID of the carrier frequency, provided by the oscillator.
         * @return Evaluated basis expansion value at time t.
         */
        virtual double evaluate(const double t, int carrier_freq_id) { (void)t; (void)carrier_freq_id; return 0.0; }

        /**
         * @brief Evaluates the derivative of the control basis at a given time.
         *
         * @param t Time at which to evaluate.
         * @param carrier_freq_id ID of the carrier frequency.
         * @param grad Pointer to the derivative coefficients.
         * @param valbar Multiplier for the derivative term.
         */
        virtual void derivative(const double t, int carrier_freq_id, double* grad, const double valbar) { (void)t; (void)carrier_freq_id; (void)grad; (void)valbar; }
};

/**
 * @brief Control parameterization using quadratic (2nd order) Bspline basis functions.
 *
 * Bspline basis functions have local support with width = 3*dtknot,
 * where dtknot = T/(nsplines -2) is the time knot vector spacing.
 */
class BSpline2nd : public ControlBasis {
    protected:
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
         * @param ncarrier Number of carrier frequencies represented by this basis object.
         * @param tstart Start time of the interval.
         * @param tstop Stop time of the interval.
         * @param enforceZeroBoundary Flag to enforce zero boundary conditions.
         */
        BSpline2nd(int nsplines, int ncarrier, double tstart, double tstop, bool enforceZeroBoundary);

        ~BSpline2nd();
        
        /**
         * @brief Sets the first and last two spline coefficients to zero for each carrier wave,
         * so that the controls start and end at zero.
         */
        void enforceBoundary();

        double evaluate(const double t, int carrier_freq_id);

        void derivative(const double t, int carrier_freq_id, double* coeff_diff, const double valbar);
};

/**
 * @brief Discretization of the Controls using piece-wise constant (0-th order) Bsplines.
 *
 * This class parameterizes the control pulse using 0-th order Bspline basis functions 
 * (hat functions), with local support of width = T/nsplines.
 */
class BSpline0 : public ControlBasis {
    protected:
        double dtknot; ///< Spacing of time knot vector.
        double width; ///< Support of each basis function (m*dtknot).

    public:
        BSpline0(int nsplines, int ncarrier, double tstart, double tstop, bool enforceZeroBoundary);
        ~BSpline0();

        /**
         * @brief Sets the first and last parameter to zero for this carrier wave, 
         * so that the controls start and end at zero.
         */
        void enforceBoundary();

        /**
         * @brief Computes total variation of the control parameters.
         *
         * Computes \f$\frac{1}{n_{splines}} \sum_{splines} (\alpha_i - \alpha_{i-1})^2\f$.
         *
         * @return double Variation value
         */
        double computeVariation();

        /**
         * @brief Computes derivative of control parameter variation.
         *
         * @param grad Pointer to gradient array to update
         * @param var_bar Adjoint of variation term
         */
        virtual void computeVariation_diff(double* grad, double var_bar);

        double evaluate(const double t, int carrier_freq_id);

        void derivative(const double t, int carrier_freq_id, double* coeff_diff, const double valbar);
};

