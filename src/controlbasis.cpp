#include "controlbasis.hpp"

namespace {
constexpr double kControlTimeTolerance = 1e-10;
}


ControlBasis::ControlBasis() {
    tstart = 0.0;
    tstop = 0.0;
    controltype = ControlType::NONE;
    enforceZeroBoundary = false;
    params.resize(0);
}

ControlBasis::ControlBasis(int nbasisfunctions, int ncarrier, double tstart_, double tstop_, bool control_zero_boundary_condition_) : ControlBasis() {
    tstart = tstart_;
    tstop = tstop_;
    enforceZeroBoundary = control_zero_boundary_condition_;

    // Allocate one vector of parameters for each carrier wave
    params.resize(ncarrier);
    for (int i=0; i<ncarrier; i++) {
        params[i].resize(nbasisfunctions, 0.0);
    }
}

ControlBasis::~ControlBasis(){}

void ControlBasis::setParams(const double* x, int carrier_freq_id) {
    if (static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return;
    }
    for (size_t i = 0; i < params[carrier_freq_id].size(); i++){
        params[carrier_freq_id][i] = x[i];
    }
}

void ControlBasis::getParams(double* x, int carrier_freq_id) {
    if (static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return;
    }
    for (size_t i = 0; i < params[carrier_freq_id].size(); i++){
        x[i] = params[carrier_freq_id][i];
    }
}

BSpline2nd::BSpline2nd(int nsplines, int ncarrier, double t0, double T, bool enforceZeroBoundary_) : ControlBasis(nsplines, ncarrier, t0, T, enforceZeroBoundary_) {
    controltype = ControlType::BSPLINE;

    dtknot = (T-t0) / (double)(nsplines - 2);
	width = 3.0*dtknot;

    /* Compute center points of the splines */
    tcenter = new double[nsplines];
    for (int i = 0; i < nsplines; i++){
        tcenter[i] = t0 + dtknot * ( (i+1) - 1.5 );
    }
}

BSpline2nd::~BSpline2nd(){
    delete [] tcenter;
}

void BSpline2nd::enforceBoundary(){

    if (!enforceZeroBoundary) {
        return;
    }

    // For each carrier wave, set first and last two splines to zero so that spline starts and ends at zero 
    for (size_t f = 0; f < params.size(); f++) {
        const size_t nsplines = params[f].size();
        if (nsplines == 0) {
            continue;
        }
        if (nsplines < 4) {
            for (size_t l = 0; l < nsplines; l++) {
                params[f][l] = 0.0;
            }
            continue;
        }
        params[f][0] = 0.0;
        params[f][1] = 0.0;
        params[f][nsplines-2] = 0.0;
        params[f][nsplines-1] = 0.0;
    }
}

double BSpline2nd::evaluate(const double t, int carrier_freq_id){

    if (carrier_freq_id < 0 || static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return 0.0;
    }
    if (params[carrier_freq_id].size() == 0) {
        return 0.0;
    }
    if (tstart - kControlTimeTolerance > t || tstop + kControlTimeTolerance < t ) {
        return 0.0;
    }

    double sum = 0.0;
    /* Sum over basis function */
    for (size_t l = 0; l < params[carrier_freq_id].size(); l++) {
        if (enforceZeroBoundary && (l <= 1 || l >= params[carrier_freq_id].size() - 2)) {
          continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero 
        }
        double alpha = params[carrier_freq_id][l];
        sum += alpha * basisfunction(l,t);
    }
    return sum;
}

void BSpline2nd::derivative(const double t, int carrier_freq_id, double* grad, const double valbar) {

    if (carrier_freq_id < 0 || static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return;
    }
    if (params[carrier_freq_id].size() == 0) {
        return;
    }
    if (tstart - kControlTimeTolerance > t || tstop + kControlTimeTolerance < t ) {
        return;
    }

    /* Iterate over basis function */
    for (size_t l = 0; l < params[carrier_freq_id].size(); l++) {
        if (enforceZeroBoundary){
            if (l <= 1 || l >= params[carrier_freq_id].size() - 2) continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero
        }
        double Blt = basisfunction(static_cast<int>(l), t);
        grad[static_cast<int>(l)] += Blt * valbar;
    }
}

double BSpline2nd::basisfunction(int id, double t){

    /* compute scaled time tau = (t-tcenter[k])  */
    double tau = (t - tcenter[id]) / width;

    /* Return 0 if tau not in local support */
    if ( tau < -1./2. || tau >= 1./2. ) return 0.0;

    /* Evaluate basis function */
    double val = 0.0;
    if       (-1./2. <= tau && tau < -1./6.) val = 9./8. + 9./2. * tau + 9./2. * pow(tau,2);
    else if  (-1./6. <= tau && tau <  1./6.) val = 3./4. - 9. * pow(tau,2);
    else if  ( 1./6. <= tau && tau <  1./2.) val = 9./8. - 9./2. * tau + 9./2. * pow(tau,2);

    return val;
}



// Zeroth order B-splines, i.e., piecewise constant
BSpline0::BSpline0(int nbasisfunctions, int ncarrier, double t0, double T, bool enforceZeroBoundary_) : ControlBasis(nbasisfunctions, ncarrier, t0, T, enforceZeroBoundary_) {
    controltype = ControlType::BSPLINE0;

    dtknot = (nbasisfunctions > 1) ? (T - t0) / (nbasisfunctions - 1.0) : 0.0;
	width = dtknot;
}

BSpline0::~BSpline0(){}


double BSpline0::evaluate(const double t, int carrier_freq_id){

    if (carrier_freq_id < 0 || static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return 0.0;
    }
    if (params[carrier_freq_id].size() == 0) {
        return 0.0;
    }
    if (tstart - kControlTimeTolerance > t || tstop + kControlTimeTolerance < t ) {
        return 0.0;
    }
    if (params[carrier_freq_id].size() == 1 || dtknot <= 0.0) {
        return params[carrier_freq_id][0];
    }

    // Figure out which basis function is active at this time point 
    int splineID = static_cast<int>(ceil((t-tstart)/dtknot - 0.5));

    // Control function is defined to be zero outside [tstart, tstop].
    if (splineID < 0 || splineID >= static_cast<int>(params[carrier_freq_id].size())){
        return 0.0;
    } else {
        return params[carrier_freq_id][splineID];
    }
}

void BSpline0::derivative(const double t, int carrier_freq_id, double* grad, const double valbar) {
    if (carrier_freq_id < 0 || static_cast<size_t>(carrier_freq_id) >= params.size()) {
        return;
    }
    if (params[carrier_freq_id].size() == 0) {
        return;
    }
    if (tstart - kControlTimeTolerance > t || tstop + kControlTimeTolerance < t ) {
        return;
    }
    if (params[carrier_freq_id].size() == 1 || dtknot <= 0.0) {
        grad[carrier_freq_id] += valbar;
        return;
    }


    // Figure out which basis function is active at this time point 
    int splineID = static_cast<int>(ceil((t-tstart)/dtknot - 0.5));

    if (splineID >= 0 && splineID < static_cast<int>(params[carrier_freq_id].size())){
        grad[splineID] += valbar;
    }
}


double BSpline0::computeVariation(){
    if (params.size() == 0) {
        return 0.0;
    }

    // Iterate over carrier frequencies. 
    double var = 0.0;
    for (size_t f=0; f < params.size(); f++) {

        for (size_t lc = 1; lc < params[f].size(); lc++){
            var += SQR(params[f][lc] - params[f][lc - 1]);
        }

        if (enforceZeroBoundary) {
            var += SQR(params[f][0]); // lc = 0
            var += SQR(params[f][params[f].size() - 1]); // lc = nsplines
        }
    }
    return var;
}


void BSpline0::computeVariation_diff(double* grad, double var_bar){

    double fact = 2.0*var_bar;

    for (size_t f=0; f < params.size(); f++) {
        int nsplines = params[f].size();
        if (nsplines == 0) {
            continue;
        }
        if (nsplines == 1) {
            if (enforceZeroBoundary) {
                grad[f*nsplines] += fact * params[f][0];
            }
            continue;
        }
        int lc = 0;
        grad[f*nsplines + lc] += fact * (params[f][lc] - params[f][lc + 1]);
        // interior lc
        for (size_t l = 1; l < params[f].size() - 1; l++){
            grad[f*nsplines + l] += fact * (2*params[f][l] - params[f][l - 1] - params[f][l + 1]);
        }
        lc = params[f].size() - 1;
        grad[f*nsplines + lc] += fact * (params[f][lc] - params[f][lc - 1]);
        
        if (enforceZeroBoundary) {
            grad[f*nsplines + 0] += fact * params[f][0];
            grad[f*nsplines + params[f].size()-1] += fact * params[f][params[f].size()-1];
        }
    }
}

void BSpline0::enforceBoundary(){

    if (!enforceZeroBoundary) {
        return;
    }

    for (size_t f = 0; f < params.size(); f++) {
        if (params[f].empty()) {
            continue;
        }
        params[f][0] = 0.0; 
        params[f][params[f].size() - 1] = 0.0;
    }
}
