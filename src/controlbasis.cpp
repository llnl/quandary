#include "controlbasis.hpp"


ControlBasis::ControlBasis() {
    nparams= 0;
    nsplines = 0;
    skip = 0;
    controltype = ControlType::NONE;
    enforceZeroBoundary = false;
}

ControlBasis::ControlBasis(int nparams_, double tstart_, double tstop_, bool enforceZeroBoundary_) : ControlBasis() {
    nparams = nparams_;
    tstart = tstart_;
    tstop = tstop_;
    enforceZeroBoundary = enforceZeroBoundary_;
}
ControlBasis::~ControlBasis(){}

double ControlBasis::computeVariation(std::vector<double>& params, int carrierfreqID){
    double var = 0.0;
    //   Re params
    for (int lc=1; lc<nsplines; lc++){
        var += SQR(params[skip + 2*carrierfreqID*nsplines + lc] - params[skip+ 2*carrierfreqID*nsplines + lc - 1]);
    }
    // Im params
    for (int lc=1; lc<nsplines; lc++){
        var += SQR(params[skip + (2*carrierfreqID+1)*nsplines + lc] - params[skip + (2*carrierfreqID+1)*nsplines + lc - 1]);
    }

    if (enforceZeroBoundary) {
        // Re
        var += SQR(params[skip + 2*carrierfreqID*nsplines + 0 ]); // lc = 0
        var += SQR(params[skip+ 2*carrierfreqID*nsplines + nsplines - 1]); // lc = nsplines
        // Im
        var += SQR(params[skip + (2*carrierfreqID+1)*nsplines ]); // lc=0
        var += SQR(params[skip + (2*carrierfreqID+1)*nsplines + nsplines- 1]); // lc=nsplines
    }
    return var;
}


void ControlBasis::computeVariation_diff(double* grad, std::vector<double>&params, double var_bar, int carrierfreqID){

    double fact = 2.0*var_bar;

    // Re params
    int lc = 0;
    grad[skip + 2*carrierfreqID*nsplines + lc] += fact * (params[skip + 2*carrierfreqID*nsplines + lc] - params[skip + 2*carrierfreqID*nsplines + lc + 1]);
    // interior lc
    for (lc=1; lc<nsplines-1; lc++){
      grad[skip+ 2*carrierfreqID*nsplines + lc] += fact * (2*params[skip + 2*carrierfreqID*nsplines + lc] - params[skip + 2*carrierfreqID*nsplines + lc - 1] - params[skip + 2*carrierfreqID*nsplines + lc + 1]);
    }
    lc = nsplines-1;
    grad[skip + 2*carrierfreqID*nsplines + lc] += fact * (params[skip + 2*carrierfreqID*nsplines + lc] - params[skip + 2*carrierfreqID*nsplines + lc - 1]);
    // Im params
    lc = 0;
    grad[skip + (2*carrierfreqID+1)*nsplines + lc] += fact * (params[skip + (2*carrierfreqID+1)*nsplines + lc] - params[skip + (2*carrierfreqID+1)*nsplines + lc + 1]);
    // interior lc
    for (int lc=1; lc<nsplines-1; lc++){
      grad[skip + (2*carrierfreqID+1)*nsplines + lc] += fact * (2*params[skip + (2*carrierfreqID+1)*nsplines + lc] - params[skip + (2*carrierfreqID+1)*nsplines + lc - 1] - params[skip + (2*carrierfreqID+1)*nsplines + lc + 1]);
    }
    lc = nsplines-1;
    grad[skip + (2*carrierfreqID+1)*nsplines + lc] += fact * (params[skip + (2*carrierfreqID+1)*nsplines + lc] - params[skip + (2*carrierfreqID+1)*nsplines + lc - 1]);


    if (enforceZeroBoundary) {
        // Re
        grad[skip + 2*carrierfreqID*nsplines ] += fact * params[skip + 2*carrierfreqID*nsplines ];
        grad[skip + 2*carrierfreqID*nsplines + nsplines-1] += fact * params[skip + 2*carrierfreqID*nsplines + nsplines-1];
        // Im
        grad[skip + 2*carrierfreqID*nsplines + nsplines] += fact * params[skip + 2*carrierfreqID*nsplines + nsplines];
        grad[skip + 2*carrierfreqID*nsplines + 2*nsplines-1] += fact * params[skip + 2*carrierfreqID*nsplines + 2*nsplines-1];
    }
}

BSpline2nd::BSpline2nd(int nsplines_, double t0, double T, bool enforceZeroBoundary_) : ControlBasis(2*nsplines_, t0, T, enforceZeroBoundary_){
    nsplines = nsplines_;
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

void BSpline2nd::enforceBoundary(double* x, int carrier_id){
    // set first and last two splines to zero so that spline starts and ends at zero 
    for (int l=0; l<nsplines; l++) {
        if (l<=1 || l >= nsplines- 2) {
            x[skip + carrier_id*nsplines*2 + l] = 0.0;
            x[skip + carrier_id*nsplines*2 + l + nsplines] = 0.0;
        }
    }
}

void BSpline2nd::evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Bl1_ptr, double* Bl2_ptr){

    double sum1 = 0.0;
    double sum2 = 0.0;
    /* Sum over basis function */
    for (int l=0; l<nsplines; l++) {
        if (enforceZeroBoundary) {
            if (l<=1 || l >= nsplines- 2) continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero 
        }
        double Blt = basisfunction(l,t);
        double alpha1 = coeff[skip + carrier_freq_id*nsplines*2 + l];
        double alpha2 = coeff[skip + carrier_freq_id*nsplines*2 + l + nsplines];
        sum1 += alpha1 * Blt;
        sum2 += alpha2 * Blt;
    }
    *Bl1_ptr = sum1;
    *Bl2_ptr = sum2;

}

void BSpline2nd::derivative(const double t, const std::vector<double>& /*coeff*/, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id) {

    /* Iterate over basis function */
    for (int l=0; l<nsplines; l++) {
        if (enforceZeroBoundary){
            if (l<=1 || l >= nsplines- 2) continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero       
        }
        double Blt = basisfunction(l, t); 
        coeff_diff[skip + carrier_freq_id*nsplines*2 + l]            += Blt * valbar1;
        coeff_diff[skip + carrier_freq_id*nsplines*2 + l + nsplines] += Blt * valbar2;
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


BSpline2ndAmplitude::BSpline2ndAmplitude(int nsplines_, double scaling_, double t0, double T, bool enforceZeroBoundary_) : ControlBasis(nsplines_ + 1, t0, T, enforceZeroBoundary_){
    nsplines = nsplines_;
    scaling = scaling_;
    controltype = ControlType::BSPLINEAMP;


    dtknot = (T-t0) / (double)(nsplines - 2);
	width = 3.0*dtknot;

    /* Compute center points of the splines */
    tcenter = new double[nsplines];
    for (int i = 0; i < nsplines; i++){
        tcenter[i] = t0 + dtknot * ( (i+1) - 1.5 );
    }
}

BSpline2ndAmplitude::~BSpline2ndAmplitude(){
    delete [] tcenter;
}

void BSpline2ndAmplitude::enforceBoundary(double* x, int carrier_id){
    // set first and last two splines to zero so that spline starts and ends at zero 
    for (int l=0; l<nsplines; l++) {
        if (l<=1 || l >= nsplines- 2) {
            x[skip + carrier_id*(nsplines+1) + l] = 0.0;
        }
    }
}

void BSpline2ndAmplitude::evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Bl1_ptr, double* Bl2_ptr){

    /* Sum over basis function for amplitudes */
    double ampsum = 0.0;
    for (int l=0; l<nsplines; l++) {
        if (enforceZeroBoundary){
            if (l<=1 || l >= nsplines- 2) continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero       
        }
        double Blt = basisfunction(l,t);
        double alpha1 = coeff[skip + carrier_freq_id*(nsplines+1) + l];
        ampsum += alpha1 * Blt;
    }
    *Bl1_ptr = ampsum;
    // last one is for the phase
    *Bl2_ptr = scaling*coeff[skip + carrier_freq_id*(nsplines+1) + nsplines];  // last one is the phase
}

void BSpline2ndAmplitude::derivative(const double t, const std::vector<double>& coeff, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id) {
    // valbar1 holds the current carrierfrequency. 
    double cos_omt = cos(valbar1*t + scaling * coeff[skip + carrier_freq_id*(nsplines+1) + nsplines]);
    double sin_omt = sin(valbar1*t + scaling * coeff[skip + carrier_freq_id*(nsplines+1) + nsplines]);

    /* Iterate over basis function */
    double ampsum = 0.0;
    for (int l=0; l<nsplines; l++) {
        if (enforceZeroBoundary){
            if (l<=1 || l >= nsplines- 2) continue; // skip first and last two splines (set to zero) so that spline starts and ends at zero       
        }
        double Blt = basisfunction(l, t); 
        double alpha1 = coeff[skip + carrier_freq_id*(nsplines+1) + l];
        ampsum += alpha1 * Blt;
        // Update derivative for the amplitude splines
        // valbar2 holds the flag whether this is p or q
        if (valbar2 > 0.0) coeff_diff[skip + carrier_freq_id*(nsplines+1) + l] += Blt * cos_omt;
        else               coeff_diff[skip + carrier_freq_id*(nsplines+1) + l] += Blt * sin_omt;
    }
    // Update derivate for phase
    if (valbar2 > 0.0) coeff_diff[skip + carrier_freq_id*(nsplines+1) + nsplines] += -ampsum * scaling * sin_omt;
    else               coeff_diff[skip + carrier_freq_id*(nsplines+1) + nsplines] +=  ampsum * scaling * cos_omt;
}

double BSpline2ndAmplitude::basisfunction(int id, double t){

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
BSpline0::BSpline0(int nsplines_, double t0, double T, bool enforceZeroBoundary_) : ControlBasis(2*nsplines_, t0, T, enforceZeroBoundary_){
    nsplines = nsplines_;
    controltype = ControlType::BSPLINE0;

    dtknot = (T-t0) / (nsplines - 1.0);
	width = dtknot;
}

BSpline0::~BSpline0(){}


void BSpline0::evaluate(const double t, const std::vector<double>& coeff, int carrier_freq_id, double* Bl1_ptr, double* Bl2_ptr){

    // Figure out which basis function is active at this time point 
    int splineID = ceil((t-tstart)/dtknot - 0.5);

    // Ctrl function defined to be zero outside [tstart, tend]
    if (splineID < 0 || splineID >= nsplines){
        *Bl1_ptr = 0.0;
        *Bl2_ptr = 0.0;
    } else {
        *Bl1_ptr = coeff[skip + carrier_freq_id*nsplines*2 + splineID];
        *Bl2_ptr = coeff[skip + carrier_freq_id*nsplines*2 + splineID + nsplines];
    }
}

void BSpline0::derivative(const double t, const std::vector<double>& /*coeff*/, double* coeff_diff, const double valbar1, const double valbar2, int carrier_freq_id) {

    // Figure out which basis function is active at this time point 
    int splineID = ceil((t-tstart)/dtknot - 0.5);

    if (splineID >= 0 && splineID < nsplines){
        coeff_diff[skip + carrier_freq_id*nsplines*2 + splineID] += valbar1;
        coeff_diff[skip + carrier_freq_id*nsplines*2 + splineID + nsplines] += valbar2;
    }
}


void BSpline0::enforceBoundary(double* x, int carrierfreqID){

        x[skip + 2*carrierfreqID*nsplines + 0 ] = 0.0; // first real
        x[skip+ 2*carrierfreqID*nsplines + nsplines - 1] = 0.0; // last real
        x[skip + (2*carrierfreqID+1)*nsplines ] = 0.0; // first imag
        x[skip + (2*carrierfreqID+1)*nsplines + nsplines- 1] = 0.0; // last imag
}
