#ifndef NEGBINOM_INCLUDED
#define NEGBINOM_INCLUDED

#include <vector>
#include "logger.h"

using namespace std;

#define NEG_BINOM_ZERO_VARIANCE GSL_EZERODIV  /* error code for momentFitNegativeBinomial */

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double logNegativeBinomial (int k, double pSuccess, double nFail);
double logNegativeBinomial (const vector<double>& kFreq, double pSuccess, double nFail);

void calcIntDistribMeanVariance (const vector<double>& kFreq, double& mean, double& variance);

/* fit methods return GSL status/error codes */
/* fitNegativeBinomial wraps other methods */
int fitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail);

/* momentFitNegativeBinomial uses method of moments */
int momentFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail);  /* method of moments */
int momentFitNegativeBinomial (double mean, double variance, double& pSuccess, double& nFail);

/* bracketFitNegativeBinomial uses first derivative */
int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail);
int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail, double nFailLowerBound, double nFailUpperBound);

/* gradientFitNegativeBinomial uses first & second derivatives, requires an initial estimate */
int gradientFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail);

#endif /* NEGBINOM_INCLUDED */
