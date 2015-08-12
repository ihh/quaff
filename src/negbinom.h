#ifndef NEGBINOM_INCLUDED
#define NEGBINOM_INCLUDED

#include <vector>
#include "logger.h"

using namespace std;

#define NEG_BINOM_ZERO_VARIANCE GSL_EZERODIV  /* error code for momentFitNegativeBinomial */

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* evaluate the log-pdf of the negative binomial distribution for a point, and for a frequency distribution */
double logNegativeBinomial (int k, double pSuccess, double nSuccess);
double logNegativeBinomial (const vector<double>& kFreq, double pSuccess, double nSuccess);

/* calculate moments of an integer frequency distribution */
void calcIntDistribMoments (const vector<double>& kFreq, double& count, double& mean, double& variance);

/* fit methods return GSL status/error codes */
/* fitNegativeBinomial wraps other methods */
int fitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nSuccess);

/* momentFitNegativeBinomial uses method of moments */
int momentFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nSuccess);  /* method of moments */
int momentFitNegativeBinomial (double mean, double variance, double& pSuccess, double& nSuccess);

/* bracketFitNegativeBinomial uses first derivative */
int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nSuccess);
int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nSuccess, double nSuccessLowerBound, double nSuccessUpperBound);

/* gradientFitNegativeBinomial uses first & second derivatives, requires an initial estimate */
int gradientFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nSuccess);

#endif /* NEGBINOM_INCLUDED */
