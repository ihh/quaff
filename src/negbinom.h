#ifndef NEGBINOM_INCLUDED
#define NEGBINOM_INCLUDED

#include "vguard.h"
#include "util.h"

using namespace std;

/* error codes for momentFitNegativeBinomial */
#define NEG_BINOM_ZERO_VARIANCE GSL_EZERODIV
#define NEG_BINOM_LOW_VARIANCE  GSL_ERANGE

/* evaluate the log-pdf of the negative binomial distribution for a point, and for a frequency distribution */
double logNegativeBinomial (int k, double pSuccess, double nSuccess);
double logNegativeBinomial (const vguard<double>& kFreq, double pSuccess, double nSuccess);

/* calculate moments of negative binomial distribution */
double negativeBinomialMean (double pSuccess, double nSuccess);
double negativeBinomialVariance (double pSuccess, double nSuccess);

/* calculate moments of an integer frequency distribution */
void calcIntDistribMoments (const vguard<double>& kFreq, double& count, double& mean, double& variance);

/* fit methods return GSL status/error codes */
/* fitNegativeBinomial wraps other methods */
int fitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess);

/* momentFitNegativeBinomial uses method of moments */
int momentFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess);  /* method of moments */
int momentFitNegativeBinomial (double mean, double variance, double& pSuccess, double& nSuccess);

/* bracketFitNegativeBinomial uses first derivative */
int bracketFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess);
int bracketFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess, double nSuccessLowerBound, double nSuccessUpperBound);

/* gradientFitNegativeBinomial uses first & second derivatives, requires an initial estimate */
int gradientFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess);

#endif /* NEGBINOM_INCLUDED */
