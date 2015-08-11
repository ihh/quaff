#ifndef NEGBINOM_INCLUDED
#define NEGBINOM_INCLUDED

#include <vector>
#include "logger.h"

using namespace std;

double logNegativeBinomial (int k, double pSuccess, double nFail);
double logNegativeBinomial (const vector<double>& kFreq, double pSuccess, double nFail);

/* fit method returns a GSL error code */
int fitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail);

#endif /* NEGBINOM_INCLUDED */
