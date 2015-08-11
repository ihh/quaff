#ifndef NEGBINOM_INCLUDED
#define NEGBINOM_INCLUDED

#include <vector>

using namespace std;

void fitNegativeBinomial (double& pSuccess, double& nFail, const vector<double>& kFreq);

#endif /* NEGBINOM_INCLUDED */
