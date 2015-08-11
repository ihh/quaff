#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include "../src/negbinom.h"

int main (int argc, char **argv) {
  if (argc != 5) {
    cout << "Usage: " << argv[0] << " pSuccess nFail nSamples relativeError\n";
    exit (EXIT_FAILURE);
  }

  const double p = atof (argv[1]);
  const double r = atof (argv[2]);
  const int N = atoi (argv[3]);
  const double eps = atof (argv[4]);

  const gsl_rng_type * rngType;
  gsl_rng * rng;

  gsl_rng_env_setup();
  rngType = gsl_rng_default;
  rng = gsl_rng_alloc (rngType);

  vector<double> kFreq;
  for (int n = 0; n < N; ++n) {
    const unsigned int k = gsl_ran_negative_binomial (rng, p, r);
    while (kFreq.size() <= k)
      kFreq.push_back (0.);
    ++kFreq[k];
  }
  
  gsl_rng_free (rng);

  cerr << "Frequencies:";
  for (auto n : kFreq)
    cerr << ' ' << n;
  cerr << endl;

  double kSum = 0, kSqSum = 0;
  for (unsigned int k = 0; k < kFreq.size(); ++k) {
    kSum += kFreq[k] * k;
    kSqSum += kFreq[k] * k * k;
  }
  const double kSampleMean = kSum / N, kSampleVariance = kSqSum / N - kSampleMean*kSampleMean;
  const double kMean = r*(1-p)/p, kVariance = r*(1-p)/(p*p);
  cerr << "Distribution: mean = " << kMean << ", variance = " << kVariance << endl;
  cerr << "      Sample: mean = " << kSampleMean << ", variance = " << kSampleVariance << endl;
  
  double pFit, rFit;
  logger.verbosity = 3;
  const int status = fitNegativeBinomial (kFreq, pFit, rFit);

  const double kFitMean = rFit*(1-pFit)/pFit, kFitVariance = rFit*(1-pFit)/(pFit*pFit);
  cerr << "         Fit: mean = " << kFitMean << ", variance = " << kFitVariance << endl;

  
  if (gsl_root_test_delta (p, pFit, 0, eps) == GSL_SUCCESS
      && gsl_root_test_delta (r, rFit, 0, eps) == GSL_SUCCESS)
    cout << "ok\n";
  else
    cout << "not ok\n";
  
  exit (EXIT_SUCCESS);
}
