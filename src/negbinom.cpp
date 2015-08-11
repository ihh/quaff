#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include "negbinom.h"
#include "util.h"

/* convergence parameters */
#define MaxNegBinomBracketMaxIterations 100
#define MaxNegBinomBracketAbsoluteError 1e-3
#define MaxNegBinomBracketRelativeError 1e-3
#define MaxNegBinomPolishMaxIterations 100
#define MaxNegBinomPolishAbsoluteError 0
#define MaxNegBinomPolishRelativeError 1e-4

/* zeroth, first and second derivatives of logNegativeBinomial(const vector<double>& kFreq...) as one-dimensional function of nFail */
/* also, function to compute optimal pSuccess for given nFail */
/* Follows the method of Crowley: http://vixra.org/pdf/1211.0113v1.pdf */

double logNegativeBinomialSingle (double nFail, void* kFreqVoid);
double logNegativeBinomialSingleDeriv1 (double nFail, void* kFreqVoid);
double logNegativeBinomialSingleDeriv2 (double nFail, void* kFreqVoid);
void logNegativeBinomialSingleDeriv12 (double nFail, void* kFreqVoid, double* f, double *df);
double optimalNegativeBinomialSuccessProb (double nFail, const vector<double>& kFreq);

/* function bodies */
double logNegativeBinomial (int k, double pSuccess, double nFail) {
  return log (gsl_ran_negative_binomial_pdf (k, pSuccess, nFail));
}

double logNegativeBinomial (const vector<double>& kFreq, double pSuccess, double nFail) {
  double lp = 0;
  for (int k = 0; k < (int) kFreq.size(); ++k)
    lp += kFreq[k] * logNegativeBinomial (pSuccess, nFail, k);
  return lp;
}

double logNegativeBinomialSingle (double nFail, void* kFreqVoid) {
  const vector<double>& kFreq = *(const vector<double>*) kFreqVoid;
  const double pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);
  return logNegativeBinomial (kFreq, pSuccess, nFail);
}

double logNegativeBinomialSingleDeriv1 (double nFail, void* kFreqVoid) {
  const vector<double>& kFreq = *(const vector<double>*) kFreqVoid;
  double freqSum = 0, kSum = 0, kDigammaSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    if (freq > 0) {  /* avoid unnecessary digamma evaluation */
      freqSum += freq;
      kSum += freq * k;
      kDigammaSum += freq * gsl_sf_psi (nFail + k);
    }
  }
  return -freqSum * log(1. + kSum / freqSum) - freqSum * gsl_sf_psi (nFail) + kDigammaSum;
}

double logNegativeBinomialSingleDeriv2 (double nFail, void* kFreqVoid) {
  const vector<double>& kFreq = *(const vector<double>*) kFreqVoid;
  double freqSum = 0, kTrigammaSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    if (freq > 0) {  /* avoid unnecessary trigamma evaluation */
      freqSum += freq;
      kTrigammaSum += freq * gsl_sf_psi_1 (nFail + k);
    }
  }
  return -freqSum * gsl_sf_psi_1 (nFail) + kTrigammaSum;
}

void logNegativeBinomialSingleDeriv12 (double nFail, void* kFreqVoid, double* f, double *df) {
  *f = logNegativeBinomialSingleDeriv1 (nFail, kFreqVoid);
  *df = logNegativeBinomialSingleDeriv2 (nFail, kFreqVoid);
}

double optimalNegativeBinomialSuccessProb (double nFail, const vector<double>& kFreq) {
  double freqSum = 0, kSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    freqSum += freq;
    kSum += freq * k;
  }
  return -log (1 + kSum / freqSum);
}

int fitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail) {
  int status, iter;
  const gsl_root_fsolver_type *bracketType;
  gsl_root_fsolver *bracketSolver;
  const gsl_root_fdfsolver_type *derivType;
  gsl_root_fdfsolver *derivSolver;
  double nFailLowerBound = 1., nFailUpperBound = kFreq.size() - 1., nFailGuess;
  gsl_function F;
  gsl_function_fdf FDF;

  // bracket
  F.function = &logNegativeBinomialSingleDeriv1;
  F.params = (void*) &kFreq;

  bracketType = gsl_root_fsolver_brent;
  bracketSolver = gsl_root_fsolver_alloc (bracketType);
  gsl_root_fsolver_set (bracketSolver, &F, nFailLowerBound, nFailUpperBound);

  if (LogThisAt(3)) {
    fprintf (stderr, "Bracketing max of negative binomial using %s method\n", 
	     gsl_root_fsolver_name (bracketSolver));

    fprintf (stderr, "%5s [%9s, %9s] %9s %10s %10s %9s\n",
	     "iter", "lower", "upper", "max",
	     "f", "df", "err(est)");
  }

  iter = 0;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (bracketSolver);
      nFailGuess = gsl_root_fsolver_root (bracketSolver);
      nFailLowerBound = gsl_root_fsolver_x_lower (bracketSolver);
      nFailUpperBound = gsl_root_fsolver_x_upper (bracketSolver);
      status = gsl_root_test_interval (nFailLowerBound, nFailUpperBound,
                                       MaxNegBinomBracketAbsoluteError, MaxNegBinomBracketRelativeError);

      if (LogThisAt(3)) {
	if (status == GSL_SUCCESS)
	  fprintf (stderr, "Converged:\n");

	fprintf (stderr, "%5d [%.7f, %.7f] %.7f %+.7f %.7f %.7f\n",
		 iter, nFailLowerBound, nFailUpperBound,
		 nFailGuess,
		 logNegativeBinomialSingle (nFailGuess, (void*) &kFreq),
		 logNegativeBinomialSingleDeriv1 (nFailGuess, (void*) &kFreq),
		 nFailUpperBound - nFailLowerBound);
      }
    }
  while (status == GSL_CONTINUE && iter < MaxNegBinomBracketMaxIterations);

  gsl_root_fsolver_free (bracketSolver);

  // solve
  FDF.f = &logNegativeBinomialSingleDeriv1;
  FDF.df = &logNegativeBinomialSingleDeriv2;
  FDF.fdf = &logNegativeBinomialSingleDeriv12;
  FDF.params = (void*) &kFreq;

  derivType = gsl_root_fdfsolver_newton;
  derivSolver = gsl_root_fdfsolver_alloc (derivType);
  gsl_root_fdfsolver_set (derivSolver, &FDF, nFailGuess);

  if (LogThisAt(3)) {
    fprintf (stderr, "Polishing max of negative binomial using %s method\n", 
	     gsl_root_fdfsolver_name (derivSolver));

    fprintf (stderr, "%-5s %10s %10s %10s %10s\n",
	     "iter", "max", "f", "df", "err(est)");
  }
  
  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (derivSolver);
      const double nFailLast = nFailGuess;
      nFailGuess = gsl_root_fdfsolver_root (derivSolver);
      status = gsl_root_test_delta (nFailGuess, nFailLast, MaxNegBinomPolishAbsoluteError, MaxNegBinomPolishRelativeError);

      if (LogThisAt(3)) {
	if (status == GSL_SUCCESS)
	  fprintf (stderr, "Converged:\n");

	fprintf (stderr, "%5d %10.7f %+10.7f %10.7f %10.7f\n",
		 iter,
		 nFailGuess,
		 logNegativeBinomialSingle (nFailGuess, (void*) &kFreq),
		 logNegativeBinomialSingleDeriv1 (nFailGuess, (void*) &kFreq),
		 nFailGuess - nFailLast);
      }
    }
  while (status == GSL_CONTINUE && iter < MaxNegBinomPolishMaxIterations);

  gsl_root_fdfsolver_free (derivSolver);

  nFail = nFailGuess;
  pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);

  if (LogThisAt(3))
    cerr << "Fit: pSuccess=" << pSuccess << ", nFail=" << nFail << endl;

  return status;
}
