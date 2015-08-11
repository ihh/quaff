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
    lp += kFreq[k] * logNegativeBinomial (k, pSuccess, nFail);
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
  return -freqSum * log(1. + kSum / (freqSum * nFail)) - freqSum * gsl_sf_psi (nFail) + kDigammaSum;
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
  return 1. / (1 + kSum / (freqSum * nFail));
}

void calcIntDistribMeanVariance (const vector<double>& kFreq, double& mean, double& variance) {
  double freqSum = 0, kSum = 0, kSqSum = 0;
  for (unsigned int k = 0; k < kFreq.size(); ++k) {
    freqSum += kFreq[k];
    kSum += kFreq[k] * k;
    kSqSum += kFreq[k] * k * k;
  }
  mean = kSum / freqSum;
  variance = kSqSum / freqSum - mean*mean;
}

int fitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail) {
  int status;
  double mean, variance;
  calcIntDistribMeanVariance (kFreq, mean, variance);
  if (variance > 0) {
    status = momentFitNegativeBinomial (mean, variance, pSuccess, nFail);
    if (status == GSL_SUCCESS)
      status = bracketFitNegativeBinomial (kFreq, pSuccess, nFail, max(1.,nFail/2), min(kFreq.size()-1.,nFail*2));
  }
  else
    status = bracketFitNegativeBinomial (kFreq, pSuccess, nFail);
  if (status == GSL_SUCCESS)
    status = gradientFitNegativeBinomial (kFreq, pSuccess, nFail);
  return status;
}

int momentFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail) {
  double kSampleMean, kSampleVariance;
  calcIntDistribMeanVariance (kFreq, kSampleMean, kSampleVariance);
  return momentFitNegativeBinomial (kSampleMean, kSampleVariance, pSuccess, nFail);
}

int momentFitNegativeBinomial (double mean, double variance, double& pSuccess, double& nFail) {
  if (variance <= 0) {
    /* guess sensible values */
    pSuccess = 1. / (mean + 1);
    nFail = 1;

    if (LogThisAt(3))
      cerr << "Method-of-moments fit failed: zero variance" << endl;
    GSL_ERROR ("Zero variance in method-of-moments fit", NEG_BINOM_ZERO_VARIANCE);
  }

  pSuccess = mean / variance;
  nFail = mean * pSuccess / (1 - pSuccess);

  if (LogThisAt(3))
    cerr << "Method-of-moments fit: pSuccess=" << pSuccess << ", nFail=" << nFail << endl;

  return GSL_SUCCESS;
}

int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail) {
  return bracketFitNegativeBinomial (kFreq, pSuccess, nFail, 1., max (1., kFreq.size() - 1.));
}

int bracketFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail, double nFailLowerBound, double nFailUpperBound) {
  int status;
  const gsl_root_fsolver_type *bracketType;
  gsl_root_fsolver *bracketSolver;
  double nFailGuess;
  gsl_function F;

  /* guess sensible values in case of failure */
  nFail = 1;
  pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);

  // bracket
  F.function = &logNegativeBinomialSingleDeriv1;
  F.params = (void*) &kFreq;

  bracketType = gsl_root_fsolver_brent;
  bracketSolver = gsl_root_fsolver_alloc (bracketType);

  const double derivAtLowerBound = GSL_FN_EVAL(&F,nFailLowerBound);
  const double derivAtUpperBound = GSL_FN_EVAL(&F,nFailUpperBound);
  if (sgn(derivAtLowerBound) == sgn(derivAtUpperBound)) {

    const double loglikeAtLowerBound = logNegativeBinomialSingle (nFailLowerBound, (void*) &kFreq);
    const double loglikeAtUpperBound = logNegativeBinomialSingle (nFailUpperBound, (void*) &kFreq);

    nFail = loglikeAtLowerBound > loglikeAtUpperBound ? nFailLowerBound : nFailUpperBound;
    pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);

    if (LogThisAt(3))
      cerr << "Bracket fit failed; derivative has same sign at both endpoints. Choosing " << (loglikeAtLowerBound > loglikeAtUpperBound ? "lower" : "upper") << " endpoint" << endl
	   << "Bracket fit (fallback): pSuccess=" << pSuccess << ", nFail=" << nFail << endl;

    status = GSL_SUCCESS;

  } else {
    status = gsl_root_fsolver_set (bracketSolver, &F, nFailLowerBound, nFailUpperBound);

    if (status == GSL_SUCCESS) {
      if (LogThisAt(3)) {
	fprintf (stderr, "Bracketing max of negative binomial using %s method\n", 
		 gsl_root_fsolver_name (bracketSolver));

	fprintf (stderr, "%5s [%9s, %9s] %9s %10s %10s %9s\n",
		 "iter", "lower", "upper", "r",
		 "f", "df/dr", "err(est)");
      }

      int iter = 0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (bracketSolver);
	  if (status == GSL_SUCCESS) {
	    nFailGuess = gsl_root_fsolver_root (bracketSolver);
	    nFailLowerBound = gsl_root_fsolver_x_lower (bracketSolver);
	    nFailUpperBound = gsl_root_fsolver_x_upper (bracketSolver);

	    status = gsl_root_test_interval (nFailLowerBound, nFailUpperBound,
					     MaxNegBinomBracketAbsoluteError, MaxNegBinomBracketRelativeError);

	    if (LogThisAt(3)) {
	      if (status == GSL_SUCCESS)
		fprintf (stderr, "Converged:\n");

	      fprintf (stderr, "%5d [%.7f, %.7f] %.7f %+.7f %+.7f %.7f\n",
		       iter, nFailLowerBound, nFailUpperBound,
		       nFailGuess,
		       logNegativeBinomialSingle (nFailGuess, (void*) &kFreq),
		       logNegativeBinomialSingleDeriv1 (nFailGuess, (void*) &kFreq),
		       nFailUpperBound - nFailLowerBound);
	    }
	  }
	}
      while (status == GSL_CONTINUE && iter < MaxNegBinomBracketMaxIterations);

      gsl_root_fsolver_free (bracketSolver);

      nFail = nFailGuess;
      pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);

      if (LogThisAt(3))
	cerr << "Bracket fit: pSuccess=" << pSuccess << ", nFail=" << nFail << endl;
    }
  }

  return status;
}

int gradientFitNegativeBinomial (const vector<double>& kFreq, double& pSuccess, double& nFail) {
  int status;
  const gsl_root_fdfsolver_type *derivType;
  gsl_root_fdfsolver *derivSolver;
  double nFailGuess = nFail;
  gsl_function F;
  gsl_function_fdf FDF;

  // solve
  FDF.f = &logNegativeBinomialSingleDeriv1;
  FDF.df = &logNegativeBinomialSingleDeriv2;
  FDF.fdf = &logNegativeBinomialSingleDeriv12;
  FDF.params = (void*) &kFreq;

  derivType = gsl_root_fdfsolver_newton;
  derivSolver = gsl_root_fdfsolver_alloc (derivType);
  status = gsl_root_fdfsolver_set (derivSolver, &FDF, nFailGuess);

  if (status == GSL_SUCCESS) {
    if (LogThisAt(3)) {
      fprintf (stderr, "Polishing max of negative binomial using %s method\n", 
	       gsl_root_fdfsolver_name (derivSolver));

      fprintf (stderr, "%-5s %10s %10s %10s %10s\n",
	       "iter", "r", "f", "df/dr", "err(est)");
    }
  
    int iter = 0;
    do
      {
	iter++;
	status = gsl_root_fdfsolver_iterate (derivSolver);
	if (status == GSL_SUCCESS) {
	  const double nFailLast = nFailGuess;
	  nFailGuess = gsl_root_fdfsolver_root (derivSolver);
	  status = gsl_root_test_delta (nFailGuess, nFailLast, MaxNegBinomPolishAbsoluteError, MaxNegBinomPolishRelativeError);

	  if (status == GSL_CONTINUE && nFailGuess > kFreq.size())
	    status = GSL_ERUNAWAY;

	  if (LogThisAt(3)) {
	    if (status == GSL_SUCCESS)
	      fprintf (stderr, "Converged:\n");
	    else if (status == GSL_ERUNAWAY)
	      fprintf (stderr, "Runaway iteration:\n");

	    fprintf (stderr, "%5d %10.7f %+10.7f %10.7f %10.7f\n",
		     iter,
		     nFailGuess,
		     logNegativeBinomialSingle (nFailGuess, (void*) &kFreq),
		     logNegativeBinomialSingleDeriv1 (nFailGuess, (void*) &kFreq),
		     nFailGuess - nFailLast);
	  }
	}
      }
    while (status == GSL_CONTINUE && iter < MaxNegBinomPolishMaxIterations);

    gsl_root_fdfsolver_free (derivSolver);

    nFail = nFailGuess;
    pSuccess = optimalNegativeBinomialSuccessProb (nFail, kFreq);

    if (LogThisAt(3))
      cerr << "Gradient fit: pSuccess=" << pSuccess << ", nFail=" << nFail << endl;
  }

  return status;
}
