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

/* zeroth, first and second derivatives of logNegativeBinomial(const vguard<double>& kFreq...) as one-dimensional function of nSuccess */
/* also, function to compute optimal pSuccess for given nSuccess */
/* Follows the method of Crowley: http://vixra.org/pdf/1211.0113v1.pdf */

double logNegativeBinomialSingle (double nSuccess, void* kFreqVoid);
double logNegativeBinomialSingleDeriv1 (double nSuccess, void* kFreqVoid);
double logNegativeBinomialSingleDeriv2 (double nSuccess, void* kFreqVoid);
void logNegativeBinomialSingleDeriv12 (double nSuccess, void* kFreqVoid, double* f, double *df);
double optimalNegativeBinomialSuccessProb (double nSuccess, const vguard<double>& kFreq);

/* function bodies */
double logNegativeBinomial (int k, double pSuccess, double nSuccess) {
  return log (gsl_ran_negative_binomial_pdf (k, pSuccess, nSuccess));
}

double logNegativeBinomial (const vguard<double>& kFreq, double pSuccess, double nSuccess) {
  double lp = 0;
  for (int k = 0; k < (int) kFreq.size(); ++k)
    lp += kFreq[k] * logNegativeBinomial (k, pSuccess, nSuccess);
  return lp;
}

double logNegativeBinomialSingle (double nSuccess, void* kFreqVoid) {
  const vguard<double>& kFreq = *(const vguard<double>*) kFreqVoid;
  const double pSuccess = optimalNegativeBinomialSuccessProb (nSuccess, kFreq);
  return logNegativeBinomial (kFreq, pSuccess, nSuccess);
}

double logNegativeBinomialSingleDeriv1 (double nSuccess, void* kFreqVoid) {
  const vguard<double>& kFreq = *(const vguard<double>*) kFreqVoid;
  double freqSum = 0, kSum = 0, kDigammaSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    if (freq > 0) {  /* avoid unnecessary digamma evaluation */
      freqSum += freq;
      kSum += freq * k;
      kDigammaSum += freq * gsl_sf_psi (nSuccess + k);
    }
  }
  return -freqSum * log(1. + kSum / (freqSum * nSuccess)) - freqSum * gsl_sf_psi (nSuccess) + kDigammaSum;
}

double logNegativeBinomialSingleDeriv2 (double nSuccess, void* kFreqVoid) {
  const vguard<double>& kFreq = *(const vguard<double>*) kFreqVoid;
  double freqSum = 0, kTrigammaSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    if (freq > 0) {  /* avoid unnecessary trigamma evaluation */
      freqSum += freq;
      kTrigammaSum += freq * gsl_sf_psi_1 (nSuccess + k);
    }
  }
  return -freqSum * gsl_sf_psi_1 (nSuccess) + kTrigammaSum;
}

void logNegativeBinomialSingleDeriv12 (double nSuccess, void* kFreqVoid, double* f, double *df) {
  *f = logNegativeBinomialSingleDeriv1 (nSuccess, kFreqVoid);
  *df = logNegativeBinomialSingleDeriv2 (nSuccess, kFreqVoid);
}

double optimalNegativeBinomialSuccessProb (double nSuccess, const vguard<double>& kFreq) {
  double freqSum = 0, kSum = 0;
  for (size_t k = 0; k < kFreq.size(); ++k) {
    const double freq = kFreq[k];
    freqSum += freq;
    kSum += freq * k;
  }
  return 1. / (1 + kSum / (freqSum * nSuccess));
}

void calcIntDistribMoments (const vguard<double>& kFreq, double& count, double& mean, double& variance) {
  count = 0;
  double kSum = 0, kSqSum = 0;
  for (unsigned int k = 0; k < kFreq.size(); ++k) {
    count += kFreq[k];
    kSum += kFreq[k] * k;
    kSqSum += kFreq[k] * k * k;
  }
  if (count > 0) {
    mean = kSum / count;
    variance = kSqSum / count - mean*mean;
  } else
    mean = variance = nan("");
}

int fitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess) {
  int status;
  double count, mean, variance;
  calcIntDistribMoments (kFreq, count, mean, variance);
  if (count <= 0) {
    pSuccess = nSuccess = nan("");
    return GSL_EZERODIV;
  }
  if (variance > 0) {
    status = momentFitNegativeBinomial (mean, variance, pSuccess, nSuccess);
    if (status == GSL_SUCCESS)
      status = bracketFitNegativeBinomial (kFreq, pSuccess, nSuccess, max(1.,nSuccess/2), min(kFreq.size()-1.,nSuccess*2));
  } else
    status = bracketFitNegativeBinomial (kFreq, pSuccess, nSuccess);
  if (status == GSL_SUCCESS)
    status = gradientFitNegativeBinomial (kFreq, pSuccess, nSuccess);
  return status;
}

int momentFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess) {
  double count, mean, variance;
  calcIntDistribMoments (kFreq, count, mean, variance);
  if (count <= 0) {
    pSuccess = nSuccess = nan("");
    return GSL_EZERODIV;
  }
  return momentFitNegativeBinomial (mean, variance, pSuccess, nSuccess);
}

int momentFitNegativeBinomial (double mean, double variance, double& pSuccess, double& nSuccess) {
  if (variance <= 0) {
    /* guess sensible values */
    pSuccess = 1 / (1 + mean);
    nSuccess = 1;

    if (LogThisAt(5))
      cerr << "Method-of-moments fit failed: zero variance" << endl;
    GSL_ERROR ("Zero variance in method-of-moments fit", NEG_BINOM_ZERO_VARIANCE);
  }

  pSuccess = mean / variance;
  nSuccess = mean * pSuccess / (1 - pSuccess);

  if (LogThisAt(5))
    cerr << "Method-of-moments fit: pSuccess=" << pSuccess << ", nSuccess=" << nSuccess << endl;

  return GSL_SUCCESS;
}

int bracketFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess) {
  return bracketFitNegativeBinomial (kFreq, pSuccess, nSuccess, 1., max (1., kFreq.size() - 1.));
}

int bracketFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess, double nSuccessLowerBound, double nSuccessUpperBound) {
  int status;
  const gsl_root_fsolver_type *bracketType;
  gsl_root_fsolver *bracketSolver;
  double nSuccessGuess = 0;
  gsl_function F;

  /* guess sensible values in case of failure */
  nSuccess = 1;
  pSuccess = optimalNegativeBinomialSuccessProb (nSuccess, kFreq);

  // bracket
  F.function = &logNegativeBinomialSingleDeriv1;
  F.params = (void*) &kFreq;

  bracketType = gsl_root_fsolver_brent;
  bracketSolver = gsl_root_fsolver_alloc (bracketType);

  const double derivAtLowerBound = GSL_FN_EVAL(&F,nSuccessLowerBound);
  const double derivAtUpperBound = GSL_FN_EVAL(&F,nSuccessUpperBound);
  if (sgn(derivAtLowerBound) == sgn(derivAtUpperBound)) {

    const double loglikeAtLowerBound = logNegativeBinomialSingle (nSuccessLowerBound, (void*) &kFreq);
    const double loglikeAtUpperBound = logNegativeBinomialSingle (nSuccessUpperBound, (void*) &kFreq);

    nSuccess = loglikeAtLowerBound > loglikeAtUpperBound ? nSuccessLowerBound : nSuccessUpperBound;
    pSuccess = optimalNegativeBinomialSuccessProb (nSuccess, kFreq);

    if (LogThisAt(5))
      cerr << "Bracket fit failed; derivative has same sign at both endpoints. Choosing " << (loglikeAtLowerBound > loglikeAtUpperBound ? "lower" : "upper") << " endpoint" << endl
	   << "Bracket fit (fallback): pSuccess=" << pSuccess << ", nSuccess=" << nSuccess << endl;

    status = GSL_SUCCESS;

  } else {
    status = gsl_root_fsolver_set (bracketSolver, &F, nSuccessLowerBound, nSuccessUpperBound);

    if (status == GSL_SUCCESS) {
      if (LogThisAt(5)) {
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
	    nSuccessGuess = gsl_root_fsolver_root (bracketSolver);
	    nSuccessLowerBound = gsl_root_fsolver_x_lower (bracketSolver);
	    nSuccessUpperBound = gsl_root_fsolver_x_upper (bracketSolver);

	    status = gsl_root_test_interval (nSuccessLowerBound, nSuccessUpperBound,
					     MaxNegBinomBracketAbsoluteError, MaxNegBinomBracketRelativeError);

	    if (LogThisAt(5)) {
	      if (status == GSL_SUCCESS)
		fprintf (stderr, "Converged:\n");

	      fprintf (stderr, "%5d [%.7f, %.7f] %.7f %+.7f %+.7f %.7f\n",
		       iter, nSuccessLowerBound, nSuccessUpperBound,
		       nSuccessGuess,
		       logNegativeBinomialSingle (nSuccessGuess, (void*) &kFreq),
		       logNegativeBinomialSingleDeriv1 (nSuccessGuess, (void*) &kFreq),
		       nSuccessUpperBound - nSuccessLowerBound);
	    }
	  }
	}
      while (status == GSL_CONTINUE && iter < MaxNegBinomBracketMaxIterations);

      gsl_root_fsolver_free (bracketSolver);

      nSuccess = nSuccessGuess;
      pSuccess = optimalNegativeBinomialSuccessProb (nSuccess, kFreq);

      if (LogThisAt(5))
	cerr << "Bracket fit: pSuccess=" << pSuccess << ", nSuccess=" << nSuccess << endl;
    }
  }

  return status;
}

int gradientFitNegativeBinomial (const vguard<double>& kFreq, double& pSuccess, double& nSuccess) {
  int status;
  const gsl_root_fdfsolver_type *derivType;
  gsl_root_fdfsolver *derivSolver;
  double nSuccessGuess = nSuccess;
  gsl_function_fdf FDF;

  // solve
  FDF.f = &logNegativeBinomialSingleDeriv1;
  FDF.df = &logNegativeBinomialSingleDeriv2;
  FDF.fdf = &logNegativeBinomialSingleDeriv12;
  FDF.params = (void*) &kFreq;

  derivType = gsl_root_fdfsolver_newton;
  derivSolver = gsl_root_fdfsolver_alloc (derivType);
  status = gsl_root_fdfsolver_set (derivSolver, &FDF, nSuccessGuess);

  if (status == GSL_SUCCESS) {
    if (LogThisAt(5)) {
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
	  const double nSuccessLast = nSuccessGuess;
	  nSuccessGuess = gsl_root_fdfsolver_root (derivSolver);
	  status = gsl_root_test_delta (nSuccessGuess, nSuccessLast, MaxNegBinomPolishAbsoluteError, MaxNegBinomPolishRelativeError);

	  if (status == GSL_CONTINUE && nSuccessGuess > kFreq.size())
	    status = GSL_ERUNAWAY;

	  if (LogThisAt(5)) {
	    if (status == GSL_SUCCESS)
	      fprintf (stderr, "Converged:\n");
	    else if (status == GSL_ERUNAWAY)
	      fprintf (stderr, "Runaway iteration:\n");

	    fprintf (stderr, "%5d %10.7f %+10.7f %10.7f %10.7f\n",
		     iter,
		     nSuccessGuess,
		     logNegativeBinomialSingle (nSuccessGuess, (void*) &kFreq),
		     logNegativeBinomialSingleDeriv1 (nSuccessGuess, (void*) &kFreq),
		     nSuccessGuess - nSuccessLast);
	  }
	}
      }
    while (status == GSL_CONTINUE && iter < MaxNegBinomPolishMaxIterations);

    gsl_root_fdfsolver_free (derivSolver);

    nSuccess = nSuccessGuess;
    pSuccess = optimalNegativeBinomialSuccessProb (nSuccess, kFreq);

    if (LogThisAt(5))
      cerr << "Gradient fit: pSuccess=" << pSuccess << ", nSuccess=" << nSuccess << endl;
  }

  return status;
}
