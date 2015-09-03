#include <iostream>
#include <cmath>
#include "logsumexp.h"

using namespace std;

#define LOG_SUM_EXP_LOOKUP_MAX 10
#define LOG_SUM_EXP_LOOKUP_PRECISION .0001

#define LOG_SUM_EXP_LOOKUP_ENTRIES (((int) (LOG_SUM_EXP_LOOKUP_MAX / LOG_SUM_EXP_LOOKUP_PRECISION)) + 1)

struct LogSumExpLookupTable {
  double *lookup;
  LogSumExpLookupTable();
  ~LogSumExpLookupTable();
};

LogSumExpLookupTable logSumExpLookupTable;

LogSumExpLookupTable::LogSumExpLookupTable() {
  lookup = new double [LOG_SUM_EXP_LOOKUP_ENTRIES];
  int n;
  double x;
  for (n = 0; n < LOG_SUM_EXP_LOOKUP_ENTRIES; ++n) {
    x = n * LOG_SUM_EXP_LOOKUP_PRECISION;
    lookup[n] = log_sum_exp_unary_slow(x);
  }
}

LogSumExpLookupTable::~LogSumExpLookupTable() {
  delete[] lookup;
}

double log_sum_exp (double a, double b) {
  double max, diff, ret;
  // Note: Infinity plus or minus a finite quantity is still Infinity,
  // but Infinity - Infinity = NaN.
  // Thus, we are susceptible to NaN errors when trying to add 0+0 in log-space.
  // To work around this, we explicitly test for a==b.
  if (a == b) { max = a; diff = 0; }
  else if (a < b) { max = b; diff = b - a; }
  else { max = a; diff = a - b; }
  ret = max + log_sum_exp_unary (diff);
#if defined(NAN_DEBUG)
  if (std::isnan(ret)) {
    cerr << "NaN error in log_sum_exp" << endl;
    throw;
  }
#endif
  return ret;
}

double log_sum_exp (double a, double b, double c) {
  return log_sum_exp (log_sum_exp (a, b), c);
}

double log_sum_exp (double a, double b, double c, double d) {
  return log_sum_exp (log_sum_exp (log_sum_exp (a, b), c), d);
}

double log_sum_exp_slow (double a, double b) {
  double min, max, diff, ret;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  ret = max + log_sum_exp_unary_slow (diff);
#if defined(NAN_DEBUG)
  if (std::isnan(ret)) {
    cerr << "NaN error in log_sum_exp" << endl;
    throw;
  }
#endif
  return ret;
}

double log_sum_exp_slow (double a, double b, double c) {
  return log_sum_exp_slow (log_sum_exp_slow (a, b), c);
}

double log_sum_exp_slow (double a, double b, double c, double d) {
  return log_sum_exp_slow (log_sum_exp_slow (log_sum_exp_slow (a, b), c), d);
}

double log_sum_exp_unary (double x) {
#ifdef LOGSUMEXP_DEBUG
  return log_sum_exp_unary_slow(x);
#else /* LOGSUMEXP_DEBUG */
  int n;
  double dx, f0, f1, df;
  if (x >= LOG_SUM_EXP_LOOKUP_MAX || std::isnan(x) || std::isinf(x))
    return 0;
  if (x < 0) {  /* really dumb approximation for x < 0. Should never be encountered, so issue a warning */
    cerr << "Called log_sum_exp_unary(x) for negative x = " << x << endl;
    return -x;
  }
  n = (int) (x / LOG_SUM_EXP_LOOKUP_PRECISION);
  dx = x - (n * LOG_SUM_EXP_LOOKUP_PRECISION);
  f0 = logSumExpLookupTable.lookup[n];
  f1 = logSumExpLookupTable.lookup[n+1];
  df = f1 - f0;
  return f0 + df * (dx / LOG_SUM_EXP_LOOKUP_PRECISION);
#endif /* LOGSUMEXP_DEBUG */
}

double log_sum_exp_unary_slow (double x) {
  return log (1. + exp(-x));
}
