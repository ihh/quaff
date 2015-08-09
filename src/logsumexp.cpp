#include <iostream>
#include <cmath>
#include "logsumexp.h"

using namespace std;

#define LOG_SUM_EXP_LOOKUP_MAX 10
#define LOG_SUM_EXP_LOOKUP_PRECISION .0001

#define LOG_SUM_EXP_LOOKUP_ENTRIES (((int) (LOG_SUM_EXP_LOOKUP_MAX / LOG_SUM_EXP_LOOKUP_PRECISION)) + 1)

struct LogSumExpLookupTable {
  long double *lookup;
  LogSumExpLookupTable();
  ~LogSumExpLookupTable();
};

LogSumExpLookupTable logSumExpLookupTable;

LogSumExpLookupTable::LogSumExpLookupTable() {
  lookup = new long double [LOG_SUM_EXP_LOOKUP_ENTRIES];
  int n;
  long double x;
  for (n = 0; n < LOG_SUM_EXP_LOOKUP_ENTRIES; ++n) {
    x = n * LOG_SUM_EXP_LOOKUP_PRECISION;
    lookup[n] = log_sum_exp_unary_slow(x);
  }
}

LogSumExpLookupTable::~LogSumExpLookupTable() {
  delete lookup;
}

long double log_sum_exp (long double a, long double b) {
  double min, max, diff, ret;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  ret = max + log_sum_exp_unary (diff);
#if defined(NAN_DEBUG)
  if (isnan(ret)) {
    cerr << "NaN error in log_sum_exp" << endl;
    throw;
  }
#endif
  return ret;
}

long double log_sum_exp_slow (long double a, long double b) {
  double min, max, diff, ret;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  ret = max + log_sum_exp_unary_slow (diff);
#if defined(NAN_DEBUG)
  if (isnan(ret)) {
    cerr << "NaN error in log_sum_exp" << endl;
    throw;
  }
#endif
  return ret;
}

long double log_sum_exp_unary (long double x) {
#ifdef LOGSUMEXP_DEBUG
  return log_sum_exp_unary_slow(x);
#else /* LOGSUMEXP_DEBUG */
  int n;
  long double dx, f0, f1, df;
  if (x >= LOG_SUM_EXP_LOOKUP_MAX || isnan(x) || isinf(x))
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

long double log_sum_exp_unary_slow (long double x) {
  return log (1. + exp(-x));
}
