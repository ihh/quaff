#ifndef LOGSUMEXP_INCLUDED
#define LOGSUMEXP_INCLUDED

/* uncomment to disable lookup table */
/*
#define LOGSUMEXP_DEBUG
*/

/* uncomment to catch NaN errors */
/*
#define NAN_DEBUG
*/

long double log_sum_exp (long double a, long double b);  /* returns log(exp(a) + exp(b)) */
long double log_sum_exp_unary (long double x);  /* returns log(1 + exp(-x)) for nonnegative x */

long double log_sum_exp_slow (long double a, long double b);  /* does not use lookup table */
long double log_sum_exp_unary_slow (long double x);  /* does not use lookup table */

#endif /* LOGSUMEXP_INCLUDED */
