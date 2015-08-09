#include <cmath>
#include <gsl/gsl_randist.h>
#include "quaff.h"

SymQualDist::SymQualDist()
  : symProb(1. / dnaAlphabetSize),
    qualTrialSuccessProb(.5),
    qualNumFailedTrials(FastSeq::qualScoreRange / 2)
{ }

SymQualScores::SymQualScores (const SymQualDist& sqd)
  : symLogProb (log (sqd.symProb)),
    logQualProb (FastSeq::qualScoreRange)
{
  for (int k = 0; k < FastSeq::qualScoreRange; ++k)
    logQualProb[k] = log (gsl_ran_negative_binomial_pdf (k, sqd.qualTrialSuccessProb, sqd.qualNumFailedTrials));
}

SymQualCounts::SymQualCounts()
  : symCount(0.),
    qualCount (FastSeq::qualScoreRange, 0.)
{ }

QuaffParams::QuaffParams()
  : beginInsert(.5),
    extendInsert(.5),
    beginDelete(.5),
    extendDelete(.5),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vector<SymQualDist> (dnaAlphabetSize))
{ }

QuaffScores::QuaffScores (const QuaffParams& qp)
  : pqp(&qp),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vector<SymQualScores> (dnaAlphabetSize))
{
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    insert[i] = SymQualScores (qp.insert[i]);
    for (int j = 0; j < dnaAlphabetSize; ++j)
      match[i][j] = SymQualScores (qp.match[i][j]);
  }

  m2m = log(1-qp.beginInsert) + log(1-qp.beginDelete);
  m2i = log(qp.beginInsert);
  m2d = log(1-qp.beginInsert) + log(qp.beginDelete);
  m2e = log(qp.beginInsert);

  d2d = log(qp.extendDelete);
  d2m = log(1-qp.extendDelete);

  i2i = log(qp.extendInsert);
  i2m = log(1-qp.extendInsert);
}
    
QuaffCounts::QuaffCounts()
  : insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vector<SymQualCounts> (dnaAlphabetSize)),
    m2m(0),
    m2i(0),
    m2d(0),
    m2e(0),
    d2d(0),
    d2m(0),
    i2i(0),
    i2m(0)
{ }

