#include <iostream>
#include <regex>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include "qmodel.h"
#include "logsumexp.h"

SymQualDist::SymQualDist()
  : symProb(1. / dnaAlphabetSize),
    qualTrialSuccessProb(.5),
    qualNumFailedTrials(FastSeq::qualScoreRange / 2)
{ }

void SymQualDist::write (ostream& out, const string& prefix) const {
  out << prefix << ": " << symProb << endl;
  out << prefix << "qp: " << qualTrialSuccessProb << endl;
  out << prefix << "qr: " << qualNumFailedTrials << endl;
}

void SymQualDist::read (map<string,double>& paramVal, const string& prefix) {
  symProb = paramVal[prefix];
  qualTrialSuccessProb = paramVal[prefix + "qp"];
  qualNumFailedTrials = paramVal[prefix + "qr"];
}

SymQualScores::SymQualScores (const SymQualDist& sqd)
  : logQualProb (FastSeq::qualScoreRange)
{
  const double symLogProb = log (sqd.symProb);
  for (int k = 0; k < FastSeq::qualScoreRange; ++k)
    logQualProb[k] = symLogProb + log (gsl_ran_negative_binomial_pdf (k, sqd.qualTrialSuccessProb, sqd.qualNumFailedTrials));
}

SymQualCounts::SymQualCounts()
  : qualCount (FastSeq::qualScoreRange, 0.)
{ }

void SymQualCounts::write (ostream& out, const string& prefix) const {
  if (qualCount.size()) {
    out << prefix << ": [" << qualCount[0];
    for (unsigned int i = 1; i < qualCount.size(); ++i)
      out << ", " << qualCount[i];
    out << ']' << endl;
  }
}

QuaffParams::QuaffParams()
  : logger (NULL),
    beginInsert(.5),
    extendInsert(.5),
    beginDelete(.5),
    extendDelete(.5),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vector<SymQualDist> (dnaAlphabetSize))
{ }

#define QuaffParamWrite(X) out << #X ": " << X << endl
void QuaffParams::write (ostream& out) const {
  QuaffParamWrite(beginInsert);
  QuaffParamWrite(extendInsert);
  QuaffParamWrite(beginDelete);
  QuaffParamWrite(extendDelete);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    for (int j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

#define QuaffParamRead(X) X = val[#X]
void QuaffParams::read (istream& in) {
  map<string,double> val;
  regex paramVal ("^\\s*(\\S+)\\s*:\\s+([\\d\\+\\-\\.eE]+)\\s*$");
  regex nonWhite ("\\S");
  while (!in.eof()) {
    string line;
    getline(in,line);
    smatch sm;
    if (regex_match (line, sm, paramVal))
      val[sm.str(1)] = atof (sm.str(2).c_str());
    else if (regex_search (line, nonWhite))
      cerr << "Ignoring line: " << line;
  }

  QuaffParamRead(beginInsert);
  QuaffParamRead(extendInsert);
  QuaffParamRead(beginDelete);
  QuaffParamRead(extendDelete);

  for (int i = 0; i < dnaAlphabetSize; ++i)
    insert[i].read (val, string("insert") + dnaAlphabet[i]);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    for (int j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].read (val, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

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

void QuaffCounts::write (ostream& out) const {
  QuaffParamWrite(m2m);
  QuaffParamWrite(m2i);
  QuaffParamWrite(m2d);
  QuaffParamWrite(m2e);
  QuaffParamWrite(d2d);
  QuaffParamWrite(d2m);
  QuaffParamWrite(i2i);
  QuaffParamWrite(i2m);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    for (int j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

void Alignment::write (ostream& out) const {
  for (auto s : gappedSeq)
    s.writeFasta (out);
}

FastSeq Alignment::getUngapped (int row) const {
  const FastSeq& g = gappedSeq[row];
  FastSeq s = gappedSeq[row];
  s.seq.clear();
  s.qual.clear();
  for (int i = 0; i < g.length(); ++i) {
    char c = g.seq[i];
    if (c != '-' && c != '.') {
      s.seq.push_back (c);
      if (g.hasQual())
	s.qual.push_back (g.qual[i]);
    }
  }
  return s;
}

QuaffDPMatrix::QuaffDPMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp)
  : logger (qp.logger),
    px (&x),
    py (&y),
    qs (qp),
    xTok (x.tokens(dnaAlphabet)),
    yTok (y.tokens(dnaAlphabet)),
    yQual (y.qualScores()),
    xLen (x.length()),
    yLen (y.length()),
    mat (x.length() + 1, vector<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    ins (x.length() + 1, vector<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    del (x.length() + 1, vector<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    start (-numeric_limits<double>::infinity()),
    end (-numeric_limits<double>::infinity()),
    result (-numeric_limits<double>::infinity())
{ }

QuaffForwardMatrix::QuaffForwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp)
  : QuaffDPMatrix (x, y, qp)
{
  if (LogThisAt(2))
    initProgress ("Forward matrix");

  start = 0;
  for (int i = 1; i <= xLen; ++i) {
    if (LogThisAt(2))
      logProgress (i / (double) xLen, "base %d/%d", i, xLen);

    for (int j = 1; j <= yLen; ++j) {
      mat[i][j] = log_sum_exp (mat[i-1][j-1] + qs.m2m,
			       del[i-1][j-1] + qs.d2m,
			       ins[i-1][j-1] + qs.i2m);
      if (j == 1)
	mat[i][j] = log_sum_exp (mat[i][j],
				 start);

      mat[i][j] += matchScore(i,j);

      ins[i][j] = insertScore(j) + log_sum_exp (ins[i][j-1] + qs.i2i,
						mat[i][j-1] + qs.m2i);

      del[i][j] = log_sum_exp (del[i-1][j] + qs.d2d,
			       mat[i-1][j] + qs.m2d);
    }

    end = log_sum_exp (end,
		       mat[i][yLen] + qs.m2e);
  }

  result = end;
}

QuaffBackwardMatrix::QuaffBackwardMatrix (const QuaffForwardMatrix& fwd)
  : QuaffDPMatrix (*fwd.px, *fwd.py, *fwd.qs.pqp),
    pfwd (&fwd),
    qc()
{
  if (LogThisAt(2))
    initProgress ("Backward matrix");

  end = 0;
  for (int i = xLen; i > 0; --i) {
    if (LogThisAt(2))
      logProgress ((xLen - i) / (double) xLen, "base %d/%d", i, xLen);

    const double m2e = transCount (mat[i][yLen],
				   fwd.mat[i][yLen],
				   qs.m2e,
				   end);
    qc.m2e += m2e;

    for (int j = yLen; j > 0; --j) {
      
      if (i < xLen && j < yLen) {
	const double matEmit = matchScore(i+1,j+1);
	const double matDest = mat[i+1][j+1];
	double& matCount = matchCount(i+1,j+1);

	const double m2m = transCount (mat[i][j],
				       fwd.mat[i][j],
				       qs.m2m + matEmit,
				       matDest);
	qc.m2m += m2m;
	matCount += m2m;

      	const double i2m = transCount (ins[i][j],
				       fwd.ins[i][j],
				       qs.i2m + matEmit,
				       matDest);
	qc.i2m += i2m;
	matCount += i2m;

      	const double d2m = transCount (del[i][j],
				       fwd.del[i][j],
				       qs.d2m + matEmit,
				       matDest);
	qc.d2m += d2m;
	matCount += d2m;
      }

      if (i < xLen) {
	const double insEmit = insertScore(j+1);
	const double insDest = ins[i+1][j];
	double& insCount = insertCount(j+1);

	const double m2i = transCount (mat[i][j],
				       fwd.mat[i][j],
				       qs.m2i + insEmit,
				       insDest);
	qc.m2i += m2i;
	insCount += m2i;

      	const double i2i = transCount (ins[i][j],
				       fwd.ins[i][j],
				       qs.i2i + insEmit,
				       insDest);
	qc.i2i += i2i;
	insCount += i2i;
      }

      if (j < yLen) {
	const double delDest = del[i][j+1];

	const double m2d = transCount (mat[i][j],
				       fwd.mat[i][j],
				       qs.m2d,
				       delDest);
	qc.m2d += m2d;

      	const double d2d = transCount (del[i][j],
				       fwd.del[i][j],
				       qs.d2d,
				       delDest);
	qc.d2d += d2d;
      }
    }

    if (i < xLen && yLen > 0) {
      const double matEmit = matchScore(i+1,1);
      const double matDest = mat[i+1][1];
      double& matCount = matchCount(i+1,1);

      const double m2s = transCount (start,
				     fwd.start,
				     matEmit,
				     matDest);
      matCount += m2s;
    }
  }

  result = start;
  if (result != fwd.result)
    cerr << "Warning: forward score (" << fwd.result << ") does not match backward score (" << result << ")" << endl;
}

double QuaffBackwardMatrix::transCount (double& backSrc, double fwdSrc, double trans, double backDest) const {
  const double transBackDest = trans + backDest;
  const double count = exp (fwdSrc + transBackDest - pfwd->result);
  backSrc = log_sum_exp (backSrc, transBackDest);
  return count;
}
