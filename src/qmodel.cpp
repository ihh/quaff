#include <iostream>
#include <regex>
#include <list>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include "qmodel.h"
#include "logsumexp.h"
#include "negbinom.h"

#define MAX_FRACTIONAL_FWDBACK_ERROR .0001

SymQualDist::SymQualDist()
  : symProb(1. / dnaAlphabetSize),
    qualTrialSuccessProb(.5),
    qualNumSuccessfulTrials(FastSeq::qualScoreRange / 2)
{ }

void SymQualDist::write (ostream& out, const string& prefix) const {
  out << prefix << ": " << symProb << endl;
  out << prefix << "qp: " << qualTrialSuccessProb << endl;
  out << prefix << "qr: " << qualNumSuccessfulTrials << endl;
}

void SymQualDist::read (map<string,double>& paramVal, const string& prefix) {
  symProb = paramVal[prefix];
  qualTrialSuccessProb = paramVal[prefix + "qp"];
  qualNumSuccessfulTrials = paramVal[prefix + "qr"];
}

double SymQualDist::logQualProb (int k) const {
  return logNegativeBinomial (k, qualTrialSuccessProb, qualNumSuccessfulTrials);
}

double SymQualDist::logQualProb (const vguard<double>& kFreq) const {
  return logNegativeBinomial (kFreq, qualTrialSuccessProb, qualNumSuccessfulTrials);
}

SymQualScores::SymQualScores (const SymQualDist& sqd)
  : logSymQualProb (FastSeq::qualScoreRange)
{
  const double symLogProb = log (sqd.symProb);
  for (int k = 0; k < FastSeq::qualScoreRange; ++k)
    logSymQualProb[k] = symLogProb + sqd.logQualProb(k);
}

SymQualCounts::SymQualCounts()
  : qualCount (FastSeq::qualScoreRange, 0.)
{ }

void SymQualCounts::write (ostream& out, const string& prefix) const {
  if (qualCount.size()) {
    out << prefix << ": {";
    int n = 0;
    for (size_t i = 0; i < qualCount.size(); ++i)
      if (qualCount[i] > 0)
	out << (n++ ? ", " : "") << i << ": " << qualCount[i];
    out << '}' << endl;
  }
}

QuaffParams::QuaffParams()
  : beginInsert(.125),
    extendInsert(.5),
    beginDelete(.125),
    extendDelete(.5),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualDist> (dnaAlphabetSize))
{
  // give the default params a match bonus
  for (size_t i = 0; i < dnaAlphabetSize; ++i)
    for (size_t j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].symProb = (i == j ? .625 : .125);
}

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
    match (dnaAlphabetSize, vguard<SymQualScores> (dnaAlphabetSize))
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
    match (dnaAlphabetSize, vguard<SymQualCounts> (dnaAlphabetSize)),
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

QuaffParamCounts::QuaffParamCounts()
  : insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualCounts> (dnaAlphabetSize))
{
  initCounts (0);
}
    
QuaffParamCounts::QuaffParamCounts (const QuaffCounts& counts)
  : insert (counts.insert),
    match (counts.match),
    beginInsertNo (counts.m2m + counts.m2d),
    beginInsertYes (counts.m2i + counts.m2e),
    extendInsertNo (counts.i2m),
    extendInsertYes (counts.i2i),
    beginDeleteNo (counts.m2m),
    beginDeleteYes (counts.m2d),
    extendDeleteNo (counts.d2m),
    extendDeleteYes (counts.d2d)
{ }

void QuaffParamCounts::initCounts (double count) {
  for (auto& ins : insert)
    for (auto& qc : ins.qualCount)
      qc = count / FastSeq::qualScoreRange;
  for (auto& mx : match)
    for (auto& mxy : mx)
      for (auto& qc : mxy.qualCount)
	qc = count / FastSeq::qualScoreRange;
  beginInsertNo = count;
  beginInsertYes = count;
  extendInsertNo = count;
  extendInsertYes = count;
  beginDeleteNo = count;
  beginDeleteYes = count;
  extendDeleteNo = count;
  extendDeleteYes = count;
}

void QuaffParamCounts::write (ostream& out) const {
  QuaffParamWrite(beginInsertNo);
  QuaffParamWrite(beginInsertYes);
  QuaffParamWrite(extendInsertNo);
  QuaffParamWrite(extendInsertYes);
  QuaffParamWrite(beginDeleteNo);
  QuaffParamWrite(beginDeleteYes);
  QuaffParamWrite(extendDeleteNo);
  QuaffParamWrite(extendDeleteYes);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    for (int j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

Alignment& Alignment::addScoreComment() {
  if (rows())
    ((FastSeq&)gappedSeq[0]).comment = string("score(") + to_string(score) + ") " + gappedSeq[0].comment;
  return *this;
}

void Alignment::writeGappedFasta (ostream& out) const {
  for (const auto& s : gappedSeq)
    s.writeFasta (out);
}

void Alignment::writeStockholm (ostream& out) const {
  out << "# STOCKHOLM 1.0" << endl;
  out << "#=GF SC " << score << endl;
  size_t nameWidth = 0;
  for (const auto& s : gappedSeq)
    nameWidth = max (s.name.size(), nameWidth);
  for (const auto& s : gappedSeq) {
    const streamsize w = out.width (nameWidth);
    out << s.name;
    out.width (w);
    out << ' ' << s.seq << endl;
  }
  out << "//" << endl;
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
  : px (&x),
    py (&y),
    qs (qp),
    xTok (x.tokens(dnaAlphabet)),
    yTok (y.tokens(dnaAlphabet)),
    yQual (y.qualScores()),
    xLen (x.length()),
    yLen (y.length()),
    mat (x.length() + 1, vguard<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    ins (x.length() + 1, vguard<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    del (x.length() + 1, vguard<double> (y.length() + 1, -numeric_limits<double>::infinity())),
    cachedInsertEmitScore (y.length() + 1, -numeric_limits<double>::infinity()),
    start (-numeric_limits<double>::infinity()),
    end (-numeric_limits<double>::infinity()),
    result (-numeric_limits<double>::infinity())
{
  Assert (y.hasQual(), "Read sequences must have quality scores (FASTQ, not FASTA)");
  for (int j = 1; j <= y.length(); ++j)
    cachedInsertEmitScore[j] = insertEmitScore(j);
}

QuaffForwardMatrix::QuaffForwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp)
  : QuaffDPMatrix (x, y, qp)
{
  if (LogThisAt(2))
    initProgress ("Forward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

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

      mat[i][j] += matchEmitScore(i,j);

      ins[i][j] = cachedInsertEmitScore[j] + log_sum_exp (ins[i][j-1] + qs.i2i,
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
    initProgress ("Backward algorithm (%s vs %s)", px->name.c_str(), py->name.c_str());

  end = 0;
  for (size_t i = xLen; i > 0; --i) {
    if (LogThisAt(2))
      logProgress ((xLen - i) / (double) xLen, "base %d/%d", i, xLen);

    const double m2e = transCount (mat[i][yLen],
				   fwd.mat[i][yLen],
				   qs.m2e,
				   end);
    qc.m2e += m2e;

    for (size_t j = yLen; j > 0; --j) {
      
      const double matEmit = matchEmitScore(i,j);
      const double matDest = mat[i][j];
      double& matCount = matchCount(i,j);

      const double m2m = transCount (mat[i-1][j-1],
				     fwd.mat[i-1][j-1],
				     qs.m2m + matEmit,
				     matDest);
      qc.m2m += m2m;
      matCount += m2m;

      const double d2m = transCount (del[i-1][j-1],
				     fwd.del[i-1][j-1],
				     qs.d2m + matEmit,
				     matDest);
      qc.d2m += d2m;
      matCount += d2m;

      const double i2m = transCount (ins[i-1][j-1],
				     fwd.ins[i-1][j-1],
				     qs.i2m + matEmit,
				     matDest);
      qc.i2m += i2m;
      matCount += i2m;

      if (j == 1) {
	const double m2s = transCount (start,
				       fwd.start,
				       matEmit,
				       matDest);
	matCount += m2s;
      }

      const double insEmit = cachedInsertEmitScore[j];
      const double insDest = ins[i][j];
      double& insCount = insertCount(j);

      const double m2i = transCount (mat[i][j-1],
				     fwd.mat[i][j-1],
				     qs.m2i + insEmit,
				     insDest);
      qc.m2i += m2i;
      insCount += m2i;

      const double i2i = transCount (ins[i][j-1],
				     fwd.ins[i][j-1],
				     qs.i2i + insEmit,
				     insDest);
      qc.i2i += i2i;
      insCount += i2i;

      const double delDest = del[i][j];

      const double m2d = transCount (mat[i-1][j],
				     fwd.mat[i-1][j],
				     qs.m2d,
				     delDest);
      qc.m2d += m2d;

      const double d2d = transCount (del[i-1][j],
				     fwd.del[i-1][j],
				     qs.d2d,
				     delDest);
      qc.d2d += d2d;
    }
  }

  result = start;
  if (gsl_root_test_delta (result, fwd.result, 0, MAX_FRACTIONAL_FWDBACK_ERROR) != GSL_SUCCESS)
    cerr << endl << endl << "Warning: forward score (" << fwd.result << ") does not match backward score (" << result << ")" << endl << endl << endl;

  if (LogThisAt(4)) {
    cerr << "Forward-backward counts, " << px->name << " vs " << py->name << ':' << endl;
    qc.write (cerr);
  }
}

double QuaffBackwardMatrix::transCount (double& backSrc, double fwdSrc, double trans, double backDest) const {
  const double transBackDest = trans + backDest;
  const double count = exp (fwdSrc + transBackDest - pfwd->result);
  backSrc = log_sum_exp (backSrc, transBackDest);
  return count;
}

QuaffViterbiMatrix::QuaffViterbiMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp)
  : QuaffDPMatrix (x, y, qp)
{
  if (LogThisAt(2))
    initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (int i = 1; i <= xLen; ++i) {
    if (LogThisAt(2))
      logProgress (i / (double) xLen, "base %d/%d", i, xLen);

    for (int j = 1; j <= yLen; ++j) {
      mat[i][j] = max (max (mat[i-1][j-1] + qs.m2m,
			    del[i-1][j-1] + qs.d2m),
		       ins[i-1][j-1] + qs.i2m);
      if (j == 1)
	mat[i][j] = max (mat[i][j],
			 start);

      mat[i][j] += matchEmitScore(i,j);

      ins[i][j] = cachedInsertEmitScore[j] + max (ins[i][j-1] + qs.i2i,
						  mat[i][j-1] + qs.m2i);

      del[i][j] = max (del[i-1][j] + qs.d2d,
		       mat[i-1][j] + qs.m2d);
    }

    end = max (end,
	       mat[i][yLen] + qs.m2e);
  }

  result = end;
}

void QuaffViterbiMatrix::updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx) {
  if (candidateMax > currentMax) {
    currentMax = candidateMax;
    currentMaxIdx = candidateMaxIdx;
  }
}

Alignment QuaffViterbiMatrix::alignment() const {
  size_t xEnd = -1;
  double bestEndSc = -numeric_limits<double>::infinity(), sc;
    for (size_t iEnd = xLen; iEnd > 0; --iEnd) {
        sc = mat[iEnd][yLen] + qs.m2e;
        if (iEnd == xLen || sc > bestEndSc) {
            bestEndSc = sc;
            xEnd = iEnd;
        }
    }
  size_t i = xEnd, j = yLen;
  list<char> xRow, yRow;
  State state = Match;
  while (state != Start) {
    double srcSc = -numeric_limits<double>::infinity();
    double emitSc = 0;
    switch (state) {
    case Match:
      emitSc = matchEmitScore(i,j);
      xRow.push_front (px->seq[--i]);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat[i][j] + qs.m2m + emitSc, Match);
      updateMax (srcSc, state, ins[i][j] + qs.i2m + emitSc, Insert);
      updateMax (srcSc, state, del[i][j] + qs.d2m + emitSc, Delete);
      if (j == 0)
	updateMax (srcSc, state, emitSc, Start);
      break;

    case Insert:
      emitSc = cachedInsertEmitScore[j];
      xRow.push_front (gapChar);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat[i][j] + qs.m2i + emitSc, Match);
      updateMax (srcSc, state, ins[i][j] + qs.i2i + emitSc, Insert);
      break;

    case Delete:
      xRow.push_front (px->seq[--i]);
      yRow.push_front (gapChar);
      updateMax (srcSc, state, mat[i][j] + qs.m2d, Match);
      updateMax (srcSc, state, ins[i][j] + qs.d2d, Delete);
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  const size_t xStart = i + 1;
  Alignment align(2);
  align.gappedSeq[0].name = "substr(" + px->name + "," + to_string(xStart) + ".." + to_string(xEnd) + ")";
  align.gappedSeq[1].name = py->name;
  align.gappedSeq[0].seq = string (xRow.begin(), xRow.end());
  align.gappedSeq[1].seq = string (yRow.begin(), yRow.end());
  align.score = result;
  return align;
}

Alignment QuaffViterbiMatrix::scoreAdjustedAlignment (const QuaffNullParams& nullModel) const {
  Alignment a = alignment();
  a.score -= nullModel.logLikelihood (*py);
  return a;
}

QuaffTrainer::QuaffTrainer()
  : maxIterations (100),
    minFractionalLoglikeIncrement (.01)
{ }

bool QuaffTrainer::parseTrainingArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const char* arg = argv[0];
    if (strcmp (arg, "-maxiter") == 0) {
      Assert (argc > 1, "%s must have an argument", arg);
      const char* val = argv[1];
      maxIterations = atoi (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (strcmp (arg, "-mininc") == 0) {
      Assert (argc > 1, "%s must have an argument", arg);
      const char* val = argv[1];
      minFractionalLoglikeIncrement = atof (val);
      argv += 2;
      argc -= 2;
      return true;
    }
  }

  return false;
}

void QuaffParamCounts::addWeighted (const QuaffParamCounts& counts, double weight) {
  for (int q = 0; q < FastSeq::qualScoreRange; ++q)
    for (int i = 0; i < dnaAlphabetSize; ++i) {
      insert[i].qualCount[q] += weight * counts.insert[i].qualCount[q];
      for (int j = 0; j < dnaAlphabetSize; ++j)
	match[i][j].qualCount[q] += weight * counts.match[i][j].qualCount[q];
    }
  beginInsertNo += weight * counts.beginInsertNo;
  beginInsertYes += weight * counts.beginInsertYes;
  extendInsertNo += weight * counts.extendInsertNo;
  extendInsertYes += weight * counts.extendInsertYes;
  beginDeleteNo += weight * counts.beginDeleteNo;
  beginDeleteYes += weight * counts.beginDeleteYes;
  extendDeleteNo += weight * counts.extendDeleteNo;
  extendDeleteYes += weight * counts.extendDeleteYes;
}

double QuaffParamCounts::logPrior (const QuaffParams& qp) const {
  double lp = 0;
  double *alpha = new double[dnaAlphabetSize], *theta = new double[dnaAlphabetSize];
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    lp += qp.insert[i].logQualProb (insert[i].qualCount);  // not normalized...
    theta[i] = qp.insert[i].symProb;
    alpha[i] = accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.);
  }
  lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    for (int j = 0; j < dnaAlphabetSize; ++j) {
      lp += qp.match[i][j].logQualProb (match[i][j].qualCount);  // not normalized...
      theta[j] = qp.match[i][j].symProb;
      alpha[j] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
    }
    lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
  }
  delete[] theta;
  delete[] alpha;
  return lp;
}

double QuaffParamCounts::expectedLogLike (const QuaffParams& qp) const {
  double ll = 0;
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    ll += qp.insert[i].logQualProb (insert[i].qualCount);
    ll += log (qp.insert[i].symProb) * accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.);
  }
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    for (int j = 0; j < dnaAlphabetSize; ++j) {
      ll += qp.match[i][j].logQualProb (match[i][j].qualCount);  // not normalized...
      ll += log (qp.match[i][j].symProb) * accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
    }
  }
  return ll;
}

QuaffParams QuaffParamCounts::fit() const {
  QuaffParams qp;
  qp.beginDelete = 1. / (1. + beginDeleteNo / beginDeleteYes);
  qp.extendDelete = 1. / (1. + extendDeleteNo / extendDeleteYes);
  qp.beginInsert = 1. / (1. + beginInsertNo / beginInsertYes);
  qp.extendInsert = 1. / (1. + extendInsertNo / extendInsertYes);

  vguard<double> insFreq;
  for (int i = 0; i < dnaAlphabetSize; ++i)
    insFreq.push_back (accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.));
  const double insNorm = accumulate (insFreq.begin(), insFreq.end(), 0.);
  for (int i = 0; i < dnaAlphabetSize; ++i) {
    qp.insert[i].symProb = insFreq[i] / insNorm;
    fitNegativeBinomial (insert[i].qualCount, qp.insert[i].qualTrialSuccessProb, qp.insert[i].qualNumSuccessfulTrials);
  }

  for (int i = 0; i < dnaAlphabetSize; ++i) {
    vguard<double> iMatFreq (dnaAlphabetSize);
    for (int j = 0; j < dnaAlphabetSize; ++j)
      iMatFreq[j] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
    const double iMatNorm = accumulate (iMatFreq.begin(), iMatFreq.end(), 0.);
    for (int j = 0; j < dnaAlphabetSize; ++j) {
      qp.match[i][j].symProb = iMatFreq[j] / iMatNorm;
      fitNegativeBinomial (match[i][j].qualCount, qp.match[i][j].qualTrialSuccessProb, qp.match[i][j].qualNumSuccessfulTrials);
    }
  }

  return qp;
}

QuaffForwardBackwardMatrix::QuaffForwardBackwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp)
  : fwd (x, y, qp),
    back (fwd)
{ }

QuaffNullParams::QuaffNullParams (const vguard<FastSeq>& seqs, double pseudocount)
  : null (dnaAlphabetSize)
{
  vguard<SymQualCounts> nullCount (dnaAlphabetSize);
  for (auto& sqc : nullCount)
    for (auto& c : sqc.qualCount)
      c += pseudocount / FastSeq::qualScoreRange;
  double nullEmitYes = pseudocount, nullEmitNo = pseudocount;
  vguard<double> symCount (dnaAlphabetSize, pseudocount);
  
  for (const auto& s : seqs) {
    ++nullEmitNo;
    nullEmitYes += s.length();
    const vguard<int> tok = s.tokens (dnaAlphabet);
    for (size_t i = 0; i < s.length(); ++i) {
      ++symCount[tok[i]];
      ++nullCount[tok[i]].qualCount[s.getQualScoreAt(i)];
    }
  }

  nullEmit = 1 / (1 + nullEmitNo / nullEmitYes);
  const double symCountNorm = accumulate (symCount.begin(), symCount.end(), 0.);
  for (size_t n = 0; n < dnaAlphabetSize; ++n) {
    null[n].symProb = symCount[n] / symCountNorm;
    fitNegativeBinomial (nullCount[n].qualCount, null[n].qualTrialSuccessProb, null[n].qualNumSuccessfulTrials);
  }
}

double QuaffNullParams::logLikelihood (const vguard<FastSeq>& seqs) const {
  double ll = 0;
  for (const auto& s : seqs)
    ll += logLikelihood (s);
  return ll;
}

double QuaffNullParams::logLikelihood (const FastSeq& s) const {
  double ll = s.length() * log(nullEmit) + log(1. - nullEmit);
  const vguard<int> tok = s.tokens (dnaAlphabet);
  for (size_t i = 0; i < s.length(); ++i)
    ll += log (null[tok[i]].symProb) + null[tok[i]].logQualProb (s.getQualScoreAt(i));
  return ll;
}

void QuaffNullParams::write (ostream& out) const {
  QuaffParamWrite(nullEmit);
  for (int i = 0; i < dnaAlphabetSize; ++i)
    null[i].write (out, string("null") + dnaAlphabet[i]);
}

QuaffParams QuaffTrainer::fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffParamCounts& pseudocounts) {
  QuaffParams qp = seed;
  const QuaffNullParams qnp (y);
  if (LogThisAt(3))
    qnp.write (cerr << "Null model:" << endl);
  double prevLogLikeWithPrior = -numeric_limits<double>::infinity();
  for (int iter = 0; iter < maxIterations; ++iter) {
    QuaffParamCounts counts;
    double logLike = 0;
    for (const auto& yfs : y) {
      const double yNullLogLike = qnp.logLikelihood (yfs);
      double yLogLike = yNullLogLike;  // this initial value allows null model to "win"
      vguard<double> xyLogLike;
      vguard<QuaffParamCounts> xyCounts;
      for (const auto& xfs : x) {
	const QuaffForwardBackwardMatrix fb (xfs, yfs, qp);
	const double ll = fb.fwd.result;
	const QuaffParamCounts qpc (fb.back.qc);
	xyLogLike.push_back (ll);
	xyCounts.push_back (qpc);
	yLogLike = log_sum_exp (yLogLike, ll);
      }
      for (size_t nx = 0; nx < x.size(); ++nx) {
	const double xyPostProb = exp (xyLogLike[nx] - yLogLike);
	if (LogThisAt(2))
	  cerr << "P(read " << yfs.name << " derived from ref " << x[nx].name << ") = " << xyPostProb << endl;
	counts.addWeighted (xyCounts[nx], xyPostProb);
      }
      if (LogThisAt(2))
	cerr << "P(read " << yfs.name << " unrelated to refs) = " << exp(yNullLogLike - yLogLike) << endl;
      logLike += yLogLike;
    }
    const double logPrior = pseudocounts.logPrior (qp);
    const double logLikeWithPrior = logLike + logPrior;
    if (LogThisAt(1))
      cerr << "EM iteration " << (iter+1) << ": log-likelihood (" << logLike << ") + log-prior (" << logPrior << ") = " << logLikeWithPrior << endl;
    if (iter > 0 && logLikeWithPrior < prevLogLikeWithPrior + abs(prevLogLikeWithPrior)*minFractionalLoglikeIncrement)
      break;
    prevLogLikeWithPrior = logLikeWithPrior;
    if (LogThisAt(3)) {
      cerr << "Parameter counts computed during E-step:" << endl;
      counts.write (cerr);
    }
    counts.addWeighted (pseudocounts, 1.);
    const double oldExpectedLogLike = counts.expectedLogLike (qp);
    qp = counts.fit();
    const double newExpectedLogLike = counts.expectedLogLike (qp);
    if (LogThisAt(2))
      cerr << "Expected log-likelihood went from " << oldExpectedLogLike << " to " << newExpectedLogLike << " during M-step" << endl;
  }
  return qp;
}

QuaffAligner::QuaffAligner()
  : format (StockholmAlignment),
    logOddsThreshold (0),
    printAllAlignments (false)
{ }

bool QuaffAligner::parseAlignmentArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const char* arg = argv[0];
    if (strcmp (arg, "-format") == 0) {
      Assert (argc > 1, "%s must have an argument", arg);
      const string fmt = argv[1];
      if (fmt == "fasta")
	format = GappedFastaAlignment;
      else if (fmt == "stockholm")
	format = StockholmAlignment;
      else if (fmt == "refseq")
	format = UngappedFastaRef;
      else
	Abort ("Unknown format: %s", fmt.c_str());
      argv += 2;
      argc -= 2;
      return true;

    } else if (strcmp (arg, "-threshold") == 0) {
      Assert (argc > 1, "%s must have an argument", arg);
      const char* val = argv[1];
      logOddsThreshold = atof (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (strcmp (arg, "-nothreshold") == 0) {
      logOddsThreshold = -numeric_limits<double>::infinity();
      argv += 1;
      argc -= 1;
      return true;

    } else if (strcmp (arg, "-printall") == 0) {
      printAllAlignments = true;
      argv += 1;
      argc -= 1;
      return true;

    }
  }

  return false;
}

void QuaffAligner::align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params) {
  const QuaffNullParams qnp (y);
  if (LogThisAt(3))
    qnp.write (cerr << "Null model:" << endl);
  for (const auto& yfs : y) {
    size_t nBestAlign = 0;
    vguard<Alignment> xyAlign;
    for (const auto& xfs : x) {
      const QuaffViterbiMatrix viterbi (xfs, yfs, params);
      const Alignment align = viterbi.scoreAdjustedAlignment(qnp).addScoreComment();
      if (xyAlign.empty() || align.score > xyAlign[nBestAlign].score)
	nBestAlign = xyAlign.size();
      xyAlign.push_back (align);
    }
    if (printAllAlignments)
      for (const auto& a : xyAlign)
	writeAlignment (out, a);
    else
      writeAlignment (out, xyAlign[nBestAlign]);
  }
}

void QuaffAligner::writeAlignment (ostream& out, const Alignment& align) const {
  if (align.score >= logOddsThreshold) {
    FastSeq ref;
    switch (format) {
    case GappedFastaAlignment:
      align.writeGappedFasta (out);
      out << endl;
      break;

    case StockholmAlignment:
      align.writeStockholm (out);
      break;

    case UngappedFastaRef:
      Assert (align.rows() == 2, "Not a pairwise alignment");
      ref = align.getUngapped(0);
      ref.comment = string("matches(") + align.gappedSeq[1].name + ") " + ref.comment;
      ref.writeFasta (out);
      break;

    default:
      Abort ("Unrecognized alignment format");
      break;
    }
  }
}

