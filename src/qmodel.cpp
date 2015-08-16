#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include "qmodel.h"
#include "logsumexp.h"
#include "negbinom.h"
#include "logger.h"

// internal #defines
#define MAX_FRACTIONAL_FWDBACK_ERROR .0001
#define MAX_TRAINING_LOG_DELTA 20

// useful helper methods
double logBetaPdf (double prob, double yesCount, double noCount) {
  return log (gsl_ran_beta_pdf (prob, yesCount + 1, noCount + 1));
}

map<string,string> readParamFile (istream& in) {
  map<string,string> val;
  regex paramVal ("^\\s*(\\S+)\\s*:\\s+(.*)$");
  regex nonWhite ("\\S");
  while (!in.eof()) {
    string line;
    getline(in,line);
    smatch sm;
    if (regex_match (line, sm, paramVal))
      val[sm.str(1)] = sm.str(2);
    else if (regex_search (line, nonWhite))
      cerr << "Ignoring line: " << line;
  }
  return val;
}

// main method bodies
SymQualDist::SymQualDist()
  : symProb(1. / dnaAlphabetSize),
    qualTrialSuccessProb(.5),
    qualNumSuccessfulTrials(FastSeq::qualScoreRange / 2)
{ }

void SymQualDist::write (ostream& out, const string& prefix) const {
  out << prefix << ": " << symProb << endl;
  out << prefix << "qp: " << qualTrialSuccessProb
      << "\t# mean " << negativeBinomialMean(qualTrialSuccessProb,qualNumSuccessfulTrials)
      << endl;
  out << prefix << "qr: " << qualNumSuccessfulTrials
      << "\t# stdev " << sqrt(negativeBinomialVariance(qualTrialSuccessProb,qualNumSuccessfulTrials))
      << endl;
}

void SymQualDist::read (map<string,string>& paramVal, const string& prefix) {
  const string qp = prefix + "qp", qr = prefix + "qr";
  Assert (paramVal.find(prefix) != paramVal.end(), "Missing parameter: %s", prefix.c_str());
  Assert (paramVal.find(qp) != paramVal.end(), "Missing parameter: %s", qp.c_str());
  Assert (paramVal.find(qr) != paramVal.end(), "Missing parameter: %s", qr.c_str());
  symProb = atof (paramVal[prefix].c_str());
  qualTrialSuccessProb = atof (paramVal[qp].c_str());
  qualNumSuccessfulTrials = atof (paramVal[qr].c_str());
}

double SymQualDist::logQualProb (QualScore k) const {
  return logNegativeBinomial (k, qualTrialSuccessProb, qualNumSuccessfulTrials);
}

double SymQualDist::logQualProb (const vguard<double>& kFreq) const {
  return logNegativeBinomial (kFreq, qualTrialSuccessProb, qualNumSuccessfulTrials);
}

SymQualScores::SymQualScores (const SymQualDist& sqd)
  : logSymQualProb (FastSeq::qualScoreRange)
{
  logSymProb = log (sqd.symProb);
  for (QualScore k = 0; k < FastSeq::qualScoreRange; ++k)
    logSymQualProb[k] = logSymProb + sqd.logQualProb(k);
}

SymQualCounts::SymQualCounts()
  : qualCount (FastSeq::qualScoreRange, 0.)
{ }

void SymQualCounts::write (ostream& out, const string& prefix) const {
  if (qualCount.size()) {
    out << prefix << ": {";
    int n = 0;
    for (QualScore i = 0; i < qualCount.size(); ++i)
      if (qualCount[i] > 0)
	out << (n++ ? ", " : "") << i << ": " << qualCount[i];
    out << '}' << endl;
  }
}

void SymQualCounts::read (const string& counts) {
  string c = counts;
  qualCount = vguard<double> (FastSeq::qualScoreRange, 0.);
  regex countRegex ("(\\d+)\\s*:\\s*([\\d\\+\\-eE\\.]+)");
  smatch sm;
  while (regex_search (c, sm, countRegex)) {
    qualCount[atoi (sm.str(1).c_str())] = atof (sm.str(2).c_str());
    c = sm.suffix().str();
  }
}

QuaffParams::QuaffParams()
  : refEmit(.5),
    refBase(dnaAlphabetSize,.25),
    beginInsert(.5),
    extendInsert(.5),
    beginDelete(.5),
    extendDelete(.5),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualDist> (dnaAlphabetSize))
{ }

#define QuaffParamWrite(X) out << #X ": " << X << endl
void QuaffParams::write (ostream& out) const {
  QuaffParamWrite(refEmit);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    out << "refBase" << dnaAlphabet[i] << ": " << refBase[i] << endl;
  QuaffParamWrite(beginInsert);
  QuaffParamWrite(extendInsert);
  QuaffParamWrite(beginDelete);
  QuaffParamWrite(extendDelete);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

#define QuaffParamRead(X) Assert(val.find(#X) != val.end(),"Couldn't read " #X), X = atof(val[#X].c_str())
void QuaffParams::read (istream& in) {
  map<string,string> val = readParamFile (in);

  QuaffParamRead(beginInsert);
  QuaffParamRead(extendInsert);
  QuaffParamRead(beginDelete);
  QuaffParamRead(extendDelete);

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].read (val, string("insert") + dnaAlphabet[i]);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].read (val, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

void QuaffParams::fitRefSeqs (const vguard<FastSeq>& refs) {
  int totalLen;
  vguard<int> baseCount (dnaAlphabetSize, 0.);
  for (const auto& fs : refs) {
    totalLen += fs.length();
    for (auto i : fs.tokens(dnaAlphabet))
      ++baseCount[i];
  }
  refEmit = 1 / (1 + refs.size() / (double) totalLen);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    refBase[i] = baseCount[i] / (double) totalLen;
}

QuaffScores::QuaffScores (const QuaffParams& qp)
  : pqp(&qp),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualScores> (dnaAlphabetSize))
{
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    insert[i] = SymQualScores (qp.insert[i]);
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
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
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

QuaffParamCounts::QuaffParamCounts()
  : insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualCounts> (dnaAlphabetSize))
{
  zeroCounts();
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

void QuaffParamCounts::zeroCounts() {
  initCounts (0, 0, 0, 0, NULL);
}

void QuaffParamCounts::initCounts (double noBeginCount, double yesExtendCount, double matchIdentCount, double otherCount, const QuaffNullParams* nullModel) {
  for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
    for (QualScore k = 0; k < FastSeq::qualScoreRange; ++k)
      if (nullModel)
	insert[j].qualCount[k] = otherCount * nullModel->null[j].symProb * dnaAlphabetSize * gsl_ran_negative_binomial_pdf (k, nullModel->null[j].qualTrialSuccessProb, nullModel->null[j].qualNumSuccessfulTrials);
      else
	insert[j].qualCount[k] = otherCount / FastSeq::qualScoreRange;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      for (QualScore k = 0; k < FastSeq::qualScoreRange; ++k)
	if (nullModel)
	  match[i][j].qualCount[k] = (i == j ? matchIdentCount : (otherCount * nullModel->null[j].symProb * dnaAlphabetSize / (1 - nullModel->null[i].symProb))) * gsl_ran_negative_binomial_pdf (k, nullModel->null[j].qualTrialSuccessProb, nullModel->null[j].qualNumSuccessfulTrials);
	else
	  match[i][j].qualCount[k] = (i == j ? matchIdentCount : otherCount) / FastSeq::qualScoreRange;
  beginInsertNo = noBeginCount;
  beginInsertYes = otherCount;
  extendInsertNo = otherCount;
  extendInsertYes = yesExtendCount;
  beginDeleteNo = noBeginCount;
  beginDeleteYes = otherCount;
  extendDeleteNo = otherCount;
  extendDeleteYes = yesExtendCount;
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
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, string("insert") + dnaAlphabet[i]);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].write (out, string("match") + dnaAlphabet[i] + dnaAlphabet[j]);
}

void QuaffParamCounts::read (istream& in) {
  map<string,string> val = readParamFile (in);

  QuaffParamRead(beginInsertYes);
  QuaffParamRead(extendInsertYes);
  QuaffParamRead(beginDeleteYes);
  QuaffParamRead(extendDeleteYes);

  QuaffParamRead(beginInsertNo);
  QuaffParamRead(extendInsertNo);
  QuaffParamRead(beginDeleteNo);
  QuaffParamRead(extendDeleteNo);

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].read (val[string("insert") + dnaAlphabet[i]]);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      match[i][j].read (val[string("match") + dnaAlphabet[i] + dnaAlphabet[j]]);
}

const char Alignment::gapChar = '-';
const char Alignment::mismatchChar = ':';

void Alignment::writeGappedFasta (ostream& out) const {
  for (const auto& s : gappedSeq)
    s.writeFasta (out);
}

void Alignment::writeStockholm (ostream& out) const {
  vguard<string> rowName, rowData;
  for (const auto& s : gappedSeq) {
    rowName.push_back (s.name);
    rowData.push_back (s.seq);
  }

  if (rows() == 2) {
    string cons;
    for (SeqIdx pos = 0; pos < columns(); ++pos) {
      const char c0 = toupper(gappedSeq[0].seq[pos]),
	c1 = toupper(gappedSeq[1].seq[pos]);
      cons.push_back ((isGapChar(c0) || isGapChar(c1))
		      ? gapChar
		      : (c0 == c1 ? c0 : mismatchChar));
    }
    rowName.insert (rowName.begin() + 1, string ("#=GC id"));
    rowData.insert (rowData.begin() + 1, cons);
  }

  size_t nameWidth = 0;
  for (const auto& s : rowName)
    nameWidth = max (s.size(), nameWidth);

  const size_t dataWidth = max (nameWidth, 79 - nameWidth);  // nice 80-column width
  
  out << "# STOCKHOLM 1.0" << endl;
  out << "#=GF Score " << score << endl;
  for (const auto& s : gappedSeq)
    if (s.comment.size())
      out << "#=GS CC " << s.name << ' ' << s.comment << endl;
  for (size_t col = 0; col < columns(); col += dataWidth) {
    if (col > 0)
      out << endl;
    for (size_t row = 0; row < rowName.size(); ++row) {
      const streamsize w = out.width (nameWidth);
      out << std::left << rowName[row];
      out.width (w);
      out << ' ' << rowData[row].substr(col,dataWidth) << endl;
    }
  }
  out << "//" << endl;
}

FastSeq Alignment::getUngapped (size_t row) const {
  const FastSeq& g = gappedSeq[row];
  FastSeq s = gappedSeq[row];
  s.seq.clear();
  s.qual.clear();
  for (SeqIdx i = 0; i < g.length(); ++i) {
    char c = g.seq[i];
    if (!isGapChar(c)) {
      s.seq.push_back (c);
      if (g.hasQual())
	s.qual.push_back (g.qual[i]);
    }
  }
  return s;
}

bool QuaffDPConfig::parseConfigArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-global") {
      local = false;
      argv += 1;
      argc -= 1;
      return true;
    }
  }

  return parseOverlapConfigArgs (argc, argv);
}

bool QuaffDPConfig::parseOverlapConfigArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-kmatchband") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      bandSize = atoi (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-kmatch") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      kmerLen = atoi (val);
      Assert (kmerLen >= 5 && kmerLen <= 16, "%s out of range (%d). Try 5 to 16", arg.c_str(), kmerLen);
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-kmatchn") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      kmerThreshold = atoi (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-kmatchsd") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      kmerStDevThreshold = atof (val);
      kmerThreshold = -1;
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-dense") {
      sparse = false;
      argv += 1;
      argc -= 1;
      return true;

    }
  }

  return false;
}

DiagonalEnvelope QuaffDPConfig::makeEnvelope (const FastSeq& x, const FastSeq& y) const {
  DiagonalEnvelope env (x, y);
  if (sparse)
    env.initSparse (kmerLen, bandSize, kmerThreshold);
  else
    env.initFull();
  return env;
}

QuaffDPCell::QuaffDPCell()
  : mat (-numeric_limits<double>::infinity()),
    ins (-numeric_limits<double>::infinity()),
    del (-numeric_limits<double>::infinity())
{ }

QuaffDPMatrixContainer::QuaffDPMatrixContainer (const DiagonalEnvelope& env)
  : penv (&env),
    px (env.px),
    py (env.py),
    xLen (px->length()),
    yLen (py->length()),
    cell (py->length() + 1),
    start (-numeric_limits<double>::infinity()),
    end (-numeric_limits<double>::infinity()),
    result (-numeric_limits<double>::infinity())
{ }

QuaffDPCell QuaffDPMatrixContainer::dummy;

double QuaffDPMatrixContainer::cellScore (SeqIdx i, SeqIdx j, State state) const {
  double cs = numeric_limits<double>::quiet_NaN();
  switch (state) {
  case Match:
    cs = mat(i,j);
    break;
  case Insert:
    cs = ins(i,j);
    break;
  case Delete:
    cs = del(i,j);
    break;
  default:
    break;
  }
  return cs;
}

const char* QuaffDPMatrixContainer::stateToString (State state) {
  const char* s = "Unknown";
  switch (state) {
  case Start:
    s = "Start";
    break;
  case Match:
    s = "Match";
    break;
  case Insert:
    s = "Insert";
    break;
  case Delete:
    s = "Delete";
    break;
  default:
    break;
  }
  return s;
}

void QuaffDPMatrixContainer::updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx) {
  if (candidateMax > currentMax) {
    currentMax = candidateMax;
    currentMaxIdx = candidateMaxIdx;
  }
}

QuaffDPMatrix::QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : QuaffDPMatrixContainer (env),
    pconfig (&config),
    qs (qp),
    xTok (px->tokens(dnaAlphabet)),
    yTok (py->tokens(dnaAlphabet)),
    yQual (py->qualScores()),
    cachedInsertEmitScore (py->length() + 1, -numeric_limits<double>::infinity())
{
  Assert (py->hasQual(), "Read sequences must have quality scores (FASTQ, not FASTA)");
  for (SeqIdx j = 1; j <= py->length(); ++j)
    cachedInsertEmitScore[j] = insertEmitScore(j);
}

void QuaffDPMatrix::write (ostream& out) const {
  for (SeqIdx j = 1; j <= yLen; ++j) {
    for (SeqIdx i : penv->forward_i(j))
      out << "i=" << i << ":" << px->seq[i-1] << " j=" << j << ":" << py->seq[j-1] << (py->hasQual() ? string(1,py->qual[j-1]) : string()) << "\tmat " << mat(i,j) << "\tins " << ins(i,j) << "\tdel " << del(i,j) << endl;
    out << endl;
  }
  out << "result " << result << endl;
}

QuaffForwardMatrix::QuaffForwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : QuaffDPMatrix (env, qp, config)
{
  const FastSeq& x (*px);
  const FastSeq& y (*py);
  if (LogThisAt(2))
    initProgress ("Forward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    if (LogThisAt(2))
      logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (SeqIdx i : env.forward_i(j)) {

      mat(i,j) = log_sum_exp (mat(i-1,j-1) + qs.m2m,
			      del(i-1,j-1) + qs.d2m,
			      ins(i-1,j-1) + qs.i2m);

      if (j == 1 && (i == 1 || config.local))
	mat(i,j) = log_sum_exp (mat(i,j),
				start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = cachedInsertEmitScore[j] + log_sum_exp (ins(i,j-1) + qs.i2i,
							 mat(i,j-1) + qs.m2i);

      del(i,j) = log_sum_exp (del(i-1,j) + qs.d2d,
			      mat(i-1,j) + qs.m2d);

      if (j == yLen && (i == xLen || config.local))
	end = log_sum_exp (end,
			   mat(i,yLen) + qs.m2e);
    }
  }

  result = end;

  if (LogThisAt(2))
    cerr << "Forward score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    write (cerr);
}

QuaffBackwardMatrix::QuaffBackwardMatrix (const QuaffForwardMatrix& fwd)
  : QuaffDPMatrix (*fwd.penv, *fwd.qs.pqp, *fwd.pconfig),
    pfwd (&fwd),
    qc()
{
  Assert (py->hasQual(), "Forward-Backward algorithm requires quality scores to fit model, but sequence %s lacks quality scores", py->name.c_str());

  if (LogThisAt(2))
    initProgress ("Backward algorithm (%s vs %s)", px->name.c_str(), py->name.c_str());

  end = 0;
  for (SeqIdx j = yLen; j > 0; --j) {

    if (LogThisAt(2))
      logProgress ((yLen - j) / (double) yLen, "base %d/%d", j, yLen);

    for (SeqIdx i : penv->reverse_i(j)) {

      if (j == yLen && (i == xLen || pconfig->local)) {
	const double m2e = transCount (mat(i,yLen),
				       fwd.mat(i,yLen),
				       qs.m2e,
				       end);
	qc.m2e += m2e;
      }

      const double matEmit = matchEmitScore(i,j);
      const double matDest = mat(i,j);
      double& matCount = matchCount(i,j);

      const double m2m = transCount (mat(i-1,j-1),
				     fwd.mat(i-1,j-1),
				     qs.m2m + matEmit,
				     matDest);
      qc.m2m += m2m;
      matCount += m2m;

      const double d2m = transCount (del(i-1,j-1),
				     fwd.del(i-1,j-1),
				     qs.d2m + matEmit,
				     matDest);
      qc.d2m += d2m;
      matCount += d2m;

      const double i2m = transCount (ins(i-1,j-1),
				     fwd.ins(i-1,j-1),
				     qs.i2m + matEmit,
				     matDest);
      qc.i2m += i2m;
      matCount += i2m;

      if (j == 1 && (i == 1 || pconfig->local)) {
	const double s2m = transCount (start,
				       fwd.start,
				       matEmit,
				       matDest);
	matCount += s2m;
      }

      const double insEmit = cachedInsertEmitScore[j];
      const double insDest = ins(i,j);
      double& insCount = insertCount(j);

      const double m2i = transCount (mat(i,j-1),
				     fwd.mat(i,j-1),
				     qs.m2i + insEmit,
				     insDest);
      qc.m2i += m2i;
      insCount += m2i;

      const double i2i = transCount (ins(i,j-1),
				     fwd.ins(i,j-1),
				     qs.i2i + insEmit,
				     insDest);
      qc.i2i += i2i;
      insCount += i2i;

      const double delDest = del(i,j);

      const double m2d = transCount (mat(i-1,j),
				     fwd.mat(i-1,j),
				     qs.m2d,
				     delDest);
      qc.m2d += m2d;

      const double d2d = transCount (del(i-1,j),
				     fwd.del(i-1,j),
				     qs.d2d,
				     delDest);
      qc.d2d += d2d;
    }
  }

  result = start;

  if (LogThisAt(2))
    cerr << "Backward score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    write (cerr);

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

QuaffViterbiMatrix::QuaffViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : QuaffDPMatrix (env, qp, config)
{
  const FastSeq& x (*px);
  const FastSeq& y (*py);
  if (LogThisAt(2))
    initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    if (LogThisAt(2))
      logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (SeqIdx i : env.forward_i(j)) {

      mat(i,j) = max (max (mat(i-1,j-1) + qs.m2m,
			   del(i-1,j-1) + qs.d2m),
		      ins(i-1,j-1) + qs.i2m);

      if (j == 1 && (i == 1 || config.local))
	mat(i,j) = max (mat(i,j),
			start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = cachedInsertEmitScore[j] + max (ins(i,j-1) + qs.i2i,
						 mat(i,j-1) + qs.m2i);

      del(i,j) = max (del(i-1,j) + qs.d2d,
		      mat(i-1,j) + qs.m2d);

      if (j == yLen && (i == xLen || config.local))
	end = max (end,
		   mat(i,j) + qs.m2e);
    }
  }

  result = end;

  if (LogThisAt(2))
    cerr << "Viterbi score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    write (cerr);
}

Alignment QuaffViterbiMatrix::alignment() const {
  Assert (resultIsFinite(), "Can't do Viterbi traceback if final score is -infinity");
  SeqIdx xEnd = xLen;
  if (pconfig->local) {
    double bestEndSc = -numeric_limits<double>::infinity();
    double sc;
    for (SeqIdx iEnd = xLen; iEnd > 0; --iEnd) {
      sc = mat(iEnd,yLen) + qs.m2e;
      if (iEnd == xLen || sc > bestEndSc) {
	bestEndSc = sc;
	xEnd = iEnd;
      }
    }
  }
  SeqIdx i = xEnd, j = yLen;
  list<char> xRow, yRow;
  State state = Match;
  while (state != Start) {
    if (LogThisAt(7))
      cerr << "Traceback: i=" << i << " j=" << j << " state=" << stateToString(state) << " score=" << cellScore(i,j,state) << endl;
    double srcSc = -numeric_limits<double>::infinity();
    double emitSc = 0;
    switch (state) {
    case Match:
      emitSc = matchEmitScore(i,j);
      xRow.push_front (px->seq[--i]);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat(i,j) + qs.m2m + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + qs.i2m + emitSc, Insert);
      updateMax (srcSc, state, del(i,j) + qs.d2m + emitSc, Delete);
      if (j == 0 && (i == 0 || pconfig->local))
	updateMax (srcSc, state, emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error");
      break;

    case Insert:
      emitSc = cachedInsertEmitScore[j];
      xRow.push_front (Alignment::gapChar);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat(i,j) + qs.m2i + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + qs.i2i + emitSc, Insert);
      Assert (srcSc == ins(i,j+1), "Traceback error");
      break;

    case Delete:
      xRow.push_front (px->seq[--i]);
      yRow.push_front (Alignment::gapChar);
      updateMax (srcSc, state, mat(i,j) + qs.m2d, Match);
      updateMax (srcSc, state, del(i,j) + qs.d2d, Delete);
      Assert (srcSc == del(i+1,j), "Traceback error");
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  const SeqIdx xStart = i + 1;
  Alignment align(2);
  align.gappedSeq[0].name = "Ref";
  if (pconfig->local)
    align.gappedSeq[0].comment = "substr(" + px->name + "," + to_string(xStart) + ".." + to_string(xEnd) + ")";
  else
    align.gappedSeq[0].comment = px->name;
  align.gappedSeq[1].name = "Read";
  align.gappedSeq[1].comment = py->name;
  align.gappedSeq[0].seq = string (xRow.begin(), xRow.end());
  align.gappedSeq[1].seq = string (yRow.begin(), yRow.end());
  align.score = result;
  return align;
}

Alignment QuaffViterbiMatrix::scoreAdjustedAlignment (const QuaffNullParams& nullModel) const {
  Alignment a = alignment();
  const double nullLogLike = nullModel.logLikelihood (*py);
  if (LogThisAt(2))
    cerr << "Null model score: " << nullLogLike << endl;
  a.score -= nullLogLike;
  return a;
}

void QuaffParamCounts::addWeighted (const QuaffParamCounts& counts, double weight) {
  for (QualScore q = 0; q < FastSeq::qualScoreRange; ++q)
    for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
      insert[i].qualCount[q] += weight * counts.insert[i].qualCount[q];
      for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
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
  lp += logBetaPdf (qp.beginInsert, beginInsertYes, beginInsertNo);
  lp += logBetaPdf (qp.extendInsert, extendInsertYes, extendInsertNo);
  lp += logBetaPdf (qp.beginDelete, beginDeleteYes, beginDeleteNo);
  lp += logBetaPdf (qp.extendDelete, extendDeleteYes, extendDeleteNo);
  double *alpha = new double[dnaAlphabetSize], *theta = new double[dnaAlphabetSize];
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    lp += qp.insert[i].logQualProb (insert[i].qualCount);  // not normalized...
    theta[i] = qp.insert[i].symProb;
    alpha[i] = accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 1.);
  }
  lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j) {
      lp += qp.match[i][j].logQualProb (match[i][j].qualCount);  // not normalized...
      theta[j] = qp.match[i][j].symProb;
      alpha[j] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 1.);
    }
    lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
  }
  delete[] theta;
  delete[] alpha;
  return lp;
}

double QuaffParamCounts::expectedLogLike (const QuaffParams& qp) const {
  double ll = 0;
  ll += log (qp.beginInsert) * beginInsertYes + log (1 - qp.beginInsert) * beginInsertNo;
  ll += log (qp.extendInsert) * extendInsertYes + log (1 - qp.extendInsert) * extendInsertNo;
  ll += log (qp.beginDelete) * beginDeleteYes + log (1 - qp.beginDelete) * beginDeleteNo;
  ll += log (qp.extendDelete) * extendDeleteYes + log (1 - qp.extendDelete) * extendDeleteNo;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    ll += qp.insert[i].logQualProb (insert[i].qualCount);
    ll += log (qp.insert[i].symProb) * accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.);
  }
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j) {
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
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insFreq.push_back (accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.));
  const double insNorm = accumulate (insFreq.begin(), insFreq.end(), 0.);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    qp.insert[i].symProb = insFreq[i] / insNorm;
    fitNegativeBinomial (insert[i].qualCount, qp.insert[i].qualTrialSuccessProb, qp.insert[i].qualNumSuccessfulTrials);
  }

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    vguard<double> iMatFreq (dnaAlphabetSize);
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j)
      iMatFreq[j] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
    const double iMatNorm = accumulate (iMatFreq.begin(), iMatFreq.end(), 0.);
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j) {
      qp.match[i][j].symProb = iMatFreq[j] / iMatNorm;
      fitNegativeBinomial (match[i][j].qualCount, qp.match[i][j].qualTrialSuccessProb, qp.match[i][j].qualNumSuccessfulTrials);
    }
  }

  return qp;
}

QuaffForwardBackwardMatrix::QuaffForwardBackwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : fwd (env, qp, config),
    back (fwd)
{
  if (LogWhen("postmatrix"))
    write (cerr);
}

double QuaffForwardBackwardMatrix::postMatch (SeqIdx i, SeqIdx j) const {
  return exp (fwd.mat(i,j) + back.mat(i,j) - fwd.result);
}

double QuaffForwardBackwardMatrix::postDelete (SeqIdx i, SeqIdx j) const {
  return exp (fwd.del(i,j) + back.del(i,j) - fwd.result);
}

double QuaffForwardBackwardMatrix::postInsert (SeqIdx i, SeqIdx j) const {
  return exp (fwd.ins(i,j) + back.ins(i,j) - fwd.result);
}

void QuaffForwardBackwardMatrix::write (ostream& out) const {
  for (SeqIdx j = 1; j <= fwd.yLen; ++j) {
    for (SeqIdx i : fwd.penv->forward_i(j))
      out << "i=" << i << ":" << fwd.px->seq[i-1] << " j=" << j << ":" << fwd.py->seq[j-1] << (fwd.py->hasQual() ? string(1,fwd.py->qual[j-1]) : string()) << "\tmat " << postMatch(i,j) << "\tins " << postInsert(i,j) << "\tdel " << postDelete(i,j) << endl;
    out << endl;
  }
}

QuaffNullParams::QuaffNullParams()
  : nullEmit (.5),
    null (dnaAlphabetSize)
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
    const vguard<AlphTok> tok = s.tokens (dnaAlphabet);
    for (SeqIdx i = 0; i < s.length(); ++i) {
      ++symCount[tok[i]];
      if (s.hasQual())
	++nullCount[tok[i]].qualCount[s.getQualScoreAt(i)];
    }
  }

  nullEmit = 1 / (1 + nullEmitNo / nullEmitYes);
  const double symCountNorm = accumulate (symCount.begin(), symCount.end(), 0.);
  for (AlphTok n = 0; n < dnaAlphabetSize; ++n) {
    null[n].symProb = symCount[n] / symCountNorm;
    fitNegativeBinomial (nullCount[n].qualCount, null[n].qualTrialSuccessProb, null[n].qualNumSuccessfulTrials);
  }

  if (LogThisAt(3))
    write (cerr << "Null model:" << endl);
}

void QuaffNullParams::read (istream& in) {
  map<string,string> val = readParamFile (in);

  QuaffParamRead(nullEmit);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    null[i].read (val, string("null") + dnaAlphabet[i]);
}

double QuaffNullParams::logLikelihood (const vguard<FastSeq>& seqs) const {
  double ll = 0;
  for (const auto& s : seqs)
    ll += logLikelihood (s);
  return ll;
}

double QuaffNullParams::logLikelihood (const FastSeq& s) const {
  double ll = s.length() * log(nullEmit) + log(1. - nullEmit);
  const vguard<AlphTok> tok = s.tokens (dnaAlphabet);
    for (SeqIdx i = 0; i < s.length(); ++i) {
      ll += log (null[tok[i]].symProb);
      if (s.hasQual())
	ll += null[tok[i]].logQualProb (s.getQualScoreAt(i));
#if defined(NAN_DEBUG)
      if (isnan(ll)) {
	cerr << "NaN error in QuaffNullParams::logLikelihood" << endl;
	throw;
      }
#endif
    }
  return ll;
}

void QuaffNullParams::write (ostream& out) const {
  QuaffParamWrite(nullEmit);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    null[i].write (out, string("null") + dnaAlphabet[i]);
}

QuaffTrainer::QuaffTrainer()
  : maxIterations (100),
    minFractionalLoglikeIncrement (.01),
    allowNullModel (true)
{ }

bool QuaffTrainer::parseTrainingArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-maxiter") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      maxIterations = atoi (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-mininc") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      minFractionalLoglikeIncrement = atof (val);
      argv += 2;
      argc -= 2;
      return true;
    
    } else if (arg == "-force") {
      allowNullModel = false;
      argv += 1;
      argc -= 1;
      return true;

    } else if (arg == "-counts") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      rawCountsFilename = argv[1];
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-countswithprior") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      countsWithPriorFilename = argv[1];
      argv += 2;
      argc -= 2;
      return true;
}
  }

  return false;
}

QuaffParams QuaffTrainer::fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, const QuaffDPConfig& config) {
  QuaffParams qp = seed;
  QuaffParamCounts counts, countsWithPrior;
  double prevLogLikeWithPrior = -numeric_limits<double>::infinity();
  vguard<size_t> initialSortOrder (x.size());
  iota (initialSortOrder.begin(), initialSortOrder.end(), (size_t) 0);
  vguard<vguard<size_t> > sortOrder (y.size(), initialSortOrder);
  for (int iter = 0; iter < maxIterations; ++iter) {
    counts.zeroCounts();
    double logLike = 0;
    for (size_t ny = 0; ny < y.size(); ++ny) {
      const auto& yfs = y[ny];
      const double yNullLogLike = allowNullModel ? nullModel.logLikelihood(yfs) : -numeric_limits<double>::infinity();
      double yLogLike = yNullLogLike;  // this initial value allows null model to "win"
      if (LogThisAt(2))
	cerr << "Null model score for " << yfs.name << " is " << yNullLogLike << endl;
      vguard<double> xyLogLike;
      vguard<QuaffParamCounts> xyCounts;
      for (auto nx : sortOrder[ny]) {
	const auto& xfs = x[nx];
	DiagonalEnvelope env = config.makeEnvelope (xfs, yfs);
	const QuaffForwardMatrix fwd (env, qp, config);
	const double ll = fwd.result;
	QuaffParamCounts qpc;
	if (ll >= yLogLike - MAX_TRAINING_LOG_DELTA) {  // don't waste time computing low-weight counts
	  const QuaffBackwardMatrix back (fwd);
	  qpc = back.qc;
	}
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
      const vguard<size_t> ascendingOrder = orderedIndices (xyLogLike);
      sortOrder[ny] = vguard<size_t> (ascendingOrder.rbegin(), ascendingOrder.rend());
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
    countsWithPrior = counts;
    countsWithPrior.addWeighted (pseudocounts, 1.);
    const double oldExpectedLogLike = countsWithPrior.expectedLogLike (qp);
    qp = countsWithPrior.fit();
    const double newExpectedLogLike = countsWithPrior.expectedLogLike (qp);
    if (LogThisAt(2))
      cerr << "Expected log-likelihood went from " << oldExpectedLogLike << " to " << newExpectedLogLike << " during M-step" << endl;
  }
  if (rawCountsFilename.size()) {
    ofstream out (rawCountsFilename);
    counts.write (out);
  }
  if (countsWithPriorFilename.size()) {
    ofstream out (countsWithPriorFilename);
    countsWithPrior.write (out);
  }
  qp.fitRefSeqs(x);
  return qp;
}

QuaffAlignmentPrinter::QuaffAlignmentPrinter()
  : format (StockholmAlignment),
    logOddsThreshold (0)
{ }

bool QuaffAlignmentPrinter::parseAlignmentPrinterArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-format") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
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

    } else if (arg == "-threshold") {
      Assert (argc > 1, "%s must have an argument", arg.c_str());
      const char* val = argv[1];
      logOddsThreshold = atof (val);
      argv += 2;
      argc -= 2;
      return true;

    } else if (arg == "-nothreshold") {
      logOddsThreshold = -numeric_limits<double>::infinity();
      argv += 1;
      argc -= 1;
      return true;

    }
  }

  return false;
}

void QuaffAlignmentPrinter::writeAlignment (ostream& out, const Alignment& align) const {
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

QuaffAligner::QuaffAligner()
  : QuaffAlignmentPrinter(),
    printAllAlignments (false)
{ }

bool QuaffAligner::parseAlignmentArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-printall") {
      printAllAlignments = true;
      argv += 1;
      argc -= 1;
      return true;
    }
  }

  return parseAlignmentPrinterArgs (argc, argv);
}

void QuaffAligner::align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config) {
  for (const auto& yfs : y) {
    size_t nBestAlign = 0;
    vguard<Alignment> xyAlign;
    for (const auto& xfs : x) {
      DiagonalEnvelope env = config.makeEnvelope (xfs, yfs);
      const QuaffViterbiMatrix viterbi (env, params, config);
      if (viterbi.resultIsFinite()) {
	const Alignment align = viterbi.scoreAdjustedAlignment(nullModel);
	if (xyAlign.empty() || align.score > xyAlign[nBestAlign].score)
	  nBestAlign = xyAlign.size();
	xyAlign.push_back (align);
      }
    }
    if (printAllAlignments)
      for (const auto& a : xyAlign)
	writeAlignment (out, a);
    else if (xyAlign.size())
      writeAlignment (out, xyAlign[nBestAlign]);
  }
}
