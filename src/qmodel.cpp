#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <random>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include "qmodel.h"
#include "logsumexp.h"
#include "negbinom.h"
#include "memsize.h"
#include "aws.h"
#include "regexmacros.h"

// internal #defines
// Forward-backward tolerance
const double MAX_FRACTIONAL_FWDBACK_ERROR = .0001;

// Threshold for dropping poorly-matching refseqs
const double MAX_TRAINING_LOG_DELTA = 20;

const regex paramValRegex (" *" RE_VARNAME_GROUP " *: *" RE_GROUP(RE_NONWHITE_CHAR_CLASS RE_DOT_STAR), regex_constants::basic);
const regex lineRegex (RE_DOT_GROUP, regex_constants::basic);
const regex orderRegex (RE_NUMERIC_GROUP, regex_constants::basic);
const regex countRegex (RE_NUMERIC_GROUP " *: *" RE_FLOAT_GROUP, regex_constants::basic);
const regex remoteUserAddrRegex (RE_DOT_GROUP "@" RE_DNS_GROUP RE_DOT_STAR, regex_constants::basic);
const regex remoteAddrRegex (RE_DNS_GROUP RE_DOT_STAR, regex_constants::basic);
const regex singleRemotePortRegex (RE_DOT_STAR ":" RE_NUMERIC_GROUP, regex_constants::basic);
const regex multiRemotePortRegex (RE_DOT_STAR ":" RE_NUMERIC_GROUP "-" RE_NUMERIC_GROUP, regex_constants::basic);

// useful helper methods
double logBetaPdf (double prob, double yesCount, double noCount);
void readQuaffParamFileLine (map<string,string>& paramVal, const string& line);
vguard<size_t> readSortOrder (const string& orderString);

// main method bodies
double logBetaPdf (double prob, double yesCount, double noCount) {
  return log (gsl_ran_beta_pdf (prob, yesCount + 1, noCount + 1));
}

void readQuaffParamFileLine (map<string,string>& val, const string& line) {
  smatch sm;
  if (regex_match (line, sm, paramValRegex))
    val[sm.str(1)] = sm.str(2);
}

string readQuaffStringFromSocket (TCPSocket* sock, int bufSize) {
  return JsonUtil::readStringFromSocket (sock, eofRegex, bufSize);
}

map<string,string> readQuaffParamFile (istream& in) {
  map<string,string> val;
  while (!in.eof()) {
    string line;
    getline(in,line);
    smatch sm;
    if (regex_search (line, sm, lineRegex))  // chomp off newline
      readQuaffParamFileLine (val, sm.str(1));
  }
  return val;
}

map<string,string> readQuaffParamFile (const string& s) {
  string msg = s;
  map<string,string> val;
  smatch sm;
  while (regex_search (msg, sm, lineRegex)) {
    readQuaffParamFileLine (val, sm.str(1));
    msg = sm.suffix().str();
  }
  return val;
}

map<string,string> readQuaffParamFile (TCPSocket* sock) {
  string msg = readQuaffStringFromSocket (sock);
  return readQuaffParamFile (msg);
}

vguard<size_t> readSortOrder (const string& orderString) {
  string s = orderString;
  smatch sm;
  vguard<size_t> sortOrder;
  while (regex_search (s, sm, orderRegex)) {
    sortOrder.push_back (atoi (sm.str(1).c_str()));
    s = sm.suffix();
  }
  return sortOrder;
}

void randomDelayBeforeRetry (unsigned int minSeconds, unsigned int maxSeconds) {
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist (minSeconds, maxSeconds);

  const int delay = dist(gen);
    
  if (LogThisAt(3))
    logger << "Retrying in " << plural(delay,"second") << "..." << endl;

  usleep (delay * 1000000);
}

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

ostream& SymQualDist::writeJson (ostream& out) const {
  return
    out << "{ \"p\": " << qualTrialSuccessProb
	<< ", \"r\": " << qualNumSuccessfulTrials
	<< ", \"m\": " << negativeBinomialMean(qualTrialSuccessProb,qualNumSuccessfulTrials)
	<< ", \"sd\": " << sqrt(negativeBinomialVariance(qualTrialSuccessProb,qualNumSuccessfulTrials))
	<< " }";
}

bool SymQualDist::readJson (const JsonValue& val) {
  const JsonMap jm (val);
  Desire (jm.contains("p") && jm.contains("r"), "Missing negative binomial parameters");
  qualTrialSuccessProb = jm["p"].toNumber();
  qualNumSuccessfulTrials = jm["r"].toNumber();
  return true;
}

bool SymQualDist::read (map<string,string>& paramVal, const string& prefix) {
  const string qp = prefix + "qp", qr = prefix + "qr";
  Desire (paramVal.find(prefix) != paramVal.end(), "Missing parameter: %s", prefix.c_str());
  Desire (paramVal.find(qp) != paramVal.end(), "Missing parameter: %s", qp.c_str());
  Desire (paramVal.find(qr) != paramVal.end(), "Missing parameter: %s", qr.c_str());
  symProb = atof (paramVal[prefix].c_str());
  qualTrialSuccessProb = atof (paramVal[qp].c_str());
  qualNumSuccessfulTrials = atof (paramVal[qr].c_str());
  return true;
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

ostream& SymQualCounts::writeJson (ostream& out) const {
  return out << "[ " << join (qualCount, ", ") << " ]";
}

bool SymQualCounts::readJson (const JsonValue& value) {
  Desire (value.getTag() == JSON_ARRAY, "Not an array");
  qualCount = JsonUtil::doubleVec (value);
  Desire (qualCount.size() == FastSeq::qualScoreRange, "FASTQ counts array is wrong size");
  return true;
}

bool SymQualCounts::read (map<string,string>& paramVal, const string& param) {
  Desire (paramVal.find(param) != paramVal.end(), "Couldn't read %s", param.c_str());
  string c = paramVal[param];
  qualCount = vguard<double> (FastSeq::qualScoreRange, 0.);
  smatch sm;
  while (regex_search (c, sm, countRegex)) {
    qualCount[atoi (sm.str(1).c_str())] = atof (sm.str(2).c_str());
    c = sm.suffix().str();
  }
  return true;
}

QuaffKmerContext::QuaffKmerContext (const char* prefix, unsigned int kmerLen)
  : prefix(prefix),
    defaultKmerLen (kmerLen)
{
  initKmerContext (defaultKmerLen);
}

void QuaffKmerContext::initKmerContext (unsigned int newKmerLen) {
  kmerLen = newKmerLen;
  numKmers = numberOfKmers(newKmerLen,dnaAlphabetSize);
}

void QuaffKmerContext::readKmerLen (map<string,string>& paramVal) {
  const string tag = string(prefix) + "Order";
  if (paramVal.find(tag) == paramVal.end())
    initKmerContext (defaultKmerLen);
  else
    initKmerContext (atoi (paramVal[tag].c_str()));
}

void QuaffKmerContext::readJsonKmerLen (const JsonMap& m) {
  const string tag = string(prefix) + "Order";
  if (m.contains(tag))
    initKmerContext ((int) m[tag].toNumber());
  else
    initKmerContext (defaultKmerLen);
}

void QuaffKmerContext::writeKmerLen (ostream& out) const {
  if (kmerLen != defaultKmerLen)
    out << prefix << "Order: " << kmerLen << endl;
}

void QuaffKmerContext::writeJsonKmerLen (ostream& out) const {
  if (kmerLen != defaultKmerLen)
    out << prefix << "  \"Order\": " << kmerLen << ',' << endl;
}

string QuaffKmerContext::kmerString (Kmer kmer) const {
  return kmerToString (kmer, kmerLen, dnaAlphabet);
}

string QuaffKmerContext::insertParamName (AlphTok i) const {
  return string("insert") + dnaAlphabet[i];
}

string QuaffMatchKmerContext::matchParamName (AlphTok i, Kmer j) const {
  string name = "match";
  if (kmerLen > 1)
    name += kmerPrefix(j) + "_";
  return name + dnaAlphabet[i] + kmerSuffix(j);
}

string QuaffMatchKmerContext::kmerPrefix (Kmer j) const {
  const string ks = kmerString(j);
  return ks.substr(0,kmerLen-1);
}

char QuaffMatchKmerContext::kmerSuffix (Kmer j) const {
  const string ks = kmerString(j);
  return ks.back();
}

string QuaffIndelKmerContext::booleanParamName (const char* tag, Kmer j) const {
  string name = string(tag);
  if (kmerLen > 0)
    name = name + kmerString(j);
  return name;
}

QuaffParams::QuaffParams (unsigned int matchKmerLen, unsigned int indelKmerLen)
  : matchContext(matchKmerLen),
    indelContext(indelKmerLen),
    refBase(dnaAlphabetSize,.25),
    extendInsert(.5),
    extendDelete(.5),
    insert (dnaAlphabetSize)
{
  resize();
}

void QuaffParams::resize() {
  match = vguard<vguard<SymQualDist> > (dnaAlphabetSize, vguard<SymQualDist> (matchContext.numKmers));
  beginInsert = vguard<double> (indelContext.numKmers, .5);
  beginDelete = vguard<double> (indelContext.numKmers, .5);
}

#define QuaffParamWrite(X) out << #X ": " << X << endl
#define QuaffParamWriteK(X,KMER) out << indelContext.booleanParamName(#X,KMER) << ": " << X[KMER] << endl
void QuaffParams::write (ostream& out) const {
  matchContext.writeKmerLen (out);
  indelContext.writeKmerLen (out);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    out << "refBase" << dnaAlphabet[i] << ": " << refBase[i] << endl;
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    QuaffParamWriteK(beginInsert,j);
    QuaffParamWriteK(beginDelete,j);
  }
  QuaffParamWrite(extendInsert);
  QuaffParamWrite(extendDelete);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, matchContext.insertParamName(i));
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (Kmer j = 0; j < matchContext.numKmers; ++j)
      match[i][j].write (out, matchContext.matchParamName(i,j));
}

#define QuaffParamWriteJson(X) out << "  \"" << #X "\": " << X
#define QuaffParamWriteJsonKmer(X,KMER) out << (KMER == 0 ? "" : ",") << " \"" << indelContext.kmerString(KMER) << "\": " << X[KMER]
#define QuaffParamWriteJsonKmers(X) out << "  \"" << #X << "\": {"; for (Kmer j = 0; j < indelContext.numKmers; ++j) QuaffParamWriteJsonKmer(X,j); out << " }"
ostream& QuaffParams::writeJson (ostream& out) const {
  out << '{' << endl;
  matchContext.writeJsonKmerLen (out);
  indelContext.writeJsonKmerLen (out);
  out << "  \"refBase\": {";
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    out << " \"" << dnaAlphabet[i] << "\": " << refBase[i]
	<< (i==dnaAlphabetSize-1 ? " },\n" : ",");
  QuaffParamWriteJsonKmers(beginInsert) << "," << endl;
  QuaffParamWriteJsonKmers(beginDelete) << "," << endl;
  QuaffParamWriteJson(extendInsert) << "," << endl;
  QuaffParamWriteJson(extendDelete) << "," << endl;
  out << "  \"insert\": {" << endl;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].writeJson(out << "    \"" << dnaAlphabet[i] << "\": ")
      << (i == dnaAlphabetSize-1 ? " }," : ",")
      << endl;
  out << "  \"match\": {" << endl;
  for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize) {
    out << "    \"" << matchContext.kmerPrefix(jPrefix) << "\": {" << endl;
    for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
      out << "    \"" << dnaAlphabet[i] << "\": {" << endl;
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix)
	match[i][jPrefix+jSuffix].writeJson(out << "      \"" << dnaAlphabet[jSuffix] << "\": ")
	  << (jSuffix == dnaAlphabetSize-1 ? " }" : ",\n");
      out << (i == dnaAlphabetSize-1 ? " }" : ",\n");
    }
    out << (jPrefix == matchContext.numKmers - dnaAlphabetSize ? " }" : ",") << endl;
  }
  out << "}" << endl;
  return out;
}

#define QuaffParamReadJson(X) do { Desire(jm.contains(#X),"Missing parameter: \"" #X "\""); X = jm[#X].toNumber(); } while(0)
#define QuaffParamReadJsonKmer(X,JMX,KMER) do { const string tmpParamName = indelContext.kmerString(KMER); Desire(JMX.contains(tmpParamName),"Missing parameter: \"%s\".\"%s\"",#X,tmpParamName.c_str()); X[KMER] = JMX[tmpParamName].toNumber(); } while(0)
#define QuaffParamReadJsonKmers(X) do { Desire(jm.contains(#X),"Missing parameter: \"" #X "\""); JsonMap jmx (jm[#X]); for (Kmer j = 0; j < indelContext.numKmers; ++j) QuaffParamReadJsonKmer(X,jmx,j); } while(0)

bool QuaffParams::readJson (const JsonValue& val) {
  const JsonMap jm (val);
  return readJson (jm);
}

bool QuaffParams::readJson (const JsonMap& jm) {
  matchContext.readJsonKmerLen (jm);
  indelContext.readJsonKmerLen (jm);
  resize();
  
  QuaffParamReadJsonKmers (beginInsert);
  QuaffParamReadJsonKmers (beginDelete);
  QuaffParamReadJson (extendInsert);
  QuaffParamReadJson (extendDelete);

  Desire (jm.contains ("insert"), "Missing parameter: \"insert\"");
  const JsonValue& ins = jm["insert"];
  const JsonMap& jmIns (ins);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    const string iKey (1, dnaAlphabet[i]);
    Desire (jmIns.contains(iKey), "Missing parameter: \"insert\".\"%s\"", iKey.c_str());
    Desire (insert[i].readJson (jmIns[iKey]), "Couldn't read \"insert\".\"%s\"", iKey.c_str());
  }

  Desire (jm.contains ("match"), "Missing parameter: \"match\"");
  const JsonValue& mat = jm["match"];
  const JsonMap jmMat (mat);
  for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize) {
    const string jPrefixKey (matchContext.kmerPrefix(jPrefix));
    Desire (jmMat.contains (jPrefixKey), "Missing parameter: \"match\".\"%s\"", jPrefixKey.c_str());
    const JsonValue& matJ = jmMat[jPrefixKey];
    const JsonMap jmMatJ (matJ);
    for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
      const string iKey (1, dnaAlphabet[i]);
      Desire (jmMatJ.contains (iKey), "Missing parameter: \"match\".\"%s\".\"%s\"", jPrefixKey.c_str(), iKey.c_str());
      const JsonValue& matJI = jmMatJ[iKey];
      const JsonMap jmMatJI (matJI);
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	const string jSuffixKey (1, dnaAlphabet[jSuffix]);
	Desire (jmMatJI.contains (jSuffixKey), "Missing parameter: \"match\".\"%s\".\"%s\".\"%s\"", jPrefixKey.c_str(), iKey.c_str(), jSuffixKey.c_str());
	Desire (match[i][jPrefix+jSuffix].readJson (jmMatJI[jSuffixKey]), "Couldn't read \"match\".\"%s\".\"%s\".\"%s\"", jPrefixKey.c_str(), iKey.c_str(), jSuffixKey.c_str());
      }
    }
  }

  return true;
}

void QuaffParams::readJson (istream& in) {
  ParsedJson pj (in);
  Require (readJson(pj), "Couldn't read parameters");
}

void QuaffParams::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
}

#define QuaffParamRead(X) do { Desire(val.find(#X) != val.end(),"Missing parameter: " #X); X = atof(val[#X].c_str()); } while(0)
#define QuaffParamReadK(X,KMER) do { const string tmpParamName = indelContext.booleanParamName(#X,KMER); Desire(val.find(tmpParamName) != val.end(),"Missing parameter: %s",tmpParamName.c_str()); X[KMER] = atof(val[tmpParamName].c_str()); } while(0)

void QuaffParams::read (istream& in) {
  map<string,string> val = readQuaffParamFile (in);
  Require (read(val), "Couldn't read parameters");
}

bool QuaffParams::read (map<string,string>& val) {
  matchContext.readKmerLen(val);
  indelContext.readKmerLen(val);
  resize();
  
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    QuaffParamReadK(beginInsert,j);
    QuaffParamReadK(beginDelete,j);
  }
  QuaffParamRead(extendInsert);
  QuaffParamRead(extendDelete);

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    Desire (insert[i].read (val, matchContext.insertParamName(i)), "");
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (Kmer j = 0; j < matchContext.numKmers; ++j)
      Desire (match[i][j].read (val, matchContext.matchParamName(i,j)), "");
  return true;
}

void QuaffParams::fitRefSeqs (const vguard<FastSeq>& refs) {
  int totalLen;
  vguard<int> baseCount (dnaAlphabetSize, 0.);
  for (const auto& fs : refs) {
    totalLen += fs.length();
    for (auto i : fs.tokens(dnaAlphabet))
      ++baseCount[i];
  }
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    refBase[i] = baseCount[i] / (double) totalLen;
}

QuaffScores::QuaffScores (const QuaffParams& qp)
  : matchContext(qp.matchContext.kmerLen),
    indelContext(qp.indelContext.kmerLen),
    pqp(&qp),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualScores> (matchContext.numKmers)),
    m2m (indelContext.numKmers),
    m2i (indelContext.numKmers),
    m2d (indelContext.numKmers),
    m2e (indelContext.numKmers)    
{
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    insert[i] = SymQualScores (qp.insert[i]);
    for (Kmer j = 0; j < matchContext.numKmers; ++j)
      match[i][j] = SymQualScores (qp.match[i][j]);
  }

  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    m2m[j] = log(1-qp.beginInsert[j]) + log(1-qp.beginDelete[j]);
    m2i[j] = log(qp.beginInsert[j]);
    m2d[j] = log(1-qp.beginInsert[j]) + log(qp.beginDelete[j]);
    m2e[j] = log(qp.beginInsert[j]);
  }

  d2d = log(qp.extendDelete);
  d2m = log(1-qp.extendDelete);

  i2i = log(qp.extendInsert);
  i2m = log(1-qp.extendInsert);
}

QuaffEmitCounts::QuaffEmitCounts (unsigned int matchKmerLen, unsigned int indelKmerLen)
  : matchContext(matchKmerLen),
    indelContext(indelKmerLen),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualCounts> (matchContext.numKmers))
{ }

QuaffEmitCounts::QuaffEmitCounts (const QuaffEmitCounts& c)
  : matchContext(c.matchContext.kmerLen),
    indelContext(c.indelContext.kmerLen),
    insert (dnaAlphabetSize),
    match (dnaAlphabetSize, vguard<SymQualCounts> (matchContext.numKmers))
{ }

void QuaffEmitCounts::write (ostream& out) const {
  matchContext.writeKmerLen (out);
  indelContext.writeKmerLen (out);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].write (out, matchContext.insertParamName(i));
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (Kmer j = 0; j < matchContext.numKmers; ++j)
      match[i][j].write (out, matchContext.matchParamName(i,j));
}

ostream& QuaffEmitCounts::writeJson (ostream& out) const {
  matchContext.writeJsonKmerLen (out);
  indelContext.writeJsonKmerLen (out);
  out << "  \"insert\": {" << endl;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    insert[i].writeJson(out << "    \"" << dnaAlphabet[i] << "\": ")
      << (i == dnaAlphabetSize-1 ? " }," : ",")
      << endl;
  out << "  \"match\": {" << endl;
  for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize) {
    out << "    \"" << matchContext.kmerPrefix(jPrefix) << "\": {" << endl;
    for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
      out << "    \"" << dnaAlphabet[i] << "\": {" << endl;
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix)
	match[i][jPrefix+jSuffix].writeJson(out << "      \"" << dnaAlphabet[jSuffix] << "\": ")
	  << (jSuffix == dnaAlphabetSize-1 ? " }" : ",\n");
      out << (i == dnaAlphabetSize-1 ? " }" : ",\n");
    }
    out << (jPrefix == matchContext.numKmers - dnaAlphabetSize ? " }" : ",") << endl;
  }
  out << "}" << endl;
  return out;
}


void QuaffEmitCounts::resize() {
  match = vguard<vguard<SymQualCounts> > (dnaAlphabetSize, vguard<SymQualCounts> (matchContext.numKmers));
}

QuaffCounts::QuaffCounts (unsigned int matchKmerLen, unsigned int indelKmerLen)
  : QuaffEmitCounts(matchKmerLen,indelKmerLen),
    m2m(indelContext.numKmers,0),
    m2i(indelContext.numKmers,0),
    m2d(indelContext.numKmers,0),
    m2e(indelContext.numKmers,0),
    d2d(0),
    d2m(0),
    i2i(0),
    i2m(0)
{ }

void QuaffCounts::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
}

void QuaffCounts::write (ostream& out) const {
  QuaffEmitCounts::write (out);
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    QuaffParamWriteK(m2m,j);
    QuaffParamWriteK(m2i,j);
    QuaffParamWriteK(m2d,j);
    QuaffParamWriteK(m2e,j);
  }
  QuaffParamWrite(d2d);
  QuaffParamWrite(d2m);
  QuaffParamWrite(i2i);
  QuaffParamWrite(i2m);
}

ostream& QuaffCounts::writeJson (ostream& out) const {
  out << "{" << endl;
  QuaffEmitCounts::writeJson (out);
  QuaffParamWriteJsonKmers(m2m) << "," << endl;
  QuaffParamWriteJsonKmers(m2i) << "," << endl;
  QuaffParamWriteJsonKmers(m2d) << "," << endl;
  QuaffParamWriteJsonKmers(m2e) << "," << endl;
  QuaffParamWriteJson(d2d) << "," << endl;
  QuaffParamWriteJson(d2m) << "," << endl;
  QuaffParamWriteJson(i2i) << "," << endl;
  QuaffParamWriteJson(i2m) << endl;
  out << "}" << endl;
  return out;
}

QuaffParamCounts::QuaffParamCounts (unsigned int matchKmerLen, unsigned int indelKmerLen)
  : QuaffEmitCounts(matchKmerLen,indelKmerLen)
{
  resize();
  zeroCounts();
}

QuaffParamCounts::QuaffParamCounts (const QuaffCounts& counts)
  : QuaffEmitCounts(counts),
    beginInsertNo (vectorSum (counts.m2m, counts.m2d)),
    beginInsertYes (vectorSum (counts.m2i, counts.m2e)),
    extendInsertNo (counts.i2m),
    extendInsertYes (counts.i2i),
    beginDeleteNo (counts.m2m),
    beginDeleteYes (counts.m2d),
    extendDeleteNo (counts.d2m),
    extendDeleteYes (counts.d2d)
{ }

void QuaffParamCounts::resize() {
  QuaffEmitCounts::resize();
  beginInsertNo = vguard<double> (indelContext.numKmers, 0);
  beginInsertYes = vguard<double> (indelContext.numKmers, 0);
  beginDeleteNo = vguard<double> (indelContext.numKmers, 0);
  beginDeleteYes = vguard<double> (indelContext.numKmers, 0);
}

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
    for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize)
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	const Kmer j = jPrefix + jSuffix;
	for (QualScore k = 0; k < FastSeq::qualScoreRange; ++k)
	  if (nullModel)
	    match[i][j].qualCount[k] = (i == j ? matchIdentCount : (otherCount * nullModel->null[jSuffix].symProb * dnaAlphabetSize / (1 - nullModel->null[jSuffix].symProb))) * gsl_ran_negative_binomial_pdf (k, nullModel->null[jSuffix].qualTrialSuccessProb, nullModel->null[jSuffix].qualNumSuccessfulTrials);
	  else
	    match[i][j].qualCount[k] = (i == j ? matchIdentCount : otherCount) / FastSeq::qualScoreRange;
      }
  beginInsertNo = vguard<double> (indelContext.numKmers, noBeginCount);
  beginInsertYes = vguard<double> (indelContext.numKmers, otherCount);
  extendInsertNo = otherCount;
  extendInsertYes = yesExtendCount;
  beginDeleteNo = vguard<double> (indelContext.numKmers, noBeginCount);
  beginDeleteYes = vguard<double> (indelContext.numKmers, otherCount);
  extendDeleteNo = otherCount;
  extendDeleteYes = yesExtendCount;
}

void QuaffParamCounts::write (ostream& out) const {
  QuaffEmitCounts::write (out);
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    QuaffParamWriteK(beginInsertNo,j);
    QuaffParamWriteK(beginInsertYes,j);
    QuaffParamWriteK(beginDeleteNo,j);
    QuaffParamWriteK(beginDeleteYes,j);
  }
  QuaffParamWrite(extendInsertNo);
  QuaffParamWrite(extendInsertYes);
  QuaffParamWrite(extendDeleteNo);
  QuaffParamWrite(extendDeleteYes);
}

ostream& QuaffParamCounts::writeJson (ostream& out) const {
  out << "{" << endl;
  QuaffEmitCounts::writeJson (out);
  QuaffParamWriteJsonKmers(beginInsertNo) << "," << endl;
  QuaffParamWriteJsonKmers(beginInsertYes) << "," << endl;
  QuaffParamWriteJsonKmers(beginDeleteNo) << "," << endl;
  QuaffParamWriteJsonKmers(beginDeleteYes) << "," << endl;
  QuaffParamWriteJson(extendInsertNo) << "," << endl;
  QuaffParamWriteJson(extendInsertYes) << "," << endl;
  QuaffParamWriteJson(extendDeleteNo) << "," << endl;
  QuaffParamWriteJson(extendDeleteYes) << "," << endl;
  out << "}" << endl;
  return out;
}

void QuaffParamCounts::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
}

void QuaffParamCounts::read (istream& in) {
  map<string,string> val = readQuaffParamFile (in);
  Require (read (val), "Couldn't read counts");
}

bool QuaffParamCounts::read (map<string,string>& val) {
  matchContext.readKmerLen(val);
  indelContext.readKmerLen(val);
  resize();
  
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    QuaffParamReadK(beginInsertYes,j);
    QuaffParamReadK(beginInsertNo,j);
    QuaffParamReadK(beginDeleteYes,j);
    QuaffParamReadK(beginDeleteNo,j);
  }

  QuaffParamRead(extendInsertYes);
  QuaffParamRead(extendDeleteYes);

  QuaffParamRead(extendInsertNo);
  QuaffParamRead(extendDeleteNo);

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    Desire (insert[i].read (val, matchContext.insertParamName(i)), "");
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (Kmer j = 0; j < matchContext.numKmers; ++j)
      Desire (match[i][j].read (val, matchContext.matchParamName(i,j)), "");
  return true;
}

const char Alignment::gapChar = '-';
const char Alignment::mismatchChar = ':';

void Alignment::writeGappedFasta (ostream& out) const {
  for (const auto& s : gappedSeq)
    s.writeFasta (out);
}

void Alignment::writeStockholm (ostream& out) const {
  vguard<string> rowName, rowData;
  vguard<size_t> rowIndex;
  for (const auto& s : gappedSeq) {
    rowIndex.push_back (rowName.size());
    rowName.push_back (s.name);
    rowData.push_back (s.seq);
    if (s.hasQual()) {
      rowName.push_back (string("#=GR ") + s.name + " QS");
      rowData.push_back (s.qual);
    }
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
    rowName.insert (rowName.begin() + rowIndex[1], string ("#=GC id"));
    rowData.insert (rowData.begin() + rowIndex[1], cons);

    if (gappedSeq[0].hasQual()) {
      swap (rowName[0], rowName[1]);
      swap (rowData[0], rowData[1]);
    }
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

void Alignment::writeSam (ostream& out) const {
  Assert (rows() == 2, "Tried to print %u-row alignment in SAM format, which is for pairwise alignments", rows());
  if (gappedSeq[0].source.rev)
    revcomp().writeSam(out);
  else {
    int flag = gappedSeq[1].source.rev ? 16 : 0;
    out << gappedSeq[1].source.name << '\t' << flag << '\t' << gappedSeq[0].source.name << '\t' << gappedSeq[0].source.start << "\t0\t" << cigarString() << "\t*\t0\t0\t*\t*\tAS:i:" << ((int) round(score)) << endl;
  }
}

void Alignment::writeSamHeader (ostream& out, const vguard<FastSeq>& seqs, const char* go_so) {
  out << "@HD\tVN:1.0\t" << go_so << endl;
  for (const auto& s : seqs)
    if (s.source.isNull())
      out << "@SQ\tSN:" << s.name << "\tLN:" << s.length() << endl;
}

string Alignment::cigarString() const {
  Assert (rows() == 2, "Tried to print %u-row alignment in CIGAR format, which is for pairwise alignments", rows());
  char lastCigarChar = 0;
  size_t count = 0;
  string cigar;
  for (SeqIdx col = 0; col < columns(); ++col) {
    const bool gap0 = isGapChar(gappedSeq[0].seq[col]), gap1 = isGapChar(gappedSeq[1].seq[col]);
    char cigarChar = 0;
    if (!gap0 && !gap1)
      cigarChar = 'M';
    else if (!gap0 && gap1)
      cigarChar = 'D';
    else if (gap0 && !gap1)
      cigarChar = 'I';
    if (cigarChar) {
      if (cigarChar == lastCigarChar)
	++count;
      else {
	if (count > 0)
	  cigar = cigar + lastCigarChar + to_string(count);
	lastCigarChar = cigarChar;
	count = 1;
      }
    }
  }
  if (count > 0)
    cigar = cigar + lastCigarChar + to_string(count);
  return cigar;
}

Alignment Alignment::revcomp() const {
  Alignment a (*this);
  for (size_t row = 0; row < rows(); ++row)
    a.gappedSeq[row] = a.gappedSeq[row].revcomp();
  return a;
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

RemoteServer::RemoteServer (string addr, unsigned int port)
  : addr(addr), port(port)
{ }

string RemoteServer::toString() const {
  return addr + ':' + to_string(port);
}

RemoteServerJob::RemoteServerJob (const string& user, const string& addr, unsigned int port, unsigned int threads)
  : user(user), addr(addr), port(port), threads(threads)
{ }

string RemoteServerJob::toString() const {
  return (user.size()
	  ? (user+'@')
	  : string())
    + addr + ':' + to_string(port)
    + (threads == 1
       ? string()
       : (string("-") + to_string(port+threads-1)));
}

bool QuaffDPConfig::parseServerConfigArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-port") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      serverPort = atoi (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

bool QuaffDPConfig::parseRefSeqConfigArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-global") {
      local = false;
      argvec.pop_front();
      return true;
    }
  }

  return parseGeneralConfigArgs (argvec);
}

bool QuaffDPConfig::parseGeneralConfigArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-kmatchband") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      bandSize = atoi (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-kmatch") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      kmerLen = atoi (val);
      Require (kmerLen >= 5 && kmerLen <= 32, "%s out of range (%d). Try 5 to 32", arg.c_str(), kmerLen);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-kmatchn") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      kmerThreshold = atoi (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-kmatchmb") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      maxSize = atoi(val) << 20;
      if (maxSize == 0) {
	maxSize = getMemorySize();
	Require (maxSize > 0, "Can't figure out available system memory; you will need to specify a size");
      }
      kmerThreshold = -1;
      autoMemSize = false;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-kmatchmax") {
      maxSize = getMemorySize();
      Require (maxSize > 0, "Can't figure out available system memory; you will need to specify a size");
      kmerThreshold = -1;
      autoMemSize = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-kmatchoff") {
      sparse = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-threads") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      threads = atoi (val);
      threadsSpecified = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-maxthreads") {
      threads = std::thread::hardware_concurrency();
      if (threads == 0) {
	Warn ("Can't detect number of cores; running in single-threaded mode");
	threads = 1;
      } else if (LogThisAt(2))
	logger << "Running in " << plural(threads,"thread") << endl;
      threadsSpecified = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-remote") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const string& remoteStr = argvec[1];
      smatch sm;
      string user, addr;
      unsigned int minPort, maxPort;
      if (regex_match (remoteStr, sm, remoteUserAddrRegex)) {
	user = sm.str(1);
	addr = sm.str(2);
      } else if (regex_match (remoteStr, sm, remoteAddrRegex))
	addr = sm.str(1);
      else
	Fail ("Can't parse server address: %s", remoteStr.c_str());
      if (regex_match (remoteStr, sm, singleRemotePortRegex))
	minPort = maxPort = atoi (sm.str(1).c_str());
      else if (regex_match (remoteStr, sm, multiRemotePortRegex)) {
	minPort = atoi (sm.str(1).c_str());
	maxPort = atoi (sm.str(2).c_str());
	Require (maxPort >= minPort, "Invalid port range (%u-%u)", minPort, maxPort);
      } else
	Fail ("Can't parse server port range: %s", remoteStr.c_str());
      addRemote (user, addr, minPort, maxPort + 1 - minPort);
      if (!threadsSpecified)
	threads = 0;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-sshpath") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      sshPath = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-sshkey") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      sshKey = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-remotepath") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      remoteQuaffPath = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-s3bucket") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      bucket = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2ami") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2Ami = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2type") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2Type = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2cores") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2Cores = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2user") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2User = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2port") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2Port = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2instances") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      ec2Instances = atoi (argvec[1].c_str());
      if (!threadsSpecified)
	threads = 0;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    }
  }

  return aws.parseAWSConfigArgs (argvec);
}

void QuaffDPConfig::setServerArgs (const char* serverType, const string& args) {
  remoteServerArgs = string(serverType)
    + args
    + logger.args()
    + (sparse
       ? (string(" -kmatchband ") + to_string(bandSize)
	  + " -kmatch " + to_string(kmerLen)
	  + (kmerThreshold >= 0
	     ? (string(" -kmatchn ") + to_string(kmerThreshold))
	     : (string(" -kmatchmb ") + to_string(maxSize >> 20))))
       : string(" -kmatchoff"))
    + (bucket.size() ? (string(" -s3bucket ") + bucket) : string());
  if (LogThisAt(7))
    logger << "Remote server arguments: " << remoteServerArgs << endl;
}

void QuaffDPConfig::addFileArg (const char* tag, const string& filename) {
  fileArgs.push_back (pair<string,string> (string(tag), filename));
}

string QuaffDPConfig::makeServerArgs() const {
  string s = remoteServerArgs;
  if (bucket.size())
    for (const auto& fa : fileArgs)
      s += ' ' + fa.first + ' ' + BucketStagingDir + '/' + AWS::basenameStr(fa.second);
  else
    for (const auto& fa : fileArgs)
      s += ' ' + fa.first + ' ' + fa.second;
  return s;
}

DiagonalEnvelope QuaffDPConfig::makeEnvelope (const FastSeq& x, const KmerIndex& yKmerIndex, size_t cellSize) const {
  DiagonalEnvelope env (x, yKmerIndex.seq);
  if (sparse)
    env.initSparse (yKmerIndex, bandSize, kmerThreshold, cellSize, effectiveMaxSize());
  else
    env.initFull();
  return env;
}

size_t QuaffDPConfig::effectiveMaxSize() const {
  return autoMemSize ? (maxSize / threads) : maxSize;
}

void QuaffDPConfig::syncFromBucket (const string& filename) const {
  if (bucket.size() && filename.size())
    aws.syncFromBucket (bucket, filename);
}

void QuaffDPConfig::syncToBucket (const string& filename) const {
  if (bucket.size() && filename.size())
    aws.syncToBucket (filename, bucket);
}

void QuaffDPConfig::addRemote (const string& user, const string& addr, unsigned int port, unsigned int threads) {
  remoteJobs.push_back (RemoteServerJob (user, addr, port, threads));
  for (unsigned int p = port; p < port + threads; ++p)
    remotes.push_back (RemoteServer (addr, p));
}

void QuaffDPConfig::startRemoteServers() {
  for (const auto& fa : fileArgs)
    syncToBucket (fa.second);
  if (ec2Instances > 0) {
    ec2InstanceIds = aws.launchInstancesWithScript (ec2Instances, ec2Type, ec2Ami, ec2StartupScript());
    ec2InstanceAddresses = aws.getInstanceAddresses (ec2InstanceIds);
    for (const auto& addr : ec2InstanceAddresses) {
      addRemote (ec2User, addr, ec2Port, ec2Cores);
      const string testCmd = makeSshCommand (string ("while ! test -e " BucketStagingDir "; do sleep 1; done"), remoteJobs.back());
      const bool bucketStagingDirExists = execWithRetries(testCmd,MaxQuaffSshAttempts);
      Assert (bucketStagingDirExists, "Cloud initialization failed");
    }
  }
  for (const auto& remoteJob : remoteJobs) {
    remoteServerThreads.push_back (thread (&startRemoteQuaffServer, this, &remoteJob));
    logger.nameLastThread (remoteServerThreads, "ssh");
  }
  if (remoteJobs.size()) {
    // give server threads time to connect...
    const int delay = 2;
    if (LogThisAt(5))
      logger << "Allowing server threads " << plural(delay,"second") << " to connect" << endl;
    usleep (delay * 1000000);
  }
}

string QuaffDPConfig::makeServerCommand (const RemoteServerJob& job) const {
  return aws.bashEnvPrefix() + remoteQuaffPath + " server " + makeServerArgs() + " -port " + to_string(job.port) + " -threads " + to_string(job.threads) + " 2>&1";
}

string QuaffDPConfig::makeSshCommand (const string& cmd, const RemoteServerJob& job) const {
  string sshCmd = sshPath + " -oBatchMode=yes -oUserKnownHostsFile=/dev/null -oStrictHostKeyChecking=no";
  if (job.user.size())
    sshCmd += " -l " + job.user;
  if (sshKey.size())
    sshCmd += " -i " + sshKey;
  sshCmd += ' ' + job.addr + " '" + cmd + "' 2>&1";
  return sshCmd;
}

string QuaffDPConfig::ec2StartupScript() const {
  return aws.bashBang
    + "yum -y update\n"
    + "yum -y install git\n"
    + "cd " DefaultQuaffInstallPrefix " && git clone " QuaffGitRepository " && cd quaff && make PREFIX=" DefaultQuaffInstallPrefix " -i aws-install\n"
    + "mkdir -p -m a=rwx " BucketStagingDir "\n";  // create this last, as it is used as a test for startup completion
}

void QuaffDPConfig::stopRemoteServers() {
  for (const auto& remote : remotes) {
    TCPSocket sock (remote.addr, remote.port);
    const string msg = string ("quit: 1\n") + SocketTerminatorString;
    sock.send (msg.c_str(), (int) msg.size());
  }
  for (auto& t : remoteServerThreads) {
    t.join();
    logger.eraseThreadName (t);
  }
  aws.terminateInstances (ec2InstanceIds);
}

// buffer size for popen
#define PIPE_BUF_SIZE 1024
bool QuaffDPConfig::execWithRetries (const string& cmd, int maxAttempts) const {
  int attempts = 0;
  while (true) {
  
    if (LogThisAt(4))
      logger << "Executing " << cmd << endl;

    FILE* pipe = popen (cmd.c_str(), "r");
    char line[PIPE_BUF_SIZE];
    deque<string> lastLines;
    const int linesToKeep = 3;
  
    while (fgets(line, PIPE_BUF_SIZE, pipe)) {
      lastLines.push_back (string (line));
      if (lastLines.size() > linesToKeep)
	lastLines.erase (lastLines.begin());
      if (LogThisAt(4))
	logger << line << flush;
    }

    const int status = pclose (pipe);
    if (status == EXIT_SUCCESS)
      break;

    Warn ("Command failed: %s\nExit code: %d\nLast output:\n%s", cmd.c_str(), status, join(lastLines).c_str());

    if (++attempts >= maxAttempts) {
      Warn ("Too many failed attempts to connect");
      return false;
    }
    
    randomDelayBeforeRetry (MinQuaffRetryDelay, MaxQuaffRetryDelay);
  }

  return true;
}

void startRemoteQuaffServer (const QuaffDPConfig* config, const RemoteServerJob* remoteJob) {
  const string sshCmd = config->makeSshCommand (config->makeServerCommand(*remoteJob), *remoteJob);
  const bool remoteServerStarted = config->execWithRetries (sshCmd, MaxQuaffSshAttempts);
  Assert (remoteServerStarted, "Failed to start remote server");
}

double QuaffDPMatrixContainer::dummy = -numeric_limits<double>::infinity();

QuaffDPMatrixContainer::QuaffDPMatrixContainer (const DiagonalEnvelope& env)
  : penv (&env),
    px (env.px),
    py (env.py),
    xLen (px->length()),
    yLen (py->length()),
    cell (env.totalStorageSize * 3, -numeric_limits<double>::infinity()),
    start (-numeric_limits<double>::infinity()),
    end (-numeric_limits<double>::infinity()),
    result (-numeric_limits<double>::infinity())
{ }

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

vguard<size_t> QuaffDPMatrixContainer::getFastSeqLengths (const vguard<FastSeq>& db) {
  vguard<size_t> len;
  for (const auto& s : db)
    len.push_back (s.length());
  return len;
}

QuaffDPMatrix::QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : QuaffDPMatrixContainer (env),
    pconfig (&config),
    qs (qp),
    xTok (px->tokens(dnaAlphabet)),
    yTok (py->tokens(dnaAlphabet)),
    yMatchKmer (py->kmers(dnaAlphabet,qp.matchContext.kmerLen)),
    yIndelKmer (py->kmers(dnaAlphabet,qp.indelContext.kmerLen)),
    yQual (py->qualScores()),
    cachedInsertEmitScore (py->length() + 1, -numeric_limits<double>::infinity())
{
  for (SeqIdx j = 1; j <= py->length(); ++j)
    cachedInsertEmitScore[j] = insertEmitScore(j);

  // pad yIndelKmer with an extra dummy entry at the beginning, to avoid having to test inside DP loop
  yIndelKmer.insert (yIndelKmer.begin(), 0);
}

void QuaffDPMatrix::write (ostream& out) const {
  for (SeqIdx j = 1; j <= yLen; ++j) {
    for (DiagonalEnvelope::iterator pi = penv->begin(j);
	 !pi.finished();
	 ++pi)
      out << "i=" << *pi << ":" << px->seq[*pi-1] << " j=" << j << ":" << py->seq[j-1] << (py->hasQual() ? string(1,py->qual[j-1]) : string()) << "\tmat " << mat(*pi,j) << "\tins " << ins(*pi,j) << "\tdel " << del(*pi,j) << endl;
    out << endl;
  }
  out << "result " << result << endl;
}

void QuaffDPMatrix::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
}

QuaffForwardMatrix::QuaffForwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : QuaffDPMatrix (env, qp, config)
{
  const FastSeq& x (*px);
  const FastSeq& y (*py);

  ProgressLog (plog, 4);
  plog.initProgress ("Forward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    plog.logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (DiagonalEnvelope::iterator pi = env.begin(j);
	 !pi.finished();
	 ++pi) {

      const SeqIdx i = *pi;

      mat(i,j) = log_sum_exp (mat(i-1,j-1) + m2mScore(j-1),
			      del(i-1,j-1) + d2mScore(),
			      ins(i-1,j-1) + i2mScore());

      if (j == 1 && (i == 1 || config.local))
	mat(i,j) = log_sum_exp (mat(i,j),
				start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = cachedInsertEmitScore[j] + log_sum_exp (ins(i,j-1) + i2iScore(),
							 mat(i,j-1) + m2iScore(j-1));

      del(i,j) = log_sum_exp (del(i-1,j) + d2dScore(),
			      mat(i-1,j) + m2dScore(j));

      if (j == yLen && (i == xLen || config.local))
	end = log_sum_exp (end,
			   mat(i,yLen) + m2eScore(yLen));
    }
  }

  result = end;

  if (LogThisAt(4))
    logger << "Forward score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    writeToLog();
}

QuaffBackwardMatrix::QuaffBackwardMatrix (const QuaffForwardMatrix& fwd)
  : QuaffDPMatrix (*fwd.penv, *fwd.qs.pqp, *fwd.pconfig),
    pfwd (&fwd),
    qc(fwd.qs.matchContext.kmerLen,fwd.qs.indelContext.kmerLen)
{
  Require (py->hasQual(), "Forward-Backward algorithm requires quality scores to fit model, but sequence %s lacks quality scores", py->name.c_str());

  ProgressLog (plog, 4);
  plog.initProgress ("Backward algorithm (%s vs %s)", px->name.c_str(), py->name.c_str());

  end = 0;
  for (SeqIdx j = yLen; j > 0; --j) {

    plog.logProgress ((yLen - j) / (double) yLen, "base %d/%d", j, yLen);

    for (DiagonalEnvelope::reverse_iterator pi = penv->rbegin(j);
	 !pi.finished();
	 ++pi) {

      const SeqIdx i = *pi;

      if (j == yLen && (i == xLen || pconfig->local)) {
	const double m2e = transCount (mat(i,yLen),
				       fwd.mat(i,yLen),
				       m2eScore(yLen),
				       end);
	m2eCount(yLen) += m2e;
      }

      const double matEmit = matchEmitScore(i,j);
      const double matDest = mat(i,j);
      double& matCount = matchCount(i,j);

      const double m2m = transCount (mat(i-1,j-1),
				     fwd.mat(i-1,j-1),
				     m2mScore(j-1) + matEmit,
				     matDest);
      m2mCount(j-1) += m2m;
      matCount += m2m;

      const double d2m = transCount (del(i-1,j-1),
				     fwd.del(i-1,j-1),
				     d2mScore() + matEmit,
				     matDest);
      d2mCount() += d2m;
      matCount += d2m;

      const double i2m = transCount (ins(i-1,j-1),
				     fwd.ins(i-1,j-1),
				     i2mScore() + matEmit,
				     matDest);
      i2mCount() += i2m;
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
				     m2iScore(j-1) + insEmit,
				     insDest);
      m2iCount(j-1) += m2i;
      insCount += m2i;

      const double i2i = transCount (ins(i,j-1),
				     fwd.ins(i,j-1),
				     i2iScore() + insEmit,
				     insDest);
      i2iCount() += i2i;
      insCount += i2i;

      const double delDest = del(i,j);

      const double m2d = transCount (mat(i-1,j),
				     fwd.mat(i-1,j),
				     m2dScore(j),
				     delDest);
      m2dCount(j) += m2d;

      const double d2d = transCount (del(i-1,j),
				     fwd.del(i-1,j),
				     d2dScore(),
				     delDest);
      d2dCount() += d2d;
    }
  }

  result = start;

  if (LogThisAt(4))
    logger << "Backward score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    writeToLog();

  if (LogThisIf (gsl_root_test_delta (result, fwd.result, 0, MAX_FRACTIONAL_FWDBACK_ERROR) != GSL_SUCCESS))
    logger << endl << endl << "Warning: forward score (" << fwd.result << ") does not match backward score (" << result << ")" << endl << endl << endl;

  if (LogThisAt(6)) {
    logger << "Forward-backward counts, " << px->name << " vs " << py->name << ':' << endl;
    qc.writeToLog();
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

  ProgressLog (plog, 4);
  plog.initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    plog.logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (DiagonalEnvelope::iterator pi = env.begin(j);
	 !pi.finished();
	 ++pi) {

      const SeqIdx i = *pi;

      mat(i,j) = max (max (mat(i-1,j-1) + m2mScore(j-1),
			   del(i-1,j-1) + d2mScore()),
		      ins(i-1,j-1) + i2mScore());

      if (j == 1 && (i == 1 || config.local))
	mat(i,j) = max (mat(i,j),
			start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = cachedInsertEmitScore[j] + max (ins(i,j-1) + i2iScore(),
						 mat(i,j-1) + m2iScore(j-1));

      del(i,j) = max (del(i-1,j) + d2dScore(),
		      mat(i-1,j) + m2dScore(j));

      if (j == yLen && (i == xLen || config.local))
	end = max (end,
		   mat(i,j) + m2eScore(j));
    }
  }

  result = end;

  if (LogThisAt(4))
    logger << "Viterbi score: " << result << endl;
  
  if (LogWhen("dpmatrix"))
    writeToLog();
}

Alignment QuaffViterbiMatrix::alignment() const {
  Require (resultIsFinite(), "Can't do Viterbi traceback if final score is -infinity");
  SeqIdx xEnd = xLen;
  if (pconfig->local) {
    double bestEndSc = -numeric_limits<double>::infinity();
    double sc;
    for (SeqIdx iEnd = xLen; iEnd > 0; --iEnd) {
      sc = mat(iEnd,yLen) + m2eScore(yLen);
      if (iEnd == xLen || sc > bestEndSc) {
	bestEndSc = sc;
	xEnd = iEnd;
      }
    }
  }
  SeqIdx i = xEnd, j = yLen;
  list<char> xRow, yRow, yQual;
  State state = Match;
  while (state != Start) {
    if (LogThisAt(9))
      logger << "Traceback: i=" << i << " j=" << j << " state=" << stateToString(state) << " score=" << cellScore(i,j,state) << endl;
    double srcSc = -numeric_limits<double>::infinity();
    double emitSc = 0;
    switch (state) {
    case Match:
      emitSc = matchEmitScore(i,j);
      xRow.push_front (px->seq[--i]);
      yRow.push_front (py->seq[--j]);
      if (py->hasQual())
	yQual.push_front (py->qual[j]);
      updateMax (srcSc, state, mat(i,j) + m2mScore(j) + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + i2mScore() + emitSc, Insert);
      updateMax (srcSc, state, del(i,j) + d2mScore() + emitSc, Delete);
      if (j == 0 && (i == 0 || pconfig->local))
	updateMax (srcSc, state, emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error");
      break;

    case Insert:
      emitSc = cachedInsertEmitScore[j];
      xRow.push_front (Alignment::gapChar);
      yRow.push_front (py->seq[--j]);
      if (py->hasQual())
	yQual.push_front (py->qual[j]);
      updateMax (srcSc, state, mat(i,j) + m2iScore(j) + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + i2iScore() + emitSc, Insert);
      Assert (srcSc == ins(i,j+1), "Traceback error");
      break;

    case Delete:
      xRow.push_front (px->seq[--i]);
      yRow.push_front (Alignment::gapChar);
      if (py->hasQual())
	yQual.push_front (FastSeq::maxQualityChar);
      updateMax (srcSc, state, mat(i,j) + m2dScore(j), Match);
      updateMax (srcSc, state, del(i,j) + d2dScore(), Delete);
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
  align.gappedSeq[1].qual = string (yQual.begin(), yQual.end());
  align.gappedSeq[0].source.name = px->name;
  align.gappedSeq[0].source.start = xStart;
  align.gappedSeq[0].source.end = xEnd;
  align.gappedSeq[1].source.name = py->name;
  align.gappedSeq[1].source.start = 1;
  align.gappedSeq[1].source.end = yLen;
  align.gappedSeq[0].source = align.gappedSeq[0].source.compose (px->source);
  align.gappedSeq[1].source = align.gappedSeq[1].source.compose (py->source);
  align.score = result;
  return align;
}

Alignment QuaffViterbiMatrix::scoreAdjustedAlignment (const QuaffNullParams& nullModel) const {
  Alignment a = alignment();
  const double nullLogLike = nullModel.logLikelihood (*py);
  if (LogThisAt(4))
    logger << "Null model score: " << nullLogLike << endl;
  a.score -= nullLogLike;
  return a;
}

void QuaffParamCounts::addWeighted (const QuaffParamCounts& counts, double weight) {
  Assert (counts.matchContext.kmerLen == matchContext.kmerLen, "Cannot add two QuaffParamCounts with different match kmer lengths");
  Assert (counts.indelContext.kmerLen == indelContext.kmerLen, "Cannot add two QuaffParamCounts with different indel kmer lengths");
  for (QualScore q = 0; q < FastSeq::qualScoreRange; ++q)
    for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
      insert[i].qualCount[q] += weight * counts.insert[i].qualCount[q];
      for (Kmer j = 0; j < matchContext.numKmers; ++j)
	match[i][j].qualCount[q] += weight * counts.match[i][j].qualCount[q];
    }
  beginInsertNo = vectorSum (beginInsertNo, vectorScale (weight, counts.beginInsertNo));
  beginInsertYes = vectorSum (beginInsertYes, vectorScale (weight, counts.beginInsertYes));
  beginDeleteNo = vectorSum (beginDeleteNo, vectorScale (weight, counts.beginDeleteNo));
  beginDeleteYes = vectorSum (beginDeleteYes, vectorScale (weight, counts.beginDeleteYes));
  extendInsertNo += weight * counts.extendInsertNo;
  extendInsertYes += weight * counts.extendInsertYes;
  extendDeleteNo += weight * counts.extendDeleteNo;
  extendDeleteYes += weight * counts.extendDeleteYes;
}

QuaffParamCounts QuaffParamCounts::operator+ (const QuaffParamCounts& counts) const {
  QuaffParamCounts qpc (*this);
  qpc.addWeighted (counts, 1);
  return qpc;
}

double QuaffParamCounts::logPrior (const QuaffParams& qp) const {
  double lp = 0;
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    lp += logBetaPdf (qp.beginInsert[j], beginInsertYes[j], beginInsertNo[j]);
    lp += logBetaPdf (qp.beginDelete[j], beginDeleteYes[j], beginDeleteNo[j]);
  }
  lp += logBetaPdf (qp.extendInsert, extendInsertYes, extendInsertNo);
  lp += logBetaPdf (qp.extendDelete, extendDeleteYes, extendDeleteNo);
  double *alpha = new double[dnaAlphabetSize], *theta = new double[dnaAlphabetSize];
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    lp += qp.insert[i].logQualProb (insert[i].qualCount);  // not normalized...
    theta[i] = qp.insert[i].symProb;
    alpha[i] = accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 1.);
  }
  lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize) {
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	const Kmer j = jPrefix + jSuffix;
	lp += qp.match[i][j].logQualProb (match[i][j].qualCount);  // not normalized...
	theta[jSuffix] = qp.match[i][j].symProb;
	alpha[jSuffix] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 1.);
      }
      lp += log (gsl_ran_dirichlet_pdf (dnaAlphabetSize, alpha, theta));
    }
  }
  delete[] theta;
  delete[] alpha;
  return lp;
}

double QuaffParamCounts::expectedLogLike (const QuaffParams& qp) const {
  double ll = 0;
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    ll += log (qp.beginInsert[j]) * beginInsertYes[j] + log (1 - qp.beginInsert[j]) * beginInsertNo[j];
    ll += log (qp.beginDelete[j]) * beginDeleteYes[j] + log (1 - qp.beginDelete[j]) * beginDeleteNo[j];
  }
  ll += log (qp.extendInsert) * extendInsertYes + log (1 - qp.extendInsert) * extendInsertNo;
  ll += log (qp.extendDelete) * extendDeleteYes + log (1 - qp.extendDelete) * extendDeleteNo;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    ll += qp.insert[i].logQualProb (insert[i].qualCount);
    ll += log (qp.insert[i].symProb) * accumulate (insert[i].qualCount.begin(), insert[i].qualCount.end(), 0.);
  }
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i) {
    for (Kmer j = 0; j < matchContext.numKmers; ++j) {
      ll += qp.match[i][j].logQualProb (match[i][j].qualCount);  // not normalized...
      ll += log (qp.match[i][j].symProb) * accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
    }
  }
  return ll;
}

QuaffParams QuaffParamCounts::fit() const {
  QuaffParams qp (matchContext.kmerLen, indelContext.kmerLen);
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    qp.beginDelete[j] = 1. / (1. + beginDeleteNo[j] / beginDeleteYes[j]);
    qp.beginInsert[j] = 1. / (1. + beginInsertNo[j] / beginInsertYes[j]);
  }
  qp.extendDelete = 1. / (1. + extendDeleteNo / extendDeleteYes);
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
    for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize) {
      vguard<double> iMatFreq (dnaAlphabetSize);
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	const Kmer j = jPrefix + jSuffix;
	iMatFreq[jSuffix] = accumulate (match[i][j].qualCount.begin(), match[i][j].qualCount.end(), 0.);
      }
      const double iMatNorm = accumulate (iMatFreq.begin(), iMatFreq.end(), 0.);
      for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	const Kmer j = jPrefix + jSuffix;
	qp.match[i][j].symProb = iMatFreq[jSuffix] / iMatNorm;
	fitNegativeBinomial (match[i][j].qualCount, qp.match[i][j].qualTrialSuccessProb, qp.match[i][j].qualNumSuccessfulTrials);
      }
    }
  }

  return qp;
}

QuaffForwardBackwardMatrix::QuaffForwardBackwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config)
  : fwd (env, qp, config),
    back (fwd)
{
  if (LogWhen("postmatrix"))
    writeToLog();
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
    for (DiagonalEnvelope::iterator pi = fwd.penv->begin(j);
	 !pi.finished();
	 ++pi)
      out << "i=" << *pi << ":" << fwd.px->seq[*pi-1] << " j=" << j << ":" << fwd.py->seq[j-1] << (fwd.py->hasQual() ? string(1,fwd.py->qual[j-1]) : string()) << "\tmat " << postMatch(*pi,j) << "\tins " << postInsert(*pi,j) << "\tdel " << postDelete(*pi,j) << endl;
    out << endl;
  }
}

void QuaffForwardBackwardMatrix::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
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

  if (LogThisAt(5)) {
    logger << "Null model:" << endl;
    writeToLog();
  }
}

void QuaffNullParams::read (istream& in) {
  map<string,string> val = readQuaffParamFile (in);
  Require (read (val), "Couldn't read null model");
}

bool QuaffNullParams::read (map<string,string>& val) {
  QuaffParamRead(nullEmit);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    Desire (null[i].read (val, string("null") + dnaAlphabet[i]), "Couldn't read null param");
  return true;
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
      if (std::isnan(ll)) {
	cerr << "NaN error in QuaffNullParams::logLikelihood" << endl;
	throw;
      }
#endif
    }
  return ll;
}

ostream& QuaffNullParams::writeJson (ostream& out) const {
  QuaffParamWriteJson(nullEmit);
  out << "  \"refBase\": {" << endl;
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    null[i].writeJson (out << " \"" << dnaAlphabet[i] << "\": ")
      << (i == dnaAlphabetSize-1 ? " }" : ",");
  out << endl;
  return out;
}

void QuaffNullParams::write (ostream& out) const {
  QuaffParamWrite(nullEmit);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    null[i].write (out, string("null") + dnaAlphabet[i]);
}

void QuaffNullParams::writeToLog() const {
  logger.lock();
  write(clog);
  logger.unlock();
}

QuaffTrainer::QuaffTrainer()
  : maxIterations (QuaffMaxEMIterations),
    minFractionalLoglikeIncrement (QuaffMinEMLogLikeInc),
    allowNullModel (true)
{ }

bool QuaffTrainer::parseTrainingArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-maxiter") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      maxIterations = atoi (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-mininc") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      minFractionalLoglikeIncrement = atof (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveparams") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      saveParamsFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savecountswithprior") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      countsWithPriorFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return parseCountingArgs (argvec);
}

bool QuaffTrainer::parseCountingArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-force") {
      allowNullModel = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-savecounts") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      rawCountsFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return parseServerArgs (argvec);
}

bool QuaffTrainer::parseServerArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-force") {
      allowNullModel = false;
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

string QuaffTrainer::serverArgs() const {
  return string (allowNullModel ? "" : " -force");
}

vguard<vguard<size_t> > QuaffTrainer::defaultSortOrder (const vguard<FastSeq>& x, const vguard<FastSeq>& y) {
  vguard<size_t> initialSortOrder (x.size());
  iota (initialSortOrder.begin(), initialSortOrder.end(), (size_t) 0);
  return vguard<vguard<size_t> > (y.size(), initialSortOrder);
}

QuaffParamCounts QuaffTrainer::getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, double& logLike, const char* banner) {
  QuaffCountingScheduler qcs (x, y, params, nullModel, allowNullModel, config, sortOrder, banner, VFUNCFILE(2));
  list<thread> yThreads;
  Require (config.threads > 0 || config.ec2Instances > 0 || !config.remotes.empty(), "Please allocate at least one thread or one remote server");
  for (unsigned int n = 0; n < config.threads; ++n) {
    yThreads.push_back (thread (&runQuaffCountingTasks, &qcs));
    logger.nameLastThread (yThreads, "DP");
  }
  for (const auto& remote : config.remotes) {
    yThreads.push_back (thread (&delegateQuaffCountingTasks, &qcs, &remote));
    logger.nameLastThread (yThreads, "DP");
  }
  for (auto& t : yThreads) {
    t.join();
    logger.eraseThreadName (t);
  }
  logLike = qcs.finalLogLike();
  const QuaffParamCounts counts = qcs.finalCounts();
  if (rawCountsFilename.size()) {
    ofstream out (rawCountsFilename);
    counts.write (out);
  }
  return counts;
}

QuaffParamCounts QuaffTrainer::getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config) {
  vguard<vguard<size_t> > sortOrder = defaultSortOrder (x, y);
  double logLike;
  config.startRemoteServers();
  const QuaffParamCounts counts = getCounts (x, y, params, nullModel, config, sortOrder, logLike, "");
  config.stopRemoteServers();
  return counts;
}

void QuaffTrainer::serveCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffDPConfig& config) {
  list<thread> serverThreads;
  Require (config.threads > 0, "Please allocate at least one thread");
  for (unsigned int n = 0; n < config.threads; ++n) {
    serverThreads.push_back (thread (&serveCountsFromThread, &x, &y, allowNullModel, &config, config.serverPort + n));
    logger.nameLastThread (serverThreads, "server");
  }
  for (auto& t : serverThreads) {
    t.join();
    logger.eraseThreadName (t);
  }
}

void QuaffTrainer::serveCountsFromThread (const vguard<FastSeq>* px, const vguard<FastSeq>* py, bool useNullModel, const QuaffDPConfig* pconfig, unsigned int port) {
  map<string,const FastSeq*> yDict;
  for (const auto& s : *py)
    yDict[s.name] = &s;

  if (LogThisAt(8)) {
    logger << "Known read names:";
    for (const auto& s : *py)
      logger << " \"" << s.name << "\"";
    logger << endl;
  }
  
  TCPServerSocket servSock (port);
  if (LogThisAt(1))
    logger << "(listening on port " << port << ')' << endl;

  while (true) {
    TCPSocket *sock = NULL;
    sock = servSock.accept();
    if (LogThisAt(1))
      logger << "Handling request from " << sock->getForeignAddress() << endl;

    auto msg = readQuaffStringFromSocket (sock);
    auto paramVal = readQuaffParamFile (msg);

    if (paramVal.find("quit") != paramVal.end()) {
      if (LogThisAt(1))
	logger << "(quit)" << endl;
      delete sock;
      break;
    }

    const string& yName = paramVal["yName"];
    const string& xOrder = paramVal["xOrder"];

    vguard<size_t> sortOrder = readSortOrder (xOrder);

    QuaffParams params;
    QuaffNullParams nullModel;

    if (yDict.find(yName) != yDict.end()
	&& sortOrder.size()
	&& params.read (paramVal)
	&& nullModel.read (paramVal)) {

      if (LogThisAt(2))
	logger << "Aligning " << plural(sortOrder.size(),"reference") << " to " << yName << endl;
    
      const unsigned int matchKmerLen = params.matchContext.kmerLen;
      const unsigned int indelKmerLen = params.indelContext.kmerLen;
      QuaffParamCounts yCounts (matchKmerLen, indelKmerLen);

      double yLogLike;
      QuaffCountingTask task (*px, *yDict[yName], params, nullModel, useNullModel, *pconfig, sortOrder, yLogLike, yCounts);
      task.run();

      ostringstream out;
      out << "xOrder: { ";
      for (size_t nx = 0; nx < sortOrder.size(); ++nx)
	out << (nx == 0 ? "" : ", ") << sortOrder[nx];
      out << " }" << endl;
      out << "loglike: " << yLogLike << endl;
      yCounts.write (out);
      out << SocketTerminatorString << endl;

      const string s = out.str();
      sock->send (s.c_str(), (int) s.size());

      if (LogThisAt(2))
	logger << "Request completed" << endl;

    } else if (LogThisAt(1))
      logger << "Bad request, ignoring" << endl << "sortOrder = (" << join(sortOrder) << ')' << endl << "yName = \"" << yName << "\"" << endl << "Request follows:" << endl << msg << endl;
	
    delete sock;
  }
}

QuaffParams QuaffTrainer::fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, QuaffDPConfig& config) {
  const unsigned int matchKmerLen = seed.matchContext.kmerLen;
  const unsigned int indelKmerLen = seed.indelContext.kmerLen;
  Assert (pseudocounts.matchContext.kmerLen == matchKmerLen, "Prior must have same match kmer-length as parameters");
  Assert (pseudocounts.indelContext.kmerLen == indelKmerLen, "Prior must have same indel kmer-length as parameters");
  QuaffParams qp = seed;
  QuaffParamCounts counts(matchKmerLen,indelKmerLen), countsWithPrior(matchKmerLen,indelKmerLen);
  double prevLogLikeWithPrior = -numeric_limits<double>::infinity();
  vguard<vguard<size_t> > sortOrder = defaultSortOrder (x, y);
  config.startRemoteServers();
  for (int iter = 0; iter < maxIterations; ++iter) {
    const string banner = string(" (E-step #") + to_string(iter+1) + ")";
    double logLike = 0;
    counts = getCounts (x, y, qp, nullModel, config, sortOrder, logLike, banner.c_str());
    if (LogThisAt(3)) {
      logger << "Parameter counts computed during E-step:" << endl;
      counts.writeToLog();
    }

    const double logPrior = pseudocounts.logPrior (qp);
    const double logLikeWithPrior = logLike + logPrior;
    if (LogThisAt(1))
      logger << "EM iteration " << (iter+1) << ": log-likelihood (" << logLike << ") + log-prior (" << logPrior << ") = " << logLikeWithPrior << endl;
    if (iter > 0 && logLikeWithPrior < prevLogLikeWithPrior + abs(prevLogLikeWithPrior)*minFractionalLoglikeIncrement)
      break;
    prevLogLikeWithPrior = logLikeWithPrior;

    countsWithPrior = counts;
    countsWithPrior.addWeighted (pseudocounts, 1.);
    if (countsWithPriorFilename.size()) {
      ofstream out (countsWithPriorFilename);
      countsWithPrior.write (out);
    }

    const double oldExpectedLogLike = countsWithPrior.expectedLogLike (qp);
    qp = countsWithPrior.fit();
    qp.fitRefSeqs(x);
    const double newExpectedLogLike = countsWithPrior.expectedLogLike (qp);
    if (LogThisAt(2))
      logger << "Expected log-likelihood went from " << oldExpectedLogLike << " to " << newExpectedLogLike << " during M-step" << endl;

    if (saveParamsFilename.size()) {
      ofstream out (saveParamsFilename);
      qp.write (out);
    }
  }
  config.stopRemoteServers();
  return qp;
}

QuaffCountingTask::QuaffCountingTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, const QuaffDPConfig& config, vguard<size_t>& sortOrder, double& yLogLike, QuaffParamCounts& yCounts)
  : QuaffTask (yfs, params, nullModel, config),
    x(x), useNullModel(useNullModel), sortOrder(sortOrder), yLogLike(yLogLike), yCounts(yCounts)
{ }

void QuaffCountingTask::run() {
  const KmerIndex yKmerIndex (yfs, dnaAlphabet, config.kmerLen);
  const unsigned int matchKmerLen = params.matchContext.kmerLen;
  const unsigned int indelKmerLen = params.indelContext.kmerLen;
  const double yNullLogLike = useNullModel ? nullModel.logLikelihood(yfs) : -numeric_limits<double>::infinity();
  yLogLike = yNullLogLike;  // this initial value allows null model to "win"
  if (LogThisAt(4))
    logger << "Null model score for " << yfs.name << " is " << yNullLogLike << endl;
  vguard<double> xyLogLike (x.size(), -numeric_limits<double>::infinity());
  vguard<QuaffParamCounts> xyCounts (x.size(), QuaffParamCounts(matchKmerLen,indelKmerLen));
  for (auto nx : sortOrder) {
    const auto& xfs = x[nx];
    const DiagonalEnvelope env = config.makeEnvelope (xfs, yKmerIndex, 2*QuaffDPMatrixContainer::cellSize());
    const QuaffForwardMatrix fwd (env, params, config);
    xyLogLike[nx] = fwd.result;
    if (xyLogLike[nx] >= yLogLike - MAX_TRAINING_LOG_DELTA) {  // don't waste time computing low-weight counts
      const QuaffBackwardMatrix back (fwd);
      xyCounts[nx] = back.qc;
    }
    yLogLike = log_sum_exp (yLogLike, xyLogLike[nx]);
  }
  for (size_t nx = 0; nx < x.size(); ++nx) {
    const double xyPostProb = exp (xyLogLike[nx] - yLogLike);
    if (LogThisAt(3))
      logger << "P(read " << yfs.name << " derived from ref " << x[nx].name << ") = " << xyPostProb << endl;
    yCounts.addWeighted (xyCounts[nx], xyPostProb);
  }
  if (LogThisAt(3))
    logger << "P(read " << yfs.name << " unrelated to refs) = " << exp(yNullLogLike - yLogLike) << endl;
  const vguard<size_t> ascendingOrder = orderedIndices (xyLogLike);
  sortOrder = vguard<size_t> (ascendingOrder.rbegin(), ascendingOrder.rend());
  auto sortOrderCutoff = sortOrder.begin();
  for (; sortOrderCutoff != sortOrder.end(); ++sortOrderCutoff)
    if (xyLogLike[*sortOrderCutoff] < yLogLike - MAX_TRAINING_LOG_DELTA)
      break;
  sortOrder.erase (sortOrderCutoff, sortOrder.end());  // bye bye unproductive refseqs
}

void QuaffCountingTask::delegate (const RemoteServer& remote) {
  while (true) {
    if (LogThisAt(3))
      logger << "Delegating " << yfs.name << " to " << remote.toString() << endl;
    ostringstream out;
    out << "yName: " << yfs.name << endl;
    out << "xOrder: { ";
    for (size_t nx = 0; nx < sortOrder.size(); ++nx)
      out << (nx == 0 ? "" : ", ") << sortOrder[nx];
    out << " }" << endl;
    nullModel.write (out);
    params.write (out);
    out << SocketTerminatorString << endl;

    const string msg = out.str();

    try {
      TCPSocket sock(remote.addr, remote.port);
      sock.send (msg.c_str(), (int) msg.size());

      map<string,string> paramVal = readQuaffParamFile (&sock);

      if (LogThisAt(3))
	logger << "Parsing results from " << remote.toString() << endl;

      const string& xOrder = paramVal["xOrder"];
      const string& yLogLikeStr = paramVal["loglike"];
      if (xOrder.size()
	  && yLogLikeStr.size()
	  && yCounts.read (paramVal)) {

	sortOrder = readSortOrder (xOrder);
	yLogLike = atof (yLogLikeStr.c_str());
	break;
      }

    } catch (SocketException& e) {
      if (LogThisAt(3))
	logger << e.what() << endl;
    }

    randomDelayBeforeRetry (MinQuaffRetryDelay, MaxQuaffRetryDelay);
  }
}

QuaffCountingScheduler::QuaffCountingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, const QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, const char* banner, int verbosity, const char* function, const char* file)
  : QuaffScheduler (x, y, params, nullModel, config, verbosity, function, file),
    useNullModel(useNullModel),
    sortOrder(sortOrder),
    banner(banner),
    matchKmerLen (params.matchContext.kmerLen),
    indelKmerLen (params.indelContext.kmerLen),
    zeroCounts (matchKmerLen, indelKmerLen),
    yCounts (y.size(), zeroCounts),
    yLogLike (y.size(), 0.)
{
  plog.initProgress ("Calculating expected counts%s", banner);
}

bool QuaffCountingScheduler::finished() const {
  return ny == y.size();
}

QuaffCountingTask QuaffCountingScheduler::nextCountingTask() {
  const size_t n = ny++;
  plog.logProgress (n / (double) y.size(), "finished %lu/%lu reads", n, y.size());
  return QuaffCountingTask (x, y[n], params, nullModel, useNullModel, config, sortOrder[n], yLogLike[n], yCounts[n]);
}

QuaffParamCounts QuaffCountingScheduler::finalCounts() const {
  return accumulate (yCounts.begin(), yCounts.end(), zeroCounts);
}

double QuaffCountingScheduler::finalLogLike() const {
  return accumulate (yLogLike.begin(), yLogLike.end(), 0.);
}

void runQuaffCountingTasks (QuaffCountingScheduler* qcs) {
  while (true) {
    qcs->lock();
    if (qcs->finished()) {
      qcs->unlock();
      break;
    }
    QuaffCountingTask task = qcs->nextCountingTask();
    qcs->unlock();
    task.run();
  }
}

void delegateQuaffCountingTasks (QuaffCountingScheduler* qcs, const RemoteServer* remote) {
  while (true) {
    qcs->lock();
    if (qcs->finished()) {
      qcs->unlock();
      break;
    }
    QuaffCountingTask task = qcs->nextCountingTask();
    qcs->unlock();
    task.delegate (*remote);
  }
}

QuaffAlignmentPrinter::QuaffAlignmentPrinter()
  : format (StockholmAlignment),
    logOddsThreshold (0)
{ }

bool QuaffAlignmentPrinter::parseAlignmentPrinterArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-format") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const string fmt = argvec[1];
      if (fmt == "fasta")
	format = GappedFastaAlignment;
      else if (fmt == "stockholm")
	format = StockholmAlignment;
      else if (fmt == "sam")
	format = SamAlignment;
      else if (fmt == "refseq")
	format = UngappedFastaRef;
      else
	Fail ("Unknown format: %s", fmt.c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-threshold") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      const char* val = argvec[1].c_str();
      logOddsThreshold = atof (val);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nothreshold") {
      logOddsThreshold = -numeric_limits<double>::infinity();
      argvec.pop_front();
      return true;

    } else if (arg == "-savealign") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      alignFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    }
  }

  return false;
}

string QuaffAlignmentPrinter::serverArgs() const {
  string args = " -format ";
  switch (format) {
  case GappedFastaAlignment:
    args += "fasta";
    break;

  case StockholmAlignment:
    args += "stockholm";
    break;

  case SamAlignment:
    args += "sam";
    break;

  case UngappedFastaRef:
    args += "refseq";
    break;

  default:
    Abort ("Unrecognized alignment format");
    break;
  }

  if (logOddsThreshold > -numeric_limits<double>::infinity())
    args += " -threshold " + to_string(logOddsThreshold);
  else
    args += " -nothreshold";

  return args;
}

void QuaffAlignmentPrinter::writeAlignmentHeader (ostream& out, const vguard<FastSeq>& refs, bool groupByQuery) {
  if (usingOutputFile())
    alignFile.open (alignFilename);
  if (format == SamAlignment)
    Alignment::writeSamHeader (usingOutputFile() ? alignFile : out, refs, groupByQuery ? "GO:query" : "SO:unknown");
}

void QuaffAlignmentPrinter::writeAlignment (ostream& out, const Alignment& align) const {
  ostream& alignOut = alignmentOutputStream (out);
  if (align.score >= logOddsThreshold) {
    FastSeq ref;
    switch (format) {
    case GappedFastaAlignment:
      align.writeGappedFasta (alignOut);
      out << endl;
      break;

    case StockholmAlignment:
      align.writeStockholm (alignOut);
      break;

    case SamAlignment:
      align.writeSam (alignOut);
      break;

    case UngappedFastaRef:
      Assert (align.rows() == 2, "Not a pairwise alignment");
      ref = align.getUngapped(0);
      ref.comment = string("matches(") + align.gappedSeq[1].name + ") " + ref.comment;
      ref.writeFasta (alignOut);
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

bool QuaffAligner::parseAlignmentArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-printall") {
      printAllAlignments = true;
      argvec.pop_front();
      return true;
    }
  }

  return parseAlignmentPrinterArgs (argvec);
}

string QuaffAligner::serverArgs() const {
  return QuaffAlignmentPrinter::serverArgs() + (printAllAlignments ? " -printall" : "");
}

void QuaffAligner::align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config) {
  QuaffAlignmentScheduler qas (x, y, params, nullModel, config, printAllAlignments, out, *this, VFUNCFILE(2));
  list<thread> yThreads;
  Require (config.threads > 0 || config.ec2Instances > 0 || !config.remotes.empty(), "Please allocate at least one thread or one remote server");
  config.startRemoteServers();
  for (unsigned int n = 0; n < config.threads; ++n) {
    yThreads.push_back (thread (&runQuaffAlignmentTasks, &qas));
    logger.nameLastThread (yThreads, "align");
  }
  for (const auto& remote : config.remotes) {
    yThreads.push_back (thread (&delegateQuaffAlignmentTasks, &qas, &remote));
    logger.nameLastThread (yThreads, "align");
  }
  for (auto& t : yThreads) {
    t.join();
    logger.eraseThreadName (t);
  }
  config.stopRemoteServers();
}

void QuaffAligner::serveAlignments (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config) {
  list<thread> serverThreads;
  Require (config.threads > 0, "Please allocate at least one thread");
  for (unsigned int n = 0; n < config.threads; ++n) {
    serverThreads.push_back (thread (&serveAlignmentsFromThread, this, &x, &y, &params, &nullModel, &config, config.serverPort + n));
    logger.nameLastThread (serverThreads, "server");
  }
  for (auto& t : serverThreads) {
    t.join();
    logger.eraseThreadName (t);
  }
}

void QuaffAligner::serveAlignmentsFromThread (QuaffAligner* paligner, const vguard<FastSeq>* px, const vguard<FastSeq>* py, const QuaffParams* pparams, const QuaffNullParams* pnullModel, const QuaffDPConfig* pconfig, unsigned int port) {
  map<string,const FastSeq*> yDict;
  for (const auto& s : *py)
    yDict[s.name] = &s;

  if (LogThisAt(8)) {
    logger << "Known read names:";
    for (const auto& s : *py)
      logger << " \"" << s.name << "\"";
    logger << endl;
  }

  TCPServerSocket servSock (port);
  if (LogThisAt(1))
    logger << "(listening on port " << port << ')' << endl;

  while (true) {
    TCPSocket *sock = NULL;
    sock = servSock.accept();
    if (LogThisAt(1))
      logger << "Handling request from " << sock->getForeignAddress() << endl;

    auto msg = readQuaffStringFromSocket (sock);
    auto paramVal = readQuaffParamFile (msg);

    if (paramVal.find("quit") != paramVal.end()) {
      if (LogThisAt(1))
	logger << "(quit)" << endl;
      delete sock;
      break;
    }

    const string& yName = paramVal["yName"];

    if (yDict.find(yName) != yDict.end()) {

      if (LogThisAt(2))
	logger << "Aligning " << yName << endl;

      QuaffAlignmentTask task (*px, *yDict[yName], *pparams, *pnullModel, *pconfig, paligner->printAllAlignments);
      task.run();

      ostringstream out;
      for (const auto& a : task.alignList)
	paligner->writeAlignment (out, a);
      out << SocketTerminatorString << endl;

      const string s = out.str();
      sock->send (s.c_str(), (int) s.size());

      if (LogThisAt(2))
	logger << "Request completed" << endl;

    } else if (LogThisAt(1))
      logger << "Bad request, ignoring" << endl << "yName = \"" << yName << "\"" << endl << "Request follows:" << endl << msg << endl;
	
    delete sock;
  }
}

QuaffAlignmentTask::QuaffAlignmentTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, bool keepAllAlignments)
  : QuaffTask(yfs,params,nullModel,config),
    x(x),
    alignList(Alignment::scoreGreaterThan),
    keepAllAlignments(keepAllAlignments)
{ }

void QuaffAlignmentTask::run() {
  const KmerIndex yKmerIndex (yfs, dnaAlphabet, config.kmerLen);
  for (size_t nx = 0; nx < x.size(); ++nx) {
    const FastSeq& xfs = x[nx];
    if (LogThisAt(3))
      logger << "Aligning " << xfs.name << " (length " << xfs.length() << ") to " << yfs.name << " (length " << yfs.length() << ")" << endl;
    const DiagonalEnvelope env = config.makeEnvelope (xfs, yKmerIndex, QuaffDPMatrixContainer::cellSize());
    const QuaffViterbiMatrix viterbi (env, params, config);
    if (viterbi.resultIsFinite()) {
      const Alignment a = viterbi.scoreAdjustedAlignment (nullModel);
      alignList.insert (a);
      if (!keepAllAlignments)
	alignList.erase (++alignList.begin(), alignList.end());
    }
  }
}

string QuaffAlignmentTask::delegate (const RemoteServer& remote) {
  if (LogThisAt(3))
    logger << "Delegating " << yfs.name << " to " << remote.toString() << endl;
  ostringstream out;
  out << "yName: " << yfs.name << endl;
  out << SocketTerminatorString << endl;
  
  const string msg = out.str();
  string response;

  while (true) {
    try {
      TCPSocket sock(remote.addr, remote.port);
      sock.send (msg.c_str(), (int) msg.size());

      response = readQuaffStringFromSocket (&sock);
      break;

    } catch (SocketException& e) {
      if (LogThisAt(3))
	logger << e.what() << endl;
    }

    randomDelayBeforeRetry (MinQuaffRetryDelay, MaxQuaffRetryDelay);
  }

  return response;
}

QuaffAlignmentPrintingScheduler::QuaffAlignmentPrintingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file)
  : QuaffScheduler (x, y, params, nullModel, config, verbosity, function, file),
    out(out),
    printer(printer)
{ }

void QuaffAlignmentPrintingScheduler::printAlignments (const QuaffAlignmentPrinter::AlignmentList& alignList) {
  outMx.lock();
  logger.lockSilently();  // to avoid interleaving clog & cout
  for (const auto& a : alignList)
    printer.writeAlignment (out, a);
  logger.unlockSilently();
  outMx.unlock();
}

void QuaffAlignmentPrintingScheduler::printAlignments (const string& alignStr) {
  outMx.lock();
  logger.lockSilently();  // to avoid interleaving clog & cout
  printer.alignmentOutputStream(out) << alignStr;
  logger.unlockSilently();
  outMx.unlock();
}

QuaffAlignmentScheduler::QuaffAlignmentScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, bool keepAllAlignments, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file)
  : QuaffAlignmentPrintingScheduler (x, y, params, nullModel, config, out, printer, verbosity, function, file),
    keepAllAlignments(keepAllAlignments)
{
  printer.writeAlignmentHeader (out, x, true);
  plog.initProgress ("Alignment");
}

bool QuaffAlignmentScheduler::finished() const {
  return ny == y.size();
}

QuaffAlignmentTask QuaffAlignmentScheduler::nextAlignmentTask() {
  const size_t n = ny++;
  plog.logProgress (n / (double) y.size(), "aligned %lu/%lu reads", n, y.size());
  return QuaffAlignmentTask (x, y[n], params, nullModel, config, keepAllAlignments);
}

void runQuaffAlignmentTasks (QuaffAlignmentScheduler* qas) {
  while (true) {
    qas->lock();
    if (qas->finished()) {
      qas->unlock();
      break;
    }
    QuaffAlignmentTask task = qas->nextAlignmentTask();
    qas->unlock();
    task.run();
    qas->printAlignments (task.alignList);
  }
}

void delegateQuaffAlignmentTasks (QuaffAlignmentScheduler* qas, const RemoteServer* remote) {
  while (true) {
    qas->lock();
    if (qas->finished()) {
      qas->unlock();
      break;
    }
    QuaffAlignmentTask task = qas->nextAlignmentTask();
    qas->unlock();
    const string alignStr = task.delegate (*remote);
    qas->printAlignments (alignStr);
  }
}
