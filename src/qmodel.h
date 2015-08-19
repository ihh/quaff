#ifndef QMODEL_INCLUDED
#define QMODEL_INCLUDED

#include <map>
#include <numeric>
#include <deque>
#include <list>
#include "fastseq.h"
#include "diagenv.h"

// EM convergence parameters
#define QuaffMaxEMIterations 100
#define QuaffMinEMLogLikeInc .01

// Default contexts
#define DefaultMatchKmerContext 1
#define DefaultIndelKmerContext 0

// struct describing the probability of a given FASTA symbol,
// and a negative binomial distribution over the associated quality score
struct SymQualDist {
  double symProb; // probability of symbol
  double qualTrialSuccessProb, qualNumSuccessfulTrials;  // parameters of neg.binom. distribution
  SymQualDist();
  double logQualProb (QualScore k) const;
  double logQualProb (const vguard<double>& kFreq) const;
  void write (ostream& out, const string& prefix) const;
  void read (map<string,string>& paramVal, const string& prefix);
};

// Memo-ized log scores for a SymQualDist
struct SymQualScores {
  double logSymProb;  // logSymProb = P(this symbol), quality score marginalized
  vguard<double> logSymQualProb;  // logSymQualProb[q] = P(this symbol wih quality score q)
  SymQualScores()
    : logSymQualProb (FastSeq::qualScoreRange, -numeric_limits<double>::infinity()),
      logSymProb (-numeric_limits<double>::infinity())
  { }
  SymQualScores (const SymQualDist& sqd);
};

// Summary statistics for a SymQualDist
struct SymQualCounts {
  vguard<double> qualCount;  // no. of times each quality score seen
  SymQualCounts();
  double symCount() const { return accumulate (qualCount.begin(), qualCount.end(), 0.); }
  void write (ostream& out, const string& prefix) const;
  void read (const string& counts);
};

// Classes to manage kmer-dependence of various parameters
struct QuaffKmerContext {
  void initKmerContext (unsigned int newKmerLen);
  const char* prefix;
  unsigned int kmerLen, defaultKmerLen;
  Kmer numKmers;
  QuaffKmerContext (const char* prefix, unsigned int kmerLen);
  string kmerString (Kmer kmer) const;
  string insertParamName (AlphTok i) const;
  void readKmerLen (map<string,string>& paramVal);
  void writeKmerLen (ostream& out) const;
};

struct QuaffMatchKmerContext : QuaffKmerContext {
  QuaffMatchKmerContext (unsigned int kmerLen)
    : QuaffKmerContext ("match", kmerLen)
  { }
  string matchParamName (AlphTok i, Kmer j) const;
};

struct QuaffIndelKmerContext : QuaffKmerContext {
  QuaffIndelKmerContext (unsigned int kmerLen)
    : QuaffKmerContext ("gap", kmerLen)
  { }
  string booleanParamName (const char* tag, Kmer j) const;
};

// Parameters of a quaff model
struct QuaffParams {
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  vguard<double> refBase;
  vguard<double> beginInsert, beginDelete;
  double extendInsert, extendDelete;
  vguard<SymQualDist> insert;  // emissions from insert state
  vguard<vguard<SymQualDist> > match;  // substitutions from match state (conditional on input)
  QuaffParams (unsigned int matchKmerLen = DefaultMatchKmerContext, unsigned int indelKmerLen = DefaultIndelKmerContext);
  void resize();  // call after changing kmerLen
  void writeToLog() const;
  void write (ostream& out) const;
  void read (istream& in);
  void fitRefSeqs (const vguard<FastSeq>& refs);
};

struct QuaffNullParams {
  double nullEmit;
  vguard<SymQualDist> null;
  QuaffNullParams();
  QuaffNullParams (const vguard<FastSeq>& seqs, double pseudocount = 1);
  double logLikelihood (const FastSeq& seq) const;
  double logLikelihood (const vguard<FastSeq>& seqs) const;
  void writeToLog() const;
  void write (ostream& out) const;
  void read (istream& in);
};

// Memo-ized scores for transitions & emissions in quaff HMM
struct QuaffScores {
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  const QuaffParams *pqp;
  vguard<SymQualScores> insert;
  vguard<vguard<SymQualScores> > match;
  vguard<double> m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffScores (const QuaffParams& qp);
};

// Summary statistics for a quaff model
struct QuaffCounts {
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  vguard<SymQualCounts> insert;
  vguard<vguard<SymQualCounts> > match;
  vguard<double> m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffCounts (unsigned int matchKmerLen, unsigned int indelKmerLen);
  void writeToLog() const;
  void write (ostream& out) const;
};

struct QuaffParamCounts {
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  vguard<SymQualCounts> insert;
  vguard<vguard<SymQualCounts> > match;
  vguard<double> beginInsertNo, beginInsertYes, beginDeleteNo, beginDeleteYes;
  double extendInsertNo, extendInsertYes, extendDeleteNo, extendDeleteYes;
  QuaffParamCounts (unsigned int matchKmerLen = DefaultMatchKmerContext, unsigned int indelKmerLen = DefaultIndelKmerContext);
  QuaffParamCounts (const QuaffCounts& counts);
  void resize();  // call after changing kmerLen
  void zeroCounts();
  void initCounts (double noBeginCount, double yesExtendCount, double matchIdentCount, double otherCount, const QuaffNullParams* nullModel = NULL);
  void writeToLog() const;
  void write (ostream& out) const;
  void read (istream& in);
  void addWeighted (const QuaffParamCounts& counts, double weight);
  QuaffParamCounts operator+ (const QuaffParamCounts& counts) const;
  QuaffParams fit() const;  // maximum-likelihood fit
  double logPrior (const QuaffParams& qp) const;  // uses counts as hyperparameters to define a prior over params
  double expectedLogLike (const QuaffParams& qp) const;  // as logPrior, but unnormalized
};

// Alignment
struct Alignment {
  static const char gapChar, mismatchChar;
  vguard<FastSeq> gappedSeq;
  double score;
  Alignment (size_t numRows = 0)
    : gappedSeq(numRows), score(-numeric_limits<double>::infinity())
  { }
  const size_t rows() const { return gappedSeq.size(); }
  const size_t columns() const { return rows() ? gappedSeq[0].length() : 0; }
  void writeGappedFasta (ostream& out) const;
  void writeStockholm (ostream& out) const;
  FastSeq getUngapped (size_t row) const;
  static bool isGapChar(char c) { return c == '-' || c == '.'; }
};

// DP config
struct QuaffDPConfig {
  bool local, sparse, autoMemSize;
  int kmerLen, kmerThreshold, bandSize;
  size_t maxSize;
  unsigned int threads;
  QuaffDPConfig()
    : local(true),
      sparse(true),
      autoMemSize(false),
      kmerLen(DEFAULT_KMER_LENGTH),
      kmerThreshold(DEFAULT_KMER_THRESHOLD),
      maxSize(0),
      bandSize(DEFAULT_BAND_SIZE),
      threads(1)
  { }
  bool parseRefSeqConfigArgs (deque<string>& argvec);
  bool parseGeneralConfigArgs (deque<string>& argvec);
  DiagonalEnvelope makeEnvelope (const FastSeq& x, const FastSeq& y, size_t cellSize) const;
  size_t effectiveMaxSize() const;  // takes threading into account
  bool fixedMemoryEnvelope() const { return kmerThreshold < 0; }
};

// DP matrices
struct QuaffDPCell {
  double mat, ins, del;
  QuaffDPCell();
};

typedef map<SeqIdx,QuaffDPCell> QuaffDPColumn;

class QuaffDPMatrixContainer {
public:
  enum State { Start, Match, Insert, Delete };
  const DiagonalEnvelope* penv;
  const FastSeq *px, *py;
  SeqIdx xLen, yLen;
  vguard<QuaffDPColumn> cell;
  double start, end, result;
  static QuaffDPCell dummy;
  QuaffDPMatrixContainer (const DiagonalEnvelope& env);
  inline const QuaffDPCell& getCell (SeqIdx i, SeqIdx j) const {
    const QuaffDPColumn& col = cell[j];
    auto c = col.find (i);
    return c == col.end() ? dummy : c->second;
  }
  inline double& mat (SeqIdx i, SeqIdx j) { return cell[j][i].mat; }
  inline double& ins (SeqIdx i, SeqIdx j) { return cell[j][i].ins; }
  inline double& del (SeqIdx i, SeqIdx j) { return cell[j][i].del; }
  inline const double& mat (SeqIdx i, SeqIdx j) const { return getCell(i,j).mat; }
  inline const double& ins (SeqIdx i, SeqIdx j) const { return getCell(i,j).ins; }
  inline const double& del (SeqIdx i, SeqIdx j) const { return getCell(i,j).del; }
  double cellScore (SeqIdx i, SeqIdx j, State state) const;
  static const char* stateToString (State state);
  static vguard<size_t> getFastSeqLengths (const vguard<FastSeq>& db);
protected:
  static void updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
};

struct QuaffDPMatrix : QuaffDPMatrixContainer {
  const QuaffDPConfig* pconfig;
  QuaffScores qs;
  vguard<AlphTok> xTok, yTok;
  vector<Kmer> yMatchKmer, yIndelKmer;
  vguard<QualScore> yQual;
  vguard<double> cachedInsertEmitScore;
  QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
  inline double m2mScore (SeqIdx j) const { return j > 0 ? qs.m2m[yIndelKmer[j-1]] : 0; }
  inline double m2iScore (SeqIdx j) const { return j > 0 ? qs.m2i[yIndelKmer[j-1]] : 0; }
  inline double m2dScore (SeqIdx j) const { return j > 0 ? qs.m2d[yIndelKmer[j-1]] : 0; }
  inline double m2eScore (SeqIdx j) const { return j > 0 ? qs.m2e[yIndelKmer[j-1]] : 0; }
  inline double i2iScore() const { return qs.i2i; }
  inline double i2mScore() const { return qs.i2m; }
  inline double d2dScore() const { return qs.d2d; }
  inline double d2mScore() const { return qs.d2m; }
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    const SymQualScores& sqs = qs.match[xTok[i-1]][yMatchKmer[j-1]];
    return yQual.size() ? sqs.logSymQualProb[yQual[j-1]] : sqs.logSymProb;
  }
  inline double insertEmitScore (SeqIdx j) const {
    const SymQualScores& sqs = qs.insert[yTok[j-1]];
    return yQual.size() ? sqs.logSymQualProb[yQual[j-1]] : sqs.logSymProb;
  }
  void writeToLog() const;
  void write (ostream& out) const;
};

struct QuaffForwardMatrix : QuaffDPMatrix {
  QuaffForwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
};

struct QuaffBackwardMatrix : QuaffDPMatrix {
  const QuaffForwardMatrix *pfwd;
  QuaffCounts qc;
  QuaffBackwardMatrix (const QuaffForwardMatrix& fwd);
  double transCount (double& backSrc, double fwdSrc, double trans, double backDest) const;
  double dummyCount;
  inline double& m2mCount (SeqIdx j) { return j > 0 ? qc.m2m[yIndelKmer[j-1]] : dummyCount; }
  inline double& m2iCount (SeqIdx j) { return j > 0 ? qc.m2i[yIndelKmer[j-1]] : dummyCount; }
  inline double& m2dCount (SeqIdx j) { return j > 0 ? qc.m2d[yIndelKmer[j-1]] : dummyCount; }
  inline double& m2eCount (SeqIdx j) { return j > 0 ? qc.m2e[yIndelKmer[j-1]] : dummyCount; }
  inline double& i2iCount() { return qc.i2i; }
  inline double& i2mCount() { return qc.i2m; }
  inline double& d2dCount() { return qc.d2d; }
  inline double& d2mCount() { return qc.d2m; }
  inline double& matchCount (SeqIdx i, SeqIdx j) {
    return qc.match[xTok[i-1]][yMatchKmer[j-1]].qualCount[yQual[j-1]];
  }
  inline double& insertCount (SeqIdx j) {
    return qc.insert[yTok[j-1]].qualCount[yQual[j-1]];
  }
};

struct QuaffForwardBackwardMatrix {
  QuaffForwardMatrix fwd;
  QuaffBackwardMatrix back;
  QuaffForwardBackwardMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
  double postMatch (SeqIdx i, SeqIdx j) const;
  double postDelete (SeqIdx i, SeqIdx j) const;
  double postInsert (SeqIdx i, SeqIdx j) const;
  void writeToLog() const;
  void write (ostream& out) const;
};
  
struct QuaffViterbiMatrix : QuaffDPMatrix {
  QuaffViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  Alignment alignment() const;
  Alignment scoreAdjustedAlignment (const QuaffNullParams& nullModel) const;
};

// config/wrapper struct for Baum-Welch style EM algorithm
struct QuaffTrainer {
  int maxIterations;
  double minFractionalLoglikeIncrement;
  bool allowNullModel;
  string rawCountsFilename, countsWithPriorFilename;
  
  QuaffTrainer();
  bool parseTrainingArgs (deque<string>& argvec);
  bool parseCountingArgs (deque<string>& argvec);
  QuaffParams fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, const QuaffDPConfig& config);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, double& logLike, const char* banner);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
  static vguard<vguard<size_t> > defaultSortOrder (const vguard<FastSeq>& x, const vguard<FastSeq>& y);
};

// struct encapsulating a single training thread
struct QuaffCountingTask {
  const vguard<FastSeq>& x;
  const FastSeq& yfs;
  const QuaffParams& params;
  const QuaffDPConfig& config;
  const QuaffNullParams* nullModel;
  vguard<size_t>& sortOrder;
  double& yLogLike;
  QuaffParamCounts& yCounts;
  QuaffCountingTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams* nullModel, const QuaffDPConfig& config, vguard<size_t>& sortOrder, double& yLogLike, QuaffParamCounts& counts);
  void run();
};
void runQuaffCountingTask (QuaffCountingTask* task);

// config/wrapper structs for Viterbi alignment
struct QuaffAlignmentPrinter {
  enum OutputFormat { GappedFastaAlignment, StockholmAlignment, UngappedFastaRef } format;
  double logOddsThreshold;

  QuaffAlignmentPrinter();
  bool parseAlignmentPrinterArgs (deque<string>& argvec);
  void writeAlignment (ostream& out, const Alignment& align) const;
};

struct QuaffAligner : QuaffAlignmentPrinter {
  bool printAllAlignments;  // print alignment to every reference sequence that matches

  QuaffAligner();
  bool parseAlignmentArgs (deque<string>& argvec);
  void align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
};

// struct encapsulating a single alignment thread
struct QuaffAlignmentTask {
  const vguard<FastSeq>& x;
  const FastSeq& yfs;
  const QuaffParams& params;
  const QuaffDPConfig& config;
  const QuaffNullParams& nullModel;
  list<Alignment>& alignList;
  bool keepAllAlignments;
  QuaffAlignmentTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, list<Alignment>& alignList, bool keepAllAlignments);
  void run();
};
void runQuaffAlignmentTask (QuaffAlignmentTask* task);

#endif /* QMODEL_INCLUDED */
