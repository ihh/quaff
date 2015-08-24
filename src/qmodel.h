#ifndef QMODEL_INCLUDED
#define QMODEL_INCLUDED

#include <list>
#include <map>
#include <numeric>
#include <deque>
#include <set>
#include <mutex>
#include <fstream>
#include "fastseq.h"
#include "diagenv.h"
#include "logger.h"
#include "PracticalSocket.h"

// EM convergence parameters
#define QuaffMaxEMIterations 100
#define QuaffMinEMLogLikeInc .01

// Default contexts
#define DefaultMatchKmerContext 1
#define DefaultIndelKmerContext 0

// Default server port
#define DefaultServerPort 8000

// Default size of receive buffer for sockets
#define RCVBUFSIZE 1024

// Terminator string for socket messages
#define SocketTerminatorString "# EOF"

// useful helper methods
string readQuaffStringFromSocket (TCPSocket* sock, int bufSize = RCVBUFSIZE);
map<string,string> readQuaffParamFile (istream& in);
map<string,string> readQuaffParamFile (TCPSocket* sock);

// struct describing the probability of a given FASTA symbol,
// and a negative binomial distribution over the associated quality score
struct SymQualDist {
  double symProb; // probability of symbol
  double qualTrialSuccessProb, qualNumSuccessfulTrials;  // parameters of neg.binom. distribution
  SymQualDist();
  double logQualProb (QualScore k) const;
  double logQualProb (const vguard<double>& kFreq) const;
  void write (ostream& out, const string& prefix) const;
  bool read (map<string,string>& paramVal, const string& prefix);
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
  bool read (map<string,string>& paramVal, const string& param);
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
  bool read (map<string,string>& paramVal);
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
  bool read (map<string,string>& paramVal);
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
  bool read (map<string,string>& paramVal);
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
  void writeSam (ostream& out) const;
  static void writeSamHeader (ostream& out, const vguard<FastSeq>& seqs, const char* go_so = "SO:unknown");
  string cigarString() const;
  FastSeq getUngapped (size_t row) const;
  Alignment revcomp() const;
  static bool isGapChar(char c) { return c == '-' || c == '.'; }
  static bool scoreGreaterThan (const Alignment& a, const Alignment& b) { return a.score > b.score; }
};

// DP config
struct QuaffDPConfig {
  bool local, sparse, autoMemSize;
  int kmerLen, kmerThreshold, bandSize;
  size_t maxSize;
  unsigned int threads;
  list<RemoteServer> remotes;
  string bucket;
  int serverPort;
  QuaffDPConfig()
    : local(true),
      sparse(true),
      autoMemSize(false),
      kmerLen(DEFAULT_KMER_LENGTH),
      kmerThreshold(DEFAULT_KMER_THRESHOLD),
      maxSize(0),
      bandSize(DEFAULT_BAND_SIZE),
      threads(1),
      serverPort(DefaultServerPort)
  { }
  bool parseRefSeqConfigArgs (deque<string>& argvec);
  bool parseGeneralConfigArgs (deque<string>& argvec);
  bool parseServerConfigArgs (deque<string>& argvec);
  DiagonalEnvelope makeEnvelope (const FastSeq& x, const KmerIndex& yKmerIndex, size_t cellSize) const;
  size_t effectiveMaxSize() const;  // takes threading into account
  void loadFromBucket (const string& filename) const;
  void saveToBucket (const string& filename) const;
};

class QuaffDPMatrixContainer {
public:
  enum State { Start, Match, Insert, Delete };
  const DiagonalEnvelope* penv;
  const FastSeq *px, *py;
  SeqIdx xLen, yLen;
  vguard<double> cell;
  double start, end, result;
  static double dummy;
  QuaffDPMatrixContainer (const DiagonalEnvelope& env);
  inline double& getCell (SeqIdx i, SeqIdx j, unsigned int offset) {
    const int storageIndex = penv->getStorageIndexUnsafe (i, j);
    return cell[storageIndex*3 + offset];
  }
  inline const double& getCell (SeqIdx i, SeqIdx j, unsigned int offset) const {
    const int storageIndex = penv->getStorageIndexSafe (i, j);
    return storageIndex < 0 ? dummy : cell[storageIndex*3 + offset];
  }
  inline double& mat (SeqIdx i, SeqIdx j) { return getCell(i,j,0); }
  inline double& ins (SeqIdx i, SeqIdx j) { return getCell(i,j,1); }
  inline double& del (SeqIdx i, SeqIdx j) { return getCell(i,j,2); }
  inline const double& mat (SeqIdx i, SeqIdx j) const { return getCell(i,j,0); }
  inline const double& ins (SeqIdx i, SeqIdx j) const { return getCell(i,j,1); }
  inline const double& del (SeqIdx i, SeqIdx j) const { return getCell(i,j,2); }
  double cellScore (SeqIdx i, SeqIdx j, State state) const;
  static const char* stateToString (State state);
  static vguard<size_t> getFastSeqLengths (const vguard<FastSeq>& db);
  static size_t cellSize() { return 3*sizeof(double); }
protected:
  static void updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
};

struct QuaffDPMatrix : QuaffDPMatrixContainer {
  const QuaffDPConfig* pconfig;
  QuaffScores qs;
  vguard<AlphTok> xTok, yTok;
  // yIndelKmer is padded with a dummy entry at the start
  // this avoids the need for a bounds test in m2*Score(j) methods
  vector<Kmer> yMatchKmer, yIndelKmer;
  vguard<QualScore> yQual;
  vguard<double> cachedInsertEmitScore;
  QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
  inline double m2mScore (SeqIdx j) const { return qs.m2m[yIndelKmer[j]]; }
  inline double m2iScore (SeqIdx j) const { return qs.m2i[yIndelKmer[j]]; }
  inline double m2dScore (SeqIdx j) const { return qs.m2d[yIndelKmer[j]]; }
  inline double m2eScore (SeqIdx j) const { return qs.m2e[yIndelKmer[j]]; }
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
  inline double& m2mCount (SeqIdx j) { return qc.m2m[yIndelKmer[j]]; }
  inline double& m2iCount (SeqIdx j) { return qc.m2i[yIndelKmer[j]]; }
  inline double& m2dCount (SeqIdx j) { return qc.m2d[yIndelKmer[j]]; }
  inline double& m2eCount (SeqIdx j) { return qc.m2e[yIndelKmer[j]]; }
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
  string rawCountsFilename, countsWithPriorFilename, saveParamsFilename;
  
  QuaffTrainer();
  bool parseTrainingArgs (deque<string>& argvec);
  bool parseCountingArgs (deque<string>& argvec);
  bool parseServerArgs (deque<string>& argvec);
  QuaffParams fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, const QuaffDPConfig& config);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, double& logLike, const char* banner);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
  void serveCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffDPConfig& config);
  static void serveCountsFromThread (const vguard<FastSeq>* px, const vguard<FastSeq>* py, bool useNullModel, const QuaffDPConfig* pconfig, unsigned int port);
  static vguard<vguard<size_t> > defaultSortOrder (const vguard<FastSeq>& x, const vguard<FastSeq>& y);
  bool usingParamOutputFile() const { return !saveParamsFilename.empty(); }
  bool usingCountsOutputFile() const { return !rawCountsFilename.empty(); }
};

// structs for scheduling training tasks
struct QuaffTask {
  const FastSeq& yfs;
  const QuaffParams& params;
  const QuaffDPConfig& config;
  const QuaffNullParams& nullModel;
  QuaffTask (const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config)
    : yfs(yfs), params(params), nullModel(nullModel), config(config)
  { }
};

struct QuaffCountingTask : QuaffTask {
  const vguard<FastSeq>& x;
  vguard<size_t>& sortOrder;
  double& yLogLike;
  QuaffParamCounts& yCounts;
  const bool useNullModel;
  QuaffCountingTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, const QuaffDPConfig& config, vguard<size_t>& sortOrder, double& yLogLike, QuaffParamCounts& counts);
  void run();
  void delegate (const RemoteServer& remote);
};

struct QuaffScheduler {
  const vguard<FastSeq>& x;
  const vguard<FastSeq>& y;
  const QuaffParams& params;
  const QuaffDPConfig& config;
  const QuaffNullParams& nullModel;
  size_t ny;
  mutex mx;
  ProgressLogger plog;
  QuaffScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, int verbosity, const char* function, const char* file)
    : x(x), y(y), params(params), nullModel(nullModel), config(config), ny(0), plog(verbosity,function,file)
  { }
  void lock() { mx.lock(); }
  void unlock() { mx.unlock(); }
};

struct QuaffCountingScheduler : QuaffScheduler {
  const bool useNullModel;
  vguard<vguard<size_t> >& sortOrder;
  vguard<double> yLogLike;
  const char* banner;
  const unsigned int matchKmerLen, indelKmerLen;
  const QuaffParamCounts zeroCounts;
  vguard<QuaffParamCounts> yCounts;
  QuaffCountingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, const QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, const char* banner, int verbosity, const char* function, const char* file);
  bool finished() const;
  QuaffCountingTask nextCountingTask();
  QuaffParamCounts finalCounts() const;
  double finalLogLike() const;
};

// thread entry points
void runQuaffCountingTasks (QuaffCountingScheduler* qcs);
void delegateQuaffCountingTasks (QuaffCountingScheduler* qcs, const RemoteServer* remote);

// config/wrapper structs for Viterbi alignment
struct QuaffAlignmentPrinter {
  typedef multiset<Alignment,bool(*)(const Alignment&,const Alignment&)> AlignmentList;
  enum OutputFormat { GappedFastaAlignment, StockholmAlignment, SamAlignment, UngappedFastaRef } format;
  double logOddsThreshold;
  string alignFilename;
  ofstream alignFile;
  
  QuaffAlignmentPrinter();
  bool parseAlignmentPrinterArgs (deque<string>& argvec);
  void writeAlignmentHeader (ostream& out, const vguard<FastSeq>& refs, bool groupByQuery);
  void writeAlignment (ostream& out, const Alignment& align) const;
  bool usingOutputFile() const { return !alignFilename.empty(); }
};

struct QuaffAligner : QuaffAlignmentPrinter {
  bool printAllAlignments;  // print alignment to every reference sequence that matches

  QuaffAligner();
  bool parseAlignmentArgs (deque<string>& argvec);
  void align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
  void serveAlignments (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
  static void serveAlignmentsFromThread (QuaffAligner* paligner, const vguard<FastSeq>* px, const vguard<FastSeq>* py, const QuaffParams* pparams, const QuaffNullParams* pnullModel, const QuaffDPConfig* pconfig, unsigned int port);
};

// structs for scheduling alignment tasks
struct QuaffAlignmentTask : QuaffTask {
  const vguard<FastSeq>& x;
  QuaffAlignmentPrinter::AlignmentList alignList;
  bool keepAllAlignments;
  QuaffAlignmentTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, bool keepAllAlignments);
  void run();
  string delegate (const RemoteServer& remote);
};

struct QuaffAlignmentPrintingScheduler : QuaffScheduler {
  ostream& out;
  QuaffAlignmentPrinter& printer;
  mutex outMx;
  QuaffAlignmentPrintingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file);
  void printAlignments (const QuaffAlignmentPrinter::AlignmentList& alignList);
  void printAlignments (const string& alignStr);
};

struct QuaffAlignmentScheduler : QuaffAlignmentPrintingScheduler {
  const bool keepAllAlignments;
  QuaffAlignmentScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config, bool keepAllAlignments, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file);
  bool finished() const;
  QuaffAlignmentTask nextAlignmentTask();
};

// thread entry point
void runQuaffAlignmentTasks (QuaffAlignmentScheduler* qas);
void delegateQuaffAlignmentTasks (QuaffAlignmentScheduler* qcs, const RemoteServer* remote);

#endif /* QMODEL_INCLUDED */
