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
#include "aws.h"
#include "jsonutil.h"

// EM convergence parameters
#define QuaffMaxEMIterations 100
#define QuaffMinEMLogLikeInc .01

// Default contexts
#define DefaultMatchKmerContext 1
#define DefaultIndelKmerContext 0

// Default server port
#define DefaultServerPort 8000

// Default place to find quaff binary
#define DefaultQuaffInstallPrefix "/usr/local"
#define DefaultQuaffInstallDir DefaultQuaffInstallPrefix "/quaff"
#define DefaultQuaffPath DefaultQuaffInstallPrefix "/bin/quaff"

// Quaff git repository
#define QuaffGitRepository "https://github.com/ihh/quaff.git"

// Server directory where 'rsync' & 'aws sync' will stash files
#define SyncStagingDir "/tmp/quaff"

// Directory whose creation signals that server startup is finished
#define ServerReadyDir SyncStagingDir "/.ready"

// Maximum number of times a general ssh job will be tried
// Applies to all commands executed by ssh (including e.g. mkdir's during setup)
// Does not apply to the actual quaff job
#define MaxGenericSshAttempts 10

// Maximum number of times a quaff ssh job will be tried
#define MaxQuaffSshAttempts 20

// If an ssh job emits this string, it is deemed to be a partial success;
// if it then fails, the server will be rebooted (or there will be a delay of QuaffTimeWaitDelay),
// and the attempt count will be reset
#define QuaffSshReadyString "# READY"

// Server keep-alive parameters for ssh
#define QuaffServerAliveInterval 15
#define QuaffServerAliveCountMax 3

// Number of consecutive failures to connect before a client quits to avoid blocking
#define MaxQuaffClientFailures 30

// Retry delay
#define MinQuaffRetryDelay 10

// Exponential backoff for retry delay
#define QuaffRetryDelayMultiplier 1.15

// Extra delay before attempting to re-run a quaff job (to account for sockets tied up in TIME_WAIT)
// If it's an EC2 server, we don't bother with this, just reboot the server instead...
#define QuaffTimeWaitDelay 300

// qsub file suffices
#define QuaffQsubScriptSuffix ".sh"
#define QuaffQsubInfoSuffix   ".info"
#define QuaffQsubResultSuffix ".result"
#define QuaffQsubDoneSuffix   ".done"
#define QuaffQsubOutSuffix    ".out"
#define QuaffQsubErrSuffix    ".err"

// max qsub attempts
#define MaxQuaffQsubAttempts 3

// useful helper methods
void randomDelayBeforeRetry (int attempts = 0, double minSeconds = MinQuaffRetryDelay, double multiplier = QuaffRetryDelayMultiplier);

// struct describing the probability of a given FASTA symbol,
// and a negative binomial distribution over the associated quality score
struct SymQualDist {
  double symProb; // probability of symbol
  double qualTrialSuccessProb, qualNumSuccessfulTrials;  // parameters of neg.binom. distribution
  SymQualDist();
  double logQualProb (QualScore k) const;
  double logQualProb (const vguard<double>& kFreq) const;
  ostream& writeJson (ostream& out) const;
  bool readJson (const JsonValue& val);
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
  ostream& writeJson (ostream& out) const;
  bool readJson (const JsonValue& val);
};

// Classes to manage kmer-dependence of various parameters
struct QuaffKmerContext {
  void initKmerContext (unsigned int newKmerLen);
  const char* prefix;
  unsigned int kmerLen, defaultKmerLen;
  Kmer numKmers;
  QuaffKmerContext (const char* prefix, unsigned int kmerLen, unsigned int defaultKmerLen);
  string kmerString (Kmer kmer) const;
  string insertParamName (AlphTok i) const;
  void writeJsonKmerLen (ostream& out) const;
  void readJsonKmerLen (const JsonMap& m);
};

struct QuaffMatchKmerContext : QuaffKmerContext {
  QuaffMatchKmerContext (unsigned int kmerLen = DefaultMatchKmerContext)
    : QuaffKmerContext ("match", kmerLen, DefaultMatchKmerContext)
  { }
  string matchParamName (AlphTok i, Kmer j) const;
  string kmerPrefix (Kmer j) const;
  char kmerSuffix (Kmer j) const;
};

struct QuaffIndelKmerContext : QuaffKmerContext {
  QuaffIndelKmerContext (unsigned int kmerLen = DefaultIndelKmerContext)
    : QuaffKmerContext ("gap", kmerLen, DefaultIndelKmerContext)
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
  void writeToLog (int v = 0) const;
  ostream& writeJson (ostream& out) const;
  bool readJson (const JsonValue& val);
  bool readJson (const JsonMap& jm);
  void readJson (istream& in);
  void fitRefSeqs (const vguard<FastSeq>& refs);
};

struct QuaffNullParams {
  double nullEmit;
  vguard<SymQualDist> null;
  QuaffNullParams();
  QuaffNullParams (const vguard<FastSeq>& seqs, double pseudocount = 1);
  double logLikelihood (const FastSeq& seq) const;
  double logLikelihood (const vguard<FastSeq>& seqs) const;
  void writeToLog (int v = 0) const;
  ostream& writeJson (ostream& out) const;
  void readJson (istream& in);
  bool readJson (const JsonMap& jm);
  bool readJson (const JsonValue& val);
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
struct QuaffEmitCounts {
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  vguard<SymQualCounts> insert;
  vguard<vguard<SymQualCounts> > match;
  QuaffEmitCounts (unsigned int matchKmerLen, unsigned int indelKmerLen);
  QuaffEmitCounts (const QuaffEmitCounts& c);
  void resize();  // call after changing kmerLen
  ostream& writeJson (ostream& out) const;
};

struct QuaffCounts : QuaffEmitCounts {
  vguard<double> m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffCounts (unsigned int matchKmerLen, unsigned int indelKmerLen);
  void writeToLog (int v = 0) const;
  ostream& writeJson (ostream& out) const;
};

struct QuaffParamCounts : QuaffEmitCounts {
  vguard<double> beginInsertNo, beginInsertYes, beginDeleteNo, beginDeleteYes;
  double extendInsertNo, extendInsertYes, extendDeleteNo, extendDeleteYes;
  QuaffParamCounts (unsigned int matchKmerLen = DefaultMatchKmerContext, unsigned int indelKmerLen = DefaultIndelKmerContext);
  QuaffParamCounts (const QuaffCounts& counts);
  void resize();  // call after changing kmerLen
  void zeroCounts();
  void initCounts (double noBeginCount, double yesExtendCount, double matchIdentCount, double otherCount, const QuaffNullParams* nullModel = NULL);
  void writeToLog (int v = 0) const;
  ostream& writeJson (ostream& out) const;
  ostream& writeJsonWithMeta (ostream& out, const string& name, const vguard<size_t>& sortOrder, const double logLike) const;
  void readJson (istream& in);
  bool readJson (const JsonMap& jm);
  bool readJson (const JsonValue& val);
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

// structs describing remote servers & jobs
struct RemoteServerJob {
  string addr, user, ec2instanceId;
  unsigned int port, threads;
  bool ready;
  RemoteServerJob (const string& user, const string& addr, unsigned int port, unsigned int threads, const string& ec2id);
  string toString() const;
};

class RemoteServer {
private:
  TCPSocket* sock;
  const bool* ready;
public:
  string addr;
  unsigned int port;
  RemoteServer (const RemoteServerJob& rsj, unsigned int port);
  ~RemoteServer() { closeSocket(); }
  string toString() const;
  TCPSocket* getSocket();
  void closeSocket();
};

// DP config
struct QuaffDPConfig {
  bool local, sparse, autoMemSize;
  int kmerLen, kmerThreshold, bandSize;
  size_t maxSize;
  unsigned int threads;
  bool threadsSpecified;
  int serverPort;
  list<RemoteServer> remotes;
  list<RemoteServerJob> remoteJobs;
  list<thread> remoteServerThreads;
  string bucket, sshPath, rsyncPath, sshKey, remoteQuaffPath;
  bool useRsync;
  string remoteServerArgs;
  vguard<tuple<string,string,string> > fileArgs;
  vguard<tuple<string,string,string> > nonReadFileArgs;  // excludes all (read) seq files that may be split by parallelization
  unsigned int ec2Instances, ec2Cores, ec2Port;
  vguard<string> ec2InstanceIds, ec2InstanceAddresses;
  string ec2Ami, ec2Type, ec2User;
  unsigned int qsubThreads;
  string qsubHeader, qsubPath, qsubOpts;
  TempDir qsubTempDir;
  string jobFilename;
  QuaffDPConfig()
    : local(true),
      sparse(true),
      autoMemSize(false),
      kmerLen(DEFAULT_KMER_LENGTH),
      kmerThreshold(DEFAULT_KMER_THRESHOLD),
      maxSize(0),
      bandSize(DEFAULT_BAND_SIZE),
      threads(1),
      threadsSpecified(false),
      serverPort(DefaultServerPort),
      sshPath("ssh"),
      rsyncPath("rsync"),
      useRsync(false),
      remoteQuaffPath(DefaultQuaffPath),
      ec2Instances(0),
      ec2Ami(AWS_DEFAULT_AMI),
      ec2Type(AWS_DEFAULT_INSTANCE_TYPE),
      ec2Cores(AWS_DEFAULT_INSTANCE_CORES),
      ec2User(AWS_DEFAULT_USER),
      ec2Port(DefaultServerPort),
      qsubThreads(0),
      qsubPath("qsub"),
      qsubHeader("#!/bin/sh\n")
  { }
  bool parseRefSeqConfigArgs (deque<string>& argvec);
  bool parseGeneralConfigArgs (deque<string>& argvec);
  bool parseServerConfigArgs (deque<string>& argvec);
  void setServerArgs (const char* serverType, const string& args);
  void addFileArg (const char* tag, const string& filename, const char* extra = NULL);
  void addReadFileArg (const char* tag, const string& filename, const char* extra = NULL);
  DiagonalEnvelope makeEnvelope (const FastSeq& x, const KmerIndex& yKmerIndex, size_t cellSize) const;
  size_t effectiveMaxSize() const;  // takes threading into account
  void syncFromBucket (const string& filename) const;
  void syncToBucket (const string& filename) const;
  void makeStagingDir (const RemoteServerJob& remote) const;
  void syncToRemote (const string& filename, const RemoteServerJob& remote) const;
  void addRemote (const string& user, const string& addr, unsigned int port, unsigned int threads, const string& ec2id = string());
  void startRemoteServers();
  void stopRemoteServers();
  string ec2StartupScript() const;
  string makeServerArgs() const;
  string makeServerCommand (const RemoteServerJob& job) const;
  string makeSshCommand() const;
  string makeSshCommand (const string& cmd, const RemoteServerJob& job) const;
  bool execWithRetries (const string& cmd, int maxAttempts, bool* foundReadyFlag = NULL, const char* ec2id = NULL) const;
  string makeQsubScript (const string& prefix, const string& args) const;
  string makeReadIndexOpt (const FastSeq& read) const;
  string makeQsubCommand (const string& scriptFilename) const;
  bool atLeastOneThread() const;
};

// thread entry point
void startRemoteQuaffServer (const QuaffDPConfig* config, RemoteServerJob* remoteJob);

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
  void writeToLog (int v = 0) const;
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
  void writeToLog (int v = 0) const;
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
  SeqIdx maxReadBases;
  bool allowNullModel;
  string rawCountsFilename, countsWithPriorFilename, saveParamsFilename;
  
  QuaffTrainer();
  bool parseTrainingArgs (deque<string>& argvec);
  bool parseCountingArgs (deque<string>& argvec);
  bool parseServerArgs (deque<string>& argvec);
  string serverArgs() const;
  QuaffParams fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, QuaffDPConfig& config);
  QuaffParams fitUnlimited (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, QuaffDPConfig& config);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, double& logLike, const char* banner);
  QuaffParamCounts getCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  void serveCounts (const vguard<FastSeq>& x, const vguard<FastSeq>& y, QuaffDPConfig& config);
  static void serveCountsFromThread (const vguard<FastSeq>* px, const vguard<FastSeq>* py, bool useNullModel, QuaffDPConfig* pconfig, unsigned int port);
  static string getCountJobResult (const vguard<FastSeq>& x, const vguard<FastSeq>& y, bool useNullModel, QuaffDPConfig& config, map<string,size_t>& yDict, const JsonMap& jobDescription, bool& validJson);
  static vguard<vguard<size_t> > defaultSortOrder (const vguard<FastSeq>& x, const vguard<FastSeq>& y);
  bool usingParamOutputFile() const { return !saveParamsFilename.empty(); }
  bool usingCountsOutputFile() const { return !rawCountsFilename.empty(); }
};

// structs for scheduling training tasks
struct QuaffTask {
  const FastSeq& yfs;
  const QuaffParams& params;
  QuaffDPConfig& config;
  const QuaffNullParams& nullModel;
  QuaffTask (const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config)
    : yfs(yfs), params(params), nullModel(nullModel), config(config)
  { }
  string qsubResult (size_t taskId, const string& jobDescription, const string& quaffArgs);
};

struct QuaffCountingTask : QuaffTask {
  const vguard<FastSeq>& x;
  vguard<size_t>& sortOrder;
  double& yLogLike;
  QuaffParamCounts& yCounts;
  const bool useNullModel;
  const size_t ny;
  QuaffCountingTask (const vguard<FastSeq>& x, const FastSeq& yfs, size_t ny, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, QuaffDPConfig& config, vguard<size_t>& sortOrder, double& yLogLike, QuaffParamCounts& counts);
  void run();
  bool remoteRun (RemoteServer& remote);
  void qsubRun (size_t taskId);
  string toJson() const;
  bool getResultFromJson (const JsonMap& result);
};

struct QuaffScheduler {
  const vguard<FastSeq>& x;
  const vguard<FastSeq>& y;
  const QuaffParams& params;
  QuaffDPConfig& config;
  const QuaffNullParams& nullModel;
  size_t ny, pending, lockCount;  // lockCount is used by queue runners as a unique temp filename for tasks
  mutex mx;
  ProgressLogger plog;
  QuaffScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, int verbosity, const char* function, const char* file, int line)
    : x(x), y(y), params(params), nullModel(nullModel), config(config), ny(0), pending(0), lockCount(0), plog(verbosity,function,file,line)
  { }
  void lock() { mx.lock(); ++lockCount; }
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
  deque<size_t> failed;
  QuaffCountingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, bool useNullModel, QuaffDPConfig& config, vguard<vguard<size_t> >& sortOrder, const char* banner, int verbosity, const char* function, const char* file, int line);
  bool noMoreTasks() const;
  bool finished() const;
  QuaffCountingTask nextCountingTask();
  void rescheduleCountingTask (const QuaffCountingTask& task);
  QuaffParamCounts finalCounts() const;
  double finalLogLike() const;
};

// thread entry points
void runQuaffCountingTasks (QuaffCountingScheduler* qcs);
void remoteRunQuaffCountingTasks (QuaffCountingScheduler* qcs, RemoteServer* remote);
void qsubRunQuaffCountingTasks (QuaffCountingScheduler* qcs);

// config/wrapper structs for Viterbi alignment
struct QuaffAlignmentPrinter {
  typedef multiset<Alignment,bool(*)(const Alignment&,const Alignment&)> AlignmentList;
  enum OutputFormat { GappedFastaAlignment, StockholmAlignment, SamAlignment, UngappedFastaRef } format;
  double logOddsThreshold;
  string alignFilename;
  ofstream alignFile;
  
  QuaffAlignmentPrinter();
  bool parseAlignmentPrinterArgs (deque<string>& argvec);
  string serverArgs() const;
  void writeAlignmentHeader (ostream& out, const vguard<FastSeq>& refs, bool groupByQuery);
  void writeAlignment (ostream& out, const Alignment& align) const;
  bool usingOutputFile() const { return !alignFilename.empty(); }
  ostream& alignmentOutputStream (ostream& defaultOut) const { return usingOutputFile() ? (ostream&) alignFile : defaultOut; }
};

struct QuaffAligner : QuaffAlignmentPrinter {
  bool printAllAlignments;  // print alignment to every reference sequence that matches

  QuaffAligner();
  bool parseAlignmentArgs (deque<string>& argvec);
  string serverArgs() const;
  void align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  void serveAlignments (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  static void serveAlignmentsFromThread (QuaffAligner* paligner, const vguard<FastSeq>* px, const vguard<FastSeq>* py, const QuaffParams* pparams, const QuaffNullParams* pnullModel, QuaffDPConfig* pconfig, unsigned int port);
  static string getAlignmentJobResult (QuaffAligner& aligner, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, map<string,size_t>& yDict, const JsonMap& jobDescription, bool& validJson);
};

// structs for scheduling alignment tasks
struct QuaffAlignmentTask : QuaffTask {
  const vguard<FastSeq>& x;
  QuaffAlignmentPrinter::AlignmentList alignList;
  bool keepAllAlignments;
  QuaffAlignmentTask (const vguard<FastSeq>& x, const FastSeq& yfs, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, bool keepAllAlignments);
  void run();
  bool remoteRun (RemoteServer& remote, string& align);
  void qsubRun (size_t taskId, string& align);
  string toJson() const;
};

struct QuaffAlignmentPrintingScheduler : QuaffScheduler {
  ostream& out;
  QuaffAlignmentPrinter& printer;
  mutex outMx;
  QuaffAlignmentPrintingScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file, int line);
  void printAlignments (const QuaffAlignmentPrinter::AlignmentList& alignList);
  void printAlignments (const string& alignStr);
};

struct QuaffAlignmentScheduler : QuaffAlignmentPrintingScheduler {
  const bool keepAllAlignments;
  deque<const FastSeq*> failed;
  QuaffAlignmentScheduler (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, bool keepAllAlignments, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file, int line);
  bool noMoreTasks() const;
  bool finished() const;
  QuaffAlignmentTask nextAlignmentTask();
  void rescheduleAlignmentTask (const QuaffAlignmentTask& task);
};

// thread entry point
void runQuaffAlignmentTasks (QuaffAlignmentScheduler* qas);
void remoteRunQuaffAlignmentTasks (QuaffAlignmentScheduler* qcs, RemoteServer* remote);
void qsubRunQuaffAlignmentTasks (QuaffAlignmentScheduler* qcs);

#endif /* QMODEL_INCLUDED */
