#ifndef QMODEL_INCLUDED
#define QMODEL_INCLUDED

#include <map>
#include <numeric>
#include "fastseq.h"
#include "diagenv.h"

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

// Parameters of a quaff model
struct QuaffParams {
  double refEmit;
  vguard<double> refBase;
  double beginInsert, extendInsert, beginDelete, extendDelete;
  vguard<SymQualDist> insert;  // emissions from insert state
  vguard<vguard<SymQualDist> > match;  // substitutions from match state (conditional on input)
  QuaffParams();
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
  void write (ostream& out) const;
  void read (istream& in);
};

// Memo-ized scores for transitions & emissions in quaff HMM
struct QuaffScores {
  const QuaffParams *pqp;
  vguard<SymQualScores> insert;
  vguard<vguard<SymQualScores> > match;
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffScores (const QuaffParams& qp);
};

// Summary statistics for a quaff model
struct QuaffCounts {
  vguard<SymQualCounts> insert;
  vguard<vguard<SymQualCounts> > match;
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffCounts();
  void write (ostream& out) const;
};

struct QuaffParamCounts {
  vguard<SymQualCounts> insert;
  vguard<vguard<SymQualCounts> > match;
  double beginInsertNo, extendInsertNo, beginDeleteNo, extendDeleteNo;
  double beginInsertYes, extendInsertYes, beginDeleteYes, extendDeleteYes;
  QuaffParamCounts();
  QuaffParamCounts (const QuaffCounts& counts);
  void zeroCounts();
  void initCounts (double noBeginCount, double yesExtendCount, double matchIdentCount, double otherCount, const QuaffNullParams* nullModel = NULL);
  void write (ostream& out) const;
  void read (istream& in);
  void addWeighted (const QuaffParamCounts& counts, double weight);
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
    : gappedSeq(numRows), score(0)
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
  bool local, sparse;
  int kmerLen, kmerThreshold, bandSize;
  double kmerStDevThreshold;
  QuaffDPConfig()
    : local(true),
      sparse(true),
      kmerLen(DEFAULT_KMER_LENGTH),
      kmerThreshold(DEFAULT_KMER_THRESHOLD),
      kmerStDevThreshold(DEFAULT_KMER_STDEV_THRESHOLD),
      bandSize(DEFAULT_BAND_SIZE)
  { }
  bool parseConfigArgs (int& argc, char**& argv);
  bool parseOverlapConfigArgs (int& argc, char**& argv);
  DiagonalEnvelope makeEnvelope (const FastSeq& x, const FastSeq& y) const;
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
protected:
  static void updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
};

struct QuaffDPMatrix : QuaffDPMatrixContainer {
  const QuaffDPConfig* pconfig;
  QuaffScores qs;
  vguard<AlphTok> xTok, yTok;
  vguard<QualScore> yQual;
  vguard<double> cachedInsertEmitScore;
  QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    const SymQualScores& sqs = qs.match[xTok[i-1]][yTok[j-1]];
    return yQual.size() ? sqs.logSymQualProb[yQual[j-1]] : sqs.logSymProb;
  }
  inline double insertEmitScore (SeqIdx j) const {
    const SymQualScores& sqs = qs.insert[yTok[j-1]];
    return yQual.size() ? sqs.logSymQualProb[yQual[j-1]] : sqs.logSymProb;
  }
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
  inline double& matchCount (SeqIdx i, SeqIdx j) {
    return qc.match[xTok[i-1]][yTok[j-1]].qualCount[yQual[j-1]];
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
  bool parseTrainingArgs (int& argc, char**& argv);
  QuaffParams fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffNullParams& nullModel, const QuaffParamCounts& pseudocounts, const QuaffDPConfig& config);
};

// config/wrapper structs for Viterbi alignment
struct QuaffAlignmentPrinter {
  enum OutputFormat { GappedFastaAlignment, StockholmAlignment, UngappedFastaRef } format;
  double logOddsThreshold;

  QuaffAlignmentPrinter();
  bool parseAlignmentPrinterArgs (int& argc, char**& argv);
  void writeAlignment (ostream& out, const Alignment& align) const;
};

struct QuaffAligner : QuaffAlignmentPrinter {
  bool printAllAlignments;  // print alignment to every reference sequence that matches

  QuaffAligner();
  bool parseAlignmentArgs (int& argc, char**& argv);
  void align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
};

#endif /* QMODEL_INCLUDED */
