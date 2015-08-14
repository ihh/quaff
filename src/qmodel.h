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
  double logQualProb (int k) const;
  double logQualProb (const vguard<double>& kFreq) const;
  void write (ostream& out, const string& prefix) const;
  void read (map<string,double>& paramVal, const string& prefix);
};

// Memo-ized log scores for a SymQualDist
struct SymQualScores {
  vguard<double> logSymQualProb;  // logSymQualProb[q] = P(this symbol wih quality score q)
  SymQualScores() { }
  SymQualScores (const SymQualDist& sqd);
};

// Summary statistics for a SymQualDist
struct SymQualCounts {
  vguard<double> qualCount;  // no. of times each quality score seen
  SymQualCounts();
  double symCount() const { return accumulate (qualCount.begin(), qualCount.end(), 0.); }
  void write (ostream& out, const string& prefix) const;
};

// Parameters of a quaff model
struct QuaffParams {
  double beginInsert, extendInsert, beginDelete, extendDelete;
  vguard<SymQualDist> insert;  // emissions from insert state
  vguard<vguard<SymQualDist> > match;  // substitutions from match state (conditional on input)
  QuaffParams();
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
  void initCounts (double noBeginCount, double yesExtendCount, double matchIdentCount, double otherCount);
  void write (ostream& out) const;
  void addWeighted (const QuaffParamCounts& counts, double weight);
  QuaffParams fit() const;  // maximum-likelihood fit
  double logPrior (const QuaffParams& qp) const;  // uses counts as hyperparameters to define a prior over params
  double expectedLogLike (const QuaffParams& qp) const;  // as logPrior, but unnormalized
};

struct QuaffNullParams {
  double nullEmit;
  vguard<SymQualDist> null;
  QuaffNullParams (const vguard<FastSeq>& seqs, double pseudocount = 1);
  double logLikelihood (const FastSeq& seq) const;
  double logLikelihood (const vguard<FastSeq>& seqs) const;
  void write (ostream& out) const;
};

// Alignment
struct Alignment {
  vguard<FastSeq> gappedSeq;
  double score;
  Alignment (int numRows = 0)
    : gappedSeq(numRows), score(0)
  { }
  const size_t rows() const { return gappedSeq.size(); }
  Alignment& addScoreComment();  // adds score to comment field of first row
  void writeGappedFasta (ostream& out) const;
  void writeStockholm (ostream& out) const;
  FastSeq getUngapped (int row) const;
};

// DP config
struct QuaffDPConfig {
  bool local, sparse;
  int kmerLen, bandSize;
  QuaffDPConfig() : local(true), sparse(true), kmerLen(7), bandSize(20) { }
  bool parseConfigArgs (int& argc, char**& argv);
  DiagonalEnvelope makeEnvelope (const FastSeq& x, const FastSeq& y) const;
};

// DP matrices
struct QuaffDPCell {
  double mat, ins, del;
  QuaffDPCell();
};

typedef map<SeqIdx,QuaffDPCell> QuaffDPColumn;

struct QuaffDPMatrix {
  enum State { Start, Match, Insert, Delete };
  const QuaffDPConfig* pconfig;
  const DiagonalEnvelope* penv;
  const FastSeq *px, *py;
  QuaffScores qs;
  vguard<unsigned int> xTok, yTok, yQual;
  SeqIdx xLen, yLen;
  vguard<QuaffDPColumn> cell;
  double start, end, result;
  static QuaffDPCell dummy;
  vguard<double> cachedInsertEmitScore;
  QuaffDPMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, const QuaffDPConfig& config);
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
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    return qs.match[xTok[i-1]][yTok[j-1]].logSymQualProb[yQual[j-1]];
  }
  inline double insertEmitScore (SeqIdx j) const {
    return qs.insert[yTok[j-1]].logSymQualProb[yQual[j-1]];
  }
  double cellScore (SeqIdx i, SeqIdx j, State state) const;
  static const char* stateToString (State state);
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
  double postMatch (int i, int j) const;
  double postDelete (int i, int j) const;
  double postInsert (int i, int j) const;
  void write (ostream& out) const;
};
  
class QuaffViterbiMatrix : public QuaffDPMatrix {
private:
  const char gapChar = '-';
  static void updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
public:
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
  
  QuaffTrainer();
  bool parseTrainingArgs (int& argc, char**& argv);
  QuaffParams fit (const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& seed, const QuaffParamCounts& pseudocounts, const QuaffDPConfig& config);
};

// config/wrapper struct for Viterbi alignment
struct QuaffAligner {
  enum OutputFormat { GappedFastaAlignment, StockholmAlignment, UngappedFastaRef } format;
  double logOddsThreshold;
  bool printAllAlignments;  // print alignment to every reference sequence that matches
  
  QuaffAligner();
  bool parseAlignmentArgs (int& argc, char**& argv);
  void align (ostream& out, const vguard<FastSeq>& x, const vguard<FastSeq>& y, const QuaffParams& params, const QuaffDPConfig& config);
  void writeAlignment (ostream& out, const Alignment& align) const;
};

#endif /* QMODEL_INCLUDED */
