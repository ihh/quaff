#ifndef QMODEL_INCLUDED
#define QMODEL_INCLUDED

#include <map>
#include <numeric>
#include "fastseq.h"
#include "logger.h"

// struct describing the probability of a given FASTA symbol,
// and a negative binomial distribution over the associated quality score
struct SymQualDist {
  double symProb; // probability of symbol
  double qualTrialSuccessProb, qualNumFailedTrials;  // parameters of neg.binom. distribution
  SymQualDist();
  double logQualProb (int k) const;
  void write (ostream& out, const string& prefix) const;
  void read (map<string,double>& paramVal, const string& prefix);
};

// Memo-ized log scores for a SymQualDist
struct SymQualScores {
  vector<double> logQualProb;  // logQualProb[q] = P(this symbol wih quality score q)
  SymQualScores() { }
  SymQualScores (const SymQualDist& sqd);
};

// Summary statistics for a SymQualDist
struct SymQualCounts {
  vector<double> qualCount;  // no. of times each quality score seen
  SymQualCounts();
  double symCount() const { return accumulate (qualCount.begin(), qualCount.end(), 0.); }
  void write (ostream& out, const string& prefix) const;
};

// Parameters of a quaff model
struct QuaffParams {
  Logger *logger;
  double beginInsert, extendInsert, beginDelete, extendDelete;
  vector<SymQualDist> insert;  // emissions from insert state
  vector<vector<SymQualDist> > match;  // substitutions from match state (conditional on input)
  QuaffParams();
  void write (ostream& out) const;
  void read (istream& in);
};

// Memo-ized scores for transitions & emissions in quaff HMM
struct QuaffScores {
  const QuaffParams *pqp;
  vector<SymQualScores> insert;
  vector<vector<SymQualScores> > match;
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffScores (const QuaffParams& qp);
};

// Summary statistics for a quaff model
struct QuaffCounts {
  vector<SymQualCounts> insert;
  vector<vector<SymQualCounts> > match;
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffCounts();
};

struct QuaffParamCounts {
  vector<SymQualCounts> insert;
  vector<vector<SymQualCounts> > match;
  double beginInsertNo, extendInsertNo, beginDeleteNo, extendDeleteNo;
  double beginInsertYes, extendInsertYes, beginDeleteYes, extendDeleteYes;
  QuaffParamCounts (const QuaffCounts& counts);
  void write (ostream& out) const;
  void addWeighted (const QuaffParamCounts& counts, double weight);
  QuaffParams fit() const;  // maximum-likelihood fit
  double logPrior (const QuaffParams& qp) const;  // uses counts as hyperparameters to define a prior over params
};

// Alignment
struct Alignment {
  vector<FastSeq> gappedSeq;
  Alignment (int numRows = 0)
    : gappedSeq(numRows)
  { }
  void write (ostream& out) const;
  FastSeq getUngapped (int row) const;
};

// DP matrices
struct QuaffDPMatrix {
  enum State { Start, Match, Insert, Delete };
  Logger *logger;
  const FastSeq *px, *py;
  QuaffScores qs;
  vector<int> xTok, yTok, yQual;
  int xLen, yLen;
  vector<vector<double> > mat, ins, del;
  double start, end, result;
  inline double matchEmitScore (int i, int j) const {
    return qs.match[xTok[i-1]][yTok[j-1]].logQualProb[yQual[j-1]];
  }
  inline double insertEmitScore (int j) const {
    return qs.insert[yTok[j-1]].logQualProb[yQual[j-1]];
  }
  QuaffDPMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp);
};

struct QuaffForwardMatrix : QuaffDPMatrix {
  QuaffForwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp);
};

struct QuaffBackwardMatrix : QuaffDPMatrix {
  const QuaffForwardMatrix *pfwd;
  QuaffCounts qc;
  QuaffBackwardMatrix (const QuaffForwardMatrix& fwd);
  double transCount (double& backSrc, double fwdSrc, double trans, double backDest) const;
  inline double& matchCount (int i, int j) {
    return qc.match[xTok[i-1]][yTok[j-1]].qualCount[yQual[j-1]];
  }
  inline double& insertCount (int j) {
    return qc.insert[yTok[j-1]].qualCount[yQual[j-1]];
  }
};

class QuaffViterbiMatrix : public QuaffDPMatrix {
private:
  const char gapChar = '-';
  static void updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
public:
  QuaffViterbiMatrix (const FastSeq& x, const FastSeq& y, const QuaffParams& qp);
  Alignment alignment() const;
};

// Baum-Welch style EM algorithm (also fits quality score distributions)
struct QuaffTrainer {
  int maxIterations;
  double minFractionalLoglikeIncrement;

  bool parseTrainingArgs (int* argcPtr, char*** argvPtr);
  QuaffParams fit (const FastSeq& x, const FastSeq& y, const QuaffParams& seed, const QuaffParamCounts& pseudocounts);
};

#endif /* QMODEL_INCLUDED */
