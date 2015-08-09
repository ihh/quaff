#ifndef QUAFF_INCLUDED
#define QUAFF_INCLUDED

#include "fastseq.h"

// struct describing the probability of a given FASTA symbol,
// and the (negative binomial) probability distribution over the associated quality score
struct SymQualDist {
  double symProb; // probability of symbol
  double qualTrialProb, qualNumTrials;  // resp.: {p,r} params of neg.binom. distrib.
};

// Memo-ized log scores for a SymQualDist
struct SymQualScores {
  double symLogProb; // log-probability of symbol
  vector<double> logQualProb;  // log-probability distribution over qual scores
  SymQualScores (const SymQualDist& sqd);
};

// Fitting the negative binomial distribution requires taking some unusual moments,
// e.g.  \sum_i digamma(r+y_i)  as per http://vixra.org/pdf/1211.0113v1.pdf
// Since a manageably small range of quality scores is encountered in practice
// (and because there are several alternative heuristics for fitting negative binomials),
// rather than trying to deal with these moments inside the DP recursion,
// we maintain a count for every observed quality score.
struct SymQualCounts {
  double symCount;
  vector<double> qualCount;
};

// Parameters of a quaff model
struct QuaffParams {
  double beginInsert, extendInsert, beginDelete, extendDelete;
  vector<SymQualDist> insert;  // emissions from insert state
  vector<vector<SymQualDist> > match;  // substitutions from match state (conditional on input)
  QuaffParams();
  void write (ostream& out) const;
  void read (istream& in);
};

// Memo-ized scores for transitions & emissions in quaff HMM
struct QuaffScores {
  const QuaffParams *qp;
  vector<SymQualScores> insert;
  vector<vector<SymQualScores> > match;
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffScores (const QuaffParams& params);
};

// Summary statistics for a quaff model
struct QuaffCounts {
  double beginInsertYes, extendInsertYes, beginDeleteYes, extendDeleteYes;
  double beginInsertNo, extendInsertNo, beginDeleteNo, extendDeleteNo;
  vector<SymQualCounts> insert;
  vector<vector<SymQualCounts> > match;
  QuaffCounts();
  void write (ostream& out) const;
  void read (istream& in);
  void addWeighted (const QuaffCounts& counts, double weight);
  QuaffParams fit() const;  // maximum-likelihood fit
};

// Alignment
struct Alignment {
  vector<string> name, gappedSeq;
  Alignment (int numRows = 0)
    : name(numRows), gappedSeq(numRows)
  { }
  void write (ostream& out) const;
  FastSeq getUngappedSeq (int row) const;
};

// DP matrices
struct QuaffDPMatrix {
  const FastSeq *x, *y;
  const QuaffScores *qs;
  vector<vector<double> > mat, ins, del;
  double start, end, result;
  QuaffDPMatrix (const FastSeq& x, const FastSeq& y, const QuaffScores& qs);
};

struct QuaffForwardMatrix : QuaffDPMatrix {
  QuaffForwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffScores& qs);
};

struct QuaffBackwardMatrix : QuaffDPMatrix {
  QuaffBackwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffScores& qs, QuaffCounts& counts);
};

struct QuaffForwardBackwardMatrix {
  QuaffCounts counts;
  QuaffForwardMatrix fwd;
  QuaffBackwardMatrix back;
  QuaffForwardBackwardMatrix (const FastSeq& x, const FastSeq& y, const QuaffScores& qs);
};

struct QuaffViterbiMatrix : QuaffDPMatrix {
  Alignment align;
  QuaffViterbiMatrix (const FastSeq& x, const FastSeq& y, const QuaffScores& qs);
};

// Baum-Welch style EM algorithm (also fitting negative binomial distributions to quality scores)
struct QuaffTrainer {
  int maxIterations;
  double minFractionalLoglikeIncrement;

  QuaffParams fit (const FastSeq& x, const FastSeq& y, const QuaffParams& seed);
};

#endif /* QUAFF_INCLUDED */
