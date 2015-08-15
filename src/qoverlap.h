#ifndef QOVERLAP_INCLUDED
#define QOVERLAP_INCLUDED

#include "qmodel.h"

struct SymQualPairScores {
  vguard<vguard<double> > logSymQualPairProb;
  SymQualPairScores()
    : logSymQualPairProb (FastSeq::qualScoreRange, vguard<double> (FastSeq::qualScoreRange))
  { }
};

struct QuaffOverlapScores {
  const QuaffParams *pqp;
  bool yComplemented;
  vguard<SymQualScores> xInsert, yInsert;
  vguard<vguard<SymQualPairScores> > matchMinusInsert;
  double gapOpenProb, gapExtendProb, gapAdjacentProb;
  double m2m, m2i, m2d, i2m, i2i, i2d, d2m, d2i, d2d;
  QuaffOverlapScores (const QuaffParams &qp, bool yComplemented = false);
};

struct QuaffOverlapViterbiMatrix : QuaffDPMatrixContainer {
  QuaffOverlapScores qos;
  vguard<unsigned int> xTok, yTok, xQual, yQual;
  double xInsertScore, yInsertScore;
  QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented);
  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    return qos.matchMinusInsert[xTok[i-1]][yTok[j-1]].logSymQualPairProb[xQual[i-1]][yQual[j-1]];
  }
  Alignment alignment() const;
  Alignment scoreAdjustedAlignment (const QuaffNullParams& nullModel) const;
};

// config/wrapper struct for Viterbi alignment
struct QuaffOverlapAligner : QuaffAlignmentPrinter {
  QuaffOverlapAligner();
  bool parseAlignmentArgs (int& argc, char**& argv);
  void align (ostream& out, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
};

#endif /* QOVERLAP_INCLUDED */
