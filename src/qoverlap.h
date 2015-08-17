#ifndef QOVERLAP_INCLUDED
#define QOVERLAP_INCLUDED

#include "qmodel.h"

struct SymQualPairScores {
  vguard<vguard<double> > logSymQualPairProb;  // if x and y both have quality scores
  vguard<double> logSymPairXQualProb, logSymPairYQualProb;  // if only x or y has quality scores
  double logSymPairProb;  // if neither has quality scores
  SymQualPairScores()
    : logSymQualPairProb (FastSeq::qualScoreRange, vguard<double> (FastSeq::qualScoreRange)),
      logSymPairXQualProb (FastSeq::qualScoreRange, -numeric_limits<double>::infinity()),
      logSymPairYQualProb (FastSeq::qualScoreRange, -numeric_limits<double>::infinity()),
      logSymPairProb (-numeric_limits<double>::infinity())
  { }
};

struct QuaffOverlapScores : QuaffKmerContext {
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
  vguard<AlphTok> xTok, yTok;
  vguard<Kmer> xKmer, yKmer;
  vguard<QualScore> xQual, yQual;
  double xInsertScore, yInsertScore;
  QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented);
  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    const SymQualPairScores& sqps = qos.matchMinusInsert[xKmer[i-1]][yKmer[j-1]];
    return xQual.size()
      ? (yQual.size()
	 ? sqps.logSymQualPairProb[xQual[i-1]][yQual[j-1]]
	 : sqps.logSymPairXQualProb[xQual[i-1]])
      : (yQual.size()
	 ? sqps.logSymPairYQualProb[yQual[j-1]]
	 : sqps.logSymPairProb);
  }
  Alignment alignment() const;
  Alignment scoreAdjustedAlignment (const QuaffNullParams& nullModel) const;
};

// config/wrapper struct for Viterbi alignment
struct QuaffOverlapAligner : QuaffAlignmentPrinter {
  QuaffOverlapAligner();
  bool parseAlignmentArgs (deque<string>& argvec);
  void align (ostream& out, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
};

#endif /* QOVERLAP_INCLUDED */
