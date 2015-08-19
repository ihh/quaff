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

struct QuaffOverlapScores {
  const QuaffParams *pqp;
  QuaffMatchKmerContext matchContext;
  QuaffIndelKmerContext indelContext;
  bool yComplemented;
  vguard<SymQualScores> xInsert, yInsert;
  vguard<vguard<SymQualPairScores> > matchMinusInsert;
  vguard<double> gapOpenProb;
  double gapExtendProb, gapAdjacentProb;
  vguard<vguard<double> > m2m, m2i, m2d;
  double i2m, i2i, i2d, d2m, d2i, d2d;
  QuaffOverlapScores (const QuaffParams &qp, bool yComplemented = false);
};

struct QuaffOverlapViterbiMatrix : QuaffDPMatrixContainer {
  QuaffOverlapScores qos;
  vguard<AlphTok> xTok, yTok;
  vguard<Kmer> xMatchKmer, yMatchKmer, xIndelKmer, yIndelKmer;
  vguard<QualScore> xQual, yQual;
  double xInsertScore, yInsertScore;
  QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented);
  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  inline double m2mScore (SeqIdx i, SeqIdx j) const { return i > 0 && j > 0 ? qos.m2m[xIndelKmer[i-1]][yIndelKmer[j-1]] : 0; }
  inline double m2iScore (SeqIdx i, SeqIdx j) const { return i > 0 && j > 0 ? qos.m2i[xIndelKmer[i-1]][yIndelKmer[j-1]] : 0; }
  inline double m2dScore (SeqIdx i, SeqIdx j) const { return i > 0 && j > 0 ? qos.m2d[xIndelKmer[i-1]][yIndelKmer[j-1]] : 0; }
  inline double i2mScore() const { return qos.i2i; }
  inline double i2iScore() const { return qos.i2m; }
  inline double i2dScore() const { return qos.i2d; }
  inline double d2mScore() const { return qos.d2i; }
  inline double d2iScore() const { return qos.d2m; }
  inline double d2dScore() const { return qos.d2d; }
  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    const SymQualPairScores& sqps = qos.matchMinusInsert[xMatchKmer[i-1]][yMatchKmer[j-1]];
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

// struct encapsulating a single overlapper thread
struct QuaffOverlapTask {
  const FastSeq& xfs;
  const FastSeq& yfs;
  const bool yComplemented;
  const QuaffParams& params;
  const QuaffNullParams& nullModel;
  const QuaffDPConfig& config;
  Alignment align;

  QuaffOverlapTask (const FastSeq& xfs, const FastSeq& yfs, const bool yComplemented, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config);
  void run();
  bool hasAlignment() const { return align.rows() > 0; }
};
void runQuaffOverlapTask (QuaffOverlapTask& task);

#endif /* QOVERLAP_INCLUDED */
