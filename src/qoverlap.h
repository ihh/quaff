#ifndef QOVERLAP_INCLUDED
#define QOVERLAP_INCLUDED

#include <tuple>
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
  // xIndelKmer and yIndelKmer are padded with a dummy entry at the start
  // this avoids the need for a bounds test in m2*Score(i,j) methods
  vguard<Kmer> xMatchKmer, yMatchKmer, xIndelKmer, yIndelKmer;
  vguard<QualScore> xQual, yQual;
  double xInsertScore, yInsertScore;
  QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented);
  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  inline double m2mScore (SeqIdx i, SeqIdx j) const { return qos.m2m[xIndelKmer[i]][yIndelKmer[j]]; }
  inline double m2iScore (SeqIdx i, SeqIdx j) const { return qos.m2i[xIndelKmer[i]][yIndelKmer[j]]; }
  inline double m2dScore (SeqIdx i, SeqIdx j) const { return qos.m2d[xIndelKmer[i]][yIndelKmer[j]]; }
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
  void align (ostream& out, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  void serveAlignments (const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  static void serveAlignmentsFromThread (QuaffOverlapAligner* paligner, const vguard<FastSeq>* pseqs, size_t nOriginals, const QuaffParams* pparams, const QuaffNullParams* pnullModel, QuaffDPConfig* pconfig, unsigned int port);
  static string getAlignmentJobResult (QuaffOverlapAligner& aligner, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, map<string,size_t>& yDict, const JsonMap& jobDescription, bool& validJson);
};

// structs for scheduling overlap tasks
struct QuaffOverlapTask : QuaffTask {
  const FastSeq& xfs;
  const bool yComplemented;
  QuaffAlignmentPrinter::AlignmentList alignList;
  QuaffOverlapTask (const FastSeq& xfs, const FastSeq& yfs, const bool yComplemented, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config);
  void run();
  bool remoteRun (RemoteServer& remote, string& align);
  void qsubRun (size_t taskId);
  string toJson() const;
};

class QuaffOverlapScheduler : public QuaffAlignmentPrintingScheduler {
private:
  void advance();
public:
  size_t nx, nOriginals;
  deque<tuple<const FastSeq*,const FastSeq*,bool> > failed;
  QuaffOverlapScheduler (const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, QuaffDPConfig& config, ostream& out, QuaffAlignmentPrinter& printer, int verbosity, const char* function, const char* file, int line);
  bool noMoreTasks() const;
  bool finished() const;
  QuaffOverlapTask nextOverlapTask();
  void rescheduleOverlapTask (const QuaffOverlapTask& task);
};

// thread entry point
void runQuaffOverlapTasks (QuaffOverlapScheduler* qos);
void remoteRunQuaffOverlapTasks (QuaffOverlapScheduler* qos, RemoteServer* remote);
void qsubRunQuaffOverlapTasks (QuaffOverlapScheduler* qos);

#endif /* QOVERLAP_INCLUDED */
