#include <list>
#include <cmath>
#include "qoverlap.h"
#include "logsumexp.h"
#include "logger.h"

QuaffOverlapScores::QuaffOverlapScores (const QuaffParams &qp, bool yComplemented)
  : pqp(&qp),
    yComplemented(yComplemented),
    xInsert (dnaAlphabetSize),
    yInsert (dnaAlphabetSize),
    matchMinusInsert (dnaAlphabetSize, vguard<SymQualPairScores> (dnaAlphabetSize))
{
  // coarse approximation to various paths through the composite transducer...
  const double readInsertProb = qp.beginInsert;
  const double readDeleteProb = (1 - qp.beginInsert) * qp.beginDelete;
  gapOpenProb = readInsertProb + readDeleteProb;
  const double pGapIsInsert = readInsertProb / gapOpenProb;
  const double meanGapLength = pGapIsInsert / qp.extendInsert + (1 - pGapIsInsert) / (qp.refSeqEmit * qp.extendDelete * (1 - gapOpenProb));
  gapExtendProb = 1 / meanGapLength;
  gapAdjacentProb = pGapIsInsert * readInsertProb + (1 - pGapIsInsert) * gapOpenProb / (1 - qp.extendDelete * (1 - gapOpenProb));
  
  m2m = log(qp.refSeqEmit) + 2*log(1-gapOpenProb);
  m2i = m2d = log(gapOpenProb);
  i2i = d2d = log(gapExtendProb);
  i2d = d2i = log(1-gapExtendProb) + log(gapAdjacentProb);
  i2m = d2m = log(qp.refSeqEmit) + log(1-gapExtendProb) + log(1-gapAdjacentProb);

  const QuaffScores qsc (qp);
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (QualScore k = 0; k < FastSeq::qualScoreRange; ++k) {
      double ll = log(pGapIsInsert) + qsc.insert[i].logSymQualProb[k];
      for (AlphTok r = 0; r < dnaAlphabetSize; ++r)
	ll = log_sum_exp (ll, log(1-pGapIsInsert) + log(qp.refBase[r]) + qsc.match[r][i].logSymQualProb[k]);
      xInsert[i].logSymQualProb[k] = ll;
    }
  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    yInsert[i] = xInsert[yComplemented ? dnaComplement(i) : i];

  for (AlphTok i = 0; i < dnaAlphabetSize; ++i)
    for (AlphTok j = 0; j < dnaAlphabetSize; ++j) {
      const AlphTok yStrand_j = yComplemented ? dnaComplement(j) : j;
      for (QualScore ik = 0; ik < FastSeq::qualScoreRange; ++ik) {
	for (QualScore jk = 0; jk < FastSeq::qualScoreRange; ++jk) {
	  double mij = -numeric_limits<double>::infinity();
	  for (AlphTok r = 0; r < dnaAlphabetSize; ++r) {
	    const AlphTok yStrand_r = yComplemented ? dnaComplement(r) : r;
	    mij = log_sum_exp (mij, log(qp.refBase[r]) + qsc.match[r][i].logSymQualProb[ik] + qsc.match[yStrand_r][yStrand_j].logSymQualProb[jk]);
	  }
	  matchMinusInsert[i][j].logSymQualPairProb[ik][jk] = mij - xInsert[i].logSymQualProb[ik] - yInsert[j].logSymQualProb[jk];
	}
      }
    }
}

QuaffOverlapViterbiMatrix::QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented)
  : QuaffDPMatrixContainer (env),
    qos (qp, yComplemented),
    xTok (px->tokens(dnaAlphabet)),
    yTok (py->tokens(dnaAlphabet)),
    xQual (px->qualScores()),
    yQual (py->qualScores()),
    xInsertScore (0),
    yInsertScore (0)
{
  const FastSeq& x (*px);
  const FastSeq& y (*py);

  for (SeqIdx i = 0; i < xLen; ++i)
    xInsertScore += qos.xInsert[xTok[i]].logSymQualProb[xQual[i]];

  for (SeqIdx i = 0; i < yLen; ++i)
    yInsertScore += qos.yInsert[yTok[i]].logSymQualProb[yQual[i]];

  if (LogThisAt(2))
    initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    if (LogThisAt(2))
      logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (SeqIdx i : env.forward_i(j)) {

      mat(i,j) = max (max (mat(i-1,j-1) + qos.m2m,
			   del(i-1,j-1) + qos.d2m),
		      ins(i-1,j-1) + qos.i2m);

      if (j == 1 || i == 1)
	mat(i,j) = max (mat(i,j),
			start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = max (log_sum_exp (ins(i,j-1) + qos.i2i,
				   del(i,j-1) + qos.d2i),
		      mat(i,j-1) + qos.m2i);

      del(i,j) = max (log_sum_exp (del(i-1,j) + qos.d2d,
				   ins(i-1,j) + qos.d2i),
		      mat(i-1,j) + qos.m2d);

      if (j == yLen || i == xLen)
	end = max (end,
		   mat(i,j));
    }
  }

  result = end + xInsertScore + yInsertScore;

  if (LogThisAt(2))
    cerr << "Viterbi score: " << result << endl;
}

Alignment QuaffOverlapViterbiMatrix::alignment() const {
  Assert (resultIsFinite(), "Can't do Viterbi traceback if final score is -infinity");
  SeqIdx xEnd = xLen, yEnd = yLen;
  double bestEndSc = mat(xLen,yLen);
  double sc;
  for (SeqIdx iEnd = xLen; iEnd > 0; --iEnd) {
    sc = mat(iEnd,yLen);
    if (sc > bestEndSc) {
      bestEndSc = sc;
      xEnd = iEnd;
      yEnd = yLen;
    }
  }
  for (SeqIdx jEnd = yLen; jEnd > 0; --jEnd) {
    sc = mat(xLen,jEnd);
    if (sc > bestEndSc) {
      bestEndSc = sc;
      xEnd = xLen;
      yEnd = jEnd;
    }
  }

  SeqIdx i = xEnd, j = yEnd;
  list<char> xRow, yRow;
  State state = Match;
  while (state != Start) {
    if (LogThisAt(7))
      cerr << "Traceback: i=" << i << " j=" << j << " state=" << stateToString(state) << " score=" << cellScore(i,j,state) << endl;
    double srcSc = -numeric_limits<double>::infinity();
    double emitSc = 0;
    switch (state) {
    case Match:
      emitSc = matchEmitScore(i,j);
      xRow.push_front (px->seq[--i]);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat(i,j) + qos.m2m + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + qos.i2m + emitSc, Insert);
      updateMax (srcSc, state, del(i,j) + qos.d2m + emitSc, Delete);
      if (j == 0 || i == 0)
	updateMax (srcSc, state, emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error");
      break;

    case Insert:
      xRow.push_front (Alignment::gapChar);
      yRow.push_front (py->seq[--j]);
      updateMax (srcSc, state, mat(i,j) + qos.m2i, Match);
      updateMax (srcSc, state, ins(i,j) + qos.i2i, Insert);
      updateMax (srcSc, state, del(i,j) + qos.d2i, Delete);
      break;

    case Delete:
      xRow.push_front (px->seq[--i]);
      yRow.push_front (Alignment::gapChar);
      updateMax (srcSc, state, mat(i,j) + qos.m2d, Match);
      updateMax (srcSc, state, ins(i,j) + qos.i2d, Insert);
      updateMax (srcSc, state, del(i,j) + qos.d2d, Delete);
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  const SeqIdx xStart = i + 1;
  const SeqIdx yStart = j + 1;
  Alignment align(2);
  align.gappedSeq[0].name = "substr(" + px->name + "," + to_string(xStart) + ".." + to_string(xEnd) + ")";
  align.gappedSeq[1].name = "substr(" + py->name + "," + to_string(yStart) + ".." + to_string(yEnd) + ")";
  align.gappedSeq[0].seq = string (xRow.begin(), xRow.end());
  align.gappedSeq[1].seq = string (yRow.begin(), yRow.end());
  align.score = result;
  return align;
}

Alignment QuaffOverlapViterbiMatrix::scoreAdjustedAlignment (const QuaffNullParams& nullModel) const {
  Alignment a = alignment();
  a.score -= nullModel.logLikelihood (*px);
  a.score -= nullModel.logLikelihood (qos.yComplemented ? py->revcomp() : *py);
  return a;
}

QuaffOverlapAligner::QuaffOverlapAligner()
  : QuaffAlignmentPrinter()
{ }

bool QuaffOverlapAligner::parseAlignmentArgs (int& argc, char**& argv) {
  return parseAlignmentPrinterArgs (argc, argv);
}

void QuaffOverlapAligner::align (ostream& out, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config) {
  for (size_t nx = 0; nx < nOriginals; ++nx) {
    const FastSeq& xfs = seqs[nx];
    for (size_t ny = nx + 1; ny < seqs.size(); ++ny) {
      const FastSeq& yfs = seqs[ny];
      DiagonalEnvelope env = config.makeEnvelope (xfs, yfs);
      const QuaffOverlapViterbiMatrix viterbi (env, params, ny >= nOriginals);
      if (viterbi.resultIsFinite()) {
	const Alignment align = viterbi.scoreAdjustedAlignment(nullModel).addScoreComment();
	writeAlignment (out, align);
      }
    }
  }
}


