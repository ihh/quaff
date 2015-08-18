#include <list>
#include <cmath>
#include "qoverlap.h"
#include "logsumexp.h"
#include "logger.h"

QuaffOverlapScores::QuaffOverlapScores (const QuaffParams &qp, bool yComplemented)
  : matchContext (qp.matchContext.kmerLen),
    indelContext (qp.indelContext.kmerLen),
    pqp(&qp),
    yComplemented(yComplemented),
    xInsert (dnaAlphabetSize),
    yInsert (dnaAlphabetSize),
    matchMinusInsert (matchContext.numKmers, vguard<SymQualPairScores> (matchContext.numKmers)),
    gapOpenProb (indelContext.numKmers),
    m2m (indelContext.numKmers, vguard<double> (indelContext.numKmers)),
    m2i (indelContext.numKmers, vguard<double> (indelContext.numKmers)),
    m2d (indelContext.numKmers, vguard<double> (indelContext.numKmers))
{
  // coarse approximation to various paths through the intersection of two Quaff transducers...
  vguard<double> pGapIsInsertK, gapAdjacentProbK;
  for (Kmer j = 0; j < indelContext.numKmers; ++j) {
    const double readInsertProb = qp.beginInsert[j];
    const double readDeleteProb = (1 - qp.beginInsert[j]) * qp.beginDelete[j];
    gapOpenProb[j] = readInsertProb + readDeleteProb;
    const double pGapIsInsert = readInsertProb / gapOpenProb[j];
    const double gapAdjacentProb = pGapIsInsert * readInsertProb + (1 - pGapIsInsert) * gapOpenProb[j] / (1 - qp.extendDelete * (1 - gapOpenProb[j]));
    pGapIsInsertK.push_back (pGapIsInsert);
    gapAdjacentProbK.push_back (gapAdjacentProb);
  }

  for (Kmer i = 0; i < indelContext.numKmers; ++i)
    for (Kmer j = 0; j < indelContext.numKmers; ++j) {
      m2m[i][j] = log(1-gapOpenProb[i]) + log(1-gapOpenProb[j]);
      m2i[i][j] = log(gapOpenProb[i]);
      m2d[i][j] = log(1-gapOpenProb[i]) + log(gapOpenProb[j]);
    }

  const double pGapIsInsert = accumulate (pGapIsInsertK.begin(), pGapIsInsertK.end(), 0.) / pGapIsInsertK.size();
  const double meanGapLength = pGapIsInsert / qp.extendInsert + (1 - pGapIsInsert) / qp.extendDelete;
  gapExtendProb = 1 / meanGapLength;
  gapAdjacentProb = accumulate (gapAdjacentProbK.begin(), gapAdjacentProbK.end(), 0.) / gapAdjacentProbK.size();
  
  i2i = d2d = log(gapExtendProb);
  i2d = d2i = log(1-gapExtendProb) + log(gapAdjacentProb);
  i2m = d2m = log(1-gapExtendProb) + log(1-gapAdjacentProb);

  const QuaffScores qsc (qp);
  xInsert = yInsert = qsc.insert;

  for (Kmer iPrefix = 0; iPrefix < matchContext.numKmers; iPrefix += dnaAlphabetSize)
    for (AlphTok iSuffix = 0; iSuffix < dnaAlphabetSize; ++iSuffix) {
      const Kmer i = iPrefix + iSuffix;
      for (Kmer jPrefix = 0; jPrefix < matchContext.numKmers; jPrefix += dnaAlphabetSize)
	for (AlphTok jSuffix = 0; jSuffix < dnaAlphabetSize; ++jSuffix) {
	  const Kmer j = jPrefix + jSuffix;
	  for (QualScore ik = 0; ik < FastSeq::qualScoreRange; ++ik) {
	    for (QualScore jk = 0; jk < FastSeq::qualScoreRange; ++jk) {
	      double mij = -numeric_limits<double>::infinity();
	      for (AlphTok r = 0; r < dnaAlphabetSize; ++r) {
		const AlphTok yStrand_r = yComplemented ? dnaComplement(r) : r;
		mij = log_sum_exp (mij, log(qp.refBase[r]) + qsc.match[r][i].logSymQualProb[ik] + qsc.match[yStrand_r][j].logSymQualProb[jk]);
	      }
	      SymQualPairScores& sqps = matchMinusInsert[i][j];
	      sqps.logSymQualPairProb[ik][jk] = mij - xInsert[iSuffix].logSymQualProb[ik] - yInsert[jSuffix].logSymQualProb[jk];
	      sqps.logSymPairXQualProb[ik] = log_sum_exp (sqps.logSymPairXQualProb[ik], mij - xInsert[iSuffix].logSymQualProb[ik] - yInsert[jSuffix].logSymProb);
	      sqps.logSymPairYQualProb[jk] = log_sum_exp (sqps.logSymPairYQualProb[jk], mij - xInsert[iSuffix].logSymProb - yInsert[jSuffix].logSymQualProb[jk]);
	      sqps.logSymPairProb = log_sum_exp (sqps.logSymPairProb, mij - xInsert[iSuffix].logSymProb - yInsert[jSuffix].logSymProb);
	    }
	  }
	}
    }
}

QuaffOverlapViterbiMatrix::QuaffOverlapViterbiMatrix (const DiagonalEnvelope& env, const QuaffParams& qp, bool yComplemented)
  : QuaffDPMatrixContainer (env),
    qos (qp, yComplemented),
    xTok (px->tokens(dnaAlphabet)),
    xMatchKmer (px->kmers(dnaAlphabet,qp.matchContext.kmerLen)),
    xIndelKmer (px->kmers(dnaAlphabet,qp.indelContext.kmerLen)),
    xQual (px->qualScores()),
    yQual (py->qualScores()),
    xInsertScore (0),
    yInsertScore (0)
{
  const FastSeq& x (*px);
  const FastSeq& y (*py);

  if (yComplemented) {
    const FastSeq yRevcomp = y.revcomp();
    const vguard<AlphTok> yRevcompTok = yRevcomp.tokens (dnaAlphabet);
    const vguard<Kmer> yRevcompMatchKmer = yRevcomp.kmers (dnaAlphabet, qp.matchContext.kmerLen);
    const vguard<Kmer> yRevcompIndelKmer = yRevcomp.kmers (dnaAlphabet, qp.indelContext.kmerLen);
    yTok = vguard<AlphTok> (yRevcompTok.rbegin(), yRevcompTok.rend());
    yMatchKmer = vguard<Kmer> (yRevcompMatchKmer.rbegin(), yRevcompMatchKmer.rend());
    yIndelKmer = vguard<Kmer> (yRevcompIndelKmer.rbegin(), yRevcompIndelKmer.rend());
  } else {
    yTok = y.tokens(dnaAlphabet);
    yMatchKmer = y.kmers(dnaAlphabet,qp.matchContext.kmerLen);
    yIndelKmer = y.kmers(dnaAlphabet,qp.indelContext.kmerLen);
  }

  for (SeqIdx j = 0; j < xLen; ++j)
    Assert (xMatchKmer[j] < qos.matchContext.numKmers, "oops. xMatchKmer %u/%u is %llu", j, xLen, xMatchKmer[j]);

  for (SeqIdx j = 0; j < xLen; ++j)
    Assert (xIndelKmer[j] < qos.indelContext.numKmers, "oops. xIndelKmer %u/%u is %llu", j, xLen, xIndelKmer[j]);

  for (SeqIdx j = 0; j < yLen; ++j)
    Assert (yMatchKmer[j] < qos.matchContext.numKmers, "oops. yMatchKmer %u/%u is %llu", j, yLen, yMatchKmer[j]);

  for (SeqIdx j = 0; j < yLen; ++j)
    Assert (yIndelKmer[j] < qos.indelContext.numKmers, "oops. yIndelKmer %u/%u is %llu", j, yLen, yIndelKmer[j]);

  for (SeqIdx i = 0; i < xLen; ++i) {
    const SymQualScores& sqs = qos.xInsert[xTok[i]];
    xInsertScore += xQual.size() ? sqs.logSymQualProb[xQual[i]] : sqs.logSymProb;
  }

  for (SeqIdx i = 0; i < yLen; ++i) {
    const SymQualScores& sqs = qos.yInsert[yTok[i]];
    yInsertScore += yQual.size() ? sqs.logSymQualProb[yQual[i]] : sqs.logSymProb;
  }

  if (LogThisAt(2))
    initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    if (LogThisAt(2))
      logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (SeqIdx i : env.forward_i(j)) {

      mat(i,j) = max (max (mat(i-1,j-1) + m2mScore(i-1,j-1),
			   del(i-1,j-1) + d2mScore()),
		      ins(i-1,j-1) + i2mScore());

      if (j == 1 || i == 1)
	mat(i,j) = max (mat(i,j),
			start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = max (log_sum_exp (ins(i,j-1) + i2iScore(),
				   del(i,j-1) + d2iScore()),
		      mat(i,j-1) + m2iScore(i,j-1));

      del(i,j) = max (log_sum_exp (del(i-1,j) + d2dScore(),
				   ins(i-1,j) + d2iScore()),
		      mat(i-1,j) + m2dScore(i-1,j));

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
  list<char> xRow, yRow, xQual, yQual;
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
      if (px->hasQual())
	xQual.push_front (px->qual[i]);
      if (py->hasQual())
	yQual.push_front (py->qual[j]);
      updateMax (srcSc, state, mat(i,j) + m2mScore(i,j) + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + i2mScore() + emitSc, Insert);
      updateMax (srcSc, state, del(i,j) + d2mScore() + emitSc, Delete);
      if (j == 0 || i == 0)
	updateMax (srcSc, state, emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error");
      break;

    case Insert:
      xRow.push_front (Alignment::gapChar);
      yRow.push_front (py->seq[--j]);
      if (px->hasQual())
	xQual.push_front (FastSeq::maxQualityChar);
      if (py->hasQual())
	yQual.push_front (py->qual[j]);
      updateMax (srcSc, state, mat(i,j) + m2iScore(i,j), Match);
      updateMax (srcSc, state, ins(i,j) + i2iScore(), Insert);
      updateMax (srcSc, state, del(i,j) + d2iScore(), Delete);
      break;

    case Delete:
      xRow.push_front (px->seq[--i]);
      yRow.push_front (FastSeq::maxQualityChar);
      if (px->hasQual())
	xQual.push_front (px->qual[i]);
      if (py->hasQual())
	yQual.push_front (Alignment::gapChar);
      updateMax (srcSc, state, mat(i,j) + m2dScore(i,j), Match);
      updateMax (srcSc, state, ins(i,j) + i2dScore(), Insert);
      updateMax (srcSc, state, del(i,j) + d2dScore(), Delete);
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  const SeqIdx xStart = i + 1;
  const SeqIdx yStart = j + 1;
  Alignment align(2);
  align.gappedSeq[0].name = "read_x";
  align.gappedSeq[0].comment = "substr(" + px->name + "," + to_string(xStart) + ".." + to_string(xEnd) + ")";
  align.gappedSeq[1].name = "read_y";
  align.gappedSeq[1].comment = "substr(" + py->name + "," + to_string(yStart) + ".." + to_string(yEnd) + ")";
  align.gappedSeq[0].seq = string (xRow.begin(), xRow.end());
  align.gappedSeq[1].seq = string (yRow.begin(), yRow.end());
  align.gappedSeq[0].qual = string (xQual.begin(), xQual.end());
  align.gappedSeq[1].qual = string (yQual.begin(), yQual.end());
  align.score = result;
  return align;
}

Alignment QuaffOverlapViterbiMatrix::scoreAdjustedAlignment (const QuaffNullParams& nullModel) const {
  Alignment a = alignment();
  const double xNullLogLike = nullModel.logLikelihood (*px);
  const double yNullLogLike = nullModel.logLikelihood (qos.yComplemented ? py->revcomp() : *py);
  if (LogThisAt(2))
    cerr << "Null model score: " << (xNullLogLike + yNullLogLike) << endl;
  if (LogThisAt(3))
    cerr << "Null model score for " << px->name << ": " << xNullLogLike << endl
	 << "Null model score for " << py->name << ": " << yNullLogLike << endl;
  a.score -= xNullLogLike;
  a.score -= yNullLogLike;
  return a;
}

QuaffOverlapAligner::QuaffOverlapAligner()
  : QuaffAlignmentPrinter()
{ }

bool QuaffOverlapAligner::parseAlignmentArgs (deque<string>& argvec) {
  return parseAlignmentPrinterArgs (argvec);
}

void QuaffOverlapAligner::align (ostream& out, const vguard<FastSeq>& seqs, size_t nOriginals, const QuaffParams& params, const QuaffNullParams& nullModel, const QuaffDPConfig& config) {
  for (size_t nx = 0; nx < nOriginals; ++nx) {
    const FastSeq& xfs = seqs[nx];
    for (size_t ny = nx + 1; ny < seqs.size(); ++ny) {
      const FastSeq& yfs = seqs[ny];
      if (LogThisAt(1))
	cerr << "Aligning " << xfs.name << " (length " << xfs.length() << ") to " << yfs.name << " (length " << yfs.length() << ")" << endl;
      DiagonalEnvelope env = config.makeEnvelope (xfs, yfs);
      const QuaffOverlapViterbiMatrix viterbi (env, params, ny >= nOriginals);
      if (viterbi.resultIsFinite()) {
	const Alignment align = viterbi.scoreAdjustedAlignment(nullModel);
	writeAlignment (out, align);
      }
    }
  }
}


