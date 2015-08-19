#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "diagenv.h"
#include "logger.h"

// require at least this many bases per sequence for sparse envelope
#define MIN_SEQLEN_FOR_SPARSE_ENVELOPE 1024

void DiagonalEnvelope::initFull() {
  diagonals.clear();
  diagonals.reserve (xLen + yLen - 1);
  for (int d = minDiagonal(); d <= maxDiagonal(); ++d)
    diagonals.push_back (d);
}

void DiagonalEnvelope::initSparse (unsigned int kmerLen, unsigned int bandSize, int kmerThreshold, size_t cellSize, size_t maxSize) {
  if (px->length() < MIN_SEQLEN_FOR_SPARSE_ENVELOPE || py->length() < MIN_SEQLEN_FOR_SPARSE_ENVELOPE) {
    initFull();
    return;
  }

  const vector<AlphTok> xTok = px->tokens (dnaAlphabet);
  const vector<AlphTok> yTok = py->tokens (dnaAlphabet);

  map<Kmer,set<SeqIdx> > yKmerIndex;
  for (SeqIdx j = 0; j <= yLen - kmerLen; ++j)
    yKmerIndex[makeKmer (kmerLen, yTok.begin() + j, dnaAlphabetSize)].insert (j);

  if (LogThisAt(8)) {
    logger << "Frequencies of " << kmerLen << "-mers in " << py->name << ':' << endl;
    for (const auto& yKmerIndexElt : yKmerIndex) {
      logger << kmerToString (yKmerIndexElt.first, kmerLen, dnaAlphabet) << ' ' << yKmerIndexElt.second.size() << endl;
    }
  }

  map<int,unsigned int> diagKmerCount;
  for (SeqIdx i = 0; i <= xLen - kmerLen; ++i) {
    const auto yKmerIndexIter = yKmerIndex.find (makeKmer (kmerLen, xTok.begin() + i, dnaAlphabetSize));
    if (yKmerIndexIter != yKmerIndex.end())
      for (auto j : yKmerIndexIter->second)
	++diagKmerCount[get_diag(i,j)];
  }

  map<unsigned int,set<unsigned int> > countDistrib;
  for (const auto& diagKmerCountElt : diagKmerCount)
    countDistrib[diagKmerCountElt.second].insert (diagKmerCountElt.first);

  if (LogThisAt(7)) {
    logger << "Distribution of " << kmerLen << "-mer matches per diagonal for " << px->name << " vs " << py->name << ':' << endl;
    for (const auto& countDistribElt : countDistrib)
      logger << countDistribElt.second.size() << " diagonals have " << countDistribElt.first << " matches" << endl;
  }

  set<int> diags;
  diags.insert(0);  // always add the zeroth diagonal to ensure at least one path exists

  const unsigned int halfBandSize = bandSize / 2;
  const size_t diagSize = min(xLen,yLen) * cellSize;
  unsigned int nPastThreshold = 0;

  unsigned int threshold;
  if (kmerThreshold >= 0)
    threshold = kmerThreshold;
  else if (LogThisAt(5))
    logger << "Automatically setting threshold based on memory limit of " << maxSize << " bytes (each diagonal takes " << diagSize << " bytes)" << endl;

  for (auto countDistribIter = countDistrib.crbegin();
       countDistribIter != countDistrib.crend();
       ++countDistribIter) {

    if (kmerThreshold >= 0 && countDistribIter->first < kmerThreshold)
      break;

    set<int> moreDiags = diags;
    for (auto seedDiag : countDistribIter->second) {
      ++nPastThreshold;
      for (int d = seedDiag - (int) halfBandSize;
	   d <= seedDiag + (int) halfBandSize;
	   ++d)
	moreDiags.insert (d);
    }

    if (kmerThreshold < 0) {
      if (moreDiags.size() * diagSize >= maxSize)
	break;
      threshold = countDistribIter->first;
    }
    swap (diags, moreDiags);
  }

  if (LogThisAt(5))
    logger << "Threshold # of " << kmerLen << "-mer matches for seeding a diagonal is " << threshold << endl;
  
  if (LogThisAt(5))
    logger << nPastThreshold << " diagonals above threshold; " << diags.size() << " in envelope (band size " << bandSize << "); estimated memory <" << (((diags.size() * diagSize) >> 20) + 1) << "MB" << endl;

  diagonals = vector<int> (diags.begin(), diags.end());
}

vector<SeqIdx> DiagonalEnvelope::forward_i (SeqIdx j) const {
  vector<SeqIdx> i_vec;
  i_vec.reserve (diagonals.size());
  for (auto d : diagonals)
    if (intersects (j, d))
      i_vec.push_back (get_i (j, d));
  return i_vec;
}

vector<SeqIdx> DiagonalEnvelope::reverse_i (SeqIdx j) const {
  const vector<SeqIdx> f = forward_i (j);
  return vector<SeqIdx> (f.rbegin(), f.rend());
}
