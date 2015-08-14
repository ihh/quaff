#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "diagenv.h"
#include "logger.h"

void DiagonalEnvelope::initFull() {
  diagonals.clear();
  diagonals.reserve (xLen + yLen - 1);
  for (int d = minDiagonal(); d <= maxDiagonal(); ++d)
    diagonals.push_back (d);
}

void DiagonalEnvelope::initSparse (unsigned int kmerLen, unsigned int bandSize) {
  // require at least 2^kmerLen bases on each axis, or fall back to full envelope
  const double minLen = pow (2, kmerLen);
  if (px->length() < minLen || py->length() < minLen) {
    initFull();
    return;
  }

  const vector<unsigned int> xTok = px->tokens (dnaAlphabet);
  const vector<unsigned int> yTok = py->tokens (dnaAlphabet);

  map<unsigned long,set<SeqIdx> > yKmerIndex;
  for (SeqIdx j = 0; j <= yLen - kmerLen; ++j)
    yKmerIndex[makeKmer (kmerLen, yTok.begin() + j, dnaAlphabetSize)].insert (j);

  if (LogThisAt(6)) {
    cerr << "kmer distribution for " << py->name << ':' << endl;
    for (const auto& yKmerIndexElt : yKmerIndex) {
      cerr << kmerToString (yKmerIndexElt.first, kmerLen, dnaAlphabet) << ' ' << yKmerIndexElt.second.size() << endl;
    }
  }

  map<int,unsigned int> diagKmerCount;
  for (SeqIdx i = 0; i <= xLen - kmerLen; ++i) {
    const auto yKmerIndexIter = yKmerIndex.find (makeKmer (kmerLen, xTok.begin() + i, dnaAlphabetSize));
    if (yKmerIndexIter != yKmerIndex.end())
      for (auto j : yKmerIndexIter->second)
	++diagKmerCount[get_diag(i,j)];
  }

  if (LogThisAt(5)) {
    cerr << "Distribution of " << kmerLen << "-mer matches per diagonal for " << px->name << " vs " << py->name << ':' << endl;
    map<unsigned int,unsigned int> countDistrib;
    for (const auto& diagKmerCountElt : diagKmerCount)
      ++countDistrib[diagKmerCountElt.second];
    for (const auto& countDistribElt : countDistrib)
      cerr << countDistribElt.second << " diagonals have " << countDistribElt.first << " matches" << endl;
  }

  double n = 0, sum = 0, sqsum = 0;
  for (int d = minDiagonal(); d <= maxDiagonal(); ++d) {
    const auto diagKmerCountIter = diagKmerCount.find(d);
    if (diagKmerCountIter != diagKmerCount.end()) {
      const unsigned int count = diagKmerCountIter->second;
      sum += count;
      sqsum += count*count;
    }
    ++n;
  }

  const double mean = sum / n,
    variance = sqsum / n - mean*mean,
    stdev = sqrt(variance);

  const unsigned int threshold = (unsigned int) (mean + 10*stdev);  // 10 standard deviations is arbitrary city!

  if (LogThisAt(3))
    cerr << "kmers per diagonal: mean " << mean << ", sd " << stdev << ". Threshold for inclusion: " << threshold << endl;

  set<int> diags;
  int nPastThreshold = 0;
  for (auto& diagKmerCountElt : diagKmerCount)
    if (diagKmerCountElt.second >= threshold) {
      ++nPastThreshold;
      for (int d = diagKmerCountElt.first - bandSize / 2;
	   d <= diagKmerCountElt.first + bandSize / 2;
	   ++d)
	diags.insert (d);
    }

  if (LogThisAt(3))
    cerr << nPastThreshold << " diagonals above threshold; " << diags.size() << " in envelope" << endl;

  diagonals = vector<int> (diags.begin(), diags.end());
}
