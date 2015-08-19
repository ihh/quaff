#ifndef DIAG_ENV_INCLUDED
#define DIAG_ENV_INCLUDED

#include <map>
#include "fastseq.h"

using namespace std;

// Somewhat arbitrary numbers:
// Default k-mer length
#define DEFAULT_KMER_LENGTH 6

// The number of kmer matches above which a diagonal is selected for a band
#define DEFAULT_KMER_THRESHOLD 14

// Default band size
#define DEFAULT_BAND_SIZE 64

struct DiagonalEnvelope {
  const FastSeq *px, *py;
  const SeqIdx xLen, yLen;
  vector<int> diagonals;   // Sorted. (i,j) is on diagonal d if i-j=d
  DiagonalEnvelope (const FastSeq& x, const FastSeq& y)
    : px(&x), py(&y), xLen(px->length()), yLen(py->length()) { }
  void initFull();
  void initSparse (unsigned int kmerLen = DEFAULT_KMER_LENGTH,
		   unsigned int bandSize = DEFAULT_BAND_SIZE,
		   int kmerThreshold = DEFAULT_KMER_THRESHOLD,  // negative => use memory guides
		   size_t cellSize = sizeof(double),
		   size_t maxSize = 0);
  inline int minDiagonal() const { return 1 - (int) yLen; }
  inline int maxDiagonal() const { return xLen - 1; }
  inline bool intersects (SeqIdx j, int diag) const {
    const int i = diag + j;
    return i > 0 && i <= xLen;
  }
  static inline SeqIdx get_i (SeqIdx j, int diag) { return (SeqIdx) (diag + j); }
  static inline int get_diag (SeqIdx i, SeqIdx j) { return (int) i - (int) j; }
  vector<SeqIdx> forward_i (SeqIdx j) const;
  vector<SeqIdx> reverse_i (SeqIdx j) const;
};


#endif /* DIAG_ENV_INCLUDED */
