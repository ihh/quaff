#ifndef DIAG_ENV_INCLUDED
#define DIAG_ENV_INCLUDED

#include <map>
#include "fastseq.h"

using namespace std;

struct DiagonalEnvelope {
  const FastSeq *px, *py;
  const SeqIdx xLen, yLen;
  vector<int> diagonals;   // Sorted. (i,j) is on diagonal d if i-j=d
  DiagonalEnvelope (const FastSeq& x, const FastSeq& y)
    : px(&x), py(&y), xLen(px->length()), yLen(py->length()) { }
  void initFull();
  void initSparse (unsigned int kmerLen, unsigned int bandSize);
  inline int minDiagonal() const { return 1 - (int) yLen; }
  inline int maxDiagonal() const { return xLen - 1; }
  inline bool intersects (SeqIdx j, int diag) const {
    const int i = diag + j;
    return i > 0 && i <= xLen;
  }
  static inline SeqIdx get_i (SeqIdx j, int diag) { return (SeqIdx) (diag + j); }
  static inline int get_diag (SeqIdx i, SeqIdx j) { return (int) i - (int) j; }
};


#endif /* DIAG_ENV_INCLUDED */
