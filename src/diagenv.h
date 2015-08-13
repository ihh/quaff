#ifndef DIAG_ENV_INCLUDED
#define DIAG_ENV_INCLUDED

#include <map>
#include "fastseq.h"

using namespace std;

struct DiagonalEnvelope {
  const FastSeq *px, *py;
  vector<int> diagonals;   // (i,j) is on diagonal d if i-j=d
  DiagonalEnvelope (const FastSeq& x, const FastSeq& y)
    : px(&x), py(&y) { }
  void initFull();
  void initSparse (int kmerLen);
};


#endif /* DIAG_ENV_INCLUDED */
