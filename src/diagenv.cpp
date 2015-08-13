#include <algorithm>
#include <cmath>
#include "diagenv.h"

void DiagonalEnvelope::initFull() {
  int xLen = (int) px->length(), yLen = (int) py->length();
  diagonals.clear();
  diagonals.reserve (xLen + yLen - 1);
  for (int d = 1 - yLen; d <= xLen - 1; ++d)
    diagonals.push_back (d);
}

void DiagonalEnvelope::initSparse (int kmerLen) {
  // require at least 2^kmerLen bases on each axis, or fall back to full envelope
  const double minLen = pow (2, kmerLen);
  if (px->length() < minLen || py->length() < minLen)
    initFull();
  else {
    Abort ("unimplemented");
  }
}
