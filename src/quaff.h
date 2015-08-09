#ifndef QUAFF_INCLUDED
#define QUAFF_INCLUDED

#include "fastseq.h"

struct QuaffParams {
  double beginInsert, extendInsert, beginDelete, extendDelete;
  QuaffModel();
  void write (ostream& out) const;
  void read (istream& in);
};

struct QuaffLogScores {
  double m2m, m2i, m2d, m2e;
  double d2d, d2m;
  double i2i, i2m;
  QuaffLogScores (const QuaffParams& params);
};

struct QuaffDPMatrix {
  
};

#endif /* QUAFF_INCLUDED */
