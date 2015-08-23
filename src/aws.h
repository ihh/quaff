#ifndef AWS_INCLUDED
#define AWS_INCLUDED

#include <string>
#include "util.h"

using namespace std;

struct AWS {
  static void syncFromBucket (const string& bucket, const string& filename);
  static void copyToBucket (const string& filename, const string& bucket);
};

#endif /* AWS_INCLUDED */
