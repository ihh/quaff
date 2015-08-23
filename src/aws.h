#ifndef AWS_INCLUDED
#define AWS_INCLUDED

#include <string>
#include "util.h"
#include "vguard.h"

using namespace std;

struct AWS {
  // S3
  static void syncFromBucket (const string& bucket, const string& filename);
  static void copyToBucket (const string& filename, const string& bucket);
  // EC2
  static vguard<string> launchInstancesWithScript (unsigned int nInstances, const string& keyPairName, const string& instanceType, const string& ami, const string& userDataScript);
  static string getInstancePublicDNS (const string& instanceID);
  static void terminateInstance (const string& instanceID);
};

#endif /* AWS_INCLUDED */
