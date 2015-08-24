#ifndef AWS_INCLUDED
#define AWS_INCLUDED

#include <string>
#include <deque>
#include <set>
#include <map>
#include "util.h"
#include "vguard.h"

using namespace std;

// some defaults
#define AWS_DEFAULT_AMI "ami-1ecae776"
#define AWS_DEFAULT_REGION "us-east-1"
#define AWS_DEFAULT_INSTANCE_TYPE "m3.medium"

class AWS {
private:
  // list of instances to terminate on abort, exit, etc
  static set<string> runningInstanceIds;
  static map<int,void (*)(int)> signalHandler;
  static int objectCount;  // enforce singleton
  static bool cleanupCalled;  // guard against multiple calls to cleanup()
public:
  // config
  string keyPair, securityGroup;
  const string bashBang;
  // constructor, destructor, signal handlers, cleanup
  AWS();
  ~AWS();
private:
  static void cleanup();
  static void handleSignal (int signal);
  void registerCleanup();
public:
  // config
  bool parseAWSConfigArgs (deque<string>& argvec);
  // S3
  void syncFromBucket (const string& bucket, const string& filename);
  void copyToBucket (const string& filename, const string& bucket);
  // EC2
  vguard<string> launchInstancesWithScript (unsigned int nInstances, const string& instanceType, const string& ami, const string& userDataScript);
  vguard<string> getInstanceAddresses (const vguard<string>& instanceIDs) const;
  void terminateInstances (const vguard<string>& instanceIDs);
  string terminateInstancesSilently (const vguard<string>& instanceIDs);
  // AWS environment helper for constructing user-data scripts
  string bashHeader() const;
};

// singleton
extern AWS aws;

#endif /* AWS_INCLUDED */
