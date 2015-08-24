#include <cstdlib>
#include "aws.h"
#include "gason.h"
#include "util.h"
#include "logger.h"

JsonValue* jsonFind (JsonValue& parent, const char* key);
JsonValue& jsonFindOrDie (JsonValue& parent, const char* key);
  
JsonValue* jsonFind (JsonValue& parent, const char* key) {
  Assert (parent.getTag() == JSON_OBJECT, "JSON value is not an object");
  for (auto i : parent)
    if (strcmp (i->key, key) == 0)
      return &i->value;
  return NULL;
}

JsonValue& jsonFindOrDie (JsonValue& parent, const char* key) {
  JsonValue* val = jsonFind (parent, key);
  Assert (val != NULL, "Couldn't find JSON tag %s", key);
  return *val;
}

int AWS::objectCount = 0;
set<string> AWS::runningInstanceIds;
map<int,void (*)(int)> AWS::signalHandler;
bool AWS::cleanupCalled = false;

AWS aws;

void AWS::syncFromBucket (const string& bucket, const string& filename) {
  const string cmd = string("aws sync s3://") + bucket + "/ . --include " + filename;
  const int status = system (cmd.c_str());
  if (status != 0)
    Warn ("Return code %d attempting to sync file %s from S3 bucket %s", status, filename.c_str(), bucket.c_str());
}

void AWS::copyToBucket (const string& filename, const string& bucket) {
  const string cmd = string("aws cp ") + filename + " s3://" + bucket + '/';
  const int status = system (cmd.c_str());
  if (status != 0)
    Warn ("Return code %d attempting to copy file %s to S3 bucket %s", status, filename.c_str(), bucket.c_str());
}

vguard<string> AWS::launchInstancesWithScript (unsigned int nInstances, const string& instanceType, const string& ami, const string& userDataScript) {
  vguard<string> ids;
  TempFile tempfile (userDataScript, "userdata");
  const string runCmd = string("aws ec2 run-instances --enable-api-termination --key-name " + keyPair + " --security-groups " + securityGroup + " --instance-type " + instanceType + " --count " + to_string(nInstances) + ':' + to_string(nInstances) + " --image-id " + ami + " --user-data " + tempfile.fullPath);

  if (LogThisAt(9))
    logger << "Executing: " << runCmd << endl;

  const string runOut = pipeToString (runCmd.c_str());
  char* runOutCStr = new char [runOut.length()+1];
  strcpy (runOutCStr, runOut.c_str());

  if (LogThisAt(9))
    logger << "Output: " << runOut << endl;

  char *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  const int parseStatus = jsonParse(runOutCStr, &endPtr, &value, allocator);
  Assert (parseStatus == JSON_OK, "JSON parsing error: %s at byte %zd of output of %s\n%s\n", jsonStrError(parseStatus), endPtr - runOutCStr, runCmd.c_str(), endPtr);

  JsonValue& instanceList = jsonFindOrDie (value, "Instances");
  Assert (instanceList.getTag() == JSON_ARRAY, "JSON parsing error: Instances is not a list");
  for (auto inst : instanceList) {
    JsonValue& instanceId = jsonFindOrDie (inst->value, "InstanceId");
    Assert (instanceId.getTag() == JSON_STRING, "JSON parsing error: Instances->InstanceId is not a string");
    ids.push_back (instanceId.toString());
  }
  delete[] runOutCStr;

  for (const auto& id : ids)
    runningInstanceIds.insert (id);

  if (LogThisAt(3))
    logger << "Waiting for EC2 instance-status-ok: " << join(ids) << endl;

  const string waitCmd = string("aws ec2 wait instance-status-ok --instance-id ") + join(ids);
  const int waitStatus = system(waitCmd.c_str());
  Assert (waitStatus == 0, "Failed: %s", waitCmd.c_str());

  if (LogThisAt(3))
    logger << "EC2 instances running: " << join(ids) << endl;

  return ids;
}

vguard<string> AWS::getInstanceAddresses (const vguard<string>& ids) const {
  const string describeCmd = string("aws ec2 describe-instances --instance-id ") + join(ids);
  const string describeOut = pipeToString (describeCmd.c_str());
  char* describeOutCStr = new char [describeOut.length()+1];
  strcpy (describeOutCStr, describeOut.c_str());

  char *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  const int parseStatus = jsonParse(describeOutCStr, &endPtr, &value, allocator);
  Assert (parseStatus == JSON_OK, "JSON parsing error: %s at byte %zd of output of %s\n%s\n", jsonStrError(parseStatus), endPtr - describeOutCStr, describeCmd.c_str(), endPtr);

  JsonValue& reservationList = jsonFindOrDie (value, "Reservations");
  Assert (reservationList.getTag() == JSON_ARRAY, "JSON parsing error: Reservations is not a list");

  vguard<string> addr;
  for (auto res : reservationList) {
    JsonValue& instanceList = jsonFindOrDie (res->value, "Instances");
    Assert (instanceList.getTag() == JSON_ARRAY, "JSON parsing error: Reservations->Instances is not a list");
    for (auto inst : instanceList) {
      JsonValue& publicIp = jsonFindOrDie (inst->value, "PublicIpAddress");
      Assert (publicIp.getTag() == JSON_STRING, "JSON parsing error: Reservations->Instances->PublicIpAddress is not a string");
      addr.push_back (publicIp.toString());
    }
  }

  delete[] describeOutCStr;
  return addr;
}

void AWS::terminateInstances (const vguard<string>& ids) {
  if (LogThisAt(1))
    logger << "Terminating instances: " << join(ids) << endl;
  const string termOut = terminateInstancesSilently (ids);
  if (LogThisAt(2))
    logger << termOut << flush;
  for (const auto& id : ids)
    runningInstanceIds.erase (id);
}

string AWS::terminateInstancesSilently (const vguard<string>& ids) {
  const string termCmd = string("aws ec2 terminate-instances --instance-id ") + join(ids);
  return pipeToString (termCmd.c_str());
}

string AWS::bashHeader() const {
  const char* keyid = getenv("AWS_ACCESS_KEY_ID");
  Assert (keyid != NULL, "AWS_ACCESS_KEY_ID undefined");

  const char* secretkey = getenv("AWS_SECRET_ACCESS_KEY");
  Assert (secretkey != NULL, "AWS_SECRET_ACCESS_KEY undefined");

  const char* region = getenv("AWS_DEFAULT_REGION");
  Assert (region != NULL, "AWS_DEFAULT_REGION undefined");

  return bashBang + "export AWS_ACCESS_KEY_ID=" + keyid
    + "\nexport AWS_SECRET_ACCESS_KEY=" + secretkey
    + "\nexport AWS_DEFAULT_REGION=" + region + "\n";
}

void AWS::cleanup() {
  if (!cleanupCalled) {
    cleanupCalled = true;
    if (!runningInstanceIds.empty()) {
      cerr << "Terminating EC2 instances..." << endl;
      const vguard<string> ids (runningInstanceIds.begin(), runningInstanceIds.end());
      (void) aws.terminateInstancesSilently (ids);
    }
  }
}

void AWS::handleSignal (int sig) {
  cleanup();
  auto iter = signalHandler.find(sig);
  if (iter != signalHandler.end()) {
    auto handler = iter->second;
    signalHandler.erase (iter);  // prevent infinite loops
    (*handler) (sig);
  }
}

void AWS::registerCleanup() {
  for (int sig : { SIGABRT, SIGFPE, SIGILL, SIGINT, SIGSEGV, SIGTERM }) {
    void (*oldHandler)(int) = signal (sig, &AWS::handleSignal);
    if (oldHandler != NULL)
      signalHandler[sig] = oldHandler;
  }
}

AWS::AWS()
  : bashBang ("#!/bin/bash\n")
{
  Assert (objectCount++ == 0, "aws object must be a singleton!");
  registerCleanup();
}

AWS::~AWS() {
  cleanup();
}

bool AWS::parseAWSConfigArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-keyname" || arg == "--key-name") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      keyPair = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-securitygroup" || arg == "--security-group") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      securityGroup = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}
