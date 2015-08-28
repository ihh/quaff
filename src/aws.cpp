#include <cstring>
#include <cstdlib>
#include <signal.h>
#include <libgen.h>
#include "aws.h"
#include "jsonutil.h"
#include "util.h"
#include "logger.h"
#include "base64.h"

int AWS::objectCount = 0;
set<string> AWS::runningInstanceIds;
map<int,void (*)(int)> AWS::signalHandler;
bool AWS::cleanupCalled = false;

AWS aws;

string AWS::basenameStr (const string& filename) {
  char *f = new char [filename.size() + 1];
  strcpy (f, filename.c_str());
  const string n = basename(f);
  delete[] f;
  return n;
}

string AWS::dirnameStr (const string& filename) {
  char *f = new char [filename.size() + 1];
  strcpy (f, filename.c_str());
  const string n = dirname(f);
  delete[] f;
  return n;
}

string AWS::runCommandAndTestStatus (const string& cmd) {
  if (LogThisAt(4))
    logger << "Executing: " << cmd << endl;

  int cmdStatus = -1;
  const string cmdOut = pipeToString (cmd.c_str(), &cmdStatus);

  if (cmdOut.size() && LogThisAt(4))
    logger << "Output:" << endl << cmdOut << endl;

  if (cmdStatus != 0)
    Warn ("Return code %d attempting %s\n%s", cmdStatus, cmd.c_str(), cmdOut.c_str());

  return cmdOut;
}

void AWS::syncFromBucket (const string& bucket, const string& filename) {
  const string cmd = string("aws s3 sync s3://") + bucket + ' ' + dirnameStr(filename) + " --exclude '*' --include " + basenameStr(filename);

  (void) runCommandAndTestStatus (cmd);
}

void AWS::syncToBucket (const string& filename, const string& bucket) {
  const string cmd = string("aws s3 sync ") + dirnameStr(filename) + " s3://" + bucket + " --exclude '*' --include " + basenameStr(filename);

  (void) runCommandAndTestStatus (cmd);
}

vguard<string> AWS::launchInstancesWithScript (unsigned int nInstances, const string& instanceType, const string& ami, const string& userDataScript) {
  vguard<string> ids;

  string runCmd = string("aws ec2 run-instances --enable-api-termination --key-name ") + keyPair + " --security-groups " + securityGroup + " --instance-type " + instanceType + " --count " + to_string(nInstances) + ':' + to_string(nInstances) + " --image-id " + ami;
  if (userDataScript.size()) {
    const string base64UserDataScript = base64_encode (userDataScript.c_str(), (unsigned int) userDataScript.size());
    runCmd += " --user-data " + base64UserDataScript;
  }

  const string runOut = runCommandAndTestStatus (runCmd);

  char* runOutCStr = new char [runOut.length()+1];
  strcpy (runOutCStr, runOut.c_str());

  char *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  const int parseStatus = jsonParse(runOutCStr, &endPtr, &value, allocator);
  Assert (parseStatus == JSON_OK, "JSON parsing error: %s at byte %zd of output of %s\n%s\n", jsonStrError(parseStatus), endPtr - runOutCStr, runCmd.c_str(), endPtr);

  JsonValue& instanceList = JsonUtil::findOrDie (value, "Instances");
  Assert (instanceList.getTag() == JSON_ARRAY, "JSON parsing error: Instances is not a list");
  for (auto inst : instanceList) {
    JsonValue& instanceId = JsonUtil::findOrDie (inst->value, "InstanceId");
    Assert (instanceId.getTag() == JSON_STRING, "JSON parsing error: Instances->InstanceId is not a string");
    ids.push_back (instanceId.toString());
  }

  delete[] runOutCStr;

  for (const auto& id : ids)
    runningInstanceIds.insert (id);

  if (LogThisAt(3))
    logger << "Waiting for EC2 instance-status-ok: " << join(ids) << endl;

  const string waitCmd = string("aws ec2 wait instance-status-ok --instance-id ") + join(ids);

  (void) runCommandAndTestStatus (waitCmd);

  if (LogThisAt(3))
    logger << "EC2 instances running: " << join(ids) << endl;

  return ids;
}

vguard<string> AWS::getInstanceAddresses (const vguard<string>& ids) const {
  const string describeCmd = string("aws ec2 describe-instances --instance-id ") + join(ids);
  const string describeOut = runCommandAndTestStatus (describeCmd);

  char* describeOutCStr = new char [describeOut.length()+1];
  strcpy (describeOutCStr, describeOut.c_str());

  char *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  const int parseStatus = jsonParse(describeOutCStr, &endPtr, &value, allocator);
  Assert (parseStatus == JSON_OK, "JSON parsing error: %s at byte %zd of output of %s\n%s\n", jsonStrError(parseStatus), endPtr - describeOutCStr, describeCmd.c_str(), endPtr);

  JsonValue& reservationList = JsonUtil::findOrDie (value, "Reservations");
  Assert (reservationList.getTag() == JSON_ARRAY, "JSON parsing error: Reservations is not a list");

  vguard<string> addr;
  for (auto res : reservationList) {
    JsonValue& instanceList = JsonUtil::findOrDie (res->value, "Instances");
    Assert (instanceList.getTag() == JSON_ARRAY, "JSON parsing error: Reservations->Instances is not a list");
    for (auto inst : instanceList) {
      JsonValue& publicIp = JsonUtil::findOrDie (inst->value, "PublicIpAddress");
      Assert (publicIp.getTag() == JSON_STRING, "JSON parsing error: Reservations->Instances->PublicIpAddress is not a string");
      addr.push_back (publicIp.toString());
    }
  }

  delete[] describeOutCStr;
  return addr;
}

void AWS::terminateInstances (const vguard<string>& ids) {
  if (ids.size()) {
    if (LogThisAt(1))
      logger << "Terminating instances: " << join(ids) << endl;
    const string termCmd = terminateCommand (ids);
    (void) runCommandAndTestStatus (termCmd);
    for (const auto& id : ids)
      runningInstanceIds.erase (id);
  }
}

void AWS::terminateInstancesSilently (const vguard<string>& ids) {
  const string termCmd = terminateCommand (ids);
  (void) pipeToString (termCmd.c_str());
}

string AWS::terminateCommand (const vguard<string>& ids) const {
  return string("aws ec2 terminate-instances --instance-id ") + join(ids);
}

string AWS::bashEnvPrefix() const {
  const char* keyid = getenv("AWS_ACCESS_KEY_ID");
  Assert (keyid != NULL, "AWS_ACCESS_KEY_ID undefined");

  const char* secretkey = getenv("AWS_SECRET_ACCESS_KEY");
  Assert (secretkey != NULL, "AWS_SECRET_ACCESS_KEY undefined");

  const char* region = getenv("AWS_DEFAULT_REGION");

  return string("env AWS_ACCESS_KEY_ID=") + keyid
    + " AWS_SECRET_ACCESS_KEY=" + secretkey
    + (region == NULL ? string() : (string(" AWS_DEFAULT_REGION=") + region))
    + " ";
}

void AWS::cleanup() {
  if (!cleanupCalled) {
    cleanupCalled = true;
    if (!runningInstanceIds.empty()) {
      cerr << "Terminating EC2 instances..." << endl;
      const vguard<string> ids (runningInstanceIds.begin(), runningInstanceIds.end());
      aws.terminateInstancesSilently (ids);
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
  cerr << endl;
  exit (EXIT_FAILURE);
}

void AWS::registerCleanup() {
  for (int sig : { SIGABRT, SIGFPE, SIGILL, SIGINT, SIGSEGV, SIGTERM }) {
    void (*oldHandler)(int) = signal (sig, &AWS::handleSignal);
    if (oldHandler != NULL)
      signalHandler[sig] = oldHandler;
  }
}

AWS::AWS()
  : bashBang ("#!/bin/bash\n"),
    keyPair (AWS_DEFAULT_KEY_PAIR),
    securityGroup (AWS_DEFAULT_SECURITY_GROUP)
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
    if (arg == "-ec2key") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      keyPair = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-ec2group") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      securityGroup = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    }
  }

  return false;
}
