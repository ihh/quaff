#include "aws.h"
#include "gason.h"
#include "util.h"

JsonValue* jsonFind (JsonValue& parent, const char* key);
JsonValue& jsonFindOrDie (JsonValue& parent, const char* key);
  
JsonValue* jsonFind (JsonValue& parent, const char* key) {
  Assert (parent.getTag() == JSON_OBJECT, "JSON value is not an object: %s", parent.toString());
  for (auto i : parent)
    if (strcmp (i->key, key) == 0)
      return &i->value;
  return NULL;
}

JsonValue& jsonFindOrDie (JsonValue& parent, const char* key) {
  JsonValue* val = jsonFind (parent, key);
  Assert (val != NULL, "Couldn't find %s in JSON %s", key, parent.toString());
  return *val;
}

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

vguard<string> AWS::launchInstancesWithScript (unsigned int nInstances, const string& keyPairName, const string& instanceType, const string& ami, const string& userDataScript) {
  vguard<string> ids;
  TempFile tempfile (userDataScript, "userdata");
  string runCmd = string("aws ec2 run-instances --enable-api-termination --key-name " + keyPairName + " --instance-type " + instanceType + " --count " + to_string(nInstances) + ':' + to_string(nInstances) + " --image-id " + ami + " --user-data " + tempfile.fullPath);
  string runOut = pipeToString (runCmd.c_str());
  char* runOutCStr = new char [runOut.length()+1];
  strcpy (runOutCStr, runOut.c_str());
  char *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  int status = jsonParse(runOutCStr, &endPtr, &value, allocator);
  Assert (status == JSON_OK, "JSON parsing error: %s at byte %zd of output of %s\n%s\n", jsonStrError(status), endPtr - runOutCStr, runCmd.c_str(), endPtr);
  JsonValue& instanceList = jsonFindOrDie (value, "Instances");
  Assert (instanceList.getTag() == JSON_ARRAY, "JSON parsing error: Instances is not a list");
  for (auto inst : instanceList) {
    JsonValue& instanceId = jsonFindOrDie (inst->value, "InstanceId");
    Assert (instanceId.getTag() == JSON_STRING, "JSON parsing error: InstanceId is not a string");
    ids.push_back (instanceId.toString());
  }
  delete[] runOutCStr;
  // TODO: wait for instances
  // aws wait instance-running --instance-ids ...
  return ids;
}

string AWS::getInstancePublicDNS (const string& instanceID) {
  string dns;
  // WRITE ME
  return dns;
}

void AWS::terminateInstance (const string& instanceID) {
  // WRITE ME
}
