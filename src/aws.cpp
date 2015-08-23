#include "aws.h"

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
