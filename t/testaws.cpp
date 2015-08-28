#include <iostream>
#include <fstream>
#include <unistd.h>
#include "../src/aws.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <key-pair>\n";
    exit (EXIT_FAILURE);
  }

  aws.keyPair = argv[1];

  const string script = aws.bashBang
    + "mkfifo /tmp/fifo\n"
    + "cat /tmp/fifo | nc -k -v -l 8000 | cat > /tmp/fifo\n";

  logger.setVerbose (100);
  const vguard<string> ids = aws.launchInstancesWithScript (1, string(AWS_DEFAULT_INSTANCE_TYPE), string(AWS_DEFAULT_AMI), script);

  if (ids.size() != 1) {
    cerr << "Failed to create instance" << endl;
    exit (EXIT_FAILURE);
  }
  
  const vguard<string> addrs = aws.getInstanceAddresses (ids);
  
  cerr << "Instance running at " << addrs[0] << endl;
  cerr << "The following should connect to an echo server: nc " << addrs[0] << " 8000" << endl;

  const string testEchoCmd = string("echo ok | nc ") + addrs[0] + " 8000";
  cerr << "If so, then `" << testEchoCmd << "` should yield 'ok':" << endl;
  const string testEchoOut = pipeToString (testEchoCmd.c_str());
  if (testEchoOut == "ok\n") {
    cout << testEchoOut << flush;
    cerr << "(neat, huh)" << endl;
  } else {
    cerr << testEchoOut << flush;
    cout << "not ok" << endl;
  }

  cerr << "Sleeping for 120 seconds (try killing this process to check termination handling)" << endl;

  usleep(120000000);
  
  exit (EXIT_SUCCESS);
}
