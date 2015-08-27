#include <iostream>
#include <fstream>
#include "../src/qmodel.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <quaffparams.yaml>\n";
    exit (EXIT_FAILURE);
  }

  QuaffParams params;
  ifstream in (argv[1]);
  params.read (in);
  params.writeJson (cout);

  exit (EXIT_SUCCESS);
}