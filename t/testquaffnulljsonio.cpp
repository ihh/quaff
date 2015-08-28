#include <iostream>
#include <fstream>
#include "../src/qmodel.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <quaffnullparams.json>\n";
    exit (EXIT_FAILURE);
  }

  QuaffNullParams params;
  ifstream in (argv[1]);
  params.readJson (in);
  params.writeJson (cout);

  exit (EXIT_SUCCESS);
}
