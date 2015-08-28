#include <iostream>
#include <fstream>
#include "../src/qmodel.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <quaffcounts.json>\n";
    exit (EXIT_FAILURE);
  }

  QuaffParamCounts counts;
  ifstream in (argv[1]);
  counts.readJson (in);
  counts.writeJson (cout);

  exit (EXIT_SUCCESS);
}
