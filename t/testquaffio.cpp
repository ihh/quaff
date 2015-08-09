#include <iostream>
#include <fstream>
#include "../src/quaff.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <quaffparams.txt>\n";
    exit (EXIT_FAILURE);
  }

  QuaffParams params;
  ifstream in (argv[1]);
  params.read (in);
  params.write (cout);

  exit (EXIT_SUCCESS);
}
