#include <iostream>
#include <fstream>
#include "../src/qmodel.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " [-params|-nullparams|-counts] file.yaml > file.json\n";
    exit (EXIT_FAILURE);
  }

  string s (argv[1]);
  if (s == "-params") {
  
    QuaffParams params;
    ifstream in (argv[2]);
    params.read (in);
    params.writeJson (cout);

  } else if (s == "-nullparams") {
  
    QuaffNullParams params;
    ifstream in (argv[2]);
    params.read (in);
    params.writeJson (cout);

  } else if (s == "-counts") {
  
    QuaffParamCounts counts;
    ifstream in (argv[2]);
    counts.read (in);
    counts.writeJson (cout);

  } else {
    Fail ("Unknown conversion: %s", argv[1]);

  }
    
  exit (EXIT_SUCCESS);
}
