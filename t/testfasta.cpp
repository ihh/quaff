#include <iostream>
#include "../src/fastseq.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <seqs.fasta>\n";
    exit (EXIT_FAILURE);
  }

  vector<FastSeq> seqs = readFastSeqs (argv[1]);
  writeFastaSeqs (cout, seqs);

  exit (EXIT_SUCCESS);
}
