#include <iostream>
#include "../src/qmodel.h"

bool parseUnknown (int argc, char** argv, const string& usage) {
  if (argc > 0) {
    cerr << usage << "Unknown option: " << argv[0] << endl;
    Abort ("Error parsing command-line options");
  }
  return false;
}

int main (int argc, char** argv) {

  const string prog (argv[0]);
  const string usage = "Usage: " + prog + " {align,train} [options]\n"
    + "\n"
    + " " + prog + " train [-params seed-params.yaml] -fasta refs.fasta -fastq reads.fastq  >trained-params.yaml\n"
    + "  (to fit a model to unaligned FASTQ reads/FASTA refs, using EM)\n"
    + "\n"
    + " " + prog + " align -params params.yaml -fasta refs.fasta -fastq reads.fastq\n"
    + "  (to align FASTQ reads to FASTA reference sequences, using Viterbi)\n"
    + "\n"
    + "Other options (some are only available in certain command modes):\n"
    + " -verbose, -vv, -vvv, -v4, etc.\n"
    + " -log <function_name>\n"
    + "                 various levels of logging\n";

  if (argc < 2) {
    cerr << usage;
    exit (EXIT_FAILURE);
  }

  QuaffTrainer trainer;

  const string command (argv[1]);
  argv += 2;
  argc -= 2;

  if (command == "align") {
    while (logger.parseLogArgs (argc, argv)
	   || parseUnknown (argc, argv, usage))
      { }

  } else if (command == "train") {
    while (logger.parseLogArgs (argc, argv)
	   || trainer.parseTrainingArgs (argc, argv)
	   || parseUnknown (argc, argv, usage))
      { }

  } else {
    cerr << usage << "Unrecognized command: " << command << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
