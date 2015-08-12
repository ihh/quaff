#include <iostream>
#include "../src/qmodel.h"

bool parseUnknown (int argc, char** argv, const string& usage) {
  if (argc > 0) {
    cerr << usage << "Unknown option: " << argv[0] << endl;
    Abort ("Error parsing command-line options");
  }
  return false;
}

bool parseParamFilename (int argc, char** argv, QuaffParams& qp) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-params") {
      Assert (argc > 1, "%s needs an argument", arg.c_str());
      qp.read (argv[1]);
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

bool parseSeqFilename (int argc, char** argv, const char* tag, vector<string>& filenames) {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == tag) {
      Assert (argc > 1, "%s needs an argument", arg.c_str());
      filenames.push_back (string (argv[1]));
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

vector<FastSeq> loadSequences (const vector<string>& filenames) {
  vector<FastSeq> seqs;
  for (auto s : filenames) {
    const vector<FastSeq> fs = readFastSeqs (s.c_str());
    seqs.insert (seqs.end(), fs.begin(), fs.end());
  }
  return seqs;
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

  QuaffParams params;
  QuaffTrainer trainer;

  const string command (argv[1]);
  argv += 2;
  argc -= 2;

  vector<string> refFilenames, readFilenames;
  
  if (command == "align") {
    while (logger.parseLogArgs (argc, argv)
	   || parseParamFilename (argc, argv, params)
	   || parseSeqFilename (argc, argv, "-fasta", refFilenames)
	   || parseSeqFilename (argc, argv, "-fastq", readFilenames)
	   || parseUnknown (argc, argv, usage))
      { }

    const vector<FastSeq> reads = loadSequences (readFilenames);
    const vector<FastSeq> refs = loadSequences (refFilenames);

    // WRITE ME
    
  } else if (command == "train") {
    while (logger.parseLogArgs (argc, argv)
	   || trainer.parseTrainingArgs (argc, argv)
	   || parseParamFilename (argc, argv, params)
	   || parseSeqFilename (argc, argv, "-fasta", refFilenames)
	   || parseSeqFilename (argc, argv, "-fastq", readFilenames)
	   || parseUnknown (argc, argv, usage))
      { }

    const vector<FastSeq> reads = loadSequences (readFilenames);
    const vector<FastSeq> refs = loadSequences (refFilenames);

    // WRITE ME
    
  } else {
    cerr << usage << "Unrecognized command: " << command << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
