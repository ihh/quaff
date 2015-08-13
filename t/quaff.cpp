#include <iostream>
#include <fstream>
#include "../src/qmodel.h"

struct QuaffUsage {
  int& argc;
  char**& argv;
  string prog, text;
  QuaffUsage (int& argc, char**& argv);
  string getCommand();
  bool parseUnknown();
};

struct QuaffParamsIn : QuaffParams {
  int& argc;
  char**& argv;
  bool initialized;

  QuaffParamsIn (int& argc, char**& argv)
    : QuaffParams(),
      argc(argc),
      argv(argv),
      initialized(false)
  { }

  bool parseParamFilename();
  void requireParams();
};

struct SeqList {
  vguard<string> filenames;
  string type, tag;
  int& argc;
  char**& argv;
  bool wantFastq, wantRevcomps;
  vguard<FastSeq> seqs;

  SeqList (int& argc, char**& argv, const char* type, const char* tag)
    : argc(argc),
      argv(argv),
      type(type),
      tag(tag),
      wantFastq(false),
      wantRevcomps(false)
  { }

  bool parseSeqFilename();
  bool parseFwdStrand();
  void loadSequences();
};

int main (int argc, char** argv) {

  QuaffUsage usage (argc, argv);
  const string command = usage.getCommand();

  QuaffParamsIn params (argc, argv);
  QuaffParamCounts prior;
  prior.initCounts (1);

  SeqList refs (argc, argv, "reference", "-fasta");
  refs.wantRevcomps = true;
  
  SeqList reads (argc, argv, "read", "-fastq");
  reads.wantFastq = true;

  if (command == "align") {
    QuaffAligner aligner;
    while (logger.parseLogArgs (argc, argv)
	   || aligner.parseAlignmentArgs (argc, argv)
	   || params.parseParamFilename()
	   || refs.parseSeqFilename()
	   || reads.parseSeqFilename()
	   || refs.parseFwdStrand()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    refs.loadSequences();
    params.requireParams();

    aligner.align (cout, refs.seqs, reads.seqs, params);

  } else if (command == "train") {
    QuaffTrainer trainer;
    while (logger.parseLogArgs (argc, argv)
	   || trainer.parseTrainingArgs (argc, argv)
	   || params.parseParamFilename()
	   || refs.parseSeqFilename()
	   || reads.parseSeqFilename()
	   || refs.parseFwdStrand()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    refs.loadSequences();

    QuaffParams newParams = trainer.fit (refs.seqs, reads.seqs, params, prior);
    newParams.write (cout);
    
  } else {
    cerr << usage.text << "Unrecognized command: " << command << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

bool QuaffParamsIn::parseParamFilename() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-params") {
      Assert (argc > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argv[1]);
      Assert (!inFile.fail(), "Couldn't open %s", argv[1]);
      read (inFile);
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

void QuaffParamsIn::requireParams() {
  Assert (initialized, "Please specify a parameter file using -params");
}

bool SeqList::parseSeqFilename() {
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

bool SeqList::parseFwdStrand() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-fwdstrand") {
      wantRevcomps = false;
      argc -= 1;
      argv += 1;
      return true;
    }
  }
  return false;
}

void SeqList::loadSequences() {
  Assert (filenames.size() > 0, "Please specify at least one %s file using %s", type.c_str(), tag.c_str());

  for (const auto& s : filenames) {
    const vguard<FastSeq> fsvec = readFastSeqs (s.c_str());
    if (wantFastq)
      for (const auto& fs: fsvec)
	Assert (fs.hasQual(), "Sequence %s in file %s does not have quality scores", fs.name.c_str(), s.c_str());
    seqs.insert (seqs.end(), fsvec.begin(), fsvec.end());
  }

  if (wantRevcomps)
    addRevcomps (seqs);
}


QuaffUsage::QuaffUsage (int& argc, char**& argv)
  : argc(argc),
    argv(argv),
    prog (argv[0])
{
  text = "Usage: " + prog + " {align,train} [options]\n"
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
}

string QuaffUsage::getCommand() {
  if (argc < 2) {
    cerr << text;
    exit (EXIT_FAILURE);
  }
  const string command (argv[1]);
  argv += 2;
  argc -= 2;
  return command;
}

bool QuaffUsage::parseUnknown() {
  if (argc > 0) {
    cerr << text << "Unknown option: " << argv[0] << endl;
    Abort ("Error parsing command-line options");
  }
  return false;
}
