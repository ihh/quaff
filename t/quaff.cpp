#include <iostream>
#include <fstream>
#include <regex>
#include "../src/qmodel.h"
#include "../src/qoverlap.h"
#include "../src/logger.h"

struct QuaffUsage {
  int& argc;
  char**& argv;
  string prog, briefText, text;
  QuaffUsage (int& argc, char**& argv);
  string getCommand();
  bool parseUnknown();
};

struct SeqList {
  vguard<string> filenames;
  string type, tag;
  regex tagRegex;
  int& argc;
  char**& argv;
  bool wantQualScores, wantRevcomps;
  vguard<FastSeq> seqs;
  size_t nOriginals;  // number of seqs that are NOT revcomps
  
  SeqList (int& argc, char**& argv, const char* type, const char* tag, const char* tagRegex)
    : argc(argc),
      argv(argv),
      type(type),
      tag(tag),
      tagRegex(tagRegex),
      wantQualScores(false),
      wantRevcomps(false),
      nOriginals(0)
  { }

  bool parseSeqFilename();
  bool parseRevcompArgs();
  bool parseQualScoreArgs();
  void loadSequences();
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
  void requireParamsOrUsePrior (const QuaffParamCounts& prior);
};

struct QuaffNullParamsIn : QuaffNullParams {
  int& argc;
  char**& argv;
  bool initialized;
  string saveFilename;
  
  QuaffNullParamsIn (int& argc, char**& argv)
    : QuaffNullParams(),
      argc(argc),
      argv(argv),
      initialized(false)
  { }

  bool parseNullModelFilename();
  void requireNullModelOrFit (const SeqList& seqList);
};

struct QuaffPriorIn : QuaffParamCounts {
  int& argc;
  char**& argv;
  bool initialized, kmerLenSpecified;
  string saveFilename;

  QuaffPriorIn (int& argc, char**& argv)
    : QuaffParamCounts(1),  // fix initial kmerLen at 1
      argc(argc),
      argv(argv),
      initialized(false),
      kmerLenSpecified(false)
  { }

  bool parsePriorArgs();
  void requirePrior();
  void requirePriorOrUseNullModel (const QuaffNullParams& nullModel, const QuaffParamsIn& params);
};

int main (int argc, char** argv) {

  QuaffUsage usage (argc, argv);
  const string command = usage.getCommand();

  QuaffParamsIn params (argc, argv);

  SeqList refs (argc, argv, "reference", "-ref", "^-(ref|refs|fasta)$");
  refs.wantRevcomps = true;

  SeqList reads (argc, argv, "read", "-read", "^-(read|reads|fastq)$");
  reads.wantQualScores = true;

  QuaffDPConfig config;

  if (command == "align") {
    QuaffAligner aligner;
    QuaffNullParamsIn nullModel (argc, argv);
    while (logger.parseLogArgs (argc, argv)
	   || aligner.parseAlignmentArgs (argc, argv)
	   || config.parseConfigArgs (argc, argv)
	   || params.parseParamFilename()
	   || nullModel.parseNullModelFilename()
	   || refs.parseSeqFilename()
	   || refs.parseRevcompArgs()
	   || reads.parseSeqFilename()
	   || reads.parseQualScoreArgs()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    refs.loadSequences();
    params.requireParams();
    nullModel.requireNullModelOrFit (reads);
    
    aligner.align (cout, refs.seqs, reads.seqs, params, nullModel, config);

  } else if (command == "train") {
    QuaffTrainer trainer;
    QuaffNullParamsIn nullModel (argc, argv);
    QuaffPriorIn prior (argc, argv);
    while (logger.parseLogArgs (argc, argv)
	   || trainer.parseTrainingArgs (argc, argv)
	   || config.parseConfigArgs (argc, argv)
	   || params.parseParamFilename()
	   || nullModel.parseNullModelFilename()
	   || prior.parsePriorArgs()
	   || refs.parseSeqFilename()
	   || refs.parseRevcompArgs()
	   || reads.parseSeqFilename()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    refs.loadSequences();

    nullModel.requireNullModelOrFit (reads);
    prior.requirePriorOrUseNullModel (nullModel, params);
    params.requireParamsOrUsePrior (prior);

    QuaffParams newParams = trainer.fit (refs.seqs, reads.seqs, params, nullModel, prior, config);
    newParams.write (cout);

  } else if (command == "overlap") {
    QuaffOverlapAligner aligner;
    QuaffNullParamsIn nullModel (argc, argv);
    reads.wantRevcomps = true;
    while (logger.parseLogArgs (argc, argv)
	   || aligner.parseAlignmentArgs (argc, argv)
	   || config.parseOverlapConfigArgs (argc, argv)
	   || params.parseParamFilename()
	   || nullModel.parseNullModelFilename()
	   || reads.parseSeqFilename()
	   || reads.parseRevcompArgs()
	   || reads.parseQualScoreArgs()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    params.requireParams();
    nullModel.requireNullModelOrFit (reads);

    aligner.align (cout, reads.seqs, reads.nOriginals, params, nullModel, config);

  } else if (command == "help" || command == "-help" || command == "--help" || command == "-h") {
    cout << usage.text;
    return EXIT_SUCCESS;
    
  } else {
    cerr << usage.briefText << "Unrecognized command: " << command << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

bool QuaffParamsIn::parseParamFilename() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-params") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argv[1]);
      Require (!inFile.fail(), "Couldn't open %s", argv[1]);
      read (inFile);
      initialized = true;
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

void QuaffParamsIn::requireParams() {
  Require (initialized, "Please specify a parameter file using -params");
}

void QuaffParamsIn::requireParamsOrUsePrior (const QuaffParamCounts& prior) {
  if (!initialized) {
    if (LogThisAt(1))
      cerr << "Auto-initializing model parameters from prior" << endl;
    (QuaffParams&) *this = prior.fit();
  }
}

bool QuaffNullParamsIn::parseNullModelFilename() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-null") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argv[1]);
      Require (!inFile.fail(), "Couldn't open %s", argv[1]);
      read (inFile);
      initialized = true;
      argc -= 2;
      argv += 2;
      return true;

    } else if (arg == "-savenull") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      saveFilename = argv[1];
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

void QuaffNullParamsIn::requireNullModelOrFit (const SeqList& seqList) {
  if (!initialized) {
    if (LogThisAt(1))
      cerr << "Auto-optimizing null model for read sequences" << endl;
    (QuaffNullParams&) *this = QuaffNullParams (seqList.seqs);
  }
  if (saveFilename.size()) {
    ofstream out (saveFilename);
    write (out);
  }
}

bool QuaffPriorIn::parsePriorArgs() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-prior") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argv[1]);
      Require (!inFile.fail(), "Couldn't open %s", argv[1]);
      read (inFile);
      initialized = true;
      argc -= 2;
      argv += 2;
      return true;

    } else if (arg == "-order") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      initKmerContext (atoi (argv[1]));
      resize();
      kmerLenSpecified = true;
      argc -= 2;
      argv += 2;
      return true;

    } else if (arg == "-saveprior") {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      saveFilename = argv[1];
      argc -= 2;
      argv += 2;
      return true;

    }
  }
  return false;
}

void QuaffPriorIn::requirePriorOrUseNullModel (const QuaffNullParams& nullModel, const QuaffParamsIn& params) {
  if (initialized) {
    if (params.initialized)
      Require (kmerLen == params.kmerLen, "Model order in prior file (%d) does not match model order in parameter file (%d)", kmerLen, params.kmerLen);
  } else {
    if (LogThisAt(1))
      cerr << "Auto-setting prior from null model" << endl;
    if (params.initialized) {
      if (kmerLenSpecified)
	Require (kmerLen == params.kmerLen, "Model order specified on command-line (%d) does not match model order in parameter file (%d)", kmerLen, params.kmerLen);
      else {
	initKmerContext (params.kmerLen);
	resize();
      }
    }
    initCounts (9, 9, 5, 1, &nullModel);
  }
  if (saveFilename.size()) {
    ofstream out (saveFilename);
    write (out);
  }
}

bool SeqList::parseSeqFilename() {
  if (argc > 0) {
    const string arg = argv[0];
    if (regex_match (arg, tagRegex)) {
      Require (argc > 1, "%s needs an argument", arg.c_str());
      filenames.push_back (string (argv[1]));
      argc -= 2;
      argv += 2;
      return true;
    }
  }
  return false;
}

bool SeqList::parseRevcompArgs() {
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

bool SeqList::parseQualScoreArgs() {
  if (argc > 0) {
    const string arg = argv[0];
    if (arg == "-noquals") {
      wantQualScores = false;
      argc -= 1;
      argv += 1;
      return true;
    }
  }
  return false;
}

void SeqList::loadSequences() {
  Require (filenames.size() > 0, "Please specify at least one %s file using %s", type.c_str(), tag.c_str());

  for (const auto& s : filenames) {
    vguard<FastSeq> fsvec = readFastSeqs (s.c_str());
    for (auto& fs: fsvec) {
      if (wantQualScores)
	Require (fs.hasQual(), "Sequence %s in file %s does not have quality scores", fs.name.c_str(), s.c_str());
      else
	fs.qual.clear();
      if (fs.length())  // skip zero-length sequences
	seqs.push_back (fs);
    }
  }

  nOriginals = seqs.size();
  if (wantRevcomps)
    addRevcomps (seqs);
}

QuaffUsage::QuaffUsage (int& argc, char**& argv)
  : argc(argc),
    argv(argv),
    prog (argv[0])
{
  briefText = "Usage: " + prog + " {help,train,align,overlap} [options]\n";
  
  text = briefText
    + "\n"
    + "Commands:\n"
    + "\n"
    + " " + prog + " train -ref refs.fasta -read reads.fastq  >params.yaml\n"
    + "  (to fit a model to unaligned sequences, using EM)\n"
    + "\n"
    + "   -params <file>  Optional initial parameters\n"
    + "   -maxiter <n>    Max number of EM iterations\n"
    + "   -mininc <n>     EM convergence threshold (relative log-likelihood increase)\n"
    + "   -force          Force each read to match a refseq, i.e. disallow null model\n"
    + "   -order <k>      Allow substitutions to depend on k-mer contexts\n"
    + "   -prior <file>, -saveprior <file>\n"
    + "                   Respectively: load/save prior pseudocounts from/to file\n"
    + "   -counts <file>  Save E-step counts to file, which can then be used as a prior\n"
    + "   -countswithprior <file>\n"
    + "                   Like -counts, but adds in prior pseudocounts as well\n"
    + "\n"
    + "\n"
    + " " + prog + " align -params params.yaml -ref refs.fasta -read reads.fastq\n"
    + "  (to align FASTQ reads to FASTA reference sequences, using Viterbi)\n"
    + "\n"
    + "   -printall       Print all pairwise alignments, not just best for each read\n"
    + "\n"
    + "\n"
    + " " + prog + " overlap -params params.yaml -read reads.fastq\n"
    + "  (to find overlaps between FASTQ reads, using Viterbi)\n"
    + "\n"
    + "\n"
    + "Alignment options (align/overlap commands):\n"
    + "   -format {fasta,stockholm,refseq}\n"
    + "                   Alignment output format\n"
    + "   -threshold <n>\n"
    + "   -nothreshold    Log-odds ratio score threshold for alignment reporting\n"
    + "   -noquals        Ignore read quality scores during alignment\n"
    + "\n"
    + "General options (all commands, except where indicated):\n"
    + "   -verbose, -vv, -vvv, -v4, etc.\n"
    // uncomment to document debug logging:
    //    + "   -log <function_name>\n"
    + "                   Various levels of logging\n"
    + "   -fwdstrand      Do not include reverse-complemented sequences\n"
    + "   -global         Force all of refseq to be aligned (align/train only)\n"
    + "   -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "   -kmatchn <n>    Threshold# of kmer matches to include a diagonal (default " + to_string(DEFAULT_KMER_THRESHOLD) + ")\n"
    // uncomment to document this uncertain, experimental option:
    //    + "   -kmatchsd <n>   Set kmer threshold to n standard deviations above background\n"
    + "   -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "   -dense          Do full DP, not just kmer-matching diagonals (memory hog!)\n"
    + "   -null <file>, -savenull <file>\n"
    + "                   Respectively: load/save null model from/to file\n"
    + "\n";
}

string QuaffUsage::getCommand() {
  if (argc < 2) {
    cerr << briefText;
    exit (EXIT_FAILURE);
  }
  const string command (argv[1]);
  argv += 2;
  argc -= 2;
  return command;
}

bool QuaffUsage::parseUnknown() {
  if (argc > 0) {
    string arg (argv[0]);
    if (arg == "-abort") {
      // test stack trace
      Abort ("abort triggered");

    } else {
      cerr << text << "Unknown option: " << arg << endl;
      cerr << "Error parsing command-line options\n";
      exit (EXIT_FAILURE);

    }
  }
  return false;
}
