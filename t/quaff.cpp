#include <iostream>
#include <fstream>
#include <regex>
#include <deque>
#include "../src/qmodel.h"
#include "../src/qoverlap.h"
#include "../src/logger.h"
#include "../src/defaultparams.h"

// allow a different kmer-matching threshold for refseq-read alignment
#define DEFAULT_REFSEQ_KMER_THRESHOLD 20

struct QuaffUsage {
  deque<string>& argvec;
  string prog, briefText, text;
  deque<string> implicitSwitches;
  bool unlimitImplicitSwitches;
  QuaffUsage (deque<string>& argvec);
  string getCommand();
  bool parseUnknown();
};

struct SeqList {
  vguard<string> filenames;
  string type, tag;
  regex tagRegex;
  deque<string>& argvec;
  bool wantQualScores, wantRevcomps;
  vguard<FastSeq> seqs;
  size_t nOriginals;  // number of seqs that are NOT revcomps
  
  SeqList (deque<string>& argvec, const char* type, const char* tag, const char* tagRegex)
    : argvec(argvec),
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
  deque<string>& argvec;
  bool initialized;

  QuaffParamsIn (deque<string>& argvec)
    : QuaffParams(),
      argvec(argvec),
      initialized(false)
  { }

  bool parseParamFilename();
  void requireParamsOrUseDefaults();
  void requireParamsOrUsePrior (const QuaffParamCounts& prior);
};

struct QuaffNullParamsIn : QuaffNullParams {
  deque<string>& argvec;
  bool initialized;
  string saveFilename;
  
  QuaffNullParamsIn (deque<string>& argvec)
    : QuaffNullParams(),
      argvec(argvec),
      initialized(false)
  { }

  bool parseNullModelFilename();
  void requireNullModelOrFit (const SeqList& seqList);
};

struct QuaffPriorIn : QuaffParamCounts {
  deque<string>& argvec;
  bool initialized, kmerLenSpecified;
  string saveFilename;

  QuaffPriorIn (deque<string>& argvec)
    : QuaffParamCounts(1,0),  // fix initial (matchKmerLen,indelKmerLen) at (1,0)
      argvec(argvec),
      initialized(false),
      kmerLenSpecified(false)
  { }

  bool parsePriorArgs();
  void requirePrior();
  void requirePriorOrUseNullModel (const QuaffNullParams& nullModel, const QuaffParamsIn& params);
};

int main (int argc, char** argv) {

  deque<string> argvec (argc);
  for (int n = 0; n < argc; ++n)
    argvec[n] = argv[n];
  
  QuaffUsage usage (argvec);
  const string command = usage.getCommand();

  QuaffParamsIn params (argvec);

  SeqList refs (argvec, "reference", "-ref", "^-(ref|refs|fasta)$");
  refs.wantRevcomps = true;

  SeqList reads (argvec, "read", "-read", "^-(read|reads|fastq)$");
  reads.wantQualScores = true;

  QuaffDPConfig config;

  if (command == "align") {
    QuaffAligner aligner;
    QuaffNullParamsIn nullModel (argvec);
    usage.implicitSwitches.push_back (string ("-ref"));
    usage.implicitSwitches.push_back (string ("-read"));
    usage.unlimitImplicitSwitches = true;
    config.kmerThreshold = DEFAULT_REFSEQ_KMER_THRESHOLD;
    while (logger.parseLogArgs (argvec)
	   || aligner.parseAlignmentArgs (argvec)
	   || config.parseRefSeqConfigArgs (argvec)
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
    params.requireParamsOrUseDefaults();
    nullModel.requireNullModelOrFit (reads);
    
    aligner.align (cout, refs.seqs, reads.seqs, params, nullModel, config);

  } else if (command == "train") {
    QuaffTrainer trainer;
    QuaffNullParamsIn nullModel (argvec);
    QuaffPriorIn prior (argvec);
    usage.implicitSwitches.push_back (string ("-ref"));
    usage.implicitSwitches.push_back (string ("-read"));
    usage.unlimitImplicitSwitches = true;
    config.kmerThreshold = DEFAULT_REFSEQ_KMER_THRESHOLD;
    while (logger.parseLogArgs (argvec)
	   || trainer.parseTrainingArgs (argvec)
	   || config.parseRefSeqConfigArgs (argvec)
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

  } else if (command == "count") {
    QuaffTrainer trainer;
    QuaffNullParamsIn nullModel (argvec);
    usage.implicitSwitches.push_back (string ("-ref"));
    usage.implicitSwitches.push_back (string ("-read"));
    usage.unlimitImplicitSwitches = true;
    config.kmerThreshold = DEFAULT_REFSEQ_KMER_THRESHOLD;
    while (logger.parseLogArgs (argvec)
	   || trainer.parseCountingArgs (argvec)
	   || config.parseRefSeqConfigArgs (argvec)
	   || params.parseParamFilename()
	   || nullModel.parseNullModelFilename()
	   || refs.parseSeqFilename()
	   || refs.parseRevcompArgs()
	   || reads.parseSeqFilename()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    refs.loadSequences();

    nullModel.requireNullModelOrFit (reads);
    params.requireParamsOrUseDefaults();

    QuaffParamCounts counts = trainer.getCounts (refs.seqs, reads.seqs, params, nullModel, config);
    counts.write (cout);

  } else if (command == "overlap") {
    QuaffOverlapAligner aligner;
    QuaffNullParamsIn nullModel (argvec);
    reads.wantRevcomps = true;
    usage.implicitSwitches.push_back (string ("-read"));
    usage.unlimitImplicitSwitches = true;
    while (logger.parseLogArgs (argvec)
	   || aligner.parseAlignmentArgs (argvec)
	   || config.parseGeneralConfigArgs (argvec)
	   || params.parseParamFilename()
	   || nullModel.parseNullModelFilename()
	   || reads.parseSeqFilename()
	   || reads.parseRevcompArgs()
	   || reads.parseQualScoreArgs()
	   || usage.parseUnknown())
      { }

    reads.loadSequences();
    params.requireParamsOrUseDefaults();
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
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-params") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argvec[1]);
      Require (!inFile.fail(), "Couldn't open %s", argvec[1].c_str());
      read (inFile);
      initialized = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

void QuaffParamsIn::requireParamsOrUseDefaults() {
  if (!initialized)
    (QuaffParams&) *this = defaultQuaffParams();
}

void QuaffParamsIn::requireParamsOrUsePrior (const QuaffParamCounts& prior) {
  if (!initialized) {
    if (LogThisAt(1))
      cerr << "Auto-initializing model parameters from prior" << endl;
    (QuaffParams&) *this = prior.fit();
  }
}

bool QuaffNullParamsIn::parseNullModelFilename() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-null") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argvec[1]);
      Require (!inFile.fail(), "Couldn't open %s", argvec[1].c_str());
      read (inFile);
      initialized = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savenull") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      saveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
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
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-prior") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      ifstream inFile (argvec[1]);
      Require (!inFile.fail(), "Couldn't open %s", argvec[1].c_str());
      read (inFile);
      initialized = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-order") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      matchContext.initKmerContext (atoi (argvec[1].c_str()));
      resize();
      kmerLenSpecified = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-gaporder") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      indelContext.initKmerContext (atoi (argvec[1].c_str()));
      resize();
      kmerLenSpecified = true;
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveprior") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      saveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    }
  }
  return false;
}

void QuaffPriorIn::requirePriorOrUseNullModel (const QuaffNullParams& nullModel, const QuaffParamsIn& params) {
  if (initialized) {
    if (params.initialized) {
      Require (matchContext.kmerLen == params.matchContext.kmerLen, "Order of match dependence in prior file (%d) does not match order in parameter file (%d)", matchContext.kmerLen, params.matchContext.kmerLen);
      Require (indelContext.kmerLen == params.indelContext.kmerLen, "Order of indel dependence in prior file (%d) does not match order in parameter file (%d)", indelContext.kmerLen, params.indelContext.kmerLen);
    }
  } else {
    if (LogThisAt(1))
      cerr << "Auto-setting prior from null model" << endl;
    if (params.initialized) {
      if (kmerLenSpecified) {
	Require (matchContext.kmerLen == params.matchContext.kmerLen, "Order of match dependence specified on command line (%d) does not match order in parameter file (%d)", matchContext.kmerLen, params.matchContext.kmerLen);
	Require (indelContext.kmerLen == params.indelContext.kmerLen, "Order of indel dependence specified on command line (%d) does not match order in parameter file (%d)", indelContext.kmerLen, params.indelContext.kmerLen);
      } else {
	matchContext.initKmerContext (params.matchContext.kmerLen);
	indelContext.initKmerContext (params.indelContext.kmerLen);
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
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (regex_match (arg, tagRegex)) {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      filenames.push_back (string (argvec[1]));
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

bool SeqList::parseRevcompArgs() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-fwdstrand") {
      wantRevcomps = false;
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

bool SeqList::parseQualScoreArgs() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-noquals") {
      wantQualScores = false;
      argvec.pop_front();
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

  Require (seqs.size() > 0, "Please specify a valid %s file using %s", type.c_str(), tag.c_str());
}

QuaffUsage::QuaffUsage (deque<string>& argvec)
  : argvec(argvec),
    prog(argvec[0]),
    unlimitImplicitSwitches(false)
{
  // leave "count" command undocumented... it's really for the paper
  //  briefText = "Usage: " + prog + " {help,train,count,align,overlap} [options]\n";
  briefText = "Usage: " + prog + " {help,train,align,overlap} [options]\n";
  
  text = briefText
    + "\n"
    + "Commands:\n"
    + "\n"
    + " " + prog + " train refs.fasta reads.fastq  >params.yaml\n"
    + "  (to fit a model to unaligned sequences, using EM/Forward-Backward)\n"
    + "\n"
    + "   -maxiter <n>    Max number of EM iterations (default is " + to_string(QuaffMaxEMIterations) + ")\n"
    + "   -mininc <n>     EM convergence threshold as relative log-likelihood increase\n"
    + "                    (default is " + TOSTRING(QuaffMinEMLogLikeInc) + ")\n"
    + "   -force          Force each read to match a refseq, i.e. disallow null model\n"
    + "   -order <k>      Allow substitutions to depend on k-mer contexts\n"
    + "   -gaporder <k>   Allow gap open probabilities to depend on k-mer contexts\n"
    + "   -prior <file>, -saveprior <file>\n"
    + "                   Respectively: load/save prior pseudocounts from/to file\n"
    + "   -counts <file>  Save E-step counts to file, which can then be used as a prior\n"
    + "   -countswithprior <file>\n"
    + "                   Like -counts, but adds in prior pseudocounts as well\n"
    + "\n"
    + "\n"
    // uncomment to document "count" command:
    //    + " " + prog + " count refs.fasta reads.fastq\n"
    //    + "  (to get summary statistics, i.e. E-step counts, using Forward-Backward)\n"
    //    + "\n"
    //    + "\n"
    + " " + prog + " align refs.fasta reads.fastq\n"
    + "  (to align FASTQ reads to FASTA reference sequences, using Viterbi)\n"
    + "\n"
    + "   -printall       Print all pairwise alignments, not just best for each read\n"
    + "\n"
    + "\n"
    + " " + prog + " overlap reads.fastq\n"
    + "  (to find overlaps between FASTQ reads, using Viterbi)\n"
    + "\n"
    + "\n"
    + "Alignment options (for align/overlap commands):\n"
    + "   -format {fasta,stockholm,refseq}\n"
    + "                   Alignment output format\n"
    + "   -threshold <n>\n"
    + "   -nothreshold    Log-odds ratio score threshold for alignment reporting\n"
    + "   -noquals        Ignore read quality scores during alignment\n"
    + "\n"
    + "General options (for all commands, except where indicated):\n"
    + "   -verbose, -vv, -vvv, -v4, etc.\n"
    // uncomment to document debug logging:
    //    + "   -log <function_name>\n"
    + "                   Various levels of logging\n"
    + "   -params <file>  Load model parameters from file\n"
    + "   -ref <file>     Load additional FASTA reference sequences\n"
    + "   -read <file>    Load additional FASTQ read sequences\n"
    + "   -fwdstrand      Do not include reverse-complemented sequences\n"
    + "   -global         Force all of refseq to be aligned (align/train only)\n"
    + "   -null <file>, -savenull <file>\n"
    + "                   Respectively: load/save null model from/to file\n"
    + "\n"
    + "   -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "   -kmatchn <n>    Threshold# of kmer matches to seed a diagonal\n"
    + "                    (default is " + to_string(DEFAULT_KMER_THRESHOLD) + " for overlap, " + to_string(DEFAULT_REFSEQ_KMER_THRESHOLD) + " for align/train)\n"
    + "   -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "   -kmatchmb <M>   Set kmer threshold to use M megabytes of memory\n"
    + "   -kmatchmax      Set kmer threshold to use all available memory (slow)\n"
    + "   -kmatchoff      No kmer threshold, do full DP (swapfile-heavy, v slow)\n"
    + "\n"
    + "   -threads <n>, -maxthreads\n"
    "                     Specify number of threads, or use all cores available\n"
    + "\n";
}

string QuaffUsage::getCommand() {
  if (argvec.size() < 2) {
    cerr << briefText;
    exit (EXIT_FAILURE);
  }
  const string command (argvec[1]);
  argvec.pop_front();
  argvec.pop_front();
  return command;
}

bool QuaffUsage::parseUnknown() {
  if (argvec.size()) {
    string arg (argvec[0]);
    if (arg == "-abort") {
      // test stack trace
      Abort ("abort triggered");

    } else {
      if (arg[0] == '-' || implicitSwitches.empty()) {
	cerr << text << "Unknown option: " << arg << endl;
	cerr << "Error parsing command-line options\n";
	exit (EXIT_FAILURE);

      } else {
	argvec.push_front (implicitSwitches.front());
	if (implicitSwitches.size() > 1 || !unlimitImplicitSwitches)
	  implicitSwitches.pop_front();
	return true;
      }
    }
  }
  return false;
}
