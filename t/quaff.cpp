#include <iostream>
#include <fstream>
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
  string getCommand (const char* error = NULL);
  bool parseUnknown();
};

struct SeqList {
  vguard<string> filenames;
  string type, tag;
  deque<string>& argvec;
  bool wantQualScores, wantRevcomps;
  vguard<FastSeq> seqs;
  size_t nOriginals;  // number of seqs that are NOT revcomps
  string serverArgs;
  
  SeqList (deque<string>& argvec, const char* type, const char* tag)
    : argvec(argvec),
      type(type),
      tag(tag),
      wantQualScores(false),
      wantRevcomps(false),
      nOriginals(0)
  { }

  bool parseSeqFilename();
  bool parseRevcompArgs();
  bool parseQualScoreArgs();
  void loadSequences (const QuaffDPConfig& config);
  void loadSequencesForAligner (const QuaffDPConfig& config, const QuaffAlignmentPrinter& aligner);
  SeqList& syncBucket (const QuaffDPConfig& config);
  SeqList& addFileArgs (QuaffDPConfig& config);
};

struct QuaffParamsIn : QuaffParams {
  deque<string>& argvec;
  string loadFilename;
  
  QuaffParamsIn (deque<string>& argvec)
    : QuaffParams(),
      argvec(argvec)
  { }

  bool parseParamFilename();
  void loadParams (const QuaffDPConfig& config);
  void requireParamsOrUseDefaults (const QuaffDPConfig& config);
  void requireParamsOrUsePrior (const QuaffDPConfig& config, const QuaffParamCounts& prior);
  bool initialized() const { return !loadFilename.empty(); }
  QuaffParamsIn& syncBucket (const QuaffDPConfig& config);
  QuaffParamsIn& addFileArgs (QuaffDPConfig& config);
};

struct QuaffNullParamsIn : QuaffNullParams {
  deque<string>& argvec;
  string loadFilename, saveFilename;
  
  QuaffNullParamsIn (deque<string>& argvec)
    : QuaffNullParams(),
      argvec(argvec)
  { }

  bool parseNullModelFilename();
  void loadNullModel (const QuaffDPConfig& config);
  void requireNullModelOrFit (const QuaffDPConfig& config, const SeqList& seqList);
  bool initialized() const { return !loadFilename.empty(); }
  QuaffNullParamsIn& syncBucket (const QuaffDPConfig& config);
  QuaffNullParamsIn& addFileArgs (QuaffDPConfig& config);
};

struct QuaffPriorIn : QuaffParamCounts {
  deque<string>& argvec;
  bool kmerLenSpecified;
  string loadFilename, saveFilename;

  QuaffPriorIn (deque<string>& argvec)
    : QuaffParamCounts(1,0),  // fix initial (matchKmerLen,indelKmerLen) at (1,0)
      argvec(argvec),
      kmerLenSpecified(false)
  { }

  bool parsePriorArgs();
  void loadPrior (const QuaffDPConfig& config);
  void requirePrior (const QuaffDPConfig& config);
  void requirePriorOrUseNullModel (const QuaffDPConfig& config, const QuaffNullParams& nullModel, const QuaffParamsIn& params);
  bool initialized() const { return !loadFilename.empty(); }
};

int main (int argc, char** argv) {

  try {
    deque<string> argvec (argc);
    for (int n = 0; n < argc; ++n)
      argvec[n] = argv[n];
  
    QuaffUsage usage (argvec);
    const string command = usage.getCommand();

    QuaffParamsIn params (argvec);

    SeqList refs (argvec, "reference", "-ref");
    refs.wantRevcomps = true;

    SeqList reads (argvec, "read", "-read");
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

      reads.addFileArgs(config).loadSequencesForAligner (config, aligner);
      refs.addFileArgs(config).loadSequencesForAligner (config, aligner);
      params.addFileArgs(config).requireParamsOrUseDefaults (config);
      nullModel.addFileArgs(config).requireNullModelOrFit (config, reads);

      config.setServerArgs ("align", aligner.serverArgs() + refs.serverArgs + reads.serverArgs);
    
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

      reads.addFileArgs(config).loadSequences (config);
      refs.addFileArgs(config).loadSequences (config);

      nullModel.requireNullModelOrFit (config, reads);
      prior.requirePriorOrUseNullModel (config, nullModel, params);
      params.requireParamsOrUsePrior (config, prior);

      config.setServerArgs ("count", trainer.serverArgs() + refs.serverArgs + reads.serverArgs);

      QuaffParams newParams = trainer.fit (refs.seqs, reads.seqs, params, nullModel, prior, config);
      if (!trainer.usingParamOutputFile())
	newParams.writeJson (cout);

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

      reads.addFileArgs(config).loadSequences (config);
      refs.addFileArgs(config).loadSequences (config);

      nullModel.requireNullModelOrFit (config, reads);
      params.requireParamsOrUseDefaults (config);

      config.setServerArgs ("count", trainer.serverArgs() + refs.serverArgs + reads.serverArgs);

      QuaffParamCounts counts = trainer.getCounts (refs.seqs, reads.seqs, params, nullModel, config);
      if (!trainer.usingCountsOutputFile())
	counts.writeJson (cout);
    
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

      reads.addFileArgs(config).loadSequencesForAligner (config, aligner);
      params.addFileArgs(config).requireParamsOrUseDefaults (config);
      nullModel.addFileArgs(config).requireNullModelOrFit (config, reads);

      config.setServerArgs ("overlap", aligner.serverArgs() + reads.serverArgs);

      aligner.align (cout, reads.seqs, reads.nOriginals, params, nullModel, config);

    } else if (command == "server") {
      const string serverCommand = usage.getCommand ("server needs a command");

      if (serverCommand == "align") {
	QuaffAligner aligner;
	QuaffNullParamsIn nullModel (argvec);
	usage.implicitSwitches.push_back (string ("-ref"));
	usage.implicitSwitches.push_back (string ("-read"));
	usage.unlimitImplicitSwitches = true;
	config.kmerThreshold = DEFAULT_REFSEQ_KMER_THRESHOLD;
	while (logger.parseLogArgs (argvec)
	       || aligner.parseAlignmentArgs (argvec)
	       || config.parseRefSeqConfigArgs (argvec)
	       || config.parseServerConfigArgs (argvec)
	       || params.parseParamFilename()
	       || nullModel.parseNullModelFilename()
	       || refs.parseSeqFilename()
	       || refs.parseRevcompArgs()
	       || reads.parseSeqFilename()
	       || reads.parseQualScoreArgs()
	       || usage.parseUnknown())
	  { }

	reads.syncBucket(config).loadSequencesForAligner (config, aligner);
	refs.syncBucket(config).loadSequencesForAligner (config, aligner);
	params.syncBucket(config).requireParamsOrUseDefaults (config);
	nullModel.syncBucket(config).requireNullModelOrFit (config, reads);
    
	aligner.serveAlignments (refs.seqs, reads.seqs, params, nullModel, config);

      } else if (serverCommand == "count") {
	QuaffTrainer trainer;
	usage.implicitSwitches.push_back (string ("-ref"));
	usage.implicitSwitches.push_back (string ("-read"));
	usage.unlimitImplicitSwitches = true;
	config.kmerThreshold = DEFAULT_REFSEQ_KMER_THRESHOLD;
	while (logger.parseLogArgs (argvec)
	       || trainer.parseServerArgs (argvec)
	       || config.parseRefSeqConfigArgs (argvec)
	       || config.parseServerConfigArgs (argvec)
	       || refs.parseSeqFilename()
	       || refs.parseRevcompArgs()
	       || reads.parseSeqFilename()
	       || usage.parseUnknown())
	  { }

	reads.syncBucket(config).loadSequences (config);
	refs.syncBucket(config).loadSequences (config);

	trainer.serveCounts (refs.seqs, reads.seqs, config);
    
      } else if (serverCommand == "overlap") {
	QuaffOverlapAligner aligner;
	QuaffNullParamsIn nullModel (argvec);
	reads.wantRevcomps = true;
	usage.implicitSwitches.push_back (string ("-read"));
	usage.unlimitImplicitSwitches = true;
	while (logger.parseLogArgs (argvec)
	       || aligner.parseAlignmentArgs (argvec)
	       || config.parseGeneralConfigArgs (argvec)
	       || config.parseServerConfigArgs (argvec)
	       || params.parseParamFilename()
	       || nullModel.parseNullModelFilename()
	       || reads.parseSeqFilename()
	       || reads.parseRevcompArgs()
	       || reads.parseQualScoreArgs()
	       || usage.parseUnknown())
	  { }

	reads.syncBucket(config).loadSequencesForAligner (config, aligner);
	params.syncBucket(config).requireParamsOrUseDefaults (config);
	nullModel.syncBucket(config).requireNullModelOrFit (config, reads);

	aligner.serveAlignments (reads.seqs, reads.nOriginals, params, nullModel, config);

      } else {
	cerr << "Unrecognized server command: " << serverCommand << endl;
	return EXIT_FAILURE;
      }
      
    } else if (command == "help" || command == "-help" || command == "--help" || command == "-h") {
      cout << usage.text;
      return EXIT_SUCCESS;
    
    } else {
      cerr << usage.briefText << "Unrecognized command: " << command << endl;
      return EXIT_FAILURE;
    }

  } catch (...) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

bool QuaffParamsIn::parseParamFilename() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-params") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      loadFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

QuaffParamsIn& QuaffParamsIn::syncBucket (const QuaffDPConfig& config) {
  if (initialized())
    config.syncFromBucket (loadFilename);
  return *this;
}

QuaffParamsIn& QuaffParamsIn::addFileArgs (QuaffDPConfig& config) {
  if (initialized())
    config.addFileArg ("-params", loadFilename);
  return *this;
}

void QuaffParamsIn::loadParams (const QuaffDPConfig& config) {
  if (initialized()) {
    ifstream inFile (loadFilename);
    Require (!inFile.fail(), "Couldn't open %s", loadFilename.c_str());
    readJson (inFile);
  }
}

void QuaffParamsIn::requireParamsOrUseDefaults (const QuaffDPConfig& config) {
  loadParams (config);
  if (!initialized()) {
    if (LogThisAt(1))
      logger << "Using default model parameters" << endl;
    (QuaffParams&) *this = defaultQuaffParams();
  }
}

void QuaffParamsIn::requireParamsOrUsePrior (const QuaffDPConfig& config, const QuaffParamCounts& prior) {
  loadParams (config);
  if (!initialized()) {
    if (LogThisAt(1))
      logger << "Auto-initializing model parameters from prior" << endl;
    (QuaffParams&) *this = prior.fit();
  }
}

bool QuaffNullParamsIn::parseNullModelFilename() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-null") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      loadFilename = argvec[1];
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

void QuaffNullParamsIn::loadNullModel (const QuaffDPConfig& config) {
  if (initialized()) {
    ifstream inFile (loadFilename);
    Require (!inFile.fail(), "Couldn't open %s", loadFilename.c_str());
    readJson (inFile);
  }
}

QuaffNullParamsIn& QuaffNullParamsIn::syncBucket (const QuaffDPConfig& config) {
  if (initialized())
    config.syncFromBucket (loadFilename);
  return *this;
}

QuaffNullParamsIn& QuaffNullParamsIn::addFileArgs (QuaffDPConfig& config) {
  if (initialized())
    config.addFileArg ("-null", loadFilename);
  return *this;
}

void QuaffNullParamsIn::requireNullModelOrFit (const QuaffDPConfig& config, const SeqList& seqList) {
  loadNullModel (config);
  if (!initialized()) {
    if (LogThisAt(1))
      logger << "Auto-optimizing null model for read sequences" << endl;
    (QuaffNullParams&) *this = QuaffNullParams (seqList.seqs);
  }
  if (saveFilename.size()) {
    ofstream out (saveFilename);
    writeJson (out);
  }
}

bool QuaffPriorIn::parsePriorArgs() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-prior") {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      loadFilename = argvec[1];
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

void QuaffPriorIn::loadPrior (const QuaffDPConfig& config) {
  if (initialized()) {
    ifstream inFile (loadFilename);
    Require (!inFile.fail(), "Couldn't open %s", loadFilename.c_str());
    readJson (inFile);
  }
}

void QuaffPriorIn::requirePriorOrUseNullModel (const QuaffDPConfig& config, const QuaffNullParams& nullModel, const QuaffParamsIn& params) {
  loadPrior (config);
  if (initialized()) {
    if (params.initialized()) {
      Require (matchContext.kmerLen == params.matchContext.kmerLen, "Order of match dependence in prior file (%d) does not match order in parameter file (%d)", matchContext.kmerLen, params.matchContext.kmerLen);
      Require (indelContext.kmerLen == params.indelContext.kmerLen, "Order of indel dependence in prior file (%d) does not match order in parameter file (%d)", indelContext.kmerLen, params.indelContext.kmerLen);
    }
  } else {
    if (LogThisAt(1))
      logger << "Auto-setting prior from null model" << endl;
    if (params.initialized()) {
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
    writeJson (out);
  }
}

bool SeqList::parseSeqFilename() {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == tag) {
      Require (argvec.size() > 1, "%s needs an argument", arg.c_str());
      filenames.push_back (argvec[1]);
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
      serverArgs += ' ' + arg;
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
      serverArgs += ' ' + arg;
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

SeqList& SeqList::syncBucket (const QuaffDPConfig& config) {
  for (const auto& s : filenames)
    config.syncFromBucket (s);
  return *this;
}

SeqList& SeqList::addFileArgs (QuaffDPConfig& config) {
  for (const auto& s : filenames)
    config.addFileArg (tag.c_str(), s);
  return *this;
}

void SeqList::loadSequencesForAligner (const QuaffDPConfig& config, const QuaffAlignmentPrinter& aligner) {
  loadSequences (config);
  
  const set<string> dups = fastSeqDuplicateNames (seqs);
  if (!dups.empty()) {
    cerr << "Duplicate names:";
    for (const auto& d : dups)
      cerr << ' ' << d;
    cerr << endl;
    Fail ("All %s sequence names are required to be unique", type.c_str());
  }
}

void SeqList::loadSequences (const QuaffDPConfig& config) {
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
  argvec.pop_front();
  // leave "count" command undocumented... it's really for the paper
  //  briefText = "Usage: " + prog + " {help,train,count,align,overlap} [options]\n";
  briefText = "Usage: " + prog + " {help,train,align,overlap} [options]\n";
  
  text = briefText
    + "\n"
    + "Commands:\n"
    + "\n"
    + " " + prog + " train refs.fasta reads.fastq  >params.json\n"
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
    + "   -saveparams <file>, -savecounts <file>, -savecountswithprior <file>\n"
    + "                   Save parameters (or E-step counts) to file, not stdout\n"
    + "                    (saved counts can subsequently be used as a prior)\n"
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
    + "   -threshold <n>, -nothreshold\n"
    + "                   Log-odds ratio score threshold for alignment reporting\n"
    + "   -noquals        Ignore read quality scores during alignment\n"
    + "   -savealign <file>\n"
    + "                   Stream alignments to file, instead of stdout\n"
    + "   -format {fasta,stockholm,sam,refseq}\n"
    + "                   Alignment output format\n"
    + "\n"
    + "General options (for all commands, except where indicated):\n"
    + "   -verbose, -vv, -vvv, -v4, -v5, etc.\n"
    // uncomment to document debug logging:
    //    + "   -log <function_name>\n"
    + "                   Various levels of logging (-nocolor for monochrome)\n"
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
    + "   -kmatchoff      No kmer threshold, do full DP (typically very slow)\n"
    + "\n"
    + "Parallel procesing options:\n"
    + "   -threads <N>, -maxthreads\n"
    + "                   Use N threads, or use all cores available\n"
    + "   -remote [user@]host[:port[-maxport]]\n"
    + "                   Start a (multithreaded) remote quaff server via SSH\n"
    + "   -sshkey <file>  SSH private key file\n"
    + "   -sshpath <p>    Path to ssh\n"
    + "   -rsyncpath <p>  Path to rsync\n"
    + "   -remotepath <p> Path to remote binary (default " DefaultQuaffPath ")\n"
    + "   -rsync          Client will rsync data to server dir " BucketStagingDir "\n"
    + "   -s3bucket <B>   Client/server will sync data files to/from bucket B\n"
    + "   -ec2instances <N>\n"
    + "                   Launch N temporary EC2 instances as servers\n"
    + "   -ec2ami <AMI>, -ec2type <type>, -ec2cores <numberOfCores>,\n"
    + "   -ec2user <user>, -ec2key <keypair>, -ec2group <group>, -ec2port <port>\n"
    + "                   Control various aspects of the launched instances\n"
    + "                    (defaults: " AWS_DEFAULT_AMI ", " AWS_DEFAULT_INSTANCE_TYPE ", " TOSTRING(AWS_DEFAULT_INSTANCE_CORES) ",\n"
    + "                               " AWS_DEFAULT_USER ", " AWS_DEFAULT_KEY_PAIR ", " AWS_DEFAULT_SECURITY_GROUP ", " TOSTRING(DefaultServerPort) ")\n"
    + "\n"
    + "By default, quaff assumes all data files are in the same place on the server.\n"
    + "You can copy them across using -rsync, or -s3bucket, or other means (eg NFS).\n"
    + "\n"
    + "For AWS, ensure aws CLI tools are installed and credentials are set\n"
    + "(i.e. AWS_ACCESS_KEY_ID & AWS_SECRET_ACCESS_KEY environment variables).\n"
    + "You must use an AMI consistent with your AWS_DEFAULT_REGION\n"
    + "(the default AMI is a standard Amazon EC2 Linux for us-east-1).\n"
    + "A standard AMI should be fine: quaff downloads prereqs and builds itself.\n"
    + "Also ensure that security group <group> allows incoming connections\n"
    + "on ports 22 (ssh) and the range from <port> .. <port> + <numberOfCores> - 1.\n"
    + "Any problems can often be diagnosed by turning up the logging to -v5 or so."
    + "\n"
    + "Quaff makes every effort to clean up rogue EC2 instances, but please check!\n"
    + "\n";
}

string QuaffUsage::getCommand (const char* error) {
  if (argvec.empty()) {
    if (error)
      cerr << error << endl;
    else
      cerr << briefText;
    exit (EXIT_FAILURE);
  }
  const string command (argvec[0]);
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
