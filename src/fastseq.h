#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <vector>
#include "../kseq/kseq.h"

using namespace std;

extern const string dnaAlphabet;
int tokenize (char c, const string& alphabet);

struct FastSeq {
  string name, comment, seq, qual;
  bool hasQual() const { return qual.size() == seq.size(); }
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
};

vector<FastSeq> readFastSeqs (const char* filename);
void writeFastaSeqs (ostream& out, const vector<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vector<FastSeq>& fastSeqs);

string revcomp (const string& dnaSeq);

#endif /* KSEQCONTAINER_INCLUDED */
