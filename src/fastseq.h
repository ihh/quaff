#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <vector>
#include "../kseq/kseq.h"

using namespace std;

#define dnaAlphabetSize 4
extern const string dnaAlphabet;
int tokenize (char c, const string& alphabet);

struct FastSeq {
  string name, comment, seq, qual;
  static const char minQualityChar = '!', maxQualityChar = '~';
  static const int qualScoreRange = 94;
  int length() const { return seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  int getQualScoreAt (int pos) const { return (int) (qual[pos] - minQualityChar); }
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
};

vector<FastSeq> readFastSeqs (const char* filename);
void writeFastaSeqs (ostream& out, const vector<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vector<FastSeq>& fastSeqs);

string revcomp (const string& dnaSeq);

#endif /* KSEQCONTAINER_INCLUDED */
