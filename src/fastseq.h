#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <algorithm>
#include "vguard.h"
#include "../kseq/kseq.h"

using namespace std;

#define dnaAlphabetSize 4
extern const string dnaAlphabet;
int tokenize (char c, const string& alphabet);

struct FastSeq {
  string name, comment, seq, qual;
  static const char minQualityChar = '!', maxQualityChar = '~';
  static const int qualScoreRange = 94;
  size_t length() const { return seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  static inline int qualScoreForChar (char c) {
    return max (0, min (qualScoreRange - 1, (int) (c - minQualityChar)));
  }
  static char charForQualScore (int q) {
    return max (minQualityChar, min (maxQualityChar, (char) (q + minQualityChar)));
  }
  inline int getQualScoreAt (size_t pos) const { return qualScoreForChar (qual[pos]); }
  vguard<int> tokens (const string& alphabet) const;
  vguard<int> qualScores() const;
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
  FastSeq revcomp() const;
};

vguard<FastSeq> readFastSeqs (const char* filename);
void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);

string revcomp (const string& dnaSeq);
void addRevcomps (vguard<FastSeq>& db);

#endif /* KSEQCONTAINER_INCLUDED */
