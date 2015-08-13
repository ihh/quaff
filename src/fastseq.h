#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <algorithm>
#include "vguard.h"
#include "../kseq/kseq.h"

using namespace std;

#define dnaAlphabetSize 4
extern const string dnaAlphabet;

typedef unsigned int SeqIdx;

int tokenize (char c, const string& alphabet);
unsigned long makeKmer (SeqIdx k, vector<unsigned int>::const_iterator tok, unsigned int alphabetSize);
string kmerToString (unsigned long kmer, SeqIdx k, const string& alphabet);

struct FastSeq {
  string name, comment, seq, qual;
  static const char minQualityChar = '!', maxQualityChar = '~';
  static const unsigned int qualScoreRange = 94;
  SeqIdx length() const { return (SeqIdx) seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  static inline unsigned int qualScoreForChar (char c) {
    return max (0, min ((int) qualScoreRange - 1, (int) (c - minQualityChar)));
  }
  static char charForQualScore (int q) {
    return max (minQualityChar, min (maxQualityChar, (char) (q + minQualityChar)));
  }
  inline unsigned int getQualScoreAt (SeqIdx pos) const { return qualScoreForChar (qual[pos]); }
  vguard<unsigned int> tokens (const string& alphabet) const;
  vguard<unsigned int> qualScores() const;
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
