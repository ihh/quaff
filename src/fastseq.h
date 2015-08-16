#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <algorithm>
#include "vguard.h"
#include "../kseq/kseq.h"

using namespace std;

#define dnaAlphabetSize 4
extern const string dnaAlphabet;
unsigned int dnaComplement (unsigned int token);
char dnaComplementChar (char c);

typedef unsigned int SeqIdx;
typedef unsigned int AlphTok;
typedef int UnvalidatedAlphTok;
typedef unsigned long Kmer;
typedef unsigned int QualScore;

UnvalidatedAlphTok tokenize (char c, const string& alphabet);
Kmer makeKmer (SeqIdx k, vector<AlphTok>::const_iterator tok, AlphTok alphabetSize);
string kmerToString (Kmer kmer, SeqIdx k, const string& alphabet);


struct FastSeq {
  string name, comment, seq, qual;
  static const char minQualityChar = '!', maxQualityChar = '~';
  static const QualScore qualScoreRange = 94;
  SeqIdx length() const { return (SeqIdx) seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  static inline QualScore qualScoreForChar (char c) {
    return max (0, min ((int) qualScoreRange - 1, (int) (c - minQualityChar)));
  }
  static char charForQualScore (QualScore q) {
    return max (minQualityChar, min (maxQualityChar, (char) (q + minQualityChar)));
  }
  inline QualScore getQualScoreAt (SeqIdx pos) const { return qualScoreForChar (qual[pos]); }
  vguard<AlphTok> tokens (const string& alphabet) const;
  vguard<AlphTok> qualScores() const;
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
