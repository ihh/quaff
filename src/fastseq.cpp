#include <zlib.h>
#include <iostream>
#include "fastseq.h"

KSEQ_INIT(gzFile, gzread)

const string dnaAlphabet ("ACGT");

int tokenize (char c, const string& alphabet) {
  int tok;
  const char* alphStr = alphabet.c_str(); 
  tok = (int) (strchr (alphStr, toupper(c)) - alphStr);
  return tok >= (int) strlen(alphStr) ? -1 : tok;
}

vector<int> FastSeq::tokens (const string& alphabet) const {
  vector<int> tok;
  tok.reserve (length());
  for (auto c : seq)
    tok.push_back (tokenize (c, alphabet));
  return tok;
}

vector<int> FastSeq::qualScores() const {
  vector<int> q;
  if (hasQual()) {
    q.reserve (length());
    for (auto c : qual)
      q.push_back (qualScoreForChar (c));
  }
  return q;
}

void FastSeq::writeFasta (ostream& out) const {
  out << '>' << name;
  if (comment.size())
    out << ' ' << comment;
  out << endl;
  out << seq << endl;
}

void FastSeq::writeFastq (ostream& out) const {
  out << '@' << name;
  if (comment.size())
    out << ' ' << comment;
  out << endl;
  out << seq << endl;
  if (hasQual())
    out << '+' << endl << qual << endl;
}

void writeFastaSeqs (ostream& out, const vector<FastSeq>& fastSeqs) {
  for (auto s : fastSeqs)
    s.writeFasta (out);
}

void writeFastqSeqs (ostream& out, const vector<FastSeq>& fastSeqs) {
  for (auto s : fastSeqs)
    s.writeFastq (out);
}

vector<FastSeq> readFastSeqs (const char* filename) {
  vector<FastSeq> seqs;

  gzFile fp = gzopen(filename, "r");
  if (fp != Z_NULL) {
    kseq_t *ks = kseq_init(fp);
    while (kseq_read(ks) != -1) {
      FastSeq seq;
      seq.name = string(ks->name.s);
      seq.seq = string(ks->seq.s);
      if (ks->comment.l)
	seq.comment = string(ks->comment.s);
      if (ks->qual.l == ks->seq.l)
	seq.qual = string(ks->qual.s);
      seqs.push_back (seq);
    }
    kseq_destroy (ks);
    gzclose (fp);
  }

  return seqs;
}

string revcomp (const string& dnaSeq) {
  string rev = dnaSeq;
  const int len = dnaSeq.size();
  for (int i = 0; i < len; ++i) {
    const int tok = tokenize (dnaSeq[i], dnaAlphabet);
    rev[len - 1 - i] = tok < 0 ? 'N' : dnaAlphabet[3 - tok];
  }
  return rev;
}
