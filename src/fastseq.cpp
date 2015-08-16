#include <zlib.h>
#include <iostream>
#include "fastseq.h"
#include "util.h"

KSEQ_INIT(gzFile, gzread)

const string dnaAlphabet ("ACGT");

UnvalidatedAlphTok tokenize (char c, const string& alphabet) {
  UnvalidatedAlphTok tok;
  const char* alphStr = alphabet.c_str(); 
  tok = (UnvalidatedAlphTok) (strchr (alphStr, toupper(c)) - alphStr);
  return tok >= (int) strlen(alphStr) ? -1 : tok;
}

AlphTok dnaComplement (AlphTok token) {
  return dnaAlphabetSize - 1 - token;
}

char dnaComplementChar (char c) {
  const int tok = tokenize (c, dnaAlphabet);
  return tok < 0 ? 'N' : dnaAlphabet[dnaComplement(tok)];
}

Kmer makeKmer (SeqIdx k, vector<unsigned int>::const_iterator tok, unsigned int alphabetSize) {
  Kmer kmer = 0, mul = 1;
  for (SeqIdx j = 0; j < k; ++j) {
    const unsigned int token = tok[k - j - 1];
    kmer += mul * token;
    mul *= alphabetSize;
  }
  return kmer;
}

string kmerToString (Kmer kmer, SeqIdx k, const string& alphabet) {
  string rev;
  for (SeqIdx j = 0; j < k; ++j, kmer = kmer / alphabet.size())
    rev += alphabet[kmer % alphabet.size()];
  return string (rev.rbegin(), rev.rend());
}

vguard<AlphTok> FastSeq::tokens (const string& alphabet) const {
  vguard<AlphTok> tok;
  tok.reserve (length());
  for (const auto& c : seq) {
    const int t = tokenize (c, alphabet);
    if (t < 0) {
      cerr << "Unknown symbol " << c << " in sequence " << name << endl;
      throw;
    }
    tok.push_back (t);
  }
  return tok;
}

vguard<QualScore> FastSeq::qualScores() const {
  vguard<QualScore> q;
  if (hasQual()) {
    q.reserve (length());
    for (const auto& c : qual)
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

void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs) {
  for (const auto& s : fastSeqs)
    s.writeFasta (out);
}

void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs) {
  for (const auto& s : fastSeqs)
    s.writeFastq (out);
}

vguard<FastSeq> readFastSeqs (const char* filename) {
  vguard<FastSeq> seqs;

  gzFile fp = gzopen(filename, "r");
  Assert (fp != Z_NULL, "Couldn't open %s", filename);

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

  return seqs;
}

string revcomp (const string& dnaSeq) {
  string rev = dnaSeq;
  const size_t len = dnaSeq.size();
  for (size_t i = 0; i < len; ++i)
    rev[len - 1 - i] = dnaComplementChar (dnaSeq[i]);
  return rev;
}

FastSeq FastSeq::revcomp() const {
  FastSeq fs;
  fs.name = "revcomp(" + name + ")";
  fs.comment = comment;
  fs.seq = ::revcomp(seq);
  fs.qual = string (qual.rbegin(), qual.rend());
  return fs;
}

void addRevcomps (vguard<FastSeq>& db) {
  vguard<FastSeq> revcomps;
  revcomps.reserve (db.size());
  for (const auto& fs : db)
    revcomps.push_back (fs.revcomp());
  db.insert (db.end(), revcomps.begin(), revcomps.end());
}
