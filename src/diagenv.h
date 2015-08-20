#ifndef DIAG_ENV_INCLUDED
#define DIAG_ENV_INCLUDED

#include <map>
#include "fastseq.h"

using namespace std;

// Somewhat arbitrary numbers:
// Default k-mer length
#define DEFAULT_KMER_LENGTH 6

// The number of kmer matches above which a diagonal is selected for a band
#define DEFAULT_KMER_THRESHOLD 14

// Default band size
#define DEFAULT_BAND_SIZE 64

struct DiagonalEnvelope {
  const FastSeq *px, *py;
  const SeqIdx xLen, yLen;
  vguard<int> diagonals;   // Sorted ascending. (i,j) is on diagonal d if i-j=d
  DiagonalEnvelope (const FastSeq& x, const FastSeq& y)
    : px(&x), py(&y), xLen(px->length()), yLen(py->length()) { }
  void initFull();
  void initSparse (unsigned int kmerLen = DEFAULT_KMER_LENGTH,
		   unsigned int bandSize = DEFAULT_BAND_SIZE,
		   int kmerThreshold = DEFAULT_KMER_THRESHOLD,  // negative => use memory guides
		   size_t cellSize = sizeof(double),
		   size_t maxSize = 0);
  inline int minDiagonal() const { return 1 - (int) yLen; }
  inline int maxDiagonal() const { return xLen - 1; }
  inline bool contains (SeqIdx i, SeqIdx j) const {
    const int diag = i - j;
    const vguard<int>::const_iterator iter = lower_bound (diagonals.begin(), diagonals.end(), diag);
    return iter != diagonals.end() && *iter == diag;
  }
  inline bool intersects (SeqIdx j, int diag) const {
    const int i = diag + j;
    return i > 0 && i <= xLen;
  }
  static inline SeqIdx get_i (SeqIdx j, int diag) { return (SeqIdx) (diag + j); }
  static inline int get_diag (SeqIdx i, SeqIdx j) { return (int) i - (int) j; }
  vguard<SeqIdx> forward_i (SeqIdx j) const;
  vguard<SeqIdx> reverse_i (SeqIdx j) const;

  inline vguard<int>::const_iterator beginIntersecting (SeqIdx j) const {
    // i = diag + j
    // want i > 0 && i <= xLen
    // so diag > -j && diag <= xLen - j
    return upper_bound (diagonals.begin(), diagonals.end(), -(int)j);
  }

  inline vguard<int>::const_iterator endIntersecting (SeqIdx j) const {
    // see comment in beginIntersecting()
    return upper_bound (diagonals.begin(), diagonals.end(), xLen - (int)j);
  }

  class DiagonalEnvelopeIterator : public iterator<forward_iterator_tag,SeqIdx> {
  private:
    const SeqIdx j;
    vguard<int>::const_iterator iter;
    const vguard<int>::const_iterator iterEnd;
  public:
    inline DiagonalEnvelopeIterator (SeqIdx j, const vguard<int>::const_iterator& iter, const vguard<int>::const_iterator& iterEnd)
      : j(j), iter(iter), iterEnd(iterEnd)
    { }
    inline bool operator== (const DiagonalEnvelopeIterator& dei) const { return iter == dei.iter; }
    inline bool operator!= (const DiagonalEnvelopeIterator& dei) const { return iter != dei.iter; }
    inline SeqIdx operator*() const { return (SeqIdx) (*iter + j); }
    inline DiagonalEnvelopeIterator& operator++() { ++iter; return *this; }
    inline DiagonalEnvelopeIterator operator++(int) { DiagonalEnvelopeIterator dei = *this; ++*this; return dei; }
    inline bool finished() const { return iter == iterEnd; }
  };

  class DiagonalEnvelopeReverseIterator : public iterator<forward_iterator_tag,SeqIdx> {
  private:
    const SeqIdx j;
    vguard<int>::const_iterator iter;
    const vguard<int>::const_iterator iterEnd;
  public:
    inline DiagonalEnvelopeReverseIterator (SeqIdx j, const vguard<int>::const_iterator& iter, const vguard<int>::const_iterator& iterEnd)
      : j(j), iter(iter), iterEnd(iterEnd)
    { }
    inline bool operator== (const DiagonalEnvelopeReverseIterator& dei) const { return iter == dei.iter; }
    inline bool operator!= (const DiagonalEnvelopeReverseIterator& dei) const { return iter != dei.iter; }
    inline SeqIdx operator*() const { return (SeqIdx) (*(iter - 1) + j); }
    inline DiagonalEnvelopeReverseIterator& operator++() { --iter; return *this; }
    inline DiagonalEnvelopeReverseIterator operator++(int) { DiagonalEnvelopeReverseIterator deri = *this; ++*this; return deri; }
    inline bool finished() const { return iter == iterEnd; }
  };

  typedef DiagonalEnvelopeIterator iterator;
  typedef const DiagonalEnvelopeIterator const_iterator;
  typedef DiagonalEnvelopeReverseIterator reverse_iterator;
  typedef const DiagonalEnvelopeReverseIterator const_reverse_iterator;

  inline iterator begin (SeqIdx j) const { return DiagonalEnvelopeIterator (j, beginIntersecting(j), endIntersecting(j)); }
  inline const_iterator end (SeqIdx j) const { return DiagonalEnvelopeIterator (j, endIntersecting(j), endIntersecting(j)); }

  inline reverse_iterator rbegin (SeqIdx j) const { return DiagonalEnvelopeReverseIterator (j, endIntersecting(j), beginIntersecting(j)); }
  inline const_reverse_iterator rend (SeqIdx j) const { return DiagonalEnvelopeReverseIterator (j, beginIntersecting(j), beginIntersecting(j)); }
};

#endif /* DIAG_ENV_INCLUDED */
