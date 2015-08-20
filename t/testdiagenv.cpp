#include <iostream>
#include "../src/diagenv.h"
#include "../src/util.h"

#define TESTVECS(V,EXPECTED) testVecs(V,EXPECTED,#V,j)

// #define DEBUG_TESTDIAGENV

void testVecs (const vguard<SeqIdx>& v, const vguard<SeqIdx>& expected, const char* vname, const SeqIdx j) {
  if (!(v == expected)) {
    cerr << "Vector " << vname << " (j=" << j << ") is different from expected!" << endl;
    cerr << "  Actual:"; for (auto i : v) cerr << ' ' << i; cerr << endl;
    cerr << "Expected:"; for (auto i : expected) cerr << ' ' << i; cerr << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  }
}

void testEnvelope (const DiagonalEnvelope& env) {
#ifdef DEBUG_TESTDIAGENV
  cerr << "diagonals:"; int d = 0; for (auto i : env.diagonals) cerr << ' ' << d++ << ':' << i; cerr << endl;
  cerr << "xLen: " << env.xLen << endl;
#endif
  for (SeqIdx j = 1; j <= env.yLen; ++j) {
    vguard<SeqIdx> dumb_i;
    for (SeqIdx i = 1; i <= env.xLen; ++i)
      if (env.contains(i,j))
	dumb_i.push_back(i);

#ifdef DEBUG_TESTDIAGENV
    vguard<int>::const_iterator bi = env.beginIntersecting(j);
    vguard<int>::const_iterator ei = env.endIntersecting(j);
    cerr << "j=" << j << " bi=" << (bi - env.diagonals.begin()) << " ei=" << (ei - env.diagonals.begin()) << endl;
#endif
    
    const vguard<SeqIdx> fwd_i = env.forward_i(j);
    const vguard<SeqIdx> iter_i (env.begin(j), env.end(j));

    const vguard<SeqIdx> dumb_rev_i (dumb_i.rbegin(), dumb_i.rend());
    const vguard<SeqIdx> rev_i = env.reverse_i(j);
    const vguard<SeqIdx> riter_i (env.rbegin(j), env.rend(j));

    TESTVECS (fwd_i, dumb_i);
    TESTVECS (iter_i, dumb_i);
    TESTVECS (rev_i, dumb_rev_i);
    TESTVECS (riter_i, dumb_rev_i);
  }
}

int main (int argc, char **argv) {
  if (argc != 6) {
    cout << "Usage: " << argv[0] << " <fastqfile1> <fastqfile2> <kmatch> <kmatchn> <bandsize>\n";
    exit (EXIT_FAILURE);
  }

  const vguard<FastSeq> fs1 = readFastSeqs (argv[1]);
  const vguard<FastSeq> fs2 = readFastSeqs (argv[2]);

  const unsigned int K = atoi(argv[3]);
  const unsigned int N = atoi(argv[4]);
  const unsigned int B = atoi(argv[5]);
  
  Assert (fs1.size() == 1 && fs2.size() == 1, "Each Fastq file must have exactly 1 sequence");

  DiagonalEnvelope env (fs1[0], fs2[0]);
  env.initSparse (K, B, N);
  testEnvelope (env);

  cout << "ok" << endl;
  
  exit (EXIT_SUCCESS);
}
