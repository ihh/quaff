Quaff is a pair HMM for aligning FASTQ reads to FASTA references,
with the following features:

- it pre-filters and bands the DP algorithms by looking for diagonals
  with lots of k-mer matches

- you can train the pair HMM on sequences using the Forward-Backward
  algorithm

- you can align reads to references using the Viterbi algorithm

- you can also align reads to reads (i.e. look for overlaps), using
  Viterbi

- it attempts to model the FASTQ quality scores (using a negative
  binomial distribution), and also uses k-mer context for modeling
  substitution and/or gap-opening probabilities

- it's threaded for speed

<pre><code>
Usage: quaff {help,train,align,overlap} [options]

Commands:

 quaff train refs.fasta reads.fastq  &gt;params.yaml
  (to fit a model to unaligned sequences, using EM/Forward-Backward)

   -maxiter &lt;n&gt;    Max number of EM iterations (default is 100)
   -mininc &lt;n&gt;     EM convergence threshold as relative log-likelihood increase
                    (default is .01)
   -force          Force each read to match a refseq, i.e. disallow null model
   -order &lt;k&gt;      Allow substitutions to depend on k-mer contexts
   -gaporder &lt;k&gt;   Allow gap open probabilities to depend on k-mer contexts
   -prior &lt;file&gt;, -saveprior &lt;file&gt;
                   Respectively: load/save prior pseudocounts from/to file
   -counts &lt;file&gt;  Save E-step counts to file, which can then be used as a prior
   -countswithprior &lt;file&gt;
                   Like -counts, but adds in prior pseudocounts as well


 quaff align refs.fasta reads.fastq
  (to align FASTQ reads to FASTA reference sequences, using Viterbi)

   -printall       Print all pairwise alignments, not just best for each read


 quaff overlap reads.fastq
  (to find overlaps between FASTQ reads, using Viterbi)


Alignment options (for align/overlap commands):
   -format {fasta,stockholm,refseq}
                   Alignment output format
   -threshold &lt;n&gt;
   -nothreshold    Log-odds ratio score threshold for alignment reporting
   -noquals        Ignore read quality scores during alignment

General options (for all commands, except where indicated):
   -verbose, -vv, -vvv, -v4, -v5, etc.
                   Various levels of logging (-nocolor for monochrome)
   -params &lt;file&gt;  Load model parameters from file
   -ref &lt;file&gt;     Load additional FASTA reference sequences
   -read &lt;file&gt;    Load additional FASTQ read sequences
   -fwdstrand      Do not include reverse-complemented sequences
   -global         Force all of refseq to be aligned (align/train only)
   -null &lt;file&gt;, -savenull &lt;file&gt;
                   Respectively: load/save null model from/to file

   -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
   -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
                    (default is 14 for overlap, 20 for align/train)
   -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
   -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
   -kmatchmax      Set kmer threshold to use all available memory (slow)
   -kmatchoff      No kmer threshold, do full DP (swapfile-heavy, v slow)

   -threads &lt;n&gt;, -maxthreads
                   Specify number of threads, or use all cores available

</code></pre>
