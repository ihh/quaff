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
  substitution probabilities, both of which are theoretically
  informative for nanopore reads


<pre>
Usage: quaff {help,train,align,overlap} [options]

Commands:

 quaff train -ref refs.fasta -read reads.fastq  >params.yaml
  (to fit a model to unaligned sequences, using EM)

   -params <file>  Optional initial parameters
   -maxiter <n>    Max number of EM iterations
   -mininc <n>     EM convergence threshold (relative log-likelihood increase)
   -force          Force each read to match a refseq, i.e. disallow null model
   -prior <file>, -saveprior <file>
                   Respectively: load/save prior pseudocounts from/to file
   -counts <file>  Save E-step counts to file, which can then be used as a prior
   -countswithprior <file>
                   Like -counts, but adds in prior pseudocounts as well


 quaff align -params params.yaml -ref refs.fasta -read reads.fastq
  (to align FASTQ reads to FASTA reference sequences, using Viterbi)

   -printall       Print all pairwise alignments, not just best for each read


 quaff overlap -params params.yaml -read reads.fastq
  (to find overlaps between FASTQ reads, using Viterbi)


Alignment options (align/overlap commands):
   -format {fasta,stockholm,refseq}
                   Alignment output format
   -threshold <n>
   -nothreshold    Log-odds ratio score threshold for alignment reporting
   -noquals        Ignore read quality scores during alignment

General options (all commands, except where indicated):
   -verbose, -vv, -vvv, -v4, etc.
                   Various levels of logging
   -fwdstrand      Do not include reverse-complemented sequences
   -global         Force all of refseq to be aligned (align/train only)
   -kmatch         Length of kmers for pre-filtering heuristic (default 6)
   -kmatchn <n>    Threshold# of kmer matches to include a diagonal (default 14)
   -kmatchband <n> Size of DP band around kmer-matching diagonals (default 64)
   -dense          Do full DP, not just kmer-matching diagonals (memory hog!)
   -null <file>, -savenull <file>
                   Respectively: load/save null model from/to file
</pre>
