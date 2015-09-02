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

- it's highly parallel: multithreaded, can run remote servers,
  or launch its own temporary Amazon EC2 cluster

<pre><code>
Usage: quaff {help,train,align,overlap} [options]

Commands:

 quaff train refs.fasta reads.fastq  &gt;params.json
  (to fit a model to unaligned sequences, using EM/Forward-Backward)

   -maxiter &lt;n&gt;    Max number of EM iterations (default is 100)
   -mininc &lt;n&gt;     EM convergence threshold as relative log-likelihood increase
                    (default is .01)
   -force          Force each read to match a refseq, i.e. disallow null model
   -suborder &lt;k&gt;   Allow substitutions to depend on k-mer contexts
   -gaporder &lt;k&gt;   Allow gap open probabilities to depend on k-mer contexts
   -order &lt;k&gt;      Shorthand for '-suborder &lt;k&gt; -gaporder &lt;k&gt;'
   -prior &lt;file&gt;, -saveprior &lt;file&gt;
                   Respectively: load/save prior pseudocounts from/to file
   -saveparams &lt;file&gt;, -savecounts &lt;file&gt;, -savecountswithprior &lt;file&gt;
                   Save parameters (or E-step counts) to file, not stdout
                    (saved counts can subsequently be used as a prior)


 quaff align refs.fasta reads.fastq
  (to align FASTQ reads to FASTA reference sequences, using Viterbi)

   -printall       Print all pairwise alignments, not just best for each read


 quaff overlap reads.fastq
  (to find overlaps between FASTQ reads, using Viterbi)


Alignment options (for align/overlap commands):
   -threshold &lt;n&gt;, -nothreshold
                   Log-odds ratio score threshold for alignment reporting
   -noquals        Ignore read quality scores during alignment
   -savealign &lt;file&gt;
                   Stream alignments to file, instead of stdout
   -format {fasta,stockholm,sam,refseq}
                   Alignment output format

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
   -kmatchoff      No kmer threshold, do full DP (typically very slow)

Parallel procesing options:
   -threads &lt;N&gt;, -maxthreads
                   Use N threads, or use all cores available
   -remote [user@]host[:port[-maxport]]
                   Start a (multithreaded) remote quaff server via SSH
   -sshkey &lt;file&gt;  SSH private key file
   -sshpath &lt;p&gt;    Path to ssh
   -rsyncpath &lt;p&gt;  Path to rsync
   -remotepath &lt;p&gt; Path to remote binary (default /usr/local/bin/quaff)
   -rsync          Client will rsync data to server dir (/tmp/quaff)
   -s3bucket &lt;B&gt;   Client/server will sync data files to/from bucket B
   -ec2instances &lt;N&gt;
                   Launch N temporary EC2 instances as servers
   -ec2ami &lt;AMI&gt;, -ec2type &lt;type&gt;, -ec2cores &lt;numberOfCores&gt;,
   -ec2user &lt;user&gt;, -ec2key &lt;keypair&gt;, -ec2group &lt;group&gt;, -ec2port &lt;port&gt;
                   Control various aspects of the launched instances
                    (defaults: ami-1ecae776, m3.medium, 1,
                               ec2-user, quaff, quaff, 8000)

By default, quaff assumes all data files are in the same place on the server.
You can copy them across using -rsync, or -s3bucket, or other means (eg NFS).

For AWS, ensure aws CLI tools are installed and credentials are set
(i.e. AWS_ACCESS_KEY_ID & AWS_SECRET_ACCESS_KEY environment variables).
You must use an AMI consistent with your AWS_DEFAULT_REGION
(the default AMI is a standard Amazon EC2 Linux for us-east-1).
A standard AMI should be fine: quaff downloads prereqs and builds itself.
Also ensure that security group &lt;group&gt; allows incoming connections
on ports 22 (ssh) and the range from &lt;port&gt; .. &lt;port&gt; + &lt;numberOfCores&gt; - 1.
Any problems can often be diagnosed by turning up the logging to -v5 or so.
Quaff makes every effort to clean up rogue EC2 instances, but please check!

</code></pre>
