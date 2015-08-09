.SECONDARY:

CPPFILES = $(wildcard src/*.cpp)

all: test

test: testfast testquaffio testlogsumexp

bin/%: $(CPPFILES) t/%.cpp
	test -e bin || mkdir bin
	g++ -std=c++11 -g -lstdc++ -lz -lgsl -o $@ t/$*.cpp $(CPPFILES)

testfast: bin/testfasta bin/testfastq
	perl/testexpect.pl bin/testfasta data/tiny.fasta data/tiny.fasta
	perl/testexpect.pl bin/testfasta data/tiny.fastq data/tiny.fasta
	perl/testexpect.pl bin/testfastq data/tiny.fastq data/tiny.fastq
	perl/testexpect.pl bin/testfastq data/tiny.fasta data/tiny.noqual.fastq
	perl/testexpect.pl bin/testfastq data/tiny.noqual.fastq data/tiny.noqual.fastq
	perl/testexpect.pl bin/testfastq data/tiny.truncated.fastq data/tiny.noqual.fastq

testquaffio: bin/testquaffio
	perl/testexpect.pl bin/testquaffio data/testquaffparams.yaml data/testquaffparams.yaml

testlogsumexp: bin/testlogsumexp
	bin/testlogsumexp -slow >data/logsumexp.txt
	perl/testexpect.pl bin/testlogsumexp -fast data/logsumexp.txt
