.SECONDARY:

GSL_PREFIX = /usr/local
GSLFLAGS = -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib

CPPFLAGS = -DUSE_VECTOR_GUARDS -std=c++11 -g
LIBFLAGS = -lstdc++ -lz -lgsl

CPPFILES = $(wildcard src/*.cpp)

all: test

test: testfast testquaffio testlogsumexp testnegbinom

bin/%: $(CPPFILES) t/%.cpp
	test -e bin || mkdir bin
	g++ $(CPPFLAGS) $(GSLFLAGS) $(LIBFLAGS) -o $@ t/$*.cpp $(CPPFILES)

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

testnegbinom: bin/testnegbinom
	bin/testnegbinom .1 5 10000 .1


quaff: bin/quaff
