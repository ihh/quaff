.SECONDARY:

CPPFILES = $(wildcard src/*.cpp)

all: testfast

bin:
	test -e $@ || mkdir $@

bin/%: $(CPPFILES) t/%.cpp bin
	g++ -std=c++11 -g -lstdc++ -lz -o $@ t/$*.cpp $(CPPFILES)


testfast: bin/testfasta bin/testfastq
	perl/testexpect.pl bin/testfasta data/tiny.fasta data/tiny.fasta
	perl/testexpect.pl bin/testfasta data/tiny.fastq data/tiny.fasta
	perl/testexpect.pl bin/testfastq data/tiny.fastq data/tiny.fastq
	perl/testexpect.pl bin/testfastq data/tiny.fasta data/tiny.noqual.fastq
