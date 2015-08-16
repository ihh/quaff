.SECONDARY:

GSL_PREFIX = /usr/local
GSLFLAGS = -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib

CPPFLAGS = -DUSE_VECTOR_GUARDS -std=c++11 -g
LIBFLAGS = -lstdc++ -lz -lgsl

CPPFILES = $(wildcard src/*.cpp)

all: quaff test

clean:
	rm -rf bin/*

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

# For updating README.md
README.md: bin/quaff
	PATH=bin:$(PATH); quaff help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){print};print"</code></pre>\n"' >temp.md
	mv temp.md $@
