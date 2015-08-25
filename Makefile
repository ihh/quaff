.SECONDARY:

GSLPREFIX ?= /usr/local
GSLFLAGS = -I$(GSLPREFIX)/include -L$(GSLPREFIX)/lib
GSLLIBS = -lgsl -lgslcblas

CPPFLAGS = -DUSE_VECTOR_GUARDS -std=c++11 -g
LIBFLAGS = -lstdc++ -lz $(GSLLIBS)

CPPFILES = $(wildcard src/*.cpp)

all: quaff

quaff: bin/quaff

install: bin/quaff
	cp $< /usr/local/bin
	chmod a+x /usr/local/bin/quaff

uninstall:
	rm /usr/local/bin/quaff

clean:
	rm -rf bin/*

# testaws is not included in the top-level 'make test' target
test: testfast testquaffio testlogsumexp testnegbinom testdiagenv

bin/%: $(CPPFILES) t/%.cpp
	test -e bin || mkdir bin
	clang++ $(CPPFLAGS) $(GSLFLAGS) $(LIBFLAGS) -o $@ t/$*.cpp $(CPPFILES)

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

testdiagenv: bin/testdiagenv
	bin/testdiagenv data/c8f30.fastq.gz data/c8f30.fastq.gz 6 14 64

# testaws should be run separately (it needs secret keys, could cost money, etc)
testaws: bin/testaws
	echo Run like this: 'bin/testaws keyPairName'

# For updating README.md
README.md: bin/quaff
	PATH=bin:$(PATH); quaff help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){print};print"</code></pre>\n"' >temp.md
	mv temp.md $@

# For updating default params
src/defaultparams.cpp: data/defaultparams.yaml
	perl -e 'open S,"<".shift();while(<S>){print;last if/defaultQuaffParamsText =/}close S;open A,"<".shift();while(<A>){chomp;print chr(34),$$_,"\\n",chr(34),"\n"}print";\n"' $@ $< >temp.cpp
	mv temp.cpp $@
