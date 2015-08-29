.SECONDARY:

# try to figure out where GSL is
# autoconf would be better but we just need a quick hack for now :)
GSLPREFIX ?= /usr
ifeq (,$(wildcard $(GSLPREFIX)/include/gsl/gsl_sf.h))
GSLPREFIX := /usr/local
endif

GSLFLAGS = -I$(GSLPREFIX)/include -L$(GSLPREFIX)/lib
GSLLIBS = -lgsl -lgslcblas

# figure out whether to use Boost
# Boost is optional -- it's only needed for regexes with gcc
BOOSTPREFIX = /usr
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX := /usr/local
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX :=
endif
endif

BOOSTFLAGS =
BOOSTLIBS =
ifneq (,$(BOOSTPREFIX))
BOOSTFLAGS := -DUSE_BOOST -I$(BOOSTPREFIX)/include -L$(BOOSTPREFIX)/lib
BOOSTLIBS := -lboost_regex
endif

# install dir
PREFIX ?= /usr/local

# other flags
CPPFLAGS = -DUSE_VECTOR_GUARDS -std=c++11 -g $(GSLFLAGS) $(BOOSTFLAGS)
LIBFLAGS = -lstdc++ -lz $(GSLLIBS) $(BOOSTLIBS)

CPPFILES = $(wildcard src/*.cpp)

# change this to g++ if using gcc
CPP = clang++


all: quaff

quaff: bin/quaff

install: bin/quaff
	cp $< $(PREFIX)/bin
	chmod a+x $(PREFIX)/bin/quaff

uninstall:
	rm $(PREFIX)/bin/quaff

clean:
	rm -rf bin/*

# AWS
# NB pseudotarget runs make again, to re-test $(BOOST)
aws-install: aws-dep
	make PREFIX=$(PREFIX) install

aws-dep:
	yum -y install gcc clang gsl-devel zlib-devel boost-devel

# OS X
# NB pseudotarget runs make again, to re-test $(BOOST)
osx-install: osx-dep
	make PREFIX=$(PREFIX) install

osx-dep:
	brew install gsl boost awscli

# Tests
# testaws is not included in the top-level 'make test' target
test: testfast testquaffjsonio testquaffcountsjsonio testquaffnulljsonio testlogsumexp testnegbinom testdiagenv testregex

bin/%: $(CPPFILES) t/%.cpp
	test -e bin || mkdir bin
	$(CPP) $(CPPFLAGS) $(LIBFLAGS) -o $@ t/$*.cpp $(CPPFILES)

testfast: bin/testfasta bin/testfastq
	perl/testexpect.pl bin/testfasta data/tiny.fasta data/tiny.fasta
	perl/testexpect.pl bin/testfasta data/tiny.fastq data/tiny.fasta
	perl/testexpect.pl bin/testfastq data/tiny.fastq data/tiny.fastq
	perl/testexpect.pl bin/testfastq data/tiny.fasta data/tiny.noqual.fastq
	perl/testexpect.pl bin/testfastq data/tiny.noqual.fastq data/tiny.noqual.fastq
	perl/testexpect.pl bin/testfastq data/tiny.truncated.fastq data/tiny.noqual.fastq

testquaffjsonio: bin/testquaffjsonio
	perl/testexpect.pl bin/testquaffjsonio data/testquaffparams.json data/testquaffparams.json
	perl/testexpect.pl bin/testquaffjsonio data/defaultparams.json data/defaultparams.json

testquaffnulljsonio: bin/testquaffnulljsonio
	perl/testexpect.pl bin/testquaffnulljsonio data/testquaffnullparams.json data/testquaffnullparams.json

testquaffcountsjsonio: bin/testquaffcountsjsonio
	perl/testexpect.pl bin/testquaffcountsjsonio data/testquaffcounts.json data/testquaffcounts.json

testlogsumexp: bin/testlogsumexp
	bin/testlogsumexp -slow >data/logsumexp.txt
	perl/testexpect.pl bin/testlogsumexp -fast data/logsumexp.txt

testnegbinom: bin/testnegbinom
	bin/testnegbinom .1 5 10000 .1

testdiagenv: bin/testdiagenv
	bin/testdiagenv data/c8f30.fastq.gz data/c8f30.fastq.gz 6 14 64

testregex: t/testregex.cpp src/regexmacros.h
	$(CPP) $(CPPFLAGS) $(LIBFLAGS) -o bin/$@ $<
	bin/$@

# testaws should be run separately (it needs secret keys, could cost money, etc)
testaws: bin/testaws
	echo Run like this: 'bin/testaws keyPairName'

# For updating README.md
README.md: bin/quaff
	PATH=bin:$(PATH); quaff help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){print};print"</code></pre>\n"' >temp.md
	mv temp.md $@

# For updating default params
src/defaultparams.cpp: data/defaultparams.json
	perl -e 'open S,"<".shift();while(<S>){print;last if/defaultQuaffParamsText =/}close S;open A,"<".shift();$$q=chr(34);while(<A>){chomp;s/$$q/\\$$q/g;print chr(34),$$_,"\\n",chr(34),"\n"}print";\n"' $@ $< >temp.cpp
	mv temp.cpp $@
