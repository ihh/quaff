# DART is at http://github.com/ihh/dart
# Assume it's installed in home directory
DARTDIR = $(HOME)/dart
PHYLOCOMPOSER = $(DARTDIR)/bin/phylocomposer

all: overlap.pdf overlap.sexpr overlap.dot

clean:
	rm overlap.pdf overlap.sexpr overlap.dot

overlap.sexpr overlap.dot: quaff.sexpr
	$(PHYLOCOMPOSER) $< --compdot overlap.dot | perl -pe 's/\(bit-value -0\) //g' >overlap.sexpr

%.pdf: %.dot
	dot -Tpdf $< >$@

.SECONDARY:
