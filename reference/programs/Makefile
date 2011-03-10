SHELL=/bin/bash
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )
TRUNC=$(shell hg root)
TXT2TAGS=$(TRUNC)/scripts/txt2tags

.PHONY: clean bootstrap

all: all.tex

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.BINS Makefile.incl)
include Makefile.incl

all.tex: $(BINSTEX)
	echo $(BINSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

%.t2t: help2t2t.sh
	[ $@ = "config.t2t" ] && cp ../$@ . || ./help2t2t.sh $* > $@	

%.tex: %.t2t config.t2t
	$(TXT2TAGS) -t tex --no-header --infile $< --outfile $@

clean: 
	rm -f $(BINST2T)
	rm -f $(BINSTEX) all.tex
	rm -f *~
	$(MAKE) -f Makefile.BINS clean

Makefile.incl:
	$(MAKE) -f Makefile.BINS Makefile.incl
