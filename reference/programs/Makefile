SHELL=/bin/bash
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )

.PHONY: clean bootstrap

all: all.tex
wiki: all.wiki
html: all.html

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.BINS Makefile.incl)
include Makefile.incl

all.tex: $(BINSTEX)
	echo $(BINSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

all.wiki: $(BINSWIKI)
	cat $(BINSWIKI) > $@

all.html: $(BINSHTML)
	cat $(BINSHTML) > $@

%.t2t: 
	./help2t2t.sh $* > $@	

%.tex: %.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@

%.wiki: %.t2t
	./txt2tags -t wiki --no-header --infile $< --outfile $@

%.html: %.t2t
	./txt2tags -t html --no-header --infile $< --outfile $@

upload: html
	googlesites_upload.py  -d "/Documentation/reference-programs"  -f all.html


clean: 
	rm -f $(BINST2T)
	rm -f $(BINSTEX) all.tex
	rm -f $(BINSWIKI) all.wiki
	rm -f $(BINSHTML) all.html
	rm -f *~
	$(MAKE) -f Makefile.BINS clean

Makefile.incl:
	$(MAKE) -f Makefile.BINS Makefile.incl
