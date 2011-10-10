SHELL=/bin/bash
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )
TRUNC=$(shell hg root)
TXT2TAGS=$(TRUNC)/scripts/txt2tags

.PHONY: clean bootstrap

all: all.tex
html: all.html
t2t: $(BINST2T)

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.BINS Makefile.incl)
include Makefile.incl

all.tex: $(BINSTEX)
	echo $(BINSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

$(BINST2T): help2t2t.sh
%.t2t: 
	./help2t2t.sh $* > $@	

config.t2t: ../config.t2t
	cp $< .

all.html: $(BINSHTML)
	cat $(BINSHTML) > $@

%.tex: %.t2t config.t2t
	$(TXT2TAGS) -t tex --no-header --infile $< --outfile $@

%.html: %.t2t
	$(TXT2TAGS) -t html --no-header --infile $< --outfile $@

upload: html
	googlesites_upload.py  -d "/Documentation/coarse-graining/reference-programs"  -f all.html


clean: 
	rm -f $(BINST2T) config.t2t
	rm -f $(BINSTEX) all.tex
	rm -f $(BINSHTML) all.html
	rm -f *~
	$(MAKE) -f Makefile.BINS clean

Makefile.incl:
	$(MAKE) -f Makefile.BINS Makefile.incl
