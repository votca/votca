SHELL=/bin/bash

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.XMLS Makefile.incl)
include Makefile.incl

all: $(XMLSTEX) 
html: $(XMLSHTML)

.PHONY: clean
.SECONDARY:


%.t2t: 
	./xml2t2t.sh $* > $@

%.tex: %.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@

%.html: %.t2t
	./txt2tags -t html --no-header --infile $< --outfile $@

clean:
	rm -f $(XMLST2T)
	rm -f $(XMLSTEX)
	rm -f $(XMLSHTML)
	rm -f *~
	$(MAKE) -f Makefile.XMLS clean

Makefile.incl:
	$(MAKE) -f Makefile.XMLS Makefile.incl
