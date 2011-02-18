SHELL=/bin/bash

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.XMLS Makefile.incl)
include Makefile.incl

all: $(XMLSTEX) 
html: $(XMLSHTML)

.PHONY: clean upload
.SECONDARY:


%.t2t: 
	./xml2t2t.sh $* > $@

%.tex: %.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@

%.html: %.t2t
	./txt2tags -t html --no-header --infile $< --outfile $@

upload: 
	googlesites_upload.py  -d "/Documentation/reference-cgoptions"  -f cgoptions.xml.html
	googlesites_upload.py  -d "/Documentation/reference-cginteraction"  -f cginteraction.xml.html
	googlesites_upload.py  -d "/Documentation/reference-mapping"  -f mapping.xml.html

clean:
	rm -f $(XMLST2T)
	rm -f $(XMLSTEX)
	rm -f $(XMLSHTML)
	rm -f *~
	$(MAKE) -f Makefile.XMLS clean

Makefile.incl:
	$(MAKE) -f Makefile.XMLS Makefile.incl
