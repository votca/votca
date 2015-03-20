SHELL=/bin/bash
TRUNC=../..
TXT2TAGS=$(TRUNC)/scripts/txt2tags

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.XMLS Makefile.incl)
include Makefile.incl

all: $(XMLSTEX) 
html: $(XMLSHTML)

.PHONY: clean upload
.SECONDARY:

$(XMLST2T): xml2t2t.sh
%.t2t: 
	./xml2t2t.sh $* > $@

config.t2t: ../config.t2t
	cp $< .

%.tex: %.t2t config.t2t
	$(TXT2TAGS) -t tex --no-header --infile $< --outfile $@

%.html: %.t2t config.t2t
	$(TXT2TAGS) -t html --no-header --infile $< --outfile $@

upload: html
	googlesites_upload.py  -d "/Documentation/reference-cgoptions"  -f cgoptions.xml.html
	googlesites_upload.py  -d "/Documentation/reference-cginteraction"  -f cginteraction.xml.html
	googlesites_upload.py  -d "/Documentation/reference-mapping"  -f mapping.xml.html

clean:
	rm -f $(XMLST2T) config.t2t
	rm -f $(XMLSTEX)
	rm -f $(XMLSHTML)
	rm -f *~
	$(MAKE) -f Makefile.XMLS clean

Makefile.incl:
	$(MAKE) -f Makefile.XMLS Makefile.incl
