SHELL=/bin/bash
TRUNC=$(shell hg root)
TXT2TAGS=$(TRUNC)/scripts/txt2tags

all: 
#this is a hack to generate Makefile.incl
HACK:=$(shell $(MAKE) -f Makefile.XMLS Makefile.incl)
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
	googlesites_upload.py  -d "/Documentation/coarse-graining/reference-cgoptions"  -f cgoptions.xml.html
	googlesites_upload.py  -d "/Documentation/coarse-graining/reference-cginteraction"  -f cginteraction.xml.html
	googlesites_upload.py  -d "/Documentation/coarse-graining/reference-mapping"  -f mapping.xml.html

clean:
	rm -f $(XMLST2T) config.t2t
	rm -f $(XMLSTEX)
	rm -f $(XMLSHTML)
	rm -f *~
	$(MAKE) -f Makefile.XMLS clean

Makefile.incl:
	$(MAKE) -f Makefile.XMLS Makefile.incl
