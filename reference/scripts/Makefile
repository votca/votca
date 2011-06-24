SHELL=/bin/bash

#get the last two dirs from pwd
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )

.PHONY: clean
.SECONDARY:
.DELETE_ON_ERROR: %.t2t

all: Makefile.incl all.tex csg_table.tex 

html: all.html

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.SCRIPTS Makefile.incl)
include Makefile.incl

all.tex: $(SCRIPTSTEX)
	echo $(SCRIPTSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

all.html: $(SCRIPTSHTML) csg_table.html
	cat csg_table.html $(SCRIPTSHTML) > $@

csg_table.t2t: call2t2t.sh
	./call2t2t.sh > $@

csg_table.tex: csg_table.t2t config.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@
	sed -i -e 's/{tabular}/{longtable}/' $@

%.t2t: script2t2t.sh
	[ $@ = "config.t2t" ] && cp ../$@ . || ./script2t2t.sh $* > $@	

%.tex: %.t2t config.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@

%.html: %.t2t config.t2t
	./txt2tags -t html --no-header --infile $< --outfile $@

upload: html
	googlesites_upload.py  -d "/Documentation/reference-scripts"  -f all.html


clean: 
	rm -f $(SCRIPTST2T) csg_table.t2t config.t2t
	rm -f $(SCRIPTSTEX) all.tex csg_table.tex
	rm -f $(SCRIPTSHTML) all.html
	rm -f *~
	$(MAKE) -f Makefile.SCRIPTS clean

Makefile.incl:
	$(MAKE) -f Makefile.SCRIPTS Makefile.incl
