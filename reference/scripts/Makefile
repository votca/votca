SHELL=/bin/bash

#get the last two dirs from pwd
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )

.PHONY: clean
.SECONDARY:
.DELETE_ON_ERROR: %.t2t

all: Makefile.incl all.tex csg_table.tex 
wiki: all.wiki

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.SCRIPTS Makefile.incl)
include Makefile.incl

all.tex: $(SCRIPTSTEX)
	echo $(SCRIPTSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

all.wiki: $(SCRIPTSWIKI)
	cat $(SCRIPTSWIKI) > $@

csg_table.t2t: call2t2t.sh
	./call2t2t.sh > $@

csg_table.tex: csg_table.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@
	sed -i -e 's/{tabular}/{longtable}/' $@

%.t2t: script2t2t.sh
	./script2t2t.sh $* > $@	

%.tex: %.t2t
	./txt2tags -t tex --no-header --infile $< --outfile $@

%.wiki: %.t2t
	./txt2tags -t wiki --no-header --infile $< --outfile $@

clean: 
	rm -f $(SCRIPTST2T) csg_table.t2t
	rm -f $(SCRIPTSTEX) all.tex csg_table.tex
	rm -f $(SCRIPTSWIKI) all.wiki
	rm -f *~
	$(MAKE) -f Makefile.SCRIPTS clean

Makefile.incl:
	$(MAKE) -f Makefile.SCRIPTS Makefile.incl
