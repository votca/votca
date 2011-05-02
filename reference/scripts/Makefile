SHELL=/bin/bash

#get the last two dirs from pwd
DIR:=$(shell pwd | sed 's@^.*/\([^/]*/[^/]*\)$$@\1@' )
TRUNC=$(shell hg root)
TXT2TAGS=$(TRUNC)/scripts/txt2tags

.PHONY: clean
.SECONDARY:
.DELETE_ON_ERROR: %.t2t

all: Makefile.incl all.tex csg_table.tex 

#this is a hack to generate Makefile.incl
HAKE:=$(shell $(MAKE) -f Makefile.SCRIPTS Makefile.incl)
include Makefile.incl

all.tex: $(SCRIPTSTEX)
	echo $(SCRIPTSTEX) | \
	sed -e 's/\.tex[[:space:]]*/\n/g' | \
	sed -e '$$d' -e 's@.*@\\input{$(DIR)/&}@' > $@

csg_table.t2t: call2t2t.sh
	./call2t2t.sh > $@

csg_table.tex: csg_table.t2t config.t2t
	$(TXT2TAGS) -t tex --no-header --infile $< --outfile $@
	sed -i -e 's/{tabular}/{longtable}/' $@

$(SCRIPTST2T): script2t2t.sh
%.t2t:
	./script2t2t.sh $* > $@

config.t2t: ../config.t2t
	cp $< .

%.tex: %.t2t config.t2t
	$(TXT2TAGS) -t tex --no-header --infile $< --outfile $@

clean: 
	rm -f $(SCRIPTST2T) csg_table.t2t config.t2t
	rm -f $(SCRIPTSTEX) all.tex csg_table.tex
	rm -f *~
	$(MAKE) -f Makefile.SCRIPTS clean

Makefile.incl:
	$(MAKE) -f Makefile.SCRIPTS Makefile.incl
