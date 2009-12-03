all: tex 

.SUFFIXES: .xml .t2t .tex .man .wiki
.PHONY: clean xmlclean

SHELL=/bin/sh

wiki: get cgoptions.wiki mapping.wiki
man: get cgoptions.man mapping.man
tex: get cgoptions.tex mapping.tex

get: xmlclean cgoptions.xml mapping.xml

xmlclean:
	rm -f cgoptions.xml mapping.xml

%.xml:
	@if [ -z "${CSGSHARE}" ]; then \
	  echo CSGSHARE not defined >&2; \
	  exit 1; \
	fi
	@if [ -f "${CSGSHARE}/xml/$@" ]; then \
	  cp -p $(CSGSHARE)/xml/$@ .;\
	else\
	  echo Error, did not find ${CSGSHARE}/xml/$@ >&2;\
	  exit 1;\
	fi

.xml.t2t: 
	./xmltotxt.sh $< $@

.t2t.man:
	./txt2tags -t man --infile $< --outfile $@

.t2t.html:
	./txt2tags -t html --infile $< --outfile $@

.t2t.wiki:
	./txt2tags -t wiki --infile $< --outfile $@

.t2t.tex:
	./txt2tags -t tex --no-header --infile $< --outfile $@

clean: xmlclean
	rm -f *.man *.tex *.wiki *~

