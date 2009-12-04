all: tex 

.SUFFIXES: .xml .t2t .tex 
.PHONY: clean xmlclean

SHELL=/bin/sh

tex: get cgoptions.tex mapping.tex

get: xmlclean cgoptions.xml mapping.xml

xmlclean:
	rm -f cgoptions.xml mapping.xml

t2tclean:
	rm -f cgoptions.t2t mapping.t2t

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
	./xmltotxt.sh $< > $@

.t2t.tex:
	./txt2tags -t tex --no-header --infile $< --outfile $@

clean: xmlclean t2tclean
	rm -f *.man *.tex *.wiki *~

