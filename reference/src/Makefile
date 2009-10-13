all: wiki tex man 

.SUFFIXES: .xml .t2t .tex .man .wiki .html .sh
.PHONY: clean cgoptions.xml

wiki: cgoptions.wiki mapping.wiki
man: cgoptions.man mapping.man
tex: cgoptions.tex mapping.tex
html: cgoptions.html mapping.html

cgoptions.xml:
	if [ -z "${CSGSHARE}" ]; then \
	  echo CSGSHARE not defined >&2; \
	  exit 1; \
	fi
	if [ -f "${CSGSHARE}/doc/cgoptions.xml" ]; then \
	  cp -p $(CSGSHARE)/doc/cgoptions.xml .;\
	else\
	  echo ${0##*/}: Error, did not find ${CSGSHARE}/doc/cgoptions.xml >&2;\
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

clean:
	rm -f *.man *.tex *.wiki

