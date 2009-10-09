all: wiki tex man 

.SUFFIXES: .xml .t2t .tex .man .wiki .html .sh

wiki: cgoptions.wiki mapping.wiki
man: cgoptions.man mapping.man
tex: cgoptions.tex mapping.tex
html: cgoptions.html mapping.html


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

