all: wiki tex man

.SUFFIXES: .xml .t2t .tex .man .wiki

wiki: cgoptions.wiki mapping.wiki
man: cgoptions.man mapping.man
tex: cgoptions.tex mapping.tex


.xml.t2t:
	./xmltotxt.sh $< $@

.t2t.man:
	./txt2tags -t man --infile $< --outfile $@

.t2t.wiki:
	./txt2tags -t wiki --infile $< --outfile $@

.t2t.tex:
	./txt2tags --no-header -t tex --infile $< --outfile $@

clean:
	rm -f *.man *.tex *.wiki

