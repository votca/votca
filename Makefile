rmflags=-f
PICDIR=abb

all: 
	$(MAKE) $(MFLAGS) -C reference 
	./latexmk.pl -pdfdvi manual.tex

.SUFFIXES: .tex .pdf

pics:
	$(MAKE) $(MFLAGS) -C $(PICDIR)

.tex.pdf:
	./latexmk.pl -pdfdvi $*

clean:
	rm -f manual.dvi manual.pdf
	rm -f manual.fdb_latexmk
	latexmk -c
	$(MAKE) $(MFLAGS) -C reference clean

