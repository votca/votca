rmflags=-f

all: 
	./bootstrap.sh
	$(MAKE) $(MFLAGS) -C reference 
	$(MAKE) $(MFLAGS) -C fig
	$(MAKE) $(MFLAGS) -C usage 
	./latexmk.pl -pdfdvi manual.tex

.SUFFIXES: .tex .pdf

.tex.pdf:
	./latexmk.pl -pdfdvi $*

clean:
	rm -f manual.dvi manual.pdf
	rm -f manual.fdb_latexmk
	latexmk -c
	$(MAKE) $(MFLAGS) -C fig clean
	$(MAKE) $(MFLAGS) -C usage clean
	$(MAKE) $(MFLAGS) -C reference clean

