rmflags=-f

all: 
	$(MAKE) $(MFLAGS) -C reference 
	$(MAKE) $(MFLAGS) -C fig
	$(MAKE) $(MFLAGS) -C usage 
	./latexmk.pl -pdfdvi manual.tex

.SUFFIXES: .tex .pdf

.tex.pdf:
	./latexmk.pl -pdfdvi $*

clean:
	rm -f manual.fdb_latexmk
	latexmk -C manual.tex
	$(MAKE) $(MFLAGS) -C fig clean
	$(MAKE) $(MFLAGS) -C usage clean
	$(MAKE) $(MFLAGS) -C reference clean
	rm -f *~

