rmflags=-f
PICDIR=abb

all: 
	$(MAKE) $(MFLAGS) -C usage
	./latexmk.pl -pdfdvi manual.tex

.SUFFIXES: .tex .pdf

pics:
	$(MAKE) $(MFLAGS) -C $(PICDIR)

.tex.pdf:
	./latexmk.pl -pdfdvi $*

clean:
	rm -f manual.dvi manual.pdf
	latexmk -c
	$(MAKE) $(MFLAGS) -C usage clean

