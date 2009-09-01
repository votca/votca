NAME=manual
rmflags=-f
PICDIR=abb

all: $(NAME).pdf
pdf: $(NAME).pdf


.SUFFIXES: .tex .pdf

pics:
	$(MAKE) $(MFLAGS) -C $(PICDIR)

.tex.pdf:
	./latexmk.pl -pdfdvi $*

clean:
	latexmk -c

