.SUFFIXES: .tex .pdf

all: fig_submake functionality_submake reference_submake usage_submake
	./latexmk.pl -pdfdvi manual.tex

.tex.pdf:
	./latexmk.pl -pdfdvi $*

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

clean: fig_subclean functionality_subclean reference_subclean usage_submake
	./latexmk.pl -C manual.tex
	rm -f manual.fdb_latexmk
	rm -f *~

