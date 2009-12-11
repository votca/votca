.SUFFIXES: .tex .pdf
HGID:=$(shell hg parents --template "{node|short}" | sed 's/.*/\\newcommand{\\hgid}{&}/')

all: hgid.tex fig_submake functionality_submake reference_submake usage_submake
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
	rm -f hgid.tex
	rm -f *~

hgid.tex: update_hgid
update_hgid:
	[ -f hgid.tex ] || touch hgid.tex
	echo '$(HGID)' | cmp -s hgid.tex - || echo '$(HGID)' > hgid.tex

